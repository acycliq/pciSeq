function json_parse(d) {
    if (d === '\r') {
        // Add the reasoning for this ASAP cause I will forget it in a month's time
        return null
    }
    try {
        return JSON.parse(d.replace(/'/g, '"'))
    } catch {
        return d
    }
}

const getRootName = (d) => d.split('/').pop().split('_')[0];
const count = (d) => {
    var out = {};
    var key;
    d.forEach(function (i) {
        key = getRootName(i);
        out[key] = (out[key] || 0) + 1;
    });
    return out
};

// a TSV parser that parses the data incrementally in chunks
const tsvChunkedParser = () => {
    const textDecoder = new TextDecoder("utf-8");
    let columnHeadings;
    let previousChunk = "";

    return {
        parseChunk(chunk) {
            // decode and split into lines
            const textData = previousChunk + textDecoder.decode(chunk);
            const lines = textData.split("\n");

            // the first line is our column headings
            if (!columnHeadings) {
                columnHeadings = lines[0].split("\t");
                lines.shift();
            }
            // the last line is probably partial - so append to the next chunk
            previousChunk = lines.pop();

            // convert each row to an object
            const items = lines
                .map(row => {
                    const cells = row.split("\t");
                    if (cells.length !== columnHeadings.length) {
                        return null;
                    }
                    let rowValue = {};
                    columnHeadings.forEach((h, i) => {
                        rowValue[h.trim()] = json_parse(cells[i]); //I am using trim() because I have no idea why I got whitespace at the end!!
                    });
                    return rowValue;
                })
                .filter(i => i);

            return items;
        }
    };
};


const fetchExternalData = (data) => {
    // the input data is actually the workpackage
    var filenames = data.map(d => d.download_url);
    return Promise.all(
        // filenames.forEach(d => d.map(el => fetch(el)))
        filenames.map(d => fetch(d))
    )
};


onmessage = async function (event) {
    //onmessage: receives messages from the UI thread

    var totalBytes = Array(event.data.length).fill(0);
    var perc = Array(event.data.length).fill(0);

    const tsvParser = [];
    for (var i=0; i< event.data.length; i++){
        // make an array where each element is a **NEW** instance of the parser
        tsvParser.push(tsvChunkedParser())
    }

    var results = await fetchExternalData(event.data);

    // if (!response.body) {
    //   throw Error("ReadableStream not yet supported in this browser.");
    // }

    const streamResponses = (responses) => {
        return Promise.all(
            responses.map((d,i) => streamedResponse(d, i).text())
        )
    };


    function streamedResponse(my_response, i) {
        return new Response(
            new ReadableStream({
                start(controller) {
                    const reader = my_response.body.getReader();

                    const read = async () => {
                        const {done, value} = await reader.read();
                        if (done) {
                            console.log('ok, ' + i + ' is done')
                            controller.close();
                            return;
                        }

                        const items = tsvParser[i].parseChunk(value);


                        totalBytes[i] += value.byteLength;
                        // console.log('File num: ' + i)
                        // var len = my_response.headers.get('content-length')
                        // perc[i] = totalBytes[i]/my_response.headers.get('content-length');
                        postMessage({i, items, url: my_response.url, bytes_streamed: +value.byteLength});

                        controller.enqueue(value);
                        read();
                    };

                    read();
                }
            })
        );
    };


    myData = await streamResponses(results);


    // call postMessage to send a message back to the UI thread
    postMessage({ items: [], totalBytes: [myData.length, null], finished: true });

};
