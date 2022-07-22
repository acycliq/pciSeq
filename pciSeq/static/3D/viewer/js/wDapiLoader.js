var dapiWorker;

function loadImage(){
    dapiWorker = new Worker('./src/dapiloaderWorker.js')

    dapiWorker.onmessage = receiverWorkerMessage;

    // receives from onmessage
    var LoadersVolume = AMI.default.Loaders.Volume;
    var HelpersStack = AMI.default.Helpers.Stack;
    dapiWorker.postMessage({
        img: 'dapi_image-0068.nii.gz',
    })
}

function receiverWorkerMessage(event){
    console.log('in receiverWorkerMessage')
    console.log(event.data)
}