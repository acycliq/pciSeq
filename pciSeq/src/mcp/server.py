"""The MCP server: the tools in tools.py, exposed over the protocol.

This is the only file in pciSeq that imports mcp. All the logic lives in tools.py,
this just registers each function and keeps track of which run is open, so the tools
stay usable from a notebook or the viewer with no mcp anywhere.

It is an optional extra, so a plain pciSeq install does not pull in a web server:

    pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"

That puts a `pciseq-mcp` command on the path. Point Claude Code, Claude Desktop or any
other MCP client at it over stdio, for example in ~/.claude.json:

    "mcpServers": {"pciSeq": {"type": "stdio", "command": "pciseq-mcp", "args": []}}

Use the full path to the command if the client does not share your shell's PATH, which
is usually the case for desktop apps. The docstrings below are what the agent reads to decide which tool to call,
so they say what the numbers mean, not how they are computed.

Why this is kept although the viewer no longer uses it (September 2026): the desktop
viewer has its own MCP server now, and its chat reads the run itself. But this is
the only way to ask about a finished run with no screen, on a remote or headless
machine such as a cluster over SSH, where the viewer cannot open. It can go once
the viewer has a no-window mode (bead cz1.9.7).
"""
import functools
import io
import json
from pathlib import Path
from typing import Annotated, Optional

from pydantic import Field

try:
    from mcp.server.mcpserver import Image, MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
except ModuleNotFoundError as e:  # pragma: no cover
    # the tools in tools.py work without any of this, only the server needs it
    raise ModuleNotFoundError(
        'the pciSeq MCP server needs the mcp library, which is an optional extra so\n'
        'that a plain install does not pull in a web server. Add it with:\n'
        '    pip install "pciSeq_3d[mcp]"\n'
        'or straight from github:\n'
        '    pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"'
    ) from e

from . import docs as _docs
from .tools import Run, calculate as _calculate, open_run as _open

from .server_instructions import INSTRUCTIONS

def _version():
    # the pciSeq version and commit, eg '0.0.66.dev0+9992e76d', so a client such as
    # the viewer can show which pciSeq it is talking to
    import pciSeq
    return '%s+%s' % (pciSeq.__version__, pciSeq.__commit__)


server = MCPServer(name='pciSeq', version=_version(), instructions=INSTRUCTIONS)

# the run the tools are answering about. One at a time, set by open_run.
_run: Optional[Run] = None


def _need_run() -> Run:
    if _run is None:
        raise RuntimeError('no run is open, call open_run with the run folder first')
    return _run


# the errors tools.py raises on purpose, with a message meant for whoever asked:
# no such cell, no such gene, an old run without containment data, no run open.
# mcp hides the message of any other exception from the client, so the agent would
# only see 'Error executing tool spots_in_cell' and never learn what to do about it.
# Re-raising as ToolError sends the message through.
_EXPECTED = (KeyError, ValueError, NotImplementedError, RuntimeError, FileNotFoundError)


def _tool(fn):
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except _EXPECTED as e:
            # KeyError repr()s its message, which puts quotes round it
            msg = e.args[0] if e.args else str(e)
            raise ToolError(str(msg)) from e
    return server.tool()(wrapped)


@_tool
def open_run(path: str) -> dict:
    """Open a finished pciSeq run so the other tools can answer questions about it.

    Give the run's output folder, or any folder above it that contains exactly one
    run. Returns a summary: how many cells, spots, genes and classes, the pciSeq
    version that made it, and whether it carries containment data (runs before
    September 2026 do not, and spots_in_cell will say so).
    """
    global _run
    _run = _open(path)
    # the run carries the documentation of the commit that fitted it; read that one
    # rather than whatever is on this machine (docs.use_run_pages). A run without
    # it gets no documentation at all, see _run_has_no_docs
    _docs.use_run_pages(_run.docs_pages())
    return _run.summary()


@_tool
def run_info() -> dict:
    """What produced the open run and how it ended: pciSeq version, commit and its
    date, when the run was made (run_date), python and package versions, the mean
    cell radius in pixels, the resolved settings (mrf_beta, rTheta, Inefficiency,
    nNeighbors, voxel_size and the rest), the number of iterations, and whether the
    loop converged, with a sentence saying so. Use it for 'when was this run made',
    'what is the mean cell radius', 'what settings did it use' and 'did it converge'.
    Older runs carry only the provenance and the answer says so."""
    return _need_run().run_info()


@_tool
def cell(label: int) -> dict:
    """The headline facts about one cell: its class probabilities, its top genes, its
    total counts and its scale factor theta.

    label is the cell's number in the segmentation, the one shown in the viewer.
    Counts are soft, each spot contributes its probability of belonging to the cell.
    """
    return _need_run().cell(label)


@_tool
def cell_shape(label: int) -> dict:
    """The fitted shape of one cell's Gaussian: round or elongated, by how much, and
    which way the long axis points, with the radii in pixels, each one against the
    mean cell radius, and the covariance behind them. Use it for 'is cell 18223
    round', 'which way is it elongated', 'how big is its Gaussian'. Runs saved
    before October 2026 do not keep the covariance and the answer says so."""
    return _need_run().cell_shape(label)


@_tool
def explain_cell(
        label: Annotated[int, Field(description='The cell label, as in the segmentation.')],
        vs_class: Annotated[Optional[str], Field(
            description='Class to compare against. Defaults to the runner up.')] = None,
        top_n: Annotated[int, Field(
            description='How many genes to list for each side, default 10.')] = 10) -> dict:
    """Why a cell was given its class, gene by gene.

    Compares the assigned class against another, the runner up unless vs_class is
    given: the gene log-likelihood, the class prior and the spatial term for each,
    the genes that pushed hardest for each side with the cell's count and the
    average count of that gene over the cells this run called each class
    (mean_in_assigned, mean_in_compared), the shared genes (the cell's biggest counts
    that both classes fit about equally, why these two were the finalists), and a
    narrative in plain words, with the evidence as odds rather than units. A gene the
    cell lacks can count against the class that expresses it. Counts are soft,
    weighted by assignment probability.

    Use this for questions like 'why is cell 2413 Ndnf Gaba' or 'why is cell 18223
    not CA1'. label is the cell number of the segmentation, the one the viewer shows.
    """
    return _need_run().explain_cell(label, vs_class=vs_class, top_n=top_n)


@_tool
def explain_spot(spot_id: Annotated[int, Field(description='The spot id.')]) -> dict:
    """Why a spot was assigned to the cell it was, term by term.

    One row per candidate cell plus the background: the Gaussian fit to the cell
    centroid, the class expression term, the cell scale, the cell-gene scale, the
    gene efficiency, the inside-cell bonus, their sum, and the resulting probability,
    and a narrative that tells the story in plain words with the evidence as odds.
    The background row carries the gene's misread density instead.

    A spot can sit outside every cell and still go to one, because the position is
    scored against the centroid, not the mask. spot_id is the id shown in the
    viewer, the same as the spot_id in geneData.
    """
    return _need_run().explain_spot(spot_id)


@_tool
def cell_counts(label: int, gene: Optional[str] = None) -> dict:
    """How many reads a cell holds, in total or for one gene.

    These are SOFT counts: each spot contributes its probability of belonging to the
    cell, so they are estimates and not whole numbers. For the number of spots
    physically inside the cell's boundary use spots_in_cell instead.
    """
    return _need_run().cell_counts(label, gene=gene)


@_tool
def spots_in_cell(label: int, gene: Optional[str] = None) -> dict:
    """How many spots physically sit inside a cell's segmentation mask.

    This is a HARD count with no probabilities: a spot either falls inside the mask
    or it does not. It is a different number from cell_counts, because the model can
    assign a spot to a cell it is not inside, and vice versa. Needs a run made after
    September 2026; older runs raise with a clear message.
    """
    return _need_run().spots_in_cell(label, gene=gene)


@_tool
def spots_of_cell(label: int, min_prob: Optional[float] = None, gene: Optional[str] = None) -> dict:
    """Which spots belong to a cell, with their probabilities, sorted highest first.

    Two definitions, and the answer says which it used. Without min_prob: the spots
    whose MOST LIKELY parent is this cell, the argmax. That is not a hard assignment,
    the lowest probability in the list can be well under 0.5. With min_prob, say
    0.0001: EVERY spot with probability above it on this cell, whose probabilities
    add up to the cell's soft counts. The second list is usually many times longer
    than the first.

    Use this for 'which spots are assigned to cell 18223' or 'list the spots of
    cell 18223 with their probabilities'. With gene, only that gene's spots, so
    'how many Plp1 spots are assigned to cell 18223' is answered by n_spots.
    """
    return _need_run().spots_of_cell(label, min_prob=min_prob, gene=gene)


@_tool
def gene(name: str) -> dict:
    """One gene across the run: its efficiency eta and inefficiency, its misread
    density, how many spots it has and how many were called misreads, its soft
    counts in cells split by class (soft, and summed over the cells called each
    class), and the ten cells holding most of it.

    Use it for 'what is the efficiency of Plp1', 'which classes express Plp1 in this
    run', 'which cells hold most Plp1'. Counts are soft.
    """
    return _need_run().gene(name)


@_tool
def theta(label: int) -> dict:
    """The cell scale factor theta of one cell: the overall value, and theta_bar
    under each of the top classes with the class probability. Theta scales a
    class's expected counts to the cell's total, with a Gamma(rTheta, rTheta) prior
    of mean 1; rTheta is returned too. label is the segmentation label.
    """
    return _need_run().theta(label)


@_tool
def gamma(label: int, gene: Optional[str] = None) -> dict:
    """The cell-gene scale factors gamma_bar of one cell under its assigned class,
    for every gene or for one, each with the cell's count of the gene and the count
    the class predicts for it (expected). Gamma is the per cell, per gene factor
    absorbing overdispersion; it is pulled toward 1, so compare counts with expected
    to say how far off a gene is. Only the assigned class is kept in diagnostics.db,
    and the answer says so. label is the segmentation label.
    """
    return _need_run().gamma(label, gene=gene)


@_tool
def spot(spot_id: int) -> dict:
    """One spot: its gene, position and plane, the cell it was assigned to with the
    probability, and every candidate cell with its class and probability. Lighter
    than explain_spot, which gives the terms behind each probability.

    Use it for 'which cell is spot 1642419 in', 'what gene is spot 1642419',
    'which cells was spot 1642419 scored against'.
    """
    return _need_run().spot(spot_id)


@_tool
def neighbours(label: int) -> dict:
    """The cells whose classes enter the spatial (mrf) term of one cell, nearest
    first, each with its class and probability and, when the viewer files are
    present, the centroid distance in xy pixels and the plane offset. mrf_beta is
    returned too. label is the segmentation label.
    """
    return _need_run().neighbours(label)


@_tool
def class_counts(min_counts: Optional[float] = None) -> dict:
    """How many cells each class has: hard (the number of cells whose most probable
    class it is) and soft (the class probability summed over the cells). Zero is
    listed first, the rest by size. min_counts leaves out cells with fewer soft
    counts in total.

    Use it for 'how many cells per class', 'how many cells are Zero', 'how many CA1
    cells with more than 40 reads'.
    """
    return _need_run().class_counts(min_counts=min_counts)


@_tool
def class_prior() -> dict:
    """The weight of every class at the end of the run: the prior probability that a
    cell is of the class, before any of its genes is looked at. Zero is listed
    first, the rest by weight. These are the weights the run saved and the class
    score uses, so quote them as they are and never rebuild them from class_counts.

    Use it for 'what weight did CA3 settle on', 'the class weights at convergence',
    'the class prior of this run'.
    """
    return _need_run().class_prior()


@_tool
def find_cells(class_name: Optional[str] = None, plane: Optional[int] = None,
               min_counts: Optional[float] = None, top_two_within: Optional[float] = None,
               n: int = 50) -> dict:
    """The cells matching the filters given: assigned class, plane of the centroid,
    minimum total counts, and top_two_within, the largest gap allowed between the
    probabilities of the top two classes (small values pick the uncertain cells).
    Returns the number matching and the first n by probability, each with its class,
    probability, runner up, margin and total counts.

    Use it for 'which cells are CA2 on plane 40', 'list the uncertain cells', 'the
    L5 ET cells with over 100 reads'. Labels are segmentation labels.
    """
    return _need_run().find_cells(class_name=class_name, plane=plane, min_counts=min_counts,
                                  top_two_within=top_two_within, n=n)


@_tool
def metadata(key: Optional[str] = None) -> dict:
    """The metadata table of diagnostics.db. Without a key: every key with the kind
    of value it holds. With a key: its value, parsed. Covers what has no tool of its
    own: sc_mean_expression (the reference expression, gene by class), log_prior,
    rho_bar, hard_misread_counts, gene_total_spots, config, run, label_map.
    """
    return _need_run().metadata(key=key)


@_tool
def cell_row(label: int) -> dict:
    """The cellData.tsv row of one cell, value for value: Cell_Num, X, Y, Z,
    X0, Y0, Z0, Genenames, CellGeneCount, spot_id, ClassName, Prob. X, Y, Z is the
    fitted centre and X0, Y0, Z0 the segmentation centroid. Use this when the
    question is about what the saved file says for a cell. The counts are soft."""
    return _need_run().cell_row(label)


@_tool
def spot_row(spot_id: int) -> dict:
    """The geneData.tsv row of one spot, value for value: gene, position, plane,
    neighbour, neighbour_array, neighbour_prob, omp_score, omp_intensity,
    is_hard_misread and, on runs from September 2026 on, inside_cell."""
    return _need_run().spot_row(spot_id)


@_tool
def cell_image(label: int, context: bool = False, plane: Optional[int] = None,
               width: int = 1200, channel: Optional[str] = None, neighbours: bool = False,
               save_as: Optional[str] = None, mbtiles: Optional[str] = None) -> list:
    """A picture of a cell on the background image (DAPI or another stain), for 'show
    me cell 18223', 'show me cell 18223 on the DAPI' or 'where is cell 18223 in the
    section'.

    context=False gives a close-up: the cell outlined in red, every other cell on
    that plane in blue, the nuclei underneath. context=True gives the whole section
    with a ring round the cell, to show where it sits in the tissue. Ask for both
    when the user wants to see a cell; that is the pair on the why-a-cell-got-its-type
    docs page.

    Which cells are outlined depends on the question. neighbours=False, every cell
    near it, is right for 'show me the cell' or 'is it segmented properly'.
    neighbours=True outlines only the cells the spatial term of the model listened
    to, which is what 'why did its neighbours make it this class' or 'show me its
    neighbours' needs; the answer lists them and says which sit on another plane.

    The image is stitched from the viewer's tile pyramid (the .mbtiles in
    viewer_data), so it is a close visual copy, not the raw pixels. One plane at a
    time: the cell's 3D neighbours on other planes are not drawn, which is why a
    cell can have fewer blue cells around it than neighbours. The plane is the
    centroid's when the run knows its voxel size, else the one where the cell is
    biggest; pass plane to choose. save_as writes the png to that path as well,
    which is how a user in a terminal gets to see it.

    A run can have more than one background image, for example DAPI and GCaMP, one
    .mbtiles each. channel picks one by its name ('DAPI', 'GCaMP'); part of the name
    is enough. If the user did not say which image and the run has several, do not
    guess: the tool refuses and lists them, so ask the user which one they want.
    If the run has only one background image it is always used, whatever channel
    says, because the file cannot tell which stain it is; background_note in the
    answer says so, pass that on to the user.
    """
    im, info = _need_run().cell_image(label, context=context, plane=plane, width=width,
                                      channel=channel, neighbours=neighbours,
                                      mbtiles=mbtiles)
    return _picture(im, info, save_as)


@_tool
def plane_image(plane: Optional[int] = None, bbox: Optional[list] = None,
                width: int = 1200, channel: Optional[str] = None,
                save_as: Optional[str] = None, mbtiles: Optional[str] = None) -> list:
    """The background image (DAPI or another stain) of one plane with nothing drawn
    on it, for 'show me the whole image', 'show me the DAPI of plane 54' or 'show me
    the region around x 5000 to 6000, y 500 to 1200'.

    With no bbox it is the whole plane, untrimmed. bbox is [x0, y0, x1, y1] in image
    pixels, the same coordinates as the cells and spots. The plane defaults to the
    middle of the stack. For a picture about one cell use cell_image instead, it
    finds the plane and the window for you and draws the outlines.

    Stitched from the viewer's tile pyramid (the .mbtiles in viewer_data), so it is
    a close visual copy, not the raw pixels. The answer says the scale, so a point
    of the image can be placed on the picture. save_as writes the png to that path
    as well, which is how a user in a terminal gets to see it.

    A run can have more than one background image, for example DAPI and GCaMP, one
    .mbtiles each. channel picks one by its name ('DAPI', 'GCaMP'); part of the name
    is enough. If the user did not say which image and the run has several, do not
    guess: the tool refuses and lists them, so ask the user which one they want.
    If the run has only one background image it is always used, whatever channel
    says, because the file cannot tell which stain it is; background_note in the
    answer says so, pass that on to the user.
    """
    im, info = _need_run().plane_image(plane=plane, bbox=bbox, width=width,
                                       channel=channel, mbtiles=mbtiles)
    return _picture(im, info, save_as)


def _picture(im, info, save_as):
    """What an image tool hands back: the png itself, so the agent can look at it,
    and the facts as text. Written to save_as too when asked."""
    buf = io.BytesIO()
    im.save(buf, format='PNG')
    if save_as:
        path = Path(save_as).expanduser()
        path.write_bytes(buf.getvalue())
        info['saved_as'] = str(path)
    return [Image(data=buf.getvalue(), format='png'), json.dumps(info, indent=2)]


@_tool
def list_source(dir: Annotated[str, Field(
        description='Folder path, for example "pciSeq/src/core"; "" for the package.')] = '') -> dict:
    """List a folder of the pciSeq source code, at the commit that made the open run.

    The code shown is the code that produced the run's numbers: the files on GitHub
    at that commit (needs internet). A run that does not record its commit is
    refused. Use it to find the file a question is about; the model code is under
    pciSeq/src/core.
    """
    from . import source as _source
    return _source.list_source_for_run(_need_run().commit(), dir)


@_tool
def read_source(path: Annotated[str, Field(
            description='For example "pciSeq/src/core/main.py". A path from the code '
                        'map page works as it is.')],
                start_line: Annotated[int, Field(description='First line to return, default 1.')] = 1) -> dict:
    """Read a file of the pciSeq source code, with line numbers, at the commit that
    made the open run, so it is the code that produced the numbers (see list_source
    for where it comes from). Reach for it only when a question needs the actual
    code, after the docs; never claim to run it. Long files come back in slices of
    up to 400 lines, give start_line to read on. When you cite it, give the file and
    the line.
    """
    from . import source as _source
    return _source.read_source_for_run(_need_run().commit(), path, start_line)


# What the agent is told for a run fitted before pciSeq saved its pages with the run.
# The pages on this machine belong to the pciSeq installed here, which may be newer
# than the one that made the run, so they are not offered in their place.
_NO_DOCS = ('this run has no documentation saved inside it (it was fitted before '
            'pciSeq saved its pages with the run), so there are no pages to read. Say '
            'so, and do not answer from memory as if from the documentation.')


def _run_has_no_docs():
    """True when a run is open and it carries no documentation pages."""
    return _run is not None and not _run.docs_pages()


@_tool
def docs(query: Annotated[str, Field(
             description='A few words, for example "rTheta" or "spatial term".')] = '',
         page: Annotated[str, Field(
             description='A page name from the list, for example '
                         '"the-model/scale-factors.md".')] = '',
         n: Annotated[int, Field(description='How many paragraphs, for a query; default 5.')] = 5) -> dict:
    """The pciSeq documentation, three ways.

    query searches it and returns the paragraphs holding those words, best first,
    each with its page and heading. page returns one whole page. Neither returns the
    list of all the pages, each with a line on what it is about.

    The search matches words literally, so when your words may not be the ones the
    docs use (an idea you can describe but cannot name, or a term from another
    method), ask for the list and pick the page by its subject, then read that page.
    Use this before answering any question about how pciSeq works, a term, or a
    setting such as rTheta, mrf_beta or Inefficiency, and name the page you took the
    answer from.
    """
    if _run_has_no_docs():
        return {'error': _NO_DOCS}

    def contents():
        return [{'page': p, 'title': _docs.page_title(_docs.read_page(p)),
                 'about': _docs.page_summary(_docs.read_page(p))} for p in _docs.list_pages()]

    if page:
        try:
            text = _docs.read_page(page)
        except KeyError:
            return {'error': 'no page %s; call docs with no query and no page for the list' % page,
                    'contents': contents()}
        return {'page': page, 'title': _docs.page_title(text), 'text': text,
                'docs_are': _docs.source()}
    if not query:
        return {'contents': contents(), 'docs_are': _docs.source(),
                'contents_are': 'every page with what it is about; read one with docs(page=...)'}
    hits = _docs.search_docs(query, n=n)
    for h in hits:
        h['resource'] = 'pciseq-docs://' + h['page']
    out = {'query': query, 'hits': hits, 'docs_are': _docs.source()}
    if not hits:
        # nothing matched the words: the list, so the page can be picked by subject
        out['contents'] = contents()
    return out


@_tool
def calculate(expression: Annotated[str, Field(
        description='Arithmetic, for example "15.29 + 9.8 + 3.1" or "exp(3.2)".')]) -> dict:
    """Do arithmetic instead of doing it in your head.

    Numbers, + - * / ** and brackets, and exp, log (natural), log10, sqrt, abs and
    round. Use it for a number the tools do not give, for example adding up a few
    gene differences, or exp(d) to turn a log-likelihood difference d into odds.
    Quote the result as calculate returned it. It makes a sum right, not meaningful:
    never use it to build a quantity the tools do not define, such as a ratio of two
    sums of log-likelihood differences, which is not odds.
    """
    return _calculate(expression)


@server.resource('pciseq-docs://index', name='pciSeq documentation index',
                 description='Every documentation page, with its title and what it is '
                             'about. Read a page as pciseq-docs://<page>. Use it when '
                             'the docs search finds nothing, since that search matches '
                             'words literally: pick the page by its subject instead.',
                 mime_type='text/markdown')
def docs_index() -> str:
    if _run_has_no_docs():
        return _NO_DOCS
    lines = ['# pciSeq documentation', '']
    for page in _docs.list_pages():
        text = _docs.read_page(page)
        lines.append('- pciseq-docs://%s  **%s**  %s'
                     % (page, _docs.page_title(text), _docs.page_summary(text)))
    return '\n'.join(lines) if len(lines) > 2 else 'no docs pages found on this machine'


# One resource per page, registered at import. A URI template such as
# pciseq-docs://{path} would be neater but matches a single path segment, so every
# page inside a folder (api/, the-model/, ...) would be unreachable. Static URIs
# also let the agent see every page with its title in list_resources.
def _page_reader(page):
    # a static resource's function must take no arguments, so the page is closed
    # over here rather than passed in
    def reader() -> str:
        if _run_has_no_docs():
            return _NO_DOCS
        return _docs.read_page(page)
    return reader


def _register_pages():
    for page in _docs.list_pages():
        title = _docs.page_title(_docs.read_page(page)) or page
        server.resource('pciseq-docs://' + page, name=title,
                        description='pciSeq documentation, %s' % page,
                        mime_type='text/markdown')(_page_reader(page))


_register_pages()


def main() -> None:
    """Run the server over stdio, which is what every MCP client does."""
    server.run(transport='stdio')


if __name__ == '__main__':
    main()
