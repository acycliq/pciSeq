# Real-time Viewer

A small, self-contained web viewer that shows pciSeq (VarBayes) progress while it runs.

The viewer consists of two parts:
- A lightweight Python server (Flask-SocketIO) that streams updates
- A static web page (HTML + JS) that renders cells using deck.gl

No external viewer is required.

## Components and Responsibilities

- app.py
  - Reads options and, if `realtime_viewer` is true, starts the viewer server.
  - Wires a callback so VarBayes can push updates after each iteration.

- RealtimeViewerServer (server.py)
  - Serves the static files `viewer.html` and `viewer.js`.
  - Streams two types of messages to the browser:
    - Geometry (sent once): centroids, radii, class names, mean cell radius (mcr)
    - Per-iteration classes: class index per cell and probability (prob)
  - Caches geometry and the last update so a new browser tab can catch up.

- viewer.html + viewer.js
  - Connects to the server via Socket.IO.
  - Receives geometry once, then receives class/prob arrays each iteration.
  - Renders cells with deck.gl; shows a minimal HUD and a compact legend.

## How It Fits Together

Sequence at a high level:

1) VarBayes runs one iteration
2) VarBayes calls the callback with: classProb, iteration, delta
3) Server extracts:
   - cell_classes = argmax(classProb, axis=1)
   - prob = max(classProb, axis=1)
   - geometry (centroids, radii) from VarBayes.cells (only once)
4) Server emits messages to the browser
5) Browser updates the visualization

ASCII diagram of the workflow:

  Python (VarBayes)               RealtimeViewerServer                 Browser (viewer.js)
  -------------------             ----------------------               --------------------
  main_loop() iter i  --->  send_update(classProb, i, delta)  --->  classes_update_* events
         |                          |                                   |
         |                    first iteration only                      |
         |------------------->  geometry_init_* events  --------------->|

## Message Protocol

Geometry (sent once to each client):
- geometry_init_begin
  - num_cells: int
  - chunk_size: int
  - class_names: list[str]
  - mcr: float (mean cell radius)
- geometry_init_chunk (repeated)
  - start, end: int
  - centroids_x: list[float]
  - centroids_y: list[float]
  - radii: list[float]
- geometry_init_end

Per-iteration updates:
- classes_update_begin
  - iteration: int
  - delta: float (convergence)
  - num_cells: int
  - chunk_size: int
- classes_update_chunk (repeated)
  - start, end: int
  - cell_classes: list[uint8]
  - prob: list[float]
- classes_update_end

### Geometry: what it is and why we send it once

The viewer needs a small amount of static information before it can draw anything. We call this the geometry. It contains:

- Centroids: x and y (and z in 3D) for every cell. These set the position of each circle.
- Radii: one radius per cell, derived from area. The client currently prefers to draw all circles with the mean cell radius (mcr) for a clean look, but the per‑cell radii are still sent for completeness.
- Class names: human‑readable labels for the legend and tooltips.
- Mean cell radius (mcr): a single value from `VarBayes.cells.mcr` that the client uses as the fixed display radius.
- Sizes: number of cells and a chunk size so the browser can preallocate arrays and stream the payload in parts.

Why send it only once:

- It does not change as the algorithm iterates. Positions stay the same.
- Keeping it out of the per‑iteration messages makes those messages small and fast to process.
- The array indices the browser uses match the original cell order. That lets the server send only two arrays per iteration: class index per cell and probability per cell.

Conventions:
- Background cell (index 0) is excluded on the server before streaming.
- Cell indices in the stream start at 0 for the first real cell.

## Configuration

Defaults:
- Port: 5001 (do not pass it unless overriding)
- Auto-open browser: off (open http://127.0.0.1:5001 manually)
- Fixed radius: the client uses `mcr` by default for circle size; server also sends per-cell radii

Common options (in `opts` passed to `fit`):
- `realtime_viewer`: bool (enable/disable)
- `realtime_viewer_port`: int (default 5001)
- `realtime_viewer_fixed_radius`: float or None

Call-site options (fit) quick reference:

- Enable viewer:
  - Pass `opts={'realtime_viewer': True}` to `fit(...)`.
  - Do not pass any callback; `fit` wires it internally.
- Override port (optional):
  - Pass `opts={'realtime_viewer': True, 'realtime_viewer_port': 5002}`.
  - Default is 5001 if omitted.
- Fixed radius (optional):
  - Pass `realtime_viewer_fixed_radius` to force a display radius.
  - The viewer uses `mcr` by default for a consistent circle size.

Example (automatic startup from app.py):

```
cellData, geneData = fit(
    spots=spots,
    coo=coo,
    scRNAseq=scRNAseq,
    opts={
        'realtime_viewer': True,
        # 'realtime_viewer_port': 5001,         # optional; defaults to 5001
        # 'realtime_viewer_fixed_radius': None, # optional
    }
)
```

Example (manual server control):

```
from pciSeq.src.realtime_viewer import RealtimeViewerServer

viewer = RealtimeViewerServer(port=5001, auto_open_browser=False)
viewer.start()

cellData, geneData = fit(
    spots=spots,
    coo=coo,
    scRNAseq=scRNAseq,
    opts={'realtime_viewer_callback': viewer.send_update}
)

viewer.stop()
```

## Rendering Details (viewer.js)

- Deck.gl ScatterplotLayer renders one circle per cell.
- Radius: uses `mcr` (mean cell radius) if available; otherwise uses per-cell radius derived from area.
- Color: class-based palette; opacity is a function of `prob` (higher prob = more opaque).
- HUD: shows Connected, Iteration, Convergence, Cells.
- Legend: class chips with counts; click to toggle visibility; filter box to search by name.

## File Structure

```
pciSeq/src/realtime_viewer/
  __init__.py            # exports RealtimeViewerServer
  server.py              # Flask-SocketIO streaming server
  README.md              # this document
  static/
    viewer.html          # UI layout and styles
    viewer.js            # Socket client + deck.gl renderer
```

## Troubleshooting

- Port already in use: start with another port, e.g. 5002.
- Multiple tabs: each tab is a client. Close stale tabs to reduce noise.
- No updates: check Python logs for exceptions and browser console for connection errors.

## Notes for Refactoring

- Server responsibilities are limited to:
  - serve static files
  - extract geometry once and cache it
  - stream per-iteration class/prob arrays
- Client owns all rendering and UI; protocol is documented above.
- The message names and payload shapes are the contract between server and client.
