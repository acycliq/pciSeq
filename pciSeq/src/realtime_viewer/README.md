# Real-Time Viewer Module

Self-contained web-based visualization for watching pciSeq algorithm convergence in real-time.

## Overview

This module provides a complete real-time visualization system for the VarBayes algorithm. It includes:
- **Flask-SocketIO server** for streaming updates
- **Web-based viewer** with interactive visualization
- **Automatic browser launch** for convenience
- **Zero dependency on external viewers**

All components are self-contained within the pciSeq package.

## Features

- ✅ **Self-contained**: Complete viewer included (HTML/CSS/JS)
- ✅ **Zero impact when disabled**: No overhead if not used
- ✅ **Auto-launch**: Opens browser automatically
- ✅ **Real-time updates**: See cell assignments change as algorithm iterates
- ✅ **Performance metrics**: Monitor iteration number, convergence delta
- ✅ **Cell class distribution**: Live legend with counts per class
- ✅ **Thread-safe**: Runs in background without blocking

## Quick Start

```python
from pciSeq.app import fit
from pciSeq.src.realtime_viewer import RealtimeViewerServer
import pandas as pd
from scipy.sparse import load_npz

# Load data
spots = pd.read_csv('spots.csv')
coo = load_npz('segmentation.coo.npz')
scRNAseq = pd.read_csv('scRNAseq.csv.gz', compression='gzip')

# Start viewer (automatically opens browser)
viewer = RealtimeViewerServer(port=5001)
viewer.start()

# Run analysis with real-time visualization
cellData, geneData = fit(
    spots=spots,
    coo=coo,
    scRNAseq=scRNAseq,
    opts={
        'realtime_viewer_callback': viewer.send_update,
        'max_iter': 100,
    }
)

# Cleanup
viewer.stop()
```

## What You'll See

The viewer opens automatically in your browser showing:
- **Grid visualization**: Each cell colored by its current class assignment
- **Status bar**: Connection status, iteration number, convergence delta, cell count
- **Live legend**: Cell classes with counts, updated each iteration
- **Real-time updates**: Colors change as algorithm reassigns cells

## Architecture

```
┌─────────────────────────────────────────┐
│  Python (pciSeq backend)                │
│                                         │
│  VarBayes.main_loop()                   │
│    ↓ each iteration                     │
│  viewer.send_update(classProb, i, δ)    │
│    ↓                                    │
│  Flask-SocketIO Server (port 5001)     │
│  - Serves static viewer HTML/JS        │
│  - Streams updates via WebSocket       │
└─────────────────┬───────────────────────┘
                  │ WebSocket
                  ↓
┌─────────────────────────────────────────┐
│  Browser (Auto-launched)                │
│                                         │
│  viewer.html (self-contained)           │
│  - Socket.IO client                     │
│  - Canvas visualization                 │
│  - Live status updates                  │
└─────────────────────────────────────────┘
```

## File Structure

```
pciSeq/src/realtime_viewer/
├── __init__.py          # Exports RealtimeViewerServer
├── server.py            # Flask-SocketIO server
├── README.md            # This file
└── static/
    ├── viewer.html      # Web viewer UI
    └── viewer.js        # Visualization logic
```

## Advanced Usage

### Disable Auto-launch

```python
viewer = RealtimeViewerServer(
    port=5001,
    auto_open_browser=False  # Don't open browser
)
viewer.start()
# Manually navigate to http://localhost:5001
```

### Context Manager

```python
with RealtimeViewerServer(port=5001) as viewer:
    cellData, geneData = fit(
        spots=spots,
        coo=coo,
        scRNAseq=scRNAseq,
        opts={'realtime_viewer_callback': viewer.send_update}
    )
# Automatically stops server on exit
```

### Custom Port

```python
viewer = RealtimeViewerServer(port=8080)
```

## Data Flow

1. **Algorithm iterates**: `VarBayes.main_loop()` runs
2. **Callback invoked**: After each iteration, calls `viewer.send_update(classProb, i, delta)`
3. **Data serialized**: Converts classProb to uint8 array (class indices)
4. **WebSocket emit**: Broadcasts to all connected clients
5. **Browser updates**: Redraws visualization with new assignments
6. **Legend refreshes**: Updates cell class counts

## Performance

- **Data size**: ~100-500 KB per iteration (for 100K cells)
- **Frequency**: Once per iteration (1-10 seconds typical)
- **Overhead**: < 0.1% of iteration time
- **Memory**: Negligible (no accumulation)
- **Network**: Local only (127.0.0.1)

## Data Format (WebSocket)

Each `iteration_update` event:

```javascript
{
    iteration: 42,                    // Current iteration
    cell_classes: [0, 1, 2, ..., 5],  // uint8: class per cell
    confidence: [0.95, 0.87, ...],    // float32: max prob per cell
    delta: 0.0234,                    // Convergence metric
    num_cells: 100000                 // Total cells
}
```

## Visualization

- **Grid layout**: Cells arranged in optimal grid (maintains aspect ratio)
- **Color mapping**: Each class gets distinct HSL color
- **Opacity**: Based on confidence (high confidence = more opaque)
- **Updates**: Smooth transitions as assignments change

## Requirements

- `flask-socketio` (installed with pciSeq)
- Modern browser with JavaScript enabled

## Troubleshooting

### Port Already in Use

```bash
# Check what's using the port
lsof -i :5001

# Use different port
viewer = RealtimeViewerServer(port=5002)
```

### Browser Doesn't Open

Set `auto_open_browser=False` and manually navigate to `http://localhost:5001`

### Connection Issues

- Check firewall settings
- Ensure no VPN/proxy blocking localhost
- Try different browser

### No Updates Appearing

- Verify `realtime_viewer_callback` is in opts
- Check Python console for errors
- Check browser console (F12) for WebSocket errors

## Example Script

See `example_realtime_viewer.py` in the repository root for a complete working example.

## Comparison with Main Viewer

| Feature | Real-Time Viewer | Main pciSeq_viewer |
|---------|-----------------|-------------------|
| Purpose | Watch algorithm run | Explore final results |
| When | During execution | After completion |
| Data | Live updates | Static files |
| Location | Inside pciSeq repo | Separate repo |
| Launch | Automatic | Manual |
| Performance | Minimal overhead | Full-featured |