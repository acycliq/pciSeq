"""
Optional real-time viewer server for algorithm visualization.

This module provides a self-contained Flask-SocketIO server that streams
cell assignment updates during VarBayes algorithm execution.
"""
from flask import Flask, send_from_directory
from flask_socketio import SocketIO
import numpy as np
import threading
import logging
import webbrowser
import os
from pathlib import Path

logger = logging.getLogger(__name__)


class RealtimeViewerServer:
    """
    Optional real-time streaming server for cell assignment visualization.

    This server allows web clients to receive live updates of cell type
    assignments as the VarBayes algorithm iterates. It's completely optional
    and has zero impact on the algorithm when not used.

    Usage:
        # In your run script
        from pciSeq.src.realtime_viewer import RealtimeViewerServer

        viewer = RealtimeViewerServer(port=5001)
        viewer.start()

        # Pass as callback to fit() via opts
        cellData, geneData = fit(
            spots=spots,
            coo=coo,
            scRNAseq=scRNAseq,
            opts={
                'realtime_viewer_callback': viewer.send_update,
                'max_iter': 100,
            }
        )

        # Cleanup when done
        viewer.stop()

    Args:
        port (int): Port number for the server (default: 5001)
        host (str): Host address (default: '127.0.0.1')
    """

    def __init__(self, port=5001, host='127.0.0.1', auto_open_browser=True,
                 max_cells: int = None, fixed_radius: float = None):
        self.port = port
        self.host = host
        self.auto_open_browser = auto_open_browser
        # Optional payload controls
        self.max_cells = max_cells  # if set, send only top-N cells (by confidence)
        self.fixed_radius = fixed_radius  # if set, send this radius for all cells
        self._varbayes_ref = None  # Will be set by app.py when callback is wired
        self._geometry_sent = False
        self._geometry_cache = None
        self._num_cells_expected = None  # Track expected number of cells from first iteration

        # Get static folder path (same directory as this file)
        self.static_folder = Path(__file__).parent / 'static'

        self.app = Flask(__name__, static_folder=str(self.static_folder))
        # Increase max_http_buffer_size to accommodate larger iteration payloads
        # and keep threading async mode for compatibility without eventlet/gevent.
        self.socketio = SocketIO(
            self.app,
            cors_allowed_origins="*",
            async_mode='threading',
            max_http_buffer_size=20_000_000,  # allow up to ~20MB per message
            ping_interval=25,
            ping_timeout=60,
        )
        self.server_thread = None
        self._is_running = False
        self._setup_routes()
        logger.info(f"RealtimeViewerServer initialized (not started yet)")

    def _setup_routes(self):
        """Setup Flask routes for serving the viewer and handling connections."""
        @self.app.route('/')
        def index():
            """Serve the main viewer page."""
            return send_from_directory(self.static_folder, 'viewer.html')

        @self.app.route('/<path:filename>')
        def serve_static(filename):
            """Serve static files (JS, CSS, etc.)."""
            return send_from_directory(self.static_folder, filename)

        @self.app.route('/health')
        def health():
            return {'status': 'running', 'port': self.port}

        @self.socketio.on('connect')
        def handle_connect():
            logger.info('Client connected to realtime viewer')
            # Send cached geometry and last classes update if available
            try:
                if self._geometry_cache is not None:
                    geom = self._geometry_cache
                    n = geom['num_cells']
                    chunk_size = geom.get('chunk_size', 5000)
                    class_names = geom.get('class_names', [])
                    # send geometry init in chunks
                    self.socketio.emit('geometry_init_begin', {
                        'num_cells': int(n),
                        'chunk_size': int(chunk_size),
                        'class_names': class_names
                    }, namespace='/')
                    for start in range(0, n, chunk_size):
                        end = min(start + chunk_size, n)
                        self.socketio.emit('geometry_init_chunk', {
                            'start': int(start),
                            'end': int(end),
                            'centroids_x': geom['centroids_x'][start:end],
                            'centroids_y': geom['centroids_y'][start:end],
                            'radii': geom['radii'][start:end],
                        }, namespace='/')
                    self.socketio.emit('geometry_init_end', {}, namespace='/')

                if hasattr(self, '_last_update') and self._last_update:
                    cached = self._last_update
                    n = cached.get('num_cells', 0)
                    chunk_size = cached.get('chunk_size', 5000)
                    self.socketio.emit('classes_update_begin', {
                        'iteration': int(cached['iteration']),
                        'delta': float(cached['delta']),
                        'num_cells': int(n),
                        'chunk_size': int(chunk_size)
                    }, namespace='/')
                    for start in range(0, n, chunk_size):
                        end = min(start + chunk_size, n)
                        self.socketio.emit('classes_update_chunk', {
                            'start': int(start),
                            'end': int(end),
                            'cell_classes': cached['cell_classes'][start:end],
                            'confidence': cached['confidence'][start:end],
                        }, namespace='/')
                    self.socketio.emit('classes_update_end', {'iteration': int(cached['iteration'])}, namespace='/')
            except Exception as e:
                logger.warning(f"Failed to send cached state: {e}")

        @self.socketio.on('disconnect')
        def handle_disconnect():
            logger.info('Client disconnected from realtime viewer')

    def start(self):
        """Start server in background thread and optionally open browser."""
        if self._is_running:
            logger.warning("Server already running")
            return

        self.server_thread = threading.Thread(
            target=lambda: self.socketio.run(
                self.app,
                host=self.host,
                port=self.port,
                debug=False,
                use_reloader=False,
                allow_unsafe_werkzeug=True  # Safe for local development
            )
        )
        self.server_thread.daemon = True
        self.server_thread.start()
        self._is_running = True

        url = f"http://{self.host}:{self.port}"
        logger.info(f"Realtime viewer server started at {url}")

        # Open browser after short delay to ensure server is ready
        if self.auto_open_browser:
            threading.Timer(0.5, lambda: webbrowser.open(url)).start()
            logger.info(f"Opening browser at {url}")

    def stop(self):
        """Stop server."""
        if self._is_running:
            # logger.info("Realtime viewer server stopped")
            self._is_running = False

    def send_update(self, cells_classProb, iteration, delta):
        """
        Callback to send updates during algorithm execution.

        This method is called from VarBayes.main_loop() on each iteration
        to broadcast cell assignment updates to connected web clients.

        Args:
            cells_classProb: np.ndarray, shape (nC, nK)
                Cell type probability matrix where element [c, k] is the
                probability that cell c belongs to type k
            iteration: int
                Current iteration number
            delta: float
                Convergence metric (mean probability change)
        """
        if not self._is_running:
            return

        try:
            # Get VarBayes instance from the callback context
            # We need to pass the VarBayes instance to access cells data
            # For now, we'll store it as an instance variable
            if not hasattr(self, '_varbayes_ref'):
                logger.warning("VarBayes reference not set - cannot send spatial data")
                return

            varbayes = self._varbayes_ref

            # Extract argmax (assigned class per cell) - most efficient format
            cell_classes = np.argmax(cells_classProb, axis=1).astype(np.uint8)

            # Optional: include confidence (max probability per cell)
            # Round to 3 decimals to reduce payload size
            confidence = np.round(np.max(cells_classProb, axis=1).astype(np.float32), 3)

            # Extract spatial data (centroids and radii)
            # Centroids might be DataFrame or ndarray, handle both cases
            centroids = varbayes.cells.centroid
            if hasattr(centroids, 'values'):
                centroids = centroids.values  # DataFrame

            # Round to 1 decimal to reduce payload size while preserving visual fidelity
            centroids_x = np.round(centroids[:, 0].astype(np.float32), 1)
            centroids_y = np.round(centroids[:, 1].astype(np.float32), 1)

            # Calculate radii from cell area (radius = sqrt(area / π))
            areas = varbayes.cells.ini_cell_props['area']
            if hasattr(areas, 'values'):
                areas = areas.values  # Series

            # Round to 1 decimal to reduce payload size
            radii = np.round(np.sqrt(areas / np.pi).astype(np.float32), 1)

            # Skip the first cell (index 0) which is the background
            cell_classes = cell_classes[1:]
            confidence = confidence[1:]
            centroids_x = centroids_x[1:]
            centroids_y = centroids_y[1:]
            radii = radii[1:]

            logger.info(f"[ITERATION {iteration}] Skipped background cell (index 0), sending {len(cell_classes)} real cells")

            # Apply fixed radius if requested
            if self.fixed_radius is not None:
                radii = np.full_like(centroids_x, float(self.fixed_radius), dtype=np.float32)

            # If limiting cells, select top-N by confidence
            if self.max_cells is not None and len(cell_classes) > self.max_cells:
                k = int(self.max_cells)
                # Use argpartition for efficiency, then sort those top-k indices by value desc
                idx_part = np.argpartition(confidence, -k)[-k:]
                idx_sorted = idx_part[np.argsort(confidence[idx_part])[::-1]]

                cell_classes = cell_classes[idx_sorted]
                confidence = confidence[idx_sorted]
                centroids_x = centroids_x[idx_sorted]
                centroids_y = centroids_y[idx_sorted]
                radii = radii[idx_sorted]

            num_cells = len(cell_classes)

            # Send geometry once per session and cache
            if not self._geometry_sent:
                chunk_size = 5000
                # Get class names from VarBayes
                class_names = varbayes.cells.class_names.tolist() if hasattr(varbayes.cells.class_names, 'tolist') else list(varbayes.cells.class_names)

                self.socketio.emit('geometry_init_begin', {
                    'num_cells': int(num_cells),
                    'chunk_size': int(chunk_size),
                    'class_names': class_names
                }, namespace='/')
                for start in range(0, num_cells, chunk_size):
                    end = min(start + chunk_size, num_cells)
                    self.socketio.emit('geometry_init_chunk', {
                        'start': int(start),
                        'end': int(end),
                        'centroids_x': centroids_x[start:end].tolist(),
                        'centroids_y': centroids_y[start:end].tolist(),
                        'radii': radii[start:end].tolist(),
                    }, namespace='/')
                self.socketio.emit('geometry_init_end', {}, namespace='/')
                self._geometry_sent = True
                self._geometry_cache = {
                    'num_cells': int(num_cells),
                    'chunk_size': int(chunk_size),
                    'centroids_x': centroids_x.tolist(),
                    'centroids_y': centroids_y.tolist(),
                    'radii': radii.tolist(),
                    'class_names': class_names,
                }
                self._num_cells_expected = num_cells
                logger.info(f"Geometry cached: {num_cells} cells")

            # CRITICAL FIX: Ensure num_cells matches geometry cache
            # If num_cells differs from the geometry we sent, we need to adjust the data
            if self._geometry_cache is not None:
                cached_num_cells = self._geometry_cache['num_cells']
                if num_cells != cached_num_cells:
                    logger.error(f"CRITICAL BUG DETECTED! Iteration {iteration}: num_cells={num_cells} but geometry_cache has {cached_num_cells} cells!")
                    logger.error(f"Attempting to fix by truncating/padding to match geometry...")

                    # Adjust arrays to match cached geometry size
                    if num_cells < cached_num_cells:
                        # cell_classes/confidence are too small - pad with zeros
                        logger.warning(f"Padding cell_classes from {num_cells} to {cached_num_cells}")
                        padded_classes = np.zeros(cached_num_cells, dtype=np.uint8)
                        padded_classes[:num_cells] = cell_classes
                        cell_classes = padded_classes

                        padded_confidence = np.zeros(cached_num_cells, dtype=np.float32)
                        padded_confidence[:num_cells] = confidence
                        confidence = padded_confidence

                        num_cells = cached_num_cells
                    elif num_cells > cached_num_cells:
                        # cell_classes/confidence are too large - truncate
                        logger.warning(f"Truncating cell_classes from {num_cells} to {cached_num_cells}")
                        cell_classes = cell_classes[:cached_num_cells]
                        confidence = confidence[:cached_num_cells]
                        num_cells = cached_num_cells

            # Stream classes/confidence each iteration (smaller payload)
            chunk_size = 5000 if num_cells > 10000 else num_cells

            self.socketio.emit('classes_update_begin', {
                'iteration': int(iteration),
                'delta': float(delta),
                'num_cells': int(num_cells),
                'chunk_size': int(chunk_size)
            }, namespace='/')
            for start in range(0, num_cells, chunk_size):
                end = min(start + chunk_size, num_cells)
                self.socketio.emit('classes_update_chunk', {
                    'start': int(start),
                    'end': int(end),
                    'cell_classes': cell_classes[start:end].tolist(),
                    'confidence': confidence[start:end].tolist(),
                }, namespace='/')
            self.socketio.emit('classes_update_end', {'iteration': int(iteration)}, namespace='/')

            # Cache last classes/confidence for late-joining clients
            self._last_update = {
                'iteration': int(iteration),
                'delta': float(delta),
                'num_cells': int(num_cells),
                'chunk_size': int(chunk_size),
                'cell_classes': cell_classes.tolist(),
                'confidence': confidence.tolist(),
            }

            logger.info(f"Sent update for iteration {iteration} to all clients (delta={delta:.6f}, cells={len(cell_classes)})")

        except Exception as e:
            logger.error(f"Failed to send update: {e}")

    def __enter__(self):
        """Context manager support."""
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager cleanup."""
        self.stop()
