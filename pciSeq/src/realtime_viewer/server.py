"""
Optional real-time viewer server for algorithm visualization.

This module provides a self-contained Flask-SocketIO server that streams
cell assignment updates during VarBayes algorithm execution.
"""

from flask import Flask, send_from_directory, request, jsonify
from flask_socketio import SocketIO
import numpy as np
import pandas as pd
import threading
import logging
import webbrowser
from pathlib import Path
import os

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

    def __init__(
        self,
        port=5001,
        host="127.0.0.1",
        auto_open_browser=True,
        max_cells: int = None,
        fixed_radius: float = None,
    ):
        self.port = port
        self.host = host
        self.auto_open_browser = auto_open_browser
        # Optional payload controls
        self.max_cells = max_cells  # if set, send only top-N cells (by confidence)
        self.fixed_radius = fixed_radius  # if set, send this radius for all cells
        self._varbayes_ref = None  # Will be set by app.py when callback is wired
        self._geometry_sent = False
        self._geometry_cache = None
        self._num_cells_expected = (
            None  # Track expected number of cells from first iteration
        )

        # Get static folder path (same directory as this file)
        self.static_folder = Path(__file__).parent / "static"

        self.app = Flask(__name__, static_folder=str(self.static_folder))
        # Increase max_http_buffer_size to accommodate larger iteration payloads
        # and keep threading async mode for compatibility without eventlet/gevent.
        self.socketio = SocketIO(
            self.app,
            cors_allowed_origins="*",
            async_mode="threading",
            max_http_buffer_size=20_000_000,  # allow up to ~20MB per message
            ping_interval=25,
            ping_timeout=60,
        )
        self.server_thread = None
        self._is_running = False
        self._setup_error_handlers()
        self._setup_routes()
        logger.info("RealtimeViewerServer initialized (not started yet)")

    def _setup_error_handlers(self):
        """Setup Flask error handlers to return JSON instead of HTML for API endpoints."""

        @self.app.errorhandler(404)
        def not_found(error):
            # Return JSON for API endpoints, HTML for others
            if request.path.startswith("/api/"):
                return jsonify({"error": "Endpoint not found"}), 404
            # For non-API routes, return normal 404
            return error

        @self.app.errorhandler(500)
        def internal_error(error):
            # Always return JSON for 500 errors to help with debugging
            logger.error(f"Internal server error: {error}", exc_info=True)
            return (
                jsonify({"error": "Internal server error", "details": str(error)}),
                500,
            )

        @self.app.errorhandler(Exception)
        def handle_exception(error):
            # Catch any unhandled exceptions and return JSON
            logger.error(f"Unhandled exception: {error}", exc_info=True)
            if request.path.startswith("/api/"):
                return (
                    jsonify({"error": "Internal server error", "details": str(error)}),
                    500,
                )
            raise error

    def _setup_routes(self):
        """Setup Flask routes for serving the viewer and handling connections."""

        @self.app.route("/")
        def index():
            """Serve the main viewer page."""
            return send_from_directory(self.static_folder, "viewer.html")

        @self.app.route("/<path:filename>")
        def serve_static(filename):
            """Serve static files (JS, CSS, etc.)."""
            return send_from_directory(self.static_folder, filename)

        @self.app.route("/health")
        def health():
            return {"status": "running", "port": self.port}

        @self.app.route("/api/start_job", methods=["POST"])
        def start_job():
            """
            API endpoint to start a pciSeq analysis job from the GUI.
            Receives configuration, validates inputs, and launches fit() in background.
            """
            try:
                config = request.json

                # Check if request has valid JSON body
                if config is None:
                    return jsonify({"error": "Request body must be valid JSON"}), 400

                logger.info(f"Received job start request with config: {config}")

                # Validate required fields
                required_fields = ["spots_path", "scrna_path", "coo_path"]
                for field in required_fields:
                    if not config.get(field):
                        return (
                            jsonify({"error": f"Missing required field: {field}"}),
                            400,
                        )

                # Get file paths from config (expand ~ but keep paths as provided)
                spots_path = os.path.expanduser(config["spots_path"])
                scrna_path = os.path.expanduser(config["scrna_path"])
                coo_path = os.path.expanduser(config["coo_path"])

                # Log the paths for debugging
                logger.info(f"Spots path: {spots_path}")
                logger.info(f"scRNA path: {scrna_path}")
                logger.info(f"Coo path: {coo_path}")

                # Check if files exist
                if not os.path.exists(spots_path):
                    return (
                        jsonify({"error": f"Spots file not found: {spots_path}"}),
                        400,
                    )
                if not os.path.exists(scrna_path):
                    return (
                        jsonify({"error": f"scRNAseq file not found: {scrna_path}"}),
                        400,
                    )
                if not os.path.exists(coo_path):
                    return (
                        jsonify({"error": f"Cell masks file not found: {coo_path}"}),
                        400,
                    )

                # Load data in background thread
                def run_job():
                    try:
                        logger.info("Starting pciSeq analysis job...")

                        # Import fit here to avoid circular imports
                        from pciSeq.app import fit
                        from scipy.sparse import coo_matrix

                        # Load data
                        logger.info(f"Loading spots from: {spots_path}")
                        spots = pd.read_csv(spots_path)

                        # ============================================================
                        # TODO: REMOVE THIS TEMPORARY FIX LATER!
                        # This renames columns to match expected schema
                        # Should be removed once data files are standardized
                        # ============================================================
                        spots = spots.rename(
                            columns={"z_stack": "z_plane", "Gene": "gene_name"}
                        )
                        # ============================================================

                        logger.info(f"Loading scRNAseq from: {scrna_path}")
                        scRNAseq = pd.read_csv(scrna_path)
                        if "Unnamed: 0" in scRNAseq.columns:
                            scRNAseq = scRNAseq.set_index("Unnamed: 0")

                        logger.info(f"Loading cell masks from: {coo_path}")
                        coo_data = np.load(coo_path, allow_pickle=True)
                        # Convert to list of sparse matrices if needed
                        if isinstance(coo_data, np.ndarray):
                            coo = [
                                coo_matrix(d) if not hasattr(d, "tocoo") else d
                                for d in coo_data
                            ]
                        else:
                            coo = coo_data

                        # Build opts dict
                        # IMPORTANT: Use realtime_viewer_callback instead of realtime_viewer: True
                        # to avoid creating a new server (we're already running one!)
                        opts = {
                            "Inefficiency": config.get("Inefficiency", 0.2),
                            "nNeighbors": config.get("nNeighbors", 6),
                            "CellCallTolerance": config.get("CellCallTolerance", 0.02),
                            "rSpot": config.get("rSpot", 2),
                            "max_iter": config.get("max_iter", 1000),
                            "MisreadDensity": config.get("MisreadDensity", 0.00001),
                            "voxel_size": config.get("voxel_size", [1, 1, 1]),
                            "output_path": config.get("output_path", "default"),
                            "save_data": config.get("save_data", True),
                            "remove_flat_cells": config.get("remove_flat_cells", True),
                            "launch_diagnostics": config.get(
                                "launch_diagnostics", False
                            ),
                            "launch_viewer": False,  # Don't launch separate viewer
                            # Use the existing server's callback (this server!)
                            "realtime_viewer_callback": self.send_update,
                        }

                        logger.info(f"Starting fit() with opts: {opts}")

                        # Run fit
                        cellData, geneData = fit(
                            spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts
                        )

                        logger.info("pciSeq analysis completed successfully")

                    except Exception as e:
                        logger.error(f"Error running pciSeq job: {e}", exc_info=True)

                # Start job in background thread
                job_thread = threading.Thread(target=run_job, daemon=True)
                job_thread.start()

                return (
                    jsonify(
                        {
                            "status": "started",
                            "message": "pciSeq analysis started successfully",
                        }
                    ),
                    200,
                )

            except Exception as e:
                logger.error(f"Error starting job: {e}", exc_info=True)
                return jsonify({"error": str(e)}), 500

        @self.socketio.on("connect")
        def handle_connect():
            logger.info("Client connected to realtime viewer")
            # Send cached geometry and last classes update if available
            try:
                if self._geometry_cache is not None:
                    geom = self._geometry_cache
                    n = geom["num_cells"]
                    chunk_size = geom.get("chunk_size", 5000)
                    class_names = geom.get("class_names", [])
                    # send geometry init in chunks
                    self.socketio.emit(
                        "geometry_init_begin",
                        {
                            "num_cells": int(n),
                            "chunk_size": int(chunk_size),
                            "class_names": class_names,
                            "mcr": geom.get("mcr"),
                            "is_3d": bool(geom.get("is_3d", False)),
                            "voxel_size": geom.get("voxel_size"),
                            "img_dim": geom.get("img_dim"),
                            "cell_call_tolerance": geom.get("cell_call_tolerance", 0.2),
                        },
                        namespace="/",
                    )
                    for start in range(0, n, chunk_size):
                        end = min(start + chunk_size, n)
                        self.socketio.emit(
                            "geometry_init_chunk",
                            {
                                "start": int(start),
                                "end": int(end),
                                "centroids_x": geom["centroids_x"][start:end],
                                "centroids_y": geom["centroids_y"][start:end],
                                "centroids_z": geom.get("centroids_z", [0] * n)[
                                    start:end
                                ],
                                "radii": geom["radii"][start:end],
                            },
                            namespace="/",
                        )
                    self.socketio.emit("geometry_init_end", {}, namespace="/")

                if hasattr(self, "_last_update") and self._last_update:
                    cached = self._last_update
                    n = cached.get("num_cells", 0)
                    chunk_size = cached.get("chunk_size", 5000)
                    self.socketio.emit(
                        "classes_update_begin",
                        {
                            "iteration": int(cached["iteration"]),
                            "delta": float(cached["delta"]),
                            "num_cells": int(n),
                            "chunk_size": int(chunk_size),
                        },
                        namespace="/",
                    )
                    for start in range(0, n, chunk_size):
                        end = min(start + chunk_size, n)
                        self.socketio.emit(
                            "classes_update_chunk",
                            {
                                "start": int(start),
                                "end": int(end),
                                "cell_classes": cached["cell_classes"][start:end],
                                "prob": cached["prob"][start:end],
                            },
                            namespace="/",
                        )
                    self.socketio.emit(
                        "classes_update_end",
                        {"iteration": int(cached["iteration"])},
                        namespace="/",
                    )
            except Exception as e:
                logger.warning(f"Failed to send cached state: {e}")

        @self.socketio.on("disconnect")
        def handle_disconnect():
            logger.info("Client disconnected from realtime viewer")

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
                allow_unsafe_werkzeug=True,  # Safe for local development
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
            if not hasattr(self, "_varbayes_ref"):
                logger.warning("VarBayes reference not set - cannot send spatial data")
                return

            varbayes = self._varbayes_ref

            # Extract argmax (assigned class per cell) - most efficient format
            cell_classes = np.argmax(cells_classProb, axis=1).astype(np.uint8)

            # Probability (max over classes per cell). Round to 3 decimals to reduce payload size
            prob = np.round(np.max(cells_classProb, axis=1).astype(np.float32), 3)

            # Extract spatial data (centroids and radii)
            # Centroids might be DataFrame or ndarray, handle both cases
            centroids = varbayes.cells.centroid
            if hasattr(centroids, "values"):
                centroids = centroids.values  # DataFrame

            # Round to 1 decimal to reduce payload size while preserving visual fidelity
            centroids_x = np.round(centroids[:, 0].astype(np.float32), 1)
            centroids_y = np.round(centroids[:, 1].astype(np.float32), 1)

            # Calculate radii from cell area (radius = sqrt(area / π))
            areas = varbayes.cells.ini_cell_props["area"]
            if hasattr(areas, "values"):
                areas = areas.values  # Series

            # Round to 1 decimal to reduce payload size
            radii = np.round(np.sqrt(areas / np.pi).astype(np.float32), 1)

            # Skip the first cell (index 0) which is the background
            cell_classes = cell_classes[1:]
            prob = prob[1:]
            centroids_x = centroids_x[1:]
            centroids_y = centroids_y[1:]
            radii = radii[1:]

            # logger.info(f"[ITERATION {iteration}] Skipped background cell (index 0), sending {len(cell_classes)} real cells")

            # Apply fixed radius if requested
            if self.fixed_radius is not None:
                radii = np.full_like(
                    centroids_x, float(self.fixed_radius), dtype=np.float32
                )

            # If limiting cells, select top-N by confidence
            if self.max_cells is not None and len(cell_classes) > self.max_cells:
                k = int(self.max_cells)
                # Use argpartition for efficiency, then sort those top-k indices by value desc
                idx_part = np.argpartition(prob, -k)[-k:]
                idx_sorted = idx_part[np.argsort(prob[idx_part])[::-1]]

                cell_classes = cell_classes[idx_sorted]
                prob = prob[idx_sorted]
                centroids_x = centroids_x[idx_sorted]
                centroids_y = centroids_y[idx_sorted]
                radii = radii[idx_sorted]

            num_cells = len(cell_classes)

            # Send geometry once per session and cache
            if not self._geometry_sent:
                chunk_size = 5000
                # Get class names from VarBayes
                class_names = (
                    varbayes.cells.class_names.tolist()
                    if hasattr(varbayes.cells.class_names, "tolist")
                    else list(varbayes.cells.class_names)
                )
                is3d = bool(varbayes.config.get("is3D", False))
                voxel_size = varbayes.config.get("voxel_size", None)

                # z centroids if present, otherwise zeros
                if centroids.shape[1] >= 3:
                    cz_full = np.round(centroids[:, 2].astype(np.float32), 3)
                    centroids_z = cz_full[1:]
                else:
                    centroids_z = np.zeros_like(centroids_x)

                self.socketio.emit(
                    "geometry_init_begin",
                    {
                        "num_cells": int(num_cells),
                        "chunk_size": int(chunk_size),
                        "class_names": class_names,
                        "mcr": float(varbayes.cells.mcr),
                        "is_3d": is3d,
                        "voxel_size": voxel_size,
                        "img_dim": varbayes.config.get("img_dim", None),
                        "cell_call_tolerance": float(
                            varbayes.config.get("CellCallTolerance", 0.2)
                        ),
                    },
                    namespace="/",
                )
                for start in range(0, num_cells, chunk_size):
                    end = min(start + chunk_size, num_cells)
                    self.socketio.emit(
                        "geometry_init_chunk",
                        {
                            "start": int(start),
                            "end": int(end),
                            "centroids_x": centroids_x[start:end].tolist(),
                            "centroids_y": centroids_y[start:end].tolist(),
                            "centroids_z": centroids_z[start:end].tolist(),
                            "radii": radii[start:end].tolist(),
                        },
                        namespace="/",
                    )
                self.socketio.emit("geometry_init_end", {}, namespace="/")
                self._geometry_sent = True
                self._geometry_cache = {
                    "num_cells": int(num_cells),
                    "chunk_size": int(chunk_size),
                    "centroids_x": centroids_x.tolist(),
                    "centroids_y": centroids_y.tolist(),
                    "centroids_z": centroids_z.tolist(),
                    "radii": radii.tolist(),
                    "class_names": class_names,
                    "mcr": float(varbayes.cells.mcr),
                    "is_3d": is3d,
                    "voxel_size": voxel_size,
                    "img_dim": varbayes.config.get("img_dim", None),
                    "cell_call_tolerance": float(
                        varbayes.config.get("CellCallTolerance", 0.2)
                    ),
                }
                self._num_cells_expected = num_cells
                logger.info(f"Geometry cached: {num_cells} cells")

            # CRITICAL FIX: Ensure num_cells matches geometry cache
            # If num_cells differs from the geometry we sent, we need to adjust the data
            if self._geometry_cache is not None:
                cached_num_cells = self._geometry_cache["num_cells"]
                if num_cells != cached_num_cells:
                    logger.error(
                        f"CRITICAL BUG DETECTED! Iteration {iteration}: num_cells={num_cells} but geometry_cache has {cached_num_cells} cells!"
                    )
                    logger.error(
                        "Attempting to fix by truncating/padding to match geometry..."
                    )

                    # Adjust arrays to match cached geometry size
                    if num_cells < cached_num_cells:
                        # cell_classes/prob are too small - pad with zeros
                        logger.warning(
                            f"Padding cell_classes from {num_cells} to {cached_num_cells}"
                        )
                        padded_classes = np.zeros(cached_num_cells, dtype=np.uint8)
                        padded_classes[:num_cells] = cell_classes
                        cell_classes = padded_classes

                        padded_prob = np.zeros(cached_num_cells, dtype=np.float32)
                        padded_prob[:num_cells] = prob
                        prob = padded_prob

                        num_cells = cached_num_cells
                    elif num_cells > cached_num_cells:
                        # cell_classes/prob are too large - truncate
                        logger.warning(
                            f"Truncating cell_classes from {num_cells} to {cached_num_cells}"
                        )
                        cell_classes = cell_classes[:cached_num_cells]
                        prob = prob[:cached_num_cells]
                        num_cells = cached_num_cells

            # Stream classes/prob each iteration (smaller payload)
            chunk_size = 5000 if num_cells > 10000 else num_cells

            self.socketio.emit(
                "classes_update_begin",
                {
                    "iteration": int(iteration),
                    "delta": float(delta),
                    "num_cells": int(num_cells),
                    "chunk_size": int(chunk_size),
                },
                namespace="/",
            )
            for start in range(0, num_cells, chunk_size):
                end = min(start + chunk_size, num_cells)
                self.socketio.emit(
                    "classes_update_chunk",
                    {
                        "start": int(start),
                        "end": int(end),
                        "cell_classes": cell_classes[start:end].tolist(),
                        "prob": prob[start:end].tolist(),
                    },
                    namespace="/",
                )
            self.socketio.emit(
                "classes_update_end", {"iteration": int(iteration)}, namespace="/"
            )

            # Cache last classes/prob for late-joining clients
            self._last_update = {
                "iteration": int(iteration),
                "delta": float(delta),
                "num_cells": int(num_cells),
                "chunk_size": int(chunk_size),
                "cell_classes": cell_classes.tolist(),
                "prob": prob.tolist(),
            }

            # logger.info(f"Sent update for iteration {iteration} to all clients (delta={delta:.6f}, cells={len(cell_classes)})")

        except Exception as e:
            logger.error(f"Failed to send update: {e}")

    def __enter__(self):
        """Context manager support."""
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager cleanup."""
        self.stop()
