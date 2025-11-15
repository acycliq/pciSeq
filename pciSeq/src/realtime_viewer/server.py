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
from pathlib import Path
from pciSeq._version import __version__

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
        self._setup_routes()
        logger.info("RealtimeViewerServer initialized (not started yet)")

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
                            "version": __version__,
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
                                "cell_ids": geom.get("cell_ids", list(range(n)))[start:end],
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

        @self.socketio.on("request_check_cell")
        def handle_check_cell_request(data):
            """Handle check_cell diagnostic request from client."""
            logger.info(f"Received check_cell request: {data}")

            try:
                cell_label = data.get("cell_label")
                comparison_class = data.get("comparison_class", "Zero")

                if cell_label is None:
                    self.socketio.emit("check_cell_result", {
                        "error": "Missing cell_label parameter"
                    }, namespace="/")
                    return

                # Get VarBayes instance
                if not self._varbayes_ref:
                    self.socketio.emit("check_cell_result", {
                        "error": "VarBayes instance not available"
                    }, namespace="/")
                    return

                # Import check_cell function
                from pciSeq.src.core.utils import ops_utils

                # IMPORTANT: The viewer now sends original_label directly as cell.id
                # (We map seq_idx -> original_label in send_update and send it to viewer)
                # So cell_label IS the original_label - no mapping needed!
                original_label = cell_label

                logger.info(f"Viewer sent original_label: {original_label}")

                # Call check_cell with the original label (it will handle the mapping internally)
                gene_data, contr_df, _ = ops_utils.check_cell(
                    self._varbayes_ref,
                    original_label,
                    comparison_class,
                    top_n=10,
                    show_plot=False
                )

                # Get pciSeq's internal index (seq_idx) for this cell to extract the assigned class
                label_map = self._varbayes_ref.config.get('label_map')
                if label_map is not None:
                    seq_idx = label_map[original_label]
                else:
                    seq_idx = original_label

                pciseq_class = self._varbayes_ref.cells.class_names[
                    self._varbayes_ref.cells.classProb[seq_idx].argmax()
                ]

                # Prepare data for JSON serialization
                top_genes = []
                bottom_genes = []

                if 'diff' in contr_df.columns:
                    # Top genes (positive diff - favor pciSeq class)
                    top_sorted = contr_df.nlargest(10, 'diff')
                    for gene_name, row in top_sorted.iterrows():
                        top_genes.append({
                            "gene": str(gene_name),
                            "value": float(row['diff'])
                        })

                    # Bottom genes (negative diff - favor user class)
                    bottom_sorted = contr_df.nsmallest(10, 'diff')
                    for gene_name, row in bottom_sorted.iterrows():
                        bottom_genes.append({
                            "gene": str(gene_name),
                            "value": float(row['diff'])
                        })

                # Calculate sums
                top_sum = sum(g['value'] for g in top_genes)
                bottom_sum = sum(g['value'] for g in bottom_genes)

                # Prepare gene expression data table
                gene_table_data = []
                if gene_data is not None and not gene_data.empty:
                    for gene_name, row in gene_data.iterrows():
                        gene_table_data.append({
                            "gene": str(gene_name),
                            "mean_expr_pciseq": float(row.get(pciseq_class, 0)) if pciseq_class in row else 0,
                            "mean_expr_user": float(row.get(comparison_class, 0)) if comparison_class in row else 0,
                            "gene_count": int(row.get('gene count', 0)) if 'gene count' in row else 0
                        })

                # Send response
                self.socketio.emit("check_cell_result", {
                    "cell_label": int(cell_label),
                    "pciseq_class": str(pciseq_class),
                    "user_class": str(comparison_class),
                    "top_genes": top_genes,
                    "bottom_genes": bottom_genes,
                    "top_sum": float(top_sum),
                    "bottom_sum": float(bottom_sum),
                    "gene_expression_data": gene_table_data
                }, namespace="/")

            except Exception as e:
                logger.error(f"Error in check_cell handler: {e}", exc_info=True)
                self.socketio.emit("check_cell_result", {
                    "error": str(e)
                }, namespace="/")

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

        IMPORTANT - Cell Label Mapping processing:
        ========================================
        The user provides a segmentation image with cell labels (original_label).
        These labels may be non-sequential (e.g., 0, 5, 12, 47, 100...).

        During preprocessing (pciSeq.fit), some cells may be removed (e.g., flat cells
        on single planes), creating gaps in the labeling. To ensure sequential indexing,
        preprocessing relabels cells to sequential indices (seq_idx: 0, 1, 2, 3...).

        If relabeling occurs:
        - label_map dict is created: {original_label: seq_idx}
        - Example: {0: 0, 5: 1, 12: 2, 47: 3, 100: 4, ...}
        - Background is always 0 in both references: {0: 0, ...}

        If no relabeling occurs (labels were already sequential):
        - label_map = None
        - original_label == seq_idx for all cells

        Internal arrays (like cells_classProb):
        - Use seq_idx for indexing
        - Row index in cells_classProb[seq_idx] corresponds to seq_idx
        - cells_classProb.shape = (nC, nK) where nC = total number of cells

        For the viewer:
        - We send original_label as cell.id (not seq_idx)
        - This way user sees the same labels as in their segmentation
        - When user Ctrl+Clicks a cell, viewer sends original_label to check_cell()
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

            # Create cell_ids array: map seq_idx -> original_label
            # Row index in cells_classProb = seq_idx (0, 1, 2, ...)
            # We map these to original_label for the viewer
            nC = cells_classProb.shape[0]  # Total number of cells including background
            label_map = varbayes.config.get('label_map')

            if label_map is not None:
                # Reverse map: seq_idx -> original_label
                reverse_map = {v: k for k, v in label_map.items()}
                cell_ids = np.array([reverse_map[seq_idx] for seq_idx in range(nC)], dtype=np.int32)
            else:
                # No relabeling occurred, original_label == seq_idx
                cell_ids = np.arange(nC, dtype=np.int32)

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

            # Get z centroids BEFORE filtering (so it can be filtered too if needed)
            if centroids.shape[1] >= 3:
                centroids_z = np.round(centroids[:, 2].astype(np.float32), 3)
            else:
                centroids_z = np.zeros_like(centroids_x)

            # Skip background (index 0) - it's not a real cell and its centroid is NaN
            # which breaks JSON serialization. Slice all arrays consistently [1:] to keep alignment.
            cell_ids = cell_ids[1:]
            cell_classes = cell_classes[1:]
            prob = prob[1:]
            centroids_x = centroids_x[1:]
            centroids_y = centroids_y[1:]
            centroids_z = centroids_z[1:]
            radii = radii[1:]

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

                cell_ids = cell_ids[idx_sorted]
                cell_classes = cell_classes[idx_sorted]
                prob = prob[idx_sorted]
                centroids_x = centroids_x[idx_sorted]
                centroids_y = centroids_y[idx_sorted]
                centroids_z = centroids_z[idx_sorted]
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
                        "version": __version__,
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
                            "cell_ids": [int(x) for x in cell_ids[start:end]],
                            "centroids_x": [float(x) for x in centroids_x[start:end]],
                            "centroids_y": [float(x) for x in centroids_y[start:end]],
                            "centroids_z": [float(x) for x in centroids_z[start:end]],
                            "radii": [float(x) for x in radii[start:end]],
                        },
                        namespace="/",
                    )
                self.socketio.emit("geometry_init_end", {}, namespace="/")
                self._geometry_sent = True
                self._geometry_cache = {
                    "num_cells": int(num_cells),
                    "chunk_size": int(chunk_size),
                    "cell_ids": cell_ids.tolist(),
                    "centroids_x": centroids_x.tolist(),
                    "centroids_y": centroids_y.tolist(),
                    "centroids_z": centroids_z.tolist(),
                    "radii": radii.tolist(),
                    "class_names": class_names,
                    "mcr": float(varbayes.cells.mcr),
                    "is_3d": is3d,
                    "voxel_size": voxel_size,
                    "img_dim": varbayes.config.get("img_dim", None),
                    "version": __version__,
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
