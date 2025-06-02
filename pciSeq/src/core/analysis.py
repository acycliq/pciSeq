"""
Analysis tools for VarBayes algorithm results.

This module provides diagnostic and analysis capabilities for examining the results
of the VarBayes spatial transcriptomics algorithm.

The CellAnalyzer class provides methods for:
    - Analyzing gene expression patterns within cells
    - Examining spatial relationships between spots and cells
    - Calculating likelihood contributions for cell type assignments
    - Visualizing results through an interactive dashboard

"""
# Standard library imports
import json
import logging
import os
import shutil
import time
import webbrowser
import http.server
import socketserver
import threading
import random
from typing import Dict, Optional, List, Any

# Third-party imports
import numpy as np
import pandas as pd
from scipy import spatial

# Local imports
from ...src.core.utils.ops_utils import negative_binomial_loglikelihood
from ...src.core.utils.io_utils import get_out_dir
from ...src.viewer.utils import get_pciSeq_install_dir

# Configure logging
analysis_logger = logging.getLogger(__name__)


# class CellExplorer:
#     """Provides analysis and diagnostic tools for VarBayes results."""
#
#     def __init__(self, var_bayes):
#         """
#         Initialize with a VarBayes instance.
#
#         Args:
#             var_bayes: VarBayes instance to analyze
#         """
#         self.vb = var_bayes
#
#     def counts_within_radius(self, r: float) -> pd.DataFrame:
#         """
#         Calculates gene counts within specified radius of cell centroids.
#
#         Args:
#             r: Radius for counting spots (same units as spot coordinates)
#
#         Returns:
#             pd.DataFrame: DataFrame containing:
#                 - cell_label: Cell identifier
#                 - cell_type: Assigned cell type
#                 - cell_label_old: Original cell label (if relabeled)
#                 - Gene counts columns for each gene
#
#         Note:
#             Excludes background (label=0) from the output.
#         """
#         assert self.vb.cells.ini_cell_props['cell_label'][0] == 0
#
#         gene_names, gene_id = np.unique(self.vb.spots.data.gene_name.values, return_inverse=True)
#         spots = self.vb.spots.data.assign(gene_id=gene_id)
#
#         xy_coords = spots[['x', 'y']].values
#         point_tree = spatial.cKDTree(xy_coords)
#         nearby_spots = point_tree.query_ball_point(self.vb.cells.centroid, r)
#
#         out = np.zeros([self.vb.cells.centroid.shape[0], len(gene_names)])
#
#         for i, d in enumerate(nearby_spots):
#             t = spots.gene_id[d]
#             b = np.bincount(t)
#             out[i, :] = np.pad(b, (0, len(gene_names) - len(b)), 'constant')
#
#         cell_type = []
#         for i, d in enumerate(self.vb.cells.classProb):
#             j = self.vb.cells.class_names[np.argmax(d)]
#             cell_type.append(j)
#
#         temp = pd.DataFrame({
#             'cell_label': self.vb.cells.ini_cell_props['cell_label'],
#             'cell_type': cell_type
#         })
#
#         if 'cell_label_old' in self.vb.cells.ini_cell_props.keys():
#             temp['cell_label_old'] = self.vb.cells.ini_cell_props['cell_label_old']
#             temp = temp[['cell_label', 'cell_label_old', 'cell_type']]
#
#         df = pd.DataFrame(out, index=self.vb.cells.centroid.index, columns=gene_names)
#         df = pd.merge(temp, df, right_index=True, left_on='cell_label', how='left')
#
#         return df.iloc[1:, :]
#
#     def gene_loglik_contributions(self, label: int, user_class: Optional[str] = None) -> Dict:
#         """
#         Calculate gene log-likelihood contributions for a specified cell.
#
#         Args:
#             cell_num: Cell number to analyze (0 to nC-1)
#             user_class: Optional cell class to compare against assigned class
#
#         Returns:
#             Dict containing analysis results:
#                 - assigned_class: Automatically assigned cell class
#                 - user_class: User-specified class for comparison
#                 - assigned_contr: Log-likelihood contributions for assigned class
#                 - cell_num: Analyzed cell number
#                 - gene_names: List of gene names
#                 - class_names: Available cell type classes
#                 - class_probs: Probability distribution over cell types
#                 - contr: Log-likelihood contributions for all classes
#
#         Raises:
#             ValueError: If cell_num invalid or user_class not recognized
#         """
#
#         # If original labels have been renumbered find the label it's been mapped to.
#         if self.vb.config['label_map']:
#             cell_num = self.vb.config['label_map'][label]
#         else:
#             cell_num = label
#
#         if cell_num < 0 or cell_num >= self.vb.nC:
#             raise ValueError(f"Invalid cell number. Must be between 0 and {self.vb.nC - 1}")
#
#         assigned_class_idx = np.argmax(self.vb.cells.classProb[cell_num])
#         assigned_class = self.vb.cellTypes.names[assigned_class_idx]
#
#         if user_class is None:
#             user_class = assigned_class
#
#         try:
#             user_class_idx = np.where(self.vb.cellTypes.names == user_class)[0][0]
#         except IndexError:
#             raise ValueError(
#                 f"Invalid user class: {user_class}. Available classes are: {', '.join(self.vb.cellTypes.names)}")
#
#         # Calculate log-likelihood contributions
#         ScaledExp = np.einsum('cgk,g->cgk', self.vb.scaled_exp.compute(), self.vb.genes.eta_bar) + self.vb.config['SpotReg']
#         pNegBin = ScaledExp / (self.vb.config['rSpot'] + ScaledExp)
#         cgc = self.vb.cells.geneCount
#         contr = negative_binomial_loglikelihood(cgc, self.vb.config['rSpot'], pNegBin)
#         # contr = np.exp(contr)  # loglikelihood -> likelihood
#
#         # Calculate contributions for all classes
#         all_class_contrs = contr[cell_num, :, :]
#
#         # Prepare results dictionary
#         user_data = {
#             class_name: all_class_contrs[:, class_idx].tolist()
#             for class_idx, class_name in enumerate(self.vb.cellTypes.names)
#         }
#
#         class_probs = dict(zip(
#             self.vb.cellTypes.names.tolist(),
#             self.vb.cells.classProb[cell_num].tolist()
#         ))
#
#         return {
#             'assigned_class': assigned_class,
#             'user_class': user_class,
#             'assigned_contr': all_class_contrs[:, assigned_class_idx].tolist(),
#             'cell_num': label,
#             'gene_names': self.vb.genes.gene_panel.tolist(),
#             'class_names': self.vb.cellTypes.names.tolist(),
#             'class_probs': class_probs,
#             'contr': user_data,
#             'gene_counts': self.vb.cells.geneCount[cell_num,:].tolist(),
#             'gene_efficiency': self.vb.genes.inefficiency.tolist(),
#             'scRNAseq_gene_counts': self.vb.single_cell.mean_expression.values.tolist(),
#             'mean_gene_reads_per_class': self.vb.cells.mean_gene_reads_per_class().tolist()
#         }
#
#     def spot_dist_and_prob(self, label) -> Dict:
#         """
#         Calculate the relationship between spot-to-cell distances and their assignment
#         probabilities for a given cell.
#
#         This function analyzes spatial relationships between spots and a target cell by:
#         1. Finding spots near the target cell
#         2. Calculating distances from spots to cell centroid
#         3. Computing assignment probabilities
#
#         Args:
#             label (int): The cell number to analyze
#
#         Returns:
#             dict: A dictionary containing plot data:
#                 - x (list): Distances from spots to cell centroid
#                 - y (list): Assignment probabilities
#                 - labels (list): Gene names for each spot
#                 - cell_num (int): The analyzed cell number
#                 - title (str): Plot title
#                 - xlabel (str): X-axis label
#                 - ylabel (str): Y-axis label
#
#         Note:
#             The returned data is structured for visualization in the cell analysis dashboard.
#         """
#
#         # If original labels have been renumbered find the label it's been mapped to.
#         if self.vb.config['label_map']:
#             cell_num = self.vb.config['label_map'][label]
#         else:
#             cell_num = label
#
#         # Get cell centroid
#         centroid_zyx = self.vb.cells.zyx_coords[cell_num]
#
#         # Find spots near target cell
#         is_spot_near_target_cell = self.vb.spots.parent_cell_id == cell_num
#         mask = np.any(is_spot_near_target_cell, axis=1)
#
#         # Get probabilities for spots assigned to this cell
#         prob = self.vb.spots.parent_cell_prob[is_spot_near_target_cell]
#
#         # Select relevant spots
#         spots = self.vb.spots.data[mask]
#
#         # Find most likely parent cell for each spot
#         max_idx = np.argmax(self.vb.spots.parent_cell_prob[mask], axis=1)
#         assigned_cell = np.choose(max_idx, self.vb.spots.parent_cell_id[mask].T)
#
#         # Calculate distances and create DataFrame
#         spots = spots.assign(
#             cell_x=centroid_zyx[2],
#             cell_y=centroid_zyx[1],
#             cell_z=centroid_zyx[0],
#             prob=prob,
#             assigned_cell=assigned_cell,
#             dist=np.sqrt(
#                 (spots.x - centroid_zyx[2]) ** 2 +
#                 (spots.y - centroid_zyx[1]) ** 2 +
#                 (spots.z - centroid_zyx[0]) ** 2
#             )
#         )
#
#         # Select and order columns
#         spots = spots[
#             ['x', 'y', 'z', 'gene_name', 'cell_x', 'cell_y', 'cell_z',
#              'prob', 'assigned_cell', 'dist']
#         ].reset_index()
#
#         # Prepare plot data
#         # Prepare plot data with cell-specific axis labels
#         gene_counts = dict(zip(self.vb.genes.gene_panel, self.vb.cells.geneCount[cell_num,:].tolist()))
#         data = {
#             'x': spots.dist.tolist(),
#             'y': spots.prob.tolist(),
#             'labels': spots.gene_name.tolist(),
#             'cell_num': label,
#             'gene_counts': gene_counts,
#             'title': f'Cell {label} - Distance vs Assignment Probability',
#             'xlabel': f'Distance (px) from cell {label} centroid',
#             'ylabel': f'Assignment probability to cell {label}'
#         }
#         return data
#
#     def view_cell(self, cell_num, output_dir=None) -> None:
#         """
#         Generates data and launches the cell analysis dashboard for a specific cell.
#
#         This function:
#         1. Generates analysis data for the specified cell
#         2. Saves data to JSON files
#         3. Sets up and launches a local server
#         4. Opens the analysis dashboard in a web browser
#
#         Args:
#             cell_num: The cell number to analyze
#             output_dir: Optional directory to save JSON files. Defaults to cell_analysis/
#
#         Note:
#             Creates a local server to serve the dashboard. Close terminal to stop server.
#         """
#         # Get default output directory if none specified
#         if output_dir is None:
#             output_dir = get_out_dir(self.vb.config['output_path'])
#             output_dir = os.path.join(output_dir, 'data', 'debug', 'cell_analysis')
#
#         # Ensure output directory exists
#         os.makedirs(output_dir, exist_ok=True)
#
#         # Get the data2
#         assigned_class_idx = np.argmax(self.vb.cells.classProb[cell_num])
#         user_class = self.vb.cellTypes.names[assigned_class_idx]
#
#         # Generate gene contribution data
#         spot_dist = self.spot_dist_and_prob(cell_num)
#
#         loglik_data = self.gene_loglik_contributions(cell_num, user_class)
#         with open(os.path.join(output_dir, 'gene_loglik_contr.json'), 'w') as fp:
#             json.dump(loglik_data, fp)
#             analysis_logger.info(f'saved at {os.path.join(output_dir, "gene_loglik_contr.json")}')
#
#         # Save the data files
#         with open(os.path.join(output_dir, "spot_dist.json"), "w") as f:
#             json.dump(spot_dist, f)
#             analysis_logger.info(f'saved at {os.path.join(output_dir, "spot_dist.json")}')
#
#         pciSeq_dir = get_pciSeq_install_dir()
#         src = os.path.join(pciSeq_dir, 'static', 'cell_analysis')
#
#         shutil.copytree(src, output_dir, dirs_exist_ok=True)
#         # viewer_utils_logger.info('viewer code (%s) copied from %s to %s' % (dim, src, dst))
#
#         # Launch the dashboard
#         # Start HTTP server
#         os.chdir(output_dir)
#         PORT = 8000 + random.randint(0, 999)
#
#         Handler = http.server.SimpleHTTPRequestHandler
#
#         def start_server():
#             with socketserver.TCPServer(("", PORT), Handler) as httpd:
#                 print(f"Serving cell analysis dashboard at http://localhost:{PORT}")
#                 httpd.serve_forever()
#
#         # Start server in a separate thread
#         server_thread = threading.Thread(target=start_server, daemon=True)
#         server_thread.start()
#
#         # Open the dashboard in the default browser.
#         # Add the timestamp as version number to prevent loading from the cache
#         webbrowser.open(f'http://localhost:{PORT}/dashboard/cell_index.html?v={time.time()}')
#
#         # try:
#         #     input("Press Enter to stop the server and close the dashboard...")
#         # except KeyboardInterrupt:
#         #     print("\nShutting down the server...")
