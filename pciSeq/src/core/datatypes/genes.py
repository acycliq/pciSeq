# Standard library imports
import logging
from typing import Tuple, Dict, Any

# Third party imports
import numpy as np
import pandas as pd
import scipy
from sklearn.preprocessing import MinMaxScaler
from shapely.geometry import MultiPoint, Polygon, mapping
import alphashape

genes_logger = logging.getLogger(__name__)


class Genes(object):
    """
    Manages gene-specific data and calculations, including initialization and
    computation of gene expression parameters.

    Attributes:
        gene_panel (np.array): Array of unique gene names.
        _eta_bar (np.array): Eta bar values: This is basically the expected Gene inefficiency.
        _logeta_bar (np.array): Log eta bar values for genes.
        nG (int): Number of genes.
    """

    def __init__(self, spots, config: Dict):
        """
        Initializes the Genes object with spot data.

        Parameters:
            spots (Spots): Spots object containing spot data.
        """
        self.gene_panel = np.unique(spots.data.gene_name.values)
        self._eta_bar = None
        self._logeta_bar = None
        self.nG = len(self.gene_panel)
        self._misread_density = None
        self.config = config

    @property
    def eta_bar(self):
        """Returns the eta bar values for genes."""
        return self._eta_bar

    @property
    def logeta_bar(self):
        """Returns the log eta bar for genes (estimated mean of the posterior)."""
        return self._logeta_bar

    @property
    def inefficiency(self):
        """
        Returns the gene inefficiency
        The actual gene inefficiency is the estimated mean of the posterior (eta_bar)
        multiplied by the inefficiency (user-defined) value that was passed in the algo
        via the configuration file
        """
        return self.eta_bar * self.config['Inefficiency']

    @property
    def misread_density(self):
        """
        Misread density expresses the noise of the signal. It is estimated
        using the number of points that are too far from the closest cell
        and are also on the background
        """
        return self._misread_density

    def init_eta(self, a, b):
        """
        Initializes eta values for genes.

        Parameters:
            a (float): Parameter a for eta calculation.
            b (float): Parameter b for eta calculation.
        """
        self._eta_bar = np.ones(self.nG, dtype=np.float32) * (a / b)
        self._logeta_bar = np.ones(self.nG, dtype=np.float32) * self._digamma(a, b)

    def calc_eta(self, a, b):
        """
        Calculates eta values for genes.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.
        """
        a = a.astype(np.float32)
        b = b.astype(np.float32)
        self._eta_bar = a / b
        self._logeta_bar = self._digamma(a, b)

    def _digamma(self, a, b):
        """
        Calculates the digamma function for eta calculation.

        Parameters:
            a (np.array): Array of parameter a values.
            b (np.array): Array of parameter b values.

        Returns:
            np.array: Digamma values.
        """
        return scipy.special.psi(a) - np.log(b)

    def get_inefficiency(self, gene=None):
        """
        Retrieve the inefficiency values for one or more genes.

        This is a convenience method that returns inefficiency values from the gene panel.
        it returns a  DataFrame with the genes and the corresponding inefficiencies.

        Parameters:
        ----------
        gene : str, list of str, or None (default: None)
            - If None, returns inefficiency values for all genes.
            - If a string, returns inefficiency for the specified gene as a DataFrame.
            - If a list of strings, returns inefficiency values for the specified genes.

        Returns:
        -------
        pandas.DataFrame
            A DataFrame with genes as the index and inefficiency values as the column.
            Missing genes will be included with NaN values.

        Raises:
        ------
        TypeError
            If `gene` is not a string, list of strings, or None.
        """
        df = pd.DataFrame(
            {'inefficiency': self.inefficiency},
            index=self.gene_panel
        )

        if gene is None:
            return df  # Return full DataFrame

        if isinstance(gene, str):
            return df.loc[[gene]] if gene in df.index else pd.DataFrame(columns=df.columns, index=[gene])

        if isinstance(gene, list):
            return df.reindex(gene)  # Handles missing genes gracefully (NaN for missing ones)

        raise TypeError("Expected gene to be a string, list, or None.")

    def calc_misread_density(self, spots, mcr):
        # find the mid plane
        mid_plane = int(spots.data.plane_id.mean())

        # mask = spots.data.plane_id == mid_plane
        # points_df = spots.data.loc[mask, ['x', 'y']]
        # points_df = spots.data[['x', 'y']].copy()
        area = self.pointcloud_shape(spots, mid_plane)
        misreads_per_gene = self.remote_spots(spots, mid_plane, mcr)

        return misreads_per_gene/area

    def pointcloud_shape(self, spots, mid_plane, alpha=7):

        # get the data around the mid_plane
        plane_mask = spots.data.plane_id == mid_plane
        points_df = spots.data.loc[plane_mask, ['x', 'y']]

        # ------------------------------
        # Step 2: Scale the Data Using MinMaxScaler
        # ------------------------------
        scaler = MinMaxScaler()
        points_scaled = scaler.fit_transform(points_df)

        # ------------------------------
        # Step 3: Compute the Surrounding Polygon in Scaled Space
        # ------------------------------
        # Create a MultiPoint geometry from the scaled data
        # multi_pt = MultiPoint(points_scaled)

        # Compute the convex hull (you may replace this with a concave hull method if needed)
        alpha_shape = alphashape.alphashape(points_scaled, alpha)
        # hull = multi_pt.convex_hull

        # Extract the hull coordinates using shapely.mapping
        mapped_hull = mapping(alpha_shape)
        # For a Polygon, the outer boundary is in the first element of the 'coordinates'
        hull_coords_scaled = np.array(mapped_hull['coordinates'][0])

        # #make sure it is closed:
        if ~np.all(hull_coords_scaled[0] == hull_coords_scaled[-1]):
            np.append(hull_coords_scaled, hull_coords_scaled[0])

        # ------------------------------
        # Step 4: Recover the Original Scale of the Hull Coordinates
        # ------------------------------
        hull_coords_original = scaler.inverse_transform(hull_coords_scaled)

        # ------------------------------
        # Step 5: Plot the Data and the Surrounding Polygon
        # ------------------------------
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 6))
        plt.scatter(points_df['x'], points_df['y'], color='blue', label="Data Points", s=2)

        # Ensure the polygon is closed by appending the first coordinate at the end
        # x_poly = np.append(hull_coords_original[:, 0], hull_coords_original[0, 0])
        # y_poly = np.append(hull_coords_original[:, 1], hull_coords_original[0, 1])
        plt.plot(hull_coords_original[:, 0], hull_coords_original[:, 1], 'r-', linewidth=2, label="Surrounding Polygon")

        # ------------------------------
        # Step 6: Calculate and Print the Area of the Polygon
        # ------------------------------
        polygon = Polygon(hull_coords_original)
        area = polygon.area
        print("Area of the shape:", area)

        return area

    def remote_spots(self, spots, mid_plane, mcr):
        mid_plane_mask = spots.data.plane_id == mid_plane
        mid_spots = spots.data[mid_plane_mask]
        dist_mask = spots.Dist[mid_plane_mask, 0] > 3 * mcr
        isolated_spots = mid_spots[dist_mask]

        # select those on the background
        # isolated_spots = isolated_spots[isolated_spots.label == 0]
        misreads_per_gene = isolated_spots[['gene_name', 'label']].groupby('gene_name').count()

        # a = spots.data.assign(z_stack=spots.data.z * 0.28 / 0.9)
        # b = np.floor(a.z_stack) == 30
        # spots_30 = a[b]
        # my_mask = spots.Dist[b, 0] > 3 * mcr
        # isolated_spots = spots_30[my_mask]
        # misreads_per_gene = isolated_spots[['gene_name', 'label']].groupby('gene_name').count()
        # misreads_per_gene.loc['Plp1'] # Should return 69

        return misreads_per_gene




