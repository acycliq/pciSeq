# Standard library imports
import logging
from typing import Tuple, Dict, Any
import os
import tempfile

# Third party imports
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
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
        mid_plane = int(spots.data.plane_id.max()/2)

        area = self.pointcloud_shape(spots, mid_plane)
        misreads_per_gene = self.remote_spots(spots, mid_plane, mcr)

        misread_df = misreads_per_gene/area

        # get the mean
        x = misread_df.mean()

        # if misread_df is empty, then the mean is NaN, in which case use the default value
        x = x if x == x else self.config['MisreadDensity']['default']

        # give missing genes a value equal to the mean (or the default value is mean is NaN)
        misread_df = misread_df.reindex(self.gene_panel, fill_value=x)

        # read the user-defined misread densities
        user_df = pd.Series(self.config['MisreadDensity'])

        # update with the user-defined values and return
        misread_df.update(user_df)  # in-place operation
        misread_df[:] = 1e-06
        return misread_df

    def pointcloud_shape(self, spots, mid_plane, alpha=7):
        """
        Compute the alpha shape (concave hull) of spot coordinates, plot the results, and save the plot.

        This method filters spots to include only those from the specified mid-plane,
        extracts their x and y coordinates, and scales them to a normalized range.
        It then computes the alpha shape of the scaled points using the provided alpha parameter.
        The hull coordinates are scaled back to the original coordinate system, the data points
        and the computed polygon are plotted, and the plot is saved as a PNG file. The output
        directory is taken from self.config['output_dir'] if available; otherwise, it defaults to
        the system's temporary directory under a folder named "pciSeq". Finally, the area of the
        polygon is computed, logged, and returned.

        Parameters:
            spots (object): An object with a DataFrame attribute 'data' containing spot information.
                            The DataFrame must include the columns 'plane_id', 'x', and 'y'.
            mid_plane (int or float): The identifier of the plane to filter the spots.
            alpha (float, optional): The alpha parameter controlling the concavity of the hull.
                                     Default is 7.

        Returns:
            float: The area of the polygon defined by the alpha shape.
        """
        # 1. Filter spots by the specified mid-plane and extract the x, y coordinates.
        plane_mask = spots.data.plane_id == mid_plane
        points_df = spots.data.loc[plane_mask, ['x', 'y']]

        # 2. Scale the coordinates to a normalized range [0, 1].
        scaler = MinMaxScaler()
        points_scaled = scaler.fit_transform(points_df)

        # 3. Compute the alpha shape (concave hull) of the scaled points.
        alpha_shape = alphashape.alphashape(points_scaled, alpha)

        # 4. Extract the hull coordinates using shapely.mapping.
        mapped_hull = mapping(alpha_shape)
        # For a Polygon, the exterior boundary is the first element of the 'coordinates' list.
        hull_coords_scaled = np.array(mapped_hull['coordinates'][0])

        # 5. Ensure the polygon is closed by appending the first coordinate at the end if necessary.
        if not np.allclose(hull_coords_scaled[0], hull_coords_scaled[-1]):
            hull_coords_scaled = np.vstack([hull_coords_scaled, hull_coords_scaled[0]])

        # 6. Convert the scaled hull coordinates back to the original coordinate system.
        hull_coords_original = scaler.inverse_transform(hull_coords_scaled)

        # 7. Plot the data points and the computed polygon.
        plt.figure(figsize=(8, 6))
        plt.scatter(points_df['x'], points_df['y'], color='blue', label="Data Points", s=2)
        plt.plot(hull_coords_original[:, 0], hull_coords_original[:, 1],
                 'r-', linewidth=2, label="Surrounding Polygon")
        plt.legend()

        # 8. Determine the output directory:
        #    Use self.config['output_dir'] if provided; otherwise, default to the system's tmp dir/pciSeq.
        output_dir = self.config.get('output_path') if hasattr(self, 'config') else None
        if output_dir == 'default':
            tmp_dir = tempfile.gettempdir()
            folder = os.path.join(tmp_dir, "pciSeq")
        else:
            folder = output_dir
        os.makedirs(folder, exist_ok=True)

        # Save the plot as a PNG file in the determined folder.
        file_path = os.path.join(folder, "pointcloud_shape.png")
        plt.savefig(file_path)
        plt.close()  # Close the figure to free up memory
        genes_logger.info(f"saved at {file_path}")

        # 9. Compute the area of the polygon using shapely.
        polygon = Polygon(hull_coords_original)
        area = polygon.area
        genes_logger.info(f"Area of the shape: {area}")

        return area

    def remote_spots(self, spots, mid_plane, mcr):
        """
        Identify remote spots from a specified plane and compute misread counts per gene.

        This function filters spots to include only those on the given mid-plane,
        then further isolates spots that are distant (based on the threshold defined
        as 3 * mcr) and are on the background (label == 0). It groups these isolated spots
        by gene name and counts the misreads per gene.

        Parameters:
            spots (object): An object with attributes:
                - data: A pandas DataFrame with columns 'plane_id', 'x', 'y', 'label', and 'gene_name'.
                - Dist: A NumPy array representing the distance of the spots from the cell centroid.
            mid_plane (int or float): Identifier for the plane to filter the spots.
            mcr (float): The mean cell radius. Used to define the distance threshold (threshold = 3 * mcr).

        Returns:
            pandas.Series: A series with gene names as the index and misread counts as the values.
        """
        # Filter spots based on the specified mid_plane.
        mid_plane_mask = spots.data.plane_id == mid_plane
        mid_spots = spots.data[mid_plane_mask]

        # Further filter spots based on the distance threshold.
        dist_mask = spots.Dist[mid_plane_mask, 0] > (3 * mcr)
        isolated_spots = mid_spots[dist_mask]

        # Select only those spots that are on the background (label == 0).
        isolated_spots = isolated_spots[isolated_spots.label == 0]

        # Group by gene name and count the number of misreads (using 'label' as a proxy for count).
        misreads_per_gene = isolated_spots[['gene_name', 'label']].groupby('gene_name').count()

        # Log the misread density for a specific gene, e.g. 'Plp1'.
        # (This could be parameterized if needed.)
        # genes_logger.info(f"Plp1 misread density: {misreads_per_gene.squeeze().get('Plp1', 'Not found')}")

        # Return the misread counts as a Series.
        return misreads_per_gene.squeeze()





