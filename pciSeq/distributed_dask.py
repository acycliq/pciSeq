import pandas as pd
import numpy as np
import dask
import os
import dask.array as da
import dask.dataframe as dd
from dask.array import overlap
from dask.diagnostics import ProgressBar, ResourceProfiler
from scipy.sparse import load_npz, coo_matrix
from pciSeq import fit as fit_chunk
from typing import Tuple, List, Iterator, Any
from scipy.sparse import coo_matrix
import logging
from concurrent.futures import ThreadPoolExecutor

# Set up logging
dd_logger = logging.getLogger(__name__)
dd_logger.setLevel(logging.INFO)

DEFAULT_OPTS = {
    'chunk_size': 2000,
    'batch_size': 4,
    'memory_limit': '4GB',
    'overlap_factor': 1.2,
    'min_chunk_size': 1000,
    'max_chunk_size': 4000,
    'scheduler': 'threads',
    'cell_radius': 8.5
}


def filter_spots(spots: pd.DataFrame, bbox: List[List[int]]) -> Tuple[pd.DataFrame, Tuple[int, int, int]]:
    """
    Filters spots based on spatial boundaries and adjusts coordinates.

    Args:
        spots: DataFrame with x, y, z coordinates
        bbox: Bounding box [[zmin, zmax], [ymin, ymax], [xmin, xmax]]

    Returns:
        Filtered DataFrame and the offset for each dimension
    """
    zmin, zmax = bbox[0]
    ymin, ymax = bbox[1]
    xmin, xmax = bbox[2]

    # Create a mask to filter spots within the bounding box
    mask = (
            (spots['x'] >= xmin) &
            (spots['x'] < xmax) &
            (spots['y'] >= ymin) &
            (spots['y'] < ymax) &
            (spots['z'] >= zmin) &
            (spots['z'] < zmax)
    )

    # Filter and adjust coordinates
    df = spots[mask].copy()
    df['x'] -= xmin
    df['y'] -= ymin
    df['z'] -= zmin

    return df, (xmin, ymin, zmin)


# Define the function to slice data with overlap
def slice_data(image: da.Array, spots_future, overlap_depth: int):
    """
    Slices the image data with overlap and filters spots for each chunk.

    Args:
        image: 3D image data, expected to be a Dask array.
        spots_future: Dask future for the DataFrame with x, y, z coordinates.
        overlap_depth: Depth of overlap for each chunk.

    Returns:
        Tuple of Dask arrays: (output_image, output_spots)
    """
    # Ensure spots is computed from the future
    spots = spots_future.result()

    # Add overlap to image array chunks
    image_with_overlap = overlap.overlap(image, depth=overlap_depth, boundary='reflect')

    def process_chunk(image_chunk, spots, block_info=None):
        # Extract chunk indices from block_info
        if block_info is not None:
            chunk_indices = block_info[0]['chunk-location']
        else:
            raise ValueError("block_info is None, cannot determine chunk indices.")

        # Ensure chunk_indices has three elements
        if len(chunk_indices) != 3:
            raise ValueError(f"Expected 3 indices, got {len(chunk_indices)}: {chunk_indices}")

        # Calculate chunk boundaries with overlap
        z_start, y_start, x_start = [idx * (size - 2 * overlap_depth) for idx, size in zip(chunk_indices, image_chunk.shape)]
        z_end, y_end, x_end = z_start + image_chunk.shape[0], y_start + image_chunk.shape[1], x_start + image_chunk.shape[2]

        # Filter spots within these chunk boundaries
        filtered_spots = spots[
            (spots['x'] >= x_start) & (spots['x'] < x_end) &
            (spots['y'] >= y_start) & (spots['y'] < y_end) &
            (spots['z'] >= z_start) & (spots['z'] < z_end)
        ].copy()  # Use copy to avoid SettingWithCopyWarning

        # Optional: Perform processing on image_chunk and filtered_spots
        output_image = np.zeros_like(image_chunk)
        output_image.fill(len(filtered_spots))  # Example: fill chunk with the count of spots

        return output_image, filtered_spots  # Return both image chunk and spots

    # Apply process_chunk to each chunk of image_with_overlap
    output_dask = da.map_blocks(
        process_chunk, image_with_overlap,
        dtype=(image.dtype, object),  # Return tuple of (image_chunk, filtered_spots)
        spots=spots,
        meta=(np.empty((0, 0, 0), dtype=image.dtype), pd.DataFrame(columns=spots.columns))  # Metadata for output
    )

    # Trim overlap after processing to match original array shape
    # output_image = output_dask.map_blocks(lambda x: x[0], dtype=image.dtype)
    # output_spots = output_dask.map_blocks(lambda x: x[1], dtype=object)

    return output_dask[0], output_dask[1]



def _slice_data(image: da.Array, spots: pd.DataFrame) -> Tuple[
    List[Tuple[da.Array, Tuple[int, int, int]]], List[Tuple[pd.DataFrame, Tuple[int, int, int]]]]:
    """
    Slices both image and spot data into corresponding chunks for parallel processing.

    Args:
        image: Dask array containing the 3D image data
        spots: DataFrame containing spot coordinates and gene information

    Returns:
        Tuple of (image_chunks, spot_chunks) where each chunk corresponds to the same spatial region
    """
    # Get chunk slices and indices
    slices = da.core.slices_from_chunks(image.chunks)
    indices = list(np.ndindex(*image.numblocks))

    # Prepare lists to hold chunks
    image_chunks = []
    spot_chunks = []

    # Iterate over each chunk
    for idx, s in zip(indices, slices):
        # Extract the image chunk
        img_chunk = image[s]
        image_chunks.append((img_chunk, idx))

        # Calculate bounding box for the current chunk
        z_slice, y_slice, x_slice = s
        bbox = [[z_slice.start, z_slice.stop], [y_slice.start, y_slice.stop], [x_slice.start, x_slice.stop]]

        # Filter spots for the current chunk
        filtered_spots, offset = filter_spots(spots, bbox)
        spot_chunks.append((filtered_spots, offset))

    return image_chunks, spot_chunks


def process_chunks(image_chunks: List[Tuple[da.Array, tuple]],
                   spot_chunks: List[Tuple[pd.DataFrame, tuple]],
                   scRNAseq: pd.DataFrame,
                   opts: dict) -> List[Any]:
    """
    Process data chunks with memory management and error handling.
    """
    lazy_data = []

    # Process in batches
    for batch in create_batches(zip(image_chunks, spot_chunks),
                                batch_size=opts.get('batch_size', 4)):
        batch_results = []
        for img_chunk, spot_chunk in batch:
            try:
                # Create delayed computation
                result = dask.delayed(process_single_chunk, pure=True)(
                    img_chunk=img_chunk[0],
                    spot_chunk=spot_chunk[0],
                    scRNAseq=scRNAseq,
                    opts=opts
                )
                batch_results.append(result)

            except Exception as e:
                dd_logger.error(f"Chunk processing error: {str(e)}")
                continue

        lazy_data.extend(batch_results)

    return lazy_data


def process_single_chunk(img_chunk: da.Array,
                         spot_chunk: pd.DataFrame,
                         scRNAseq: pd.DataFrame,
                         opts: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process a single data chunk.
    """
    try:
        # Convert to COO matrix for efficient processing
        coo = coo_matrix(img_chunk)

        # Your existing fit_chunk function call here
        cell_data, gene_data = fit_chunk(
            spots=spot_chunk,
            coo=coo,
            scRNAseq=scRNAseq,
            opts=opts
        )

        return cell_data, gene_data

    except Exception as e:
        dd_logger.error(f"Single chunk processing error: {str(e)}")
        raise


def create_batches(iterable: Iterator, batch_size: int) -> Iterator[List]:
    """
    Creates batches from an iterable.
    """
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def fit(spots: Any,
        image: da.Array,  # Now explicitly expecting a Dask array
        scRNAseq: Any,
        opts: dict = None) -> List[Any]:
    """
    Main processing function that distributes computation across chunks.

    Args:
        spots: Future or DataFrame with spot coordinates
        image: Future or ndarray with image data
        image_shape: Tuple of image dimensions
        scRNAseq: Future or DataFrame with scRNA-seq data
        opts: Optional dictionary of parameters
    """
    # Merge options
    # opts = {**DEFAULT_OPTS, **user_opts} if user_opts else DEFAULT_OPTS

    try:
        # No need to convert image anymore, it's already a Dask array
        overlap = int(opts['cell_radius'] * opts['overlap_factor'])
        # daskimage = da.overlap.overlap(image, depth=overlap, boundary='reflect')

        # Slice data into chunks
        image_chunks, spot_chunks = slice_data(image, spots, overlap_depth=overlap)

        # Process chunks with progress monitoring
        with ProgressBar(), ResourceProfiler() as rprof:
            res = process_chunks(image_chunks, spot_chunks, scRNAseq, opts)
            return res

    except Exception as e:
        dd_logger.error(f"Processing failed: {str(e)}")
        raise


def calculate_optimal_chunk_size(shape: Tuple[int, int, int], memory_limit: str) -> Tuple[int, int, int]:
    """
    Calculates optimal chunk size based on image shape and memory constraints.

    Args:
        shape: Tuple of image dimensions (z, y, x).
        memory_limit: Memory limit as a string, e.g., '4GB' or '500MB'.

    Returns:
        Tuple of optimal chunk sizes for each dimension (z_chunk, y_chunk, x_chunk).
    """
    # Convert memory limit to bytes
    if isinstance(memory_limit, str):
        if memory_limit.endswith('GB'):
            memory_bytes = float(memory_limit[:-2]) * 1e9
        elif memory_limit.endswith('MB'):
            memory_bytes = float(memory_limit[:-2]) * 1e6
        else:
            raise ValueError("Memory limit must be specified in GB or MB")
    else:
        raise ValueError("Memory limit must be a string ending with 'GB' or 'MB'")

    # Calculate the number of bytes per element (assuming float64, which is 8 bytes)
    bytes_per_element = 8

    # Calculate the total number of elements that can fit in the memory limit
    elements_per_chunk = memory_bytes / bytes_per_element

    # Calculate the chunk size for each dimension
    # Here, we assume a cubic chunk for simplicity, but you can adjust this based on your needs
    z_chunk = min(shape[0], int(elements_per_chunk ** (1/3)))
    y_chunk = min(shape[1], int(elements_per_chunk ** (1/3)))
    x_chunk = min(shape[2], int(elements_per_chunk ** (1/3)))

    return z_chunk, y_chunk, x_chunk


if __name__ == "__main__":
    from dask.distributed import Client, LocalCluster

    # Set up cluster
    cluster = LocalCluster(
        n_workers=4,
        threads_per_worker=2,
        memory_limit='4GB'
    )

    ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

    # Configure options
    DEFAULT_OPTS = {
        'memory_limit': '8GB',  # Adjust based on your machine
        'batch_size': 8,  # Number of chunks to process simultaneously
        'cell_radius': 8.5,  # Typical cell radius in pixels
        'overlap_factor': 1.2,  # Overlap between chunks
        'scheduler': 'threads'  # or 'processes' for CPU-bound tasks
    }

    user_opts = {'save_data': True,
                 'launch_viewer': False,
                 'launch_diagnostics': False,
                 # 'InsideCellBonus': 0,
                 # 'cell_radius': 8.5
                 }

    my_opts = {**DEFAULT_OPTS, **user_opts}


    scRNAseq = pd.read_csv('/media/dimitris/New Volume/data/Mathieu/WT94_alpha072/yao_73g_cx.csv')
    scRNAseq = scRNAseq.set_index('cluster').T

    # # read 3D some demo data
    ca1_spots = pd.read_csv(r"/media/dimitris/New Volume/data/Mathieu/WT94_alpha072/pciseq_anchor.csv")
    ca1_spots = ca1_spots[['y', 'x', 'z_stack', 'Gene']]
    ca1_spots = ca1_spots.rename(columns={'z_stack': 'z_plane',
                                                  'Gene': 'gene_name'})

    label_image = np.load(r"/media/dimitris/New Volume/data/Mathieu/WT94_alpha072/dapi_segmented.npy")
    # coo = [coo_matrix(d) for d in coo]
    # label_image = coo.toarray()


    with Client(cluster) as client:
        try:
            # Convert to dask array FIRST with optimal chunking
            chunk_size = calculate_optimal_chunk_size(label_image.shape, my_opts['memory_limit'])
            dask_image = da.from_array(label_image, chunks=chunk_size)

            # Now scatter the data
            scattered_spots = client.scatter(ca1_spots)
            scattered_scRNAseq = client.scatter(scRNAseq)

            # Get lazy results using dask array and scattered data
            lazy_results = fit(
                spots=scattered_spots,
                image=dask_image,  # Pass the dask array directly
                scRNAseq=scattered_scRNAseq,
                opts=my_opts
            )

            # Compute results
            results = []
            total = len(lazy_results)
            for i, lr in enumerate(lazy_results):
                try:
                    result = client.compute(lr).result()
                    results.append(result)
                    dd_logger.info(f"Processed chunk {i + 1}/{total}")
                except Exception as e:
                    dd_logger.error(f"Failed to compute chunk {i}: {str(e)}")
                    continue

            # Save results
            for i, (cell_data, gene_data) in enumerate(results):
                try:
                    cell_data.to_csv(f'cell_data_{i}.csv')
                    gene_data.to_csv(f'gene_data_{i}.csv')
                except Exception as e:
                    dd_logger.error(f"Failed to save result {i}: {str(e)}")

        except Exception as e:
            dd_logger.error(f"Pipeline failed: {str(e)}")
            raise
