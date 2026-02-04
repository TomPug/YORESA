import os
import sys
import numpy as np
import rasterio
from scipy.signal import periodogram
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

def compute_row_psd(args):
    """
    Compute the PSD for a single row of pixels.
    Args:
        args: tuple (row_index, row_data, nf)
    Returns:
        tuple: (row_index, psd_row) where psd_row has shape (width, nf)
    """
    idx, row_data, nf = args
    w = row_data.shape[0]
    row_psd = np.zeros((w, nf), dtype=np.float32)
    for j in range(w):
        _, psd = periodogram(row_data[j, :])
        row_psd[j, :] = psd
    return idx, row_psd


def periodogram_parallel(input_path, output_path, periods_path, num_workers=None):
    """
    Calculate the periodogram of a multiband raster and save the PSD and periods.

    Args:
        input_path (str): Path to input raster (bands=time axis).
        output_path (str): Path to output raster for PSD results.
        periods_path (str): Path to binary file for periods (1/frequency).
        num_workers (int, optional): Number of parallel workers. Defaults to CPU count.
    """
    # Read input raster
    with rasterio.open(input_path) as src:
        img = src.read()[:-2]  # shape: (bands, height, width)
        profile = src.profile

    bands, height, width = img.shape
    # Reorder to (height, width, bands)
    img_t = np.transpose(img, (1, 2, 0))

    # Compute frequency bins
    freq, psd_example = periodogram(img_t[0, 0, :])
    nf = len(freq)

    # Compute periods, avoiding divide-by-zero at freq[0]
    (1 / freq).tofile(periods_path, sep='\n')

    # Prepare result container
    psd_full = np.zeros((height, width, nf), dtype=np.float32)

    # Prepare row tasks including nf for pickling
    rows = [(i, img_t[i, :, :], nf) for i in range(height)]

    # Parallel processing
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for idx, row_psd in tqdm(executor.map(compute_row_psd, rows), total=height, desc="Rows processed"):
            psd_full[idx, :, :] = row_psd

    # Reorder back to (bands, height, width)
    psd_full_t = np.transpose(psd_full, (2, 0, 1))

    # Update profile for output raster
    profile.update(
        dtype=rasterio.float32,
        count=psd_full_t.shape[0]
    )

    # Write PSD raster
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(psd_full_t)


def main():
    if len(sys.argv) != 2:
        print("Usage: python Periodogram.py <path_to_raster_image>")
        sys.exit(1)

    input_path = sys.argv[1]
    if not os.path.exists(input_path):
        print(f"Error: The file {input_path} does not exist.")
        sys.exit(1)

    base, _ = os.path.splitext(input_path)
    output_psd = f"{base}_periodogram.dat"
    output_periods = f"{base}_periods.txt"
    num_workers = os.cpu_count()//2 # Use all available CPU cores

    periodogram_parallel(input_path, output_psd, output_periods,num_workers)
    print(f"Saved PSD raster to {output_psd}")
    print(f"Saved periods to {output_periods}")

if __name__ == '__main__':
    main()