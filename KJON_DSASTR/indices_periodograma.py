import rasterio as rio
import numpy as np
import pandas as pd
import os
import sys
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")

def read_raster(file_path):
    """
    Reads a raster file and returns its data as a numpy array.

    Parameters:
    file_path (str): The path to the raster file.

    Returns:
    numpy.ndarray: The data from the raster file.
    """
    with rio.open(file_path) as src:
        data = src.read()
        profile = src.profile
    return data,profile


def open_periodos(txt_path, freq_muestreo=73):
    """
    Opens a text file containing periods and returns them as a list of floats.

    Parameters:
    txt_path (str): The path to the text file containing periods.

    Returns:
    list: A list of periods as floats.
    """
    with open(txt_path, 'r') as archivo:
        content = archivo.read().strip()
        # Replace '/n' with '\n' to standardize line breaks, then split
        content = content.replace('/n', '\n')
        periodos = [float(val) for val in content.split('\n') if val.strip() != '']
        
    ciclos = [freq_muestreo/3, freq_muestreo/2, freq_muestreo/1]
    print("Ciclos:", ciclos)
    # pos_ciclos = [np.where(np.isclose(np.arrayperiodos, ciclo, atol=0.01))[0] for ciclo in ciclos]
    pos_ciclos = [np.argmin(np.abs(np.array(periodos) - ciclo)) for ciclo in ciclos]
    print("Posiciones de ciclos:", pos_ciclos)
    
    return pos_ciclos, ciclos

# Save the seasonality mode as a raster file
def save_raster(data, file_path,input_profile):
    """
    Saves a numpy array as a raster file.

    Parameters:
    data (numpy.ndarray): The data to save.
    file_path (str): The path where the raster file will be saved.
    profile (dict): The metadata profile for the raster file.
    """
    profile = input_profile.copy()  # Use the input profile to ensure correct metadata
    profile.update({
        'count': 5,
        'dtype':'float32'
    })
    with rio.open(file_path, 'w', **profile) as dst:
        dst.write(data)
        dst.descriptions = ['Seasonality Mode', 'Fisher Kappa', 'Seasonality Stability', 'Plurianual Cycles', 'Seasonality Amplitude']


def seasonality_Mode(ras, pos_ciclos):
    """
    Extracts the seasonality mode from the raster data based on the specified cycle positions.

    Parameters:
    ras (numpy.ndarray): The raster data.
    pos_ciclos (list): List of indices corresponding to the cycles.

    Returns:
    numpy.ndarray: The seasonality mode extracted from the raster data.
    """

    bands = ras[pos_ciclos, :, :]
    print("Shape of bands:", bands.shape)
    seasonalityMode = np.argmax(bands, axis=0)
    # set the value in ciclos acording to the position of the cycle
    seasonalityMode = np.where(seasonalityMode == 0, ciclos[0],
                            np.where(seasonalityMode == 1, ciclos[1],
                                    np.where(seasonalityMode == 2, ciclos[2], seasonalityMode)))
    seasonalityMode = seasonalityMode.astype(np.float32)
    
    return seasonalityMode


def fisher_kappa(ras, pos_ciclos):
    """
    Calculates the Fisher's kappa statistic for the raster data based on the specified cycle positions.

    Parameters:
    ras (numpy.ndarray): The raster data.
    pos_ciclos (list): List of indices corresponding to the cycles.

    Returns:
    numpy.ndarray: The Fisher's kappa statistic extracted from the raster data.
    """
    
    bands = ras[pos_ciclos, :, :]
    
    mean_period = np.mean(bands, axis=0)
    sa = np.max(bands, axis=0)
    fk = sa/ mean_period
    fisherKappa = np.where(np.isfinite(fk), fk, 0)  # Replace NaN with 0
    
    return fisherKappa


def seasonality_stability(ras, pos_ciclos):
    """
    Calculates the seasonality stability for the raster data based on the specified cycle positions.

    Parameters:
    ras (numpy.ndarray): The raster data.
    pos_ciclos (list): List of indices corresponding to the cycles.

    Returns:
    numpy.ndarray: The seasonality stability extracted from the raster data.
    """
    
    bands = ras[pos_ciclos, :, :]
    # bands_2 = ras[pos_ciclos[2][0]:, :, :]
    bands_2 = ras[pos_ciclos[2]:, :, :]
    sa = np.max(bands, axis=0)
    ta = np.sum(bands_2, axis=0)
    ss = sa / ta    
    return ss


def plurianual_cycles(ras, pos_ciclos):
    """
    Calculates the plurianual cycles for the raster data based on the specified cycle positions.

    Parameters:
    ras (numpy.ndarray): The raster data.
    pos_ciclos (list): List of indices corresponding to the cycles.

    Returns:
    numpy.ndarray: The plurianual cycles extracted from the raster data.
    """
    
    # high_cycle = np.sum(ras[:pos_ciclos[2][0]-1, :, :], axis=0)
    high_cycle = np.sum(ras[:pos_ciclos[2]-1, :, :], axis=0)
    all_ordenates = np.sum(ras, axis=0)
    pc = high_cycle / all_ordenates
    pc = pc.reshape((1, pc.shape[0], pc.shape[1]))
    return pc


def seasonality_amplitude(ras, pos_ciclos):
    """
    Calculates the seasonality amplitude for the raster data based on the specified cycle positions.

    Parameters:
    ras (numpy.ndarray): The raster data.
    pos_ciclos (list): List of indices corresponding to the cycles.

    Returns:
    numpy.ndarray: The seasonality amplitude extracted from the raster data.
    """
    
    bands = ras[pos_ciclos, :, :]
    sa = np.max(bands, axis=0)
    return sa

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python indices_periodograma.py <output_directory> <raster_file> <periods_file> <freq_muestreo>")
        print("Example: python indices_periodograma.py /path/to/output raster.dat periods.txt 73")
        sys.exit(1)

    output_dir = sys.argv[1]
    raster_file = sys.argv[2]
    periods_file = sys.argv[3]
    freq_muestreo = int(sys.argv[4])
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(raster_file):
        print(f"Error: The raster file {raster_file} does not exist.")
        sys.exit(1)
    if not os.path.exists(periods_file):
        print(f"Error: The periods file {periods_file} does not exist.")
        sys.exit(1)
        
    ras,profile = read_raster(raster_file)
    pos_ciclos, ciclos = open_periodos(periods_file, freq_muestreo=freq_muestreo)
    
    indices = np.zeros((5, ras.shape[1], ras.shape[2]), dtype=np.float32)
    
    indices[0, :, :] = seasonality_Mode(ras, pos_ciclos)
    indices[1, :, :] = fisher_kappa(ras, pos_ciclos)
    indices[2, :, :] = seasonality_stability(ras, pos_ciclos)
    indices[3, :, :] = plurianual_cycles(ras, pos_ciclos)
    indices[4, :, :] = seasonality_amplitude(ras, pos_ciclos)
    print("Indices calculated successfully.")
    fname = Path(raster_file).stem
    parts = fname.split('_')
    print(parts)
    tile = parts[1]
    index = parts[3]
    output_file = Path(output_dir) / f"indices_periodograma_{tile}_{index}.dat"
    save_raster(indices, output_file,profile)
    print(f"Indices saved to {output_file}")
    print("Done.")
    
    
    
    