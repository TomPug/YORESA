#### Script to calculate periodogram of time series extracted from raster data at given point locations
import numpy as np
from scipy.signal import periodogram
import xarray as xr

def read_raster(input_path):
    """Read multiband raster data with temporal coordinates.

    Args:
        input_path (str): Path to input raster file.

    Returns:
        xr.DataArray: DataArray with dimensions (time, y, x).
    """
    # Leer con xarray usando chunks para no cargar todo en memoria
    # ds = xr.open_dataset(input_path, chunks={'time': -1, 'y': 500, 'x': 500})
    ds = xr.open_dataset(input_path)

    
    print(f"Variables disponibles: {list(ds.data_vars)}")
    print(f"Dimensiones del dataset: {dict(ds.sizes)}")
    print(f"Coordenadas: {list(ds.coords)}")
    
    # Buscar la variable de datos principal (excluir variables de referencia espacial)
    skip_vars = {'spatial_ref', 'crs', 'transverse_mercator', 'grid_mapping'}
    data_vars = [v for v in ds.data_vars if v not in skip_vars]
    
    if not data_vars:
        raise ValueError(f"No se encontró variable de datos. Variables disponibles: {list(ds.data_vars)}")
    
    var_name = data_vars[0]
    raster = (ds[var_name]*10000).astype(np.int16)  # Reducir memoria usando int16
    
    print(f"Variable leída: {var_name}")
    print(f"Dimensiones: {raster.dims}")
    print(f"Shape: {raster.shape}")
    
    # Identificar dimensión temporal
    time_dim = None
    for dim in raster.dims:
        if dim in ('time', 'date', 't', 'band'):
            time_dim = dim
            break
        # Verificar si tiene coordenadas datetime
        if dim in raster.coords:
            if np.issubdtype(raster.coords[dim].dtype, np.datetime64):
                time_dim = dim
                break
    
    if time_dim is None:
        print("No se encontró dimensión temporal explícita, usando primera dimensión")
        time_dim = raster.dims[0]
    
    print(f"Dimensión temporal identificada: {time_dim}")
    
    # Solo resamplear si tiene coordenadas datetime
    if time_dim in raster.coords and np.issubdtype(raster.coords[time_dim].dtype, np.datetime64):
        print(f"Rango temporal: {raster[time_dim].values[0]} a {raster[time_dim].values[-1]}")
        print(f"Número de timesteps: {len(raster[time_dim])}")
        
        # Resamplear a 6 días usando dask
        print("Resampling a 6 días (esto puede tardar)...")
        raster_resampled = raster.resample({time_dim: '6D'}).mean()
        
        # Computar el resample por chunks para evitar error de memoria
        raster_resampled = raster_resampled.compute()
        
        print(f"Timesteps después de resample: {len(raster_resampled[time_dim])}")
        return raster_resampled
    else:
        print(f"No hay coordenadas datetime, retornando sin resample")
        print(f"Número de capas: {raster.shape[0]}")
        return raster.compute()


def calculate_periodogram_xarray(da):
    """Calculate periodogram for each pixel's time series using xarray.

    Args:
        da (xr.DataArray): DataArray with dimensions (time, y, x).
    
    Returns:
        xr.DataArray: DataArray with dimensions (y, x, frequency) containing PSD values.
        np.ndarray: 1D array of frequency bins.
    """
    spatial_dims = ("y", "x")
    time_dim = "time"
    
    # Calculate frequency bins (same for all pixels)
    sample_ts = da.isel({spatial_dims[0]: 0, spatial_dims[1]: 0}).values
    freq, _ = periodogram(sample_ts)
    n_freq = len(freq)
    
    # Stack spatial dimensions for vectorized processing
    stacked = da.stack(pixel=spatial_dims)
    
    # Prepare output array
    psd_values = np.zeros((stacked.sizes['pixel'], n_freq), dtype=np.float32)
    
    # Calculate periodogram for each pixel
    print(f"Calculando periodograma para {stacked.sizes['pixel']} píxeles...")
    for i in range(stacked.sizes['pixel']):
        if i % 10000 == 0:
            print(f"  Progreso: {i}/{stacked.sizes['pixel']}")
        _, psd = periodogram(stacked.isel(pixel=i).values)
        psd_values[i, :] = psd
    
    # Create DataArray with PSD results
    psd_da = xr.DataArray(
        psd_values,
        dims=['pixel', 'frequency'],
        coords={
            'pixel': stacked.coords['pixel'],
            'frequency': freq
        }
    )
    
    # Unstack back to spatial dimensions
    psd_unstacked = psd_da.unstack("pixel")
    
    # Transpose to (y, x, frequency) for consistency
    psd_unstacked = psd_unstacked.transpose(spatial_dims[0], spatial_dims[1], "frequency")
    
    # Copy spatial reference information from original if available
    if hasattr(da, 'rio') and da.rio.crs is not None:
        psd_unstacked = psd_unstacked.rio.write_crs(da.rio.crs)
        if da.rio.transform() is not None:
            psd_unstacked = psd_unstacked.rio.write_transform(da.rio.transform())
    
    return psd_unstacked, freq


def save_raster_xarray(output_path, psd_da):
    """Save PSD data to raster file.

    Args:
        output_path (str): Path to output raster file.
        psd_da (xr.DataArray): DataArray with dimensions (y, x, frequency) containing PSD values.
    """
    import rioxarray as rxr
    
    # Transpose to (frequency, y, x) for raster band structure
    psd_transposed = psd_da.transpose('frequency', 'y', 'x')
    
    # Rename dimension 'frequency' to 'band' for compatibility
    psd_transposed = psd_transposed.rename({'frequency': 'band'})
    
    # Ensure rio accessor is available
    if not hasattr(psd_transposed, 'rio'):
        psd_transposed = psd_transposed.rio.write_crs("EPSG:32630")  # Tu CRS
    
    # Save to GeoTiff
    psd_transposed.rio.to_raster(
        output_path,
        driver='GTiff',
        dtype='float32'
    )


def main():
    input_path = r"E:\SAR_UVa\coherence_stack.nc"
    output_path = r"E:\SAR_UVa\output_psd_raster.tif"
    periods_path = r"E:\SAR_UVa\output_periods.txt"

    # Read raster data
    print("Reading raster data...")
    data = read_raster(input_path)
    
    # Calculate periodogram
    print("\nCalculating periodogram...")
    psd_data, freq = calculate_periodogram_xarray(data)
    
    # Save PSD raster
    print("\nSaving PSD raster...")
    save_raster_xarray(output_path, psd_data)
    
    # Save periods (1/frequency) and frequencies
    print("Saving periods and frequencies...")
    periods = 1 / freq[1:]  # Exclude first element (DC component, freq=0)
    
    # Save with headers for clarity
    with open(periods_path, 'w') as f:
        f.write("# Frequency\tPeriod (time steps)\n")
        for i, (fr, per) in enumerate(zip(freq[1:], periods)):
            f.write(f"{fr:.6f}\t{per:.6f}\n")
    
    print(f"\nDone!")
    print(f"PSD raster saved to: {output_path}")
    print(f"Periods saved to: {periods_path}")


if __name__ == "__main__":
    main()