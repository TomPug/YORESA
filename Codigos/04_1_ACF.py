#### Script to calculate ACF (Autocorrelation Function) of time series extracted from raster data at given point locations
import numpy as np
from statsmodels.tsa.stattools import acf
import xarray as xr


def read_raster(input_path):
    """Read multiband raster data with temporal coordinates.

    Args:
        input_path (str): Path to input raster file.

    Returns:
        xr.DataArray: DataArray with dimensions (time, y, x).
    """
    # Leer con xarray usando chunks para no cargar todo en memoria
    ds = xr.open_dataset(input_path, chunks={'time': -1, 'y': 2000, 'x': 2000})

    print(f"Variables disponibles: {list(ds.data_vars)}")
    print(f"Dimensiones del dataset: {dict(ds.sizes)}")
    print(f"Coordenadas: {list(ds.coords)}")
    
    # Buscar la variable de datos principal (excluir variables de referencia espacial)
    skip_vars = {'spatial_ref', 'crs', 'transverse_mercator', 'grid_mapping'}
    data_vars = [v for v in ds.data_vars if v not in skip_vars]
    
    if not data_vars:
        raise ValueError(f"No se encontró variable de datos. Variables disponibles: {list(ds.data_vars)}")
    
    var_name = data_vars[0]
    raster = ds[var_name].astype(np.float32)  # Reducir memoria usando float32
    
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


def calculate_acf_xarray(da, nlags=None):
    """Calculate ACF (Autocorrelation Function) for each pixel's time series using xarray.

    Args:
        da (xr.DataArray): DataArray with dimensions (time, y, x).
        nlags (int, optional): Number of lags to compute. If None, uses n_timesteps - 1.
    
    Returns:
        xr.DataArray: DataArray with dimensions (y, x, lag) containing ACF values.
        np.ndarray: 1D array of lag indices.
    """
    spatial_dims = ("y", "x")
    time_dim = "time"
    
    # Determinar número de lags
    n_timesteps = da.sizes[time_dim] if time_dim in da.dims else da.shape[0]
    if nlags is None:
        nlags = min(n_timesteps - 1, 40)  # Por defecto, máximo 40 lags o n-1
    
    lags = np.arange(nlags + 1)  # Incluye lag 0
    n_lags = len(lags)
    
    print(f"Número de timesteps: {n_timesteps}")
    print(f"Número de lags a calcular: {nlags}")
    
    # Stack spatial dimensions for vectorized processing
    stacked = da.stack(pixel=spatial_dims)
    
    # Prepare output array
    acf_values = np.zeros((stacked.sizes['pixel'], n_lags), dtype=np.float32)
    
    # Calculate ACF for each pixel
    print(f"Calculando ACF para {stacked.sizes['pixel']} píxeles...")
    for i in range(stacked.sizes['pixel']):
        if i % 10000 == 0:
            print(f"  Progreso: {i}/{stacked.sizes['pixel']}")
        
        ts = stacked.isel(pixel=i).values
        
        # Manejar NaN values
        if np.isnan(ts).all():
            acf_values[i, :] = np.nan
        else:
            # Calcular ACF usando statsmodels
            try:
                acf_result = acf(ts, nlags=nlags, fft=True, missing='conservative')
                acf_values[i, :] = acf_result
            except Exception:
                acf_values[i, :] = np.nan
    
    # Create DataArray with ACF results
    acf_da = xr.DataArray(
        acf_values,
        dims=['pixel', 'lag'],
        coords={
            'pixel': stacked.coords['pixel'],
            'lag': lags
        }
    )
    
    # Unstack back to spatial dimensions
    acf_unstacked = acf_da.unstack("pixel")
    
    # Transpose to (y, x, lag) for consistency
    acf_unstacked = acf_unstacked.transpose(spatial_dims[0], spatial_dims[1], "lag")
    
    # Copy spatial reference information from original if available
    if hasattr(da, 'rio') and da.rio.crs is not None:
        acf_unstacked = acf_unstacked.rio.write_crs(da.rio.crs)
        if da.rio.transform() is not None:
            acf_unstacked = acf_unstacked.rio.write_transform(da.rio.transform())
    
    return acf_unstacked, lags


def save_raster_xarray(output_path, acf_da):
    """Save ACF data to raster file.

    Args:
        output_path (str): Path to output raster file.
        acf_da (xr.DataArray): DataArray with dimensions (y, x, lag) containing ACF values.
    """
    import rioxarray as rxr
    
    # Transpose to (lag, y, x) for raster band structure
    acf_transposed = acf_da.transpose('lag', 'y', 'x')
    
    # Rename dimension 'lag' to 'band' for compatibility
    acf_transposed = acf_transposed.rename({'lag': 'band'})
    
    # Ensure rio accessor is available
    if not hasattr(acf_transposed, 'rio'):
        acf_transposed = acf_transposed.rio.write_crs("EPSG:32630")  # Tu CRS
    
    # Save to GeoTiff
    acf_transposed.rio.to_raster(
        output_path,
        driver='GTiff',
        dtype='float32'
    )


def save_acf_stats(stats_path, acf_da, lags):
    """Save ACF statistics summary.

    Args:
        stats_path (str): Path to output statistics file.
        acf_da (xr.DataArray): DataArray with ACF values.
        lags (np.ndarray): Array of lag indices.
    """
    # Calcular estadísticas por lag
    mean_acf = acf_da.mean(dim=['y', 'x']).values
    std_acf = acf_da.std(dim=['y', 'x']).values
    
    with open(stats_path, 'w') as f:
        f.write("# Lag\tMean_ACF\tStd_ACF\n")
        for lag, mean, std in zip(lags, mean_acf, std_acf):
            f.write(f"{lag}\t{mean:.6f}\t{std:.6f}\n")


def main():
    input_path = r"E:\SAR_UVa\coherence_stack.nc"
    output_path = r"E:\SAR_UVa\output_acf_raster.tif"
    stats_path = r"E:\SAR_UVa\output_acf_stats.txt"
    nlags = 40  # Número de lags a calcular (ajustar según necesidad)

    # Read raster data
    print("Reading raster data...")
    data = read_raster(input_path)
    
    # Calculate ACF
    print("\nCalculating ACF...")
    acf_data, lags = calculate_acf_xarray(data, nlags=nlags)
    
    # Save ACF raster
    print("\nSaving ACF raster...")
    save_raster_xarray(output_path, acf_data)
    
    # Save ACF statistics
    print("Saving ACF statistics...")
    save_acf_stats(stats_path, acf_data, lags)
    
    print(f"\nDone!")
    print(f"ACF raster saved to: {output_path}")
    print(f"ACF statistics saved to: {stats_path}")
    print(f"Each band in the raster corresponds to a lag value (0 to {nlags})")


if __name__ == "__main__":
    main()
