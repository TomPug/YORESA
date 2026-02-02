import os
import rioxarray
import xarray as xr
import zipfile
import shutil
import pandas as pd
import time
from contextlib import contextmanager

def ensure_dir_exists(directory):
    """Asegura que el directorio existe, si no, lo crea."""
    if not os.path.exists(directory):
        os.makedirs(directory)

@contextmanager
def safe_open_file(path, mode='wb'):
    """Manejador de contexto para abrir archivos de manera segura."""
    attempts = 0
    max_attempts = 3
    while attempts < max_attempts:
        try:
            file = open(path, mode)
            yield file
            file.close()
            break
        except PermissionError:
            attempts += 1
            if attempts == max_attempts:
                raise
            time.sleep(1)  # Esperar 1 segundo antes de reintentar

def process_zip_files(input_dir, output_dir, temp_dir):
    """Procesa los archivos zip y crea el stack de datos."""
    # Asegurar que los directorios existen
    ensure_dir_exists(output_dir)
    ensure_dir_exists(temp_dir)
    
    # Listar y ordenar archivos
    lis_dir = os.listdir(input_dir)
    pos_dates = [(i, name.split('_')[1].split('T')[0]) for i, name in enumerate(lis_dir)]
    sorted_dates = sorted(pos_dates, key=lambda x: x[1])
    
    # Listas para almacenar resultados
    unw_arrays = []
    dates = []
    temp_files = []  # Lista para rastrear archivos temporales
    
    for pos, date in sorted_dates:
        zip_file = lis_dir[pos]
        zip_path = os.path.join(input_dir, zip_file)
        print(f"Procesando {zip_file}")
        
        try:
            # Crear un nombre único para el archivo temporal
            temp_tif = os.path.join(temp_dir, f"{date}_{int(time.time())}.tif")
            
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                for file in zip_ref.namelist():
                    if file.endswith('_corr.tif'):
                        # Extraer al directorio temporal
                        with zip_ref.open(file) as source:
                            with safe_open_file(temp_tif) as target:
                                shutil.copyfileobj(source, target)
                        
                        # Esperar un momento para asegurarse de que el archivo se ha escrito completamente
                        time.sleep(0.5)
                        
                        # Leer el archivo con rioxarray
                        try:
                            da = rioxarray.open_rasterio(temp_tif)
                            da_squeezed = da.squeeze()
                            unw_arrays.append(da_squeezed)
                            dates.append(date)
                            temp_files.append(temp_tif)  # Guardar la ruta del archivo temporal
                        except Exception as e:
                            print(f"Error al leer el archivo {temp_tif}: {e}")
                            # Si hay error, intentar eliminar el archivo temporal
                            if os.path.exists(temp_tif):
                                try:
                                    os.remove(temp_tif)
                                except Exception as cleanup_error:
                                    print(f"No se pudo eliminar el archivo temporal {temp_tif}: {cleanup_error}")
                        break
                        
        except Exception as e:
            print(f"Error al procesar {zip_file}: {e}")
            continue
    
    # Crear y guardar el stack
    if unw_arrays:
        try:
            stack = xr.concat(unw_arrays, dim=pd.DatetimeIndex(dates, name='time'))
            output_file = os.path.join(output_dir, 'coherence_stack.nc')
            
            # Intentar guardar el archivo con reintentos
            save_attempts = 0
            while save_attempts < 3:
                try:
                    stack.to_netcdf(output_file)
                    print(f"Stack guardado exitosamente en {output_file}")
                    break
                except Exception as e:
                    save_attempts += 1
                    if save_attempts == 3:
                        print(f"Error al guardar el stack después de 3 intentos: {e}")
                        break
                    time.sleep(2)
        finally:
            # Cerrar todos los DataArrays
            for da in unw_arrays:
                try:
                    da.close()
                except:
                    pass
            
            # Limpiar archivos temporales después de crear el stack
            for temp_file in temp_files:
                try:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                except Exception as e:
                    print(f"No se pudo eliminar el archivo temporal {temp_file}: {e}")
    else:
        print("No se encontraron archivos unw_phase.tif en los archivos zip.")

# Uso del código
if __name__ == "__main__":
    input_dir = r"\\tierra\L\SAR\INSAR_EXT"
    output_dir = r'E:/'
    temp_dir = r'E:/temp_processing'  # Directorio temporal para procesamiento
    
    try:
        process_zip_files(input_dir, output_dir, temp_dir)
    finally:
        # Limpiar directorio temporal
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)