from osgeo import gdal
import numpy as np
import sys, os
from tqdm import tqdm
import json

def find_envi_files(input_dir):
    envi_files = []
    for file in os.listdir(input_dir):
        if file.endswith(".hdr"):
            base_name = file[:-4]  # Remove .hdr extension
            for ext in [".dat", ".img", ".bin"]:  # Agrega aquí las extensiones posibles
                data_file = os.path.join(input_dir, base_name + ext)
                if os.path.exists(data_file):
                    envi_files.append(data_file)
    return envi_files

def save_layout(input_dir, stacks, stack_name):
    layout = {}
    band_count = 1
    for i, ds in enumerate(stacks):
        file_name = os.path.basename(ds.GetDescription())
        for j in range(ds.RasterCount):
            band_name = f"Band_{band_count}"
            layout[band_name] = {
                "file": file_name,
                "band_index": j + 1,
                "stack_index": band_count
            }
            band_count += 1
    
    layout_file = os.path.join(input_dir, f"{stack_name}_layout.json")
    with open(layout_file, 'w') as f:
        json.dump(layout, f, indent=2)
    print(f"Layout saved to {layout_file}")
    
    
def stackdestack(input_dir, stack_name):
    global tipo_dato
    envi_files = find_envi_files(input_dir)
    if not envi_files:
        print("No valid ENVI files found.")
        return

    stacks = []
    size = 0

    # Cargamos datasets y chequeamos si al menos hay uno válido
    for file in envi_files:
        ds = gdal.Open(file)
        ## extraemos el tipo de dato de la imagen
        tipo_dato = ds.GetRasterBand(1).DataType
        if ds is not None:
            stacks.append(ds)
        else:
            print(f"Couldn't open {file}")

    if not stacks:
        print("No valid ENVI datasets could be opened.")
        return

    # Usamos el primer dataset como referencia de tamaño
    ref_xsize = stacks[0].RasterXSize
    ref_ysize = stacks[0].RasterYSize

    # Solo usamos los que coincidan en tamaño
    valid_stacks = []
    for ds in stacks:
        if ds.RasterXSize == ref_xsize and ds.RasterYSize == ref_ysize:
            valid_stacks.append(ds)
        else:
            print(f"Skipping {ds.GetDescription()} due to incompatible size ({ds.RasterXSize}x{ds.RasterYSize} vs reference {ref_xsize}x{ref_ysize})")

    if not valid_stacks:
        print("No ENVI files with matching dimensions found.")
        return

    # Contamos las bandas totales de los datasets válidos
    size = sum([ds.RasterCount for ds in valid_stacks])

    driver = gdal.GetDriverByName("ENVI")
    stack_path = os.path.join(input_dir, stack_name)
    stack = driver.Create(stack_path, ref_xsize, ref_ysize, size, tipo_dato)
    pbar = tqdm(total=size)
    count = 1

    for ds in valid_stacks:
        for i in range(ds.RasterCount):
            band = ds.GetRasterBand(i + 1)
            arr = band.ReadAsArray()
            if arr.shape != (ref_ysize, ref_xsize):
                print(f"Error: Band shape {arr.shape} doesn't match stack size {(ref_ysize, ref_xsize)}. Skipping band.")
                continue
            stack.GetRasterBand(count).WriteArray(arr)
            count += 1
            pbar.update(1)

    stack.SetProjection(valid_stacks[0].GetProjection())
    stack.SetGeoTransform(valid_stacks[0].GetGeoTransform())
    stack.FlushCache()
    
    # Save the layout
    save_layout(input_dir, valid_stacks, stack_name)
    
    stack = None
    for ds in valid_stacks:
        ds = None
    pbar.close()
    
    # (Opcional: borra archivos si no son el stack)
    files = os.listdir(input_dir)
    for file in files:
        if os.path.splitext(os.path.basename(file))[0] != os.path.splitext(stack_name)[0] and not file.endswith('.csv') :
            os.remove(os.path.join(input_dir, file))
            
    print(f"Stack '{stack_path}' created successfully.")
    

def main():
    if len(sys.argv) != 3:
        print("Usage: python envi_stackdestack.py <stacks_dir> <stack_name>")
        sys.exit(1)
    input_dir = str(sys.argv[1])
    stack_name = str(sys.argv[2])
    stackdestack(input_dir, stack_name)

if __name__ == "__main__":
    main()