import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

def xls_to_gdf(ruta_xls):
    """Funcione que convierte nuestro xls a geodatabase

    Args:
        ruta_xls (str): Ruta del archivo xls

    Returns:
        gpd.GeoDataFrame: GeoDataFrame con la geometría creada
    """
    # Leer el archivo Excel
    df = pd.read_excel(ruta_xls)

    # Crear geometría de puntos
    geometry = [Point(xy) for xy in zip(df['XETRS89_H30'], df['YETRS89_H30'])]

    # Crear GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:25830")

    return gdf

def guardar_gdf(ruta_xls,gdf,tipo):
    """Guarda el geodata en el tipo de archivo que se desee

    Args:
        ruta_xls (str): Ruta del archivo xls
        gdf (gpd.GeoDataFrame): GeoDataFrame a guardar
        tipo (str): Tipo de archivo de salida ('shp', 'geojson', 'gpkg', 'gdb')

    Raises:
        ValueError: Si el tipo no es soportado
    """
    if tipo == 'shp':
        ruta_salida = ruta_xls.replace('.xlsx', '.shp')
        gdf.to_file(ruta_salida, driver='ESRI Shapefile', encoding='utf-8')
    elif tipo == 'geojson':
        ruta_salida = ruta_xls.replace('.xlsx', '.geojson')
        gdf.to_file(ruta_salida, driver='GeoJSON', encoding='utf-8')
    elif tipo == 'gpkg':
        ruta_salida = ruta_xls.replace('.xlsx', '.gpkg')
        gdf.to_file(ruta_salida, driver='GPKG', encoding='utf-8')
    elif tipo == 'gdb':
        ruta_salida = ruta_xls.replace('.xlsx', '.gdb')
        gdf.to_file(ruta_salida, driver='FileGDB', encoding='utf-8')
    else:
        raise ValueError("Tipo no soportado. Use 'shp', 'geojson', 'gpkg' o 'gdb'.")

def main():
    ruta_xls = 'ruta/al/archivo.xlsx' # Reemplazar con la ruta real del archivo Excel que nos paso Cristina
    tipo = 'geojson'  # Reemplazar con la ruta real del archivo Excel
    gdf = xls_to_gdf(ruta_xls)
    
    # Guardar en diferentes formatos
    guardar_gdf(ruta_xls, gdf, tipo)

if __name__ == "__main__":
    main()
