# =============================================================================
# Script: Crear Stack de Rasters y Exportar como NetCDF
# Descripción: Lee imágenes de un directorio, las ordena por fecha,
#              crea un stack, escala por 10000, convierte a int16 y
#              exporta como NetCDF con tiempo como coordenada.
# =============================================================================

# Cargar librerías necesarias
library(terra)


# =============================================================================
# FUNCIÓN PRINCIPAL
# =============================================================================

crear_stack_netcdf <- function(directorio_imagenes, 
                                patron_fecha = "\\d{8}",
                                formato_fecha = "%Y%m%d",
                                archivo_salida = "stack_output.nc",
                                factor_escala = 10000,
                                patron_archivos = "\\.tif$") {
  #' @param directorio_imagenes Ruta al directorio con las imágenes
  #' @param patron_fecha Patrón regex para extraer la fecha del nombre del archivo

  #' @param formato_fecha Formato de la fecha (por defecto YYYYMMDD)
  #' @param archivo_salida Nombre del archivo NetCDF de salida
  #' @param factor_escala Factor de escala (por defecto 10000)
  #' @param patron_archivos Patrón para filtrar archivos (por defecto .tif)
  
  cat("=== Iniciando procesamiento ===\n")
  
  # -----------------------------------------------------------------------------
  # 1. Listar y ordenar archivos por fecha
  # -----------------------------------------------------------------------------
  
  # Obtener lista de archivos
  archivos <- list.files(directorio_imagenes, 
                         pattern = patron_archivos, 
                         full.names = TRUE)
  
  if (length(archivos) == 0) {
    stop("No se encontraron archivos en el directorio especificado")
  }
  
  cat(sprintf("Encontrados %d archivos\n", length(archivos)))
  
  # Extraer fechas de los nombres de archivo
  nombres_archivos <- basename(archivos)
  fechas_str <- regmatches(nombres_archivos, regexpr(patron_fecha, nombres_archivos))
  
  if (length(fechas_str) != length(archivos)) {
    stop("No se pudieron extraer fechas de todos los archivos. Verificar el patrón de fecha.")
  }
  
  # Convertir a formato Date
  fechas <- as.Date(fechas_str, format = formato_fecha)
  
  # Crear data.frame y ordenar por fecha
  df_archivos <- data.frame(
    archivo = archivos,
    fecha = fechas,
    stringsAsFactors = FALSE
  )
  df_archivos <- df_archivos[order(df_archivos$fecha), ]
  
  cat("Archivos ordenados por fecha:\n")
  print(df_archivos)
  
  # -----------------------------------------------------------------------------
  # 2. Crear el stack de rasters
  # -----------------------------------------------------------------------------
  
  cat("\nCreando stack de rasters...\n")
  
  # Cargar todos los rasters en orden
  raster_stack <- rast(df_archivos$archivo)
  
  cat(sprintf("Stack creado con %d capas\n", nlyr(raster_stack)))
  cat(sprintf("Dimensiones: %d filas x %d columnas\n", nrow(raster_stack), ncol(raster_stack)))
  
  # -----------------------------------------------------------------------------
  # 3. Escalar y convertir a INT16
  # -----------------------------------------------------------------------------
  
  cat(sprintf("\nAplicando factor de escala: %d\n", factor_escala))
  
  # Multiplicar por el factor de escala
  raster_scaled <- raster_stack * factor_escala
  
  # Redondear y convertir a entero
  raster_int16 <- round(raster_scaled)
  
  # Establecer valores fuera de rango de INT16 como NA
  # INT16 rango: -32768 a 32767
  raster_int16[raster_int16 < -32768] <- NA
  raster_int16[raster_int16 > 32767] <- NA
  
  cat("Conversión a INT16 completada\n")
  
  # -----------------------------------------------------------------------------
  # 4. Asignar nombres de tiempo a las capas
  # -----------------------------------------------------------------------------
  
  # Asignar fechas como nombres de las capas
  names(raster_int16) <- as.character(df_archivos$fecha)
  
  # Establecer el tiempo en el raster
  time(raster_int16) <- df_archivos$fecha
  
  cat("\nFechas asignadas a las capas:\n")
  print(time(raster_int16))
  
  # -----------------------------------------------------------------------------
  # 5. Exportar como NetCDF
  # -----------------------------------------------------------------------------
  
  cat(sprintf("\nExportando a NetCDF: %s\n", archivo_salida))
  
  # Exportar como NetCDF con tipo de dato INT16
  writeCDF(raster_int16, 
           filename = archivo_salida,
           varname = "data",
           longname = "Datos escalados",
           unit = "scaled_units",
           missval = -32768,
           overwrite = TRUE)
  
  cat("=== Procesamiento completado exitosamente ===\n")
  cat(sprintf("Archivo guardado: %s\n", archivo_salida))
  
  # Retornar información del resultado
  return(list(
    stack = raster_int16,
    fechas = df_archivos$fecha,
    archivo_salida = archivo_salida,
    n_capas = nlyr(raster_int16)
  ))
}

# =============================================================================
# EJEMPLO DE USO
# =============================================================================

# Descomentar y modificar según sea necesario:

resultado <- crear_stack_netcdf(
  directorio_imagenes = "C:/ruta/a/tus/imagenes",
  patron_fecha = "\\d{8}",        # Busca 8 dígitos consecutivos (YYYYMMDD)
  formato_fecha = "%Y%m%d",       # Formato de fecha
  archivo_salida = "mi_stack.nc",
  factor_escala = 10000,
  patron_archivos = "\\.tif$"     # Archivos .tif
)

# =============================================================================
# EJEMPLOS DE PATRONES DE FECHA COMUNES
# =============================================================================
# 
# Patrón YYYYMMDD (ej: 20240115):
#   patron_fecha = "\\d{8}"
#   formato_fecha = "%Y%m%d"
#
# Patrón YYYY-MM-DD (ej: 2024-01-15):
#   patron_fecha = "\\d{4}-\\d{2}-\\d{2}"
#   formato_fecha = "%Y-%m-%d"
#
# Patrón YYYY_MM_DD (ej: 2024_01_15):
#   patron_fecha = "\\d{4}_\\d{2}_\\d{2}"
#   formato_fecha = "%Y_%m_%d"
#
# Patrón DD-MM-YYYY (ej: 15-01-2024):
#   patron_fecha = "\\d{2}-\\d{2}-\\d{4}"
#   formato_fecha = "%d-%m-%Y"
# =============================================================================
