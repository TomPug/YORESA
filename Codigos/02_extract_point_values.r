#' Extraer y preparar series temporales de coherencia SAR
#'
#' @param raster_path Ruta al archivo raster de coherencia (stack temporal)
#' @param puntos_path Ruta al archivo shapefile/geojson con los puntos
#' @param output_dir Directorio donde guardar los CSV (por defecto: mismo directorio que los puntos)
#' @param output_prefix Prefijo para los nombres de los archivos CSV (por defecto: "coherencia")
#' @param band_name Nombre de la banda (por defecto "Coherence")
#' @param remove_na Lógico, si TRUE elimina valores NA en formato largo (por defecto TRUE)
#' @param plot_map Lógico, si TRUE genera un mapa de los puntos (por defecto FALSE)
#' @param verbose Lógico, si TRUE muestra información del proceso (por defecto TRUE)
#'
#' @return Lista con dos elementos: valores_ancho (formato original) y valores_largo (formato largo)
#' 
#' @examples
#' resultado <- extraer_coherencia_temporal(
#'   raster_path = "E:/SAR_UVa/coherence_stack.nc",
#'   puntos_path = "E:/SAR_UVa/shps/puntos_ifn_extremadura.geojson",
#'   output_dir = "E:/SAR_UVa/resultados",
#'   plot_map = TRUE
#' )
#'
#' Extraer y preparar series temporales de coherencia SAR
#'
#' @param raster_path Ruta al archivo raster de coherencia (stack temporal)
#' @param puntos_path Ruta al archivo shapefile/geojson con los puntos
#' @param parcela_field Nombre del campo en el GeoJSON que contiene el ID de parcela (por defecto "parcela")
#' @param output_dir Directorio donde guardar los CSV (por defecto: mismo directorio que los puntos)
#' @param output_prefix Prefijo para los nombres de los archivos CSV (por defecto: "coherencia")
#' @param band_name Nombre de la banda (por defecto "Coherence")
#' @param remove_na Lógico, si TRUE elimina valores NA en formato largo (por defecto TRUE)
#' @param plot_map Lógico, si TRUE genera un mapa de los puntos (por defecto FALSE)
#' @param verbose Lógico, si TRUE muestra información del proceso (por defecto TRUE)
#'
#' @return Lista con dos elementos: valores_ancho (formato original) y valores_largo (formato largo)
#' 
#' @examples
#' resultado <- extraer_coherencia_temporal(
#'   raster_path = "E:/SAR_UVa/coherence_stack.nc",
#'   puntos_path = "E:/SAR_UVa/shps/puntos_ifn_extremadura.geojson",
#'   parcela_field = "PARCELA",
#'   output_dir = "E:/SAR_UVa/resultados",
#'   plot_map = TRUE
#' )
#'
extraer_coherencia_temporal <- function(raster_path, 
                                        puntos_path,
                                        parcela_field = "PARCELA",
                                        output_dir = NULL,
                                        output_prefix = "coherencia",
                                        band_name = "Coherence",
                                        remove_na = TRUE,
                                        plot_map = FALSE,
                                        verbose = TRUE) {
  
  # Verificar que las librerías necesarias están cargadas
  required_packages <- c("terra", "sf", "dplyr", "tidyr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if(length(missing_packages) > 0) {
    stop(paste("Instala los siguientes paquetes:", paste(missing_packages, collapse = ", ")))
  }
  
  # Cargar librerías
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  
  if(verbose) cat("=== Extrayendo coherencia temporal ===\n\n")
  
  # Definir directorio de salida si no se especifica
  if(is.null(output_dir)) {
    output_dir <- dirname(puntos_path)
  }
  
  # Crear directorio si no existe
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if(verbose) cat("Directorio creado:", output_dir, "\n\n")
  }
  
  # 1. Cargar datos
  if(verbose) cat("1. Cargando raster...\n")
  if(!file.exists(raster_path)) stop("El archivo raster no existe: ", raster_path)
  coh_stack <- rast(raster_path)
  
  if(verbose) cat("2. Cargando puntos...\n")
  if(!file.exists(puntos_path)) stop("El archivo de puntos no existe: ", puntos_path)
  puntos <- st_read(puntos_path, quiet = !verbose)
  
  # Verificar que existe el campo parcela
  if(!parcela_field %in% names(puntos)) {
    stop(paste0("El campo '", parcela_field, "' no existe en el GeoJSON. Campos disponibles: ", 
                paste(names(puntos), collapse = ", ")))
  }
  
  if(verbose) {
    cat("   - Capas en el raster:", nlyr(coh_stack), "\n")
    cat("   - Puntos cargados:", nrow(puntos), "\n")
    cat("   - Campo de parcela usado:", parcela_field, "\n\n")
  }
  
  # 2. Visualización opcional
  if(plot_map) {
    if(verbose) cat("3. Generando mapa...\n")
    plot(coh_stack[[1]], main = "Coherencia - Primera capa")
    plot(st_geometry(puntos), add = TRUE, col = "red", pch = 16, cex = 0.8)
  }
  
  # 3. Extraer valores
  if(verbose) cat("4. Extrayendo valores en los puntos...\n")
  valores <- terra::extract(coh_stack, puntos)
  
  # 4. Obtener fechas
  if(verbose) cat("5. Obteniendo fechas...\n")
  fechas <- time(coh_stack)
  
  # Verificar fechas
  if(all(is.na(fechas))) {
    warning("El raster no tiene fechas asignadas. Usando índices numéricos.")
    fechas <- 1:nlyr(coh_stack)
  }
  
  if(verbose) {
    cat("   - Primera fecha:", as.character(fechas[1]), "\n")
    cat("   - Última fecha:", as.character(fechas[length(fechas)]), "\n\n")
  }
  
  # 5. Preparar valores en formato ancho
  if(verbose) cat("6. Preparando datos en formato ancho...\n")
  
  # Añadir el campo PARCELA del GeoJSON (no el ID automático de extract)
  valores_ancho <- valores %>%
    mutate(PARCELA = puntos[[parcela_field]]) %>%
    select(-ID)  # Eliminar el ID automático que genera extract
  
  # Renombrar columnas con fechas
  col_names <- names(valores_ancho)
  col_coherencia <- grep("__xarray_dataarray_variable__", col_names, value = TRUE)
  
  if(length(col_coherencia) > 0) {
    # Crear nombres de columnas con fechas
    nuevos_nombres <- paste0("coh_", format(fechas, "%Y%m%d"))
    names(valores_ancho)[names(valores_ancho) %in% col_coherencia] <- nuevos_nombres
  }
  
  # Reordenar para que PARCELA sea la primera columna
  valores_ancho <- valores_ancho %>%
    select(PARCELA, everything())
  
  # 6. Guardar formato ancho
  archivo_ancho <- file.path(output_dir, paste0(output_prefix, "_ancho.csv"))
  if(verbose) cat("7. Guardando formato ancho:", archivo_ancho, "\n")
  write.csv(valores_ancho, archivo_ancho, row.names = FALSE)
  
  # 7. Transformar a formato largo
  if(verbose) cat("8. Transformando a formato largo...\n")
  
  # Si las columnas ya fueron renombradas
  if(length(col_coherencia) > 0) {
    valores_largo <- valores_ancho %>%
      pivot_longer(
        cols = starts_with("coh_"),
        names_to = "fecha_str",
        values_to = "coherencia"
      ) %>%
      mutate(
        fecha = as.Date(gsub("coh_", "", fecha_str), format = "%Y%m%d"),
        band = band_name
      ) %>%
      select(PARCELA, fecha, coherencia, band) %>%
      arrange(PARCELA, fecha)
  } else {
    # Método original si no hay columnas de coherencia
    valores_largo <- valores_ancho %>%
      pivot_longer(
        cols = -PARCELA,
        names_to = "banda",
        values_to = "coherencia"
      ) %>%
      group_by(PARCELA) %>%
      mutate(
        banda_num = row_number(),
        fecha = fechas[banda_num],
        band = band_name
      ) %>%
      ungroup() %>%
      select(PARCELA, fecha, coherencia, band) %>%
      arrange(PARCELA, fecha)
  }
  
  # 8. Eliminar NAs si se solicita (solo en formato largo)
  n_antes <- nrow(valores_largo)
  if(remove_na) {
    if(verbose) cat("9. Eliminando valores NA del formato largo...\n")
    valores_largo <- valores_largo %>%
      filter(!is.na(coherencia))
    n_despues <- nrow(valores_largo)
    
    if(verbose) cat("   - Valores eliminados:", n_antes - n_despues, "\n\n")
  } else {
    n_despues <- n_antes
  }
  
  # 9. Guardar formato largo
  archivo_largo <- file.path(output_dir, paste0(output_prefix, "_largo.csv"))
  if(verbose) cat("10. Guardando formato largo:", archivo_largo, "\n")
  write.csv(valores_largo, archivo_largo, row.names = FALSE)
  
  # 10. Resumen final
  if(verbose) {
    cat("\n=== Resumen de resultados ===\n")
    cat("Archivos guardados:\n")
    cat("  - Formato ancho:", archivo_ancho, "\n")
    cat("  - Formato largo:", archivo_largo, "\n\n")
    cat("Estadísticas formato largo:\n")
    cat("  Total de observaciones:", n_despues, "\n")
    cat("  Parcelas únicas:", length(unique(valores_largo$PARCELA)), "\n")
    cat("  Fechas únicas:", length(unique(valores_largo$fecha)), "\n")
    cat("  Rango de coherencia:", 
        paste(round(range(valores_largo$coherencia, na.rm = TRUE), 3), collapse = " - "), "\n")
    cat("  Media de coherencia:", round(mean(valores_largo$coherencia, na.rm = TRUE), 3), "\n\n")
  }
  
  # 11. Retornar ambos formatos
  return(invisible(list(
    valores_ancho = valores_ancho,
    valores_largo = valores_largo,
    archivos = list(
      ancho = archivo_ancho,
      largo = archivo_largo
    )
  )))
}
## tiempo de ejecucion

resultado <- extraer_coherencia_temporal(
  raster_path = "E:/SAR_UVa/coherence_stack.nc",
  puntos_path = "E:/SAR_UVa/shps/puntos_ifn_extremadura.geojson",
  output_dir = "E:/SAR_UVa/CSV'S",
  output_prefix = "valores_coherencia_extremadura",
  plot_map = TRUE
)
