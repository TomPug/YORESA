library(terra)
library(sf)


coh_stack <- rast("E:/SAR_UVa/coherence_stack.nc")
puntos <- st_read("E:/SAR_UVa/shps/puntos_ifn_extremadura.geojson")
plot(puntos)


valores <- extract(coh_stack,puntos)


valores <- na.omit(valores)
valores$ID_PARCELA <- valores$ID

library(terra)
library(sf)
library(tidyr)
library(dplyr)

# Ya tienes cargados coh_stack, puntos y valores

# 1. Obtener las fechas del raster
fechas <- time(coh_stack)

# Verificar que las fechas no sean NA
if(all(is.na(fechas))) {
  stop("El raster no tiene fechas. Necesitas asignarlas primero.")
}

# 2. Transformar de formato ancho a largo
valores_largo <- valores %>%
  # Pivotar todas las columnas excepto ID y ID_PARCELA
  pivot_longer(
    cols = starts_with("__xarray_dataarray_variable__"),
    names_to = "banda",
    values_to = "coherencia"
  ) %>%
  # Extraer el número de banda
  mutate(
    banda_num = as.integer(gsub(".*_(\\d+)$", "\\1", banda))
  ) %>%
  # Agregar las fechas
  mutate(
    fecha = fechas[banda_num]
  ) %>%
  # Seleccionar y ordenar columnas
  select(ID_PARCELA, fecha, coherencia) %>%
  # Ordenar por ID y fecha
  arrange(ID_PARCELA, fecha)

# Ver el resultado
head(valores_largo, 20)

valores_largo$band <- "Coherence"
valores_largo$SP_PPAL <- 1

# Resumen
cat("Dimensiones:", nrow(valores_largo), "filas x", ncol(valores_largo), "columnas\n")
cat("Número de parcelas:", n_distinct(valores_largo$ID_PARCELA), "\n")
cat("Rango de fechas:", range(valores_largo$fecha), "\n")

write.csv(valores,"E:/SAR_UVa/CSV'S/valores_coherencia_extremadura")
write.csv(valores_largo,"E:/SAR_UVa/CSV'S/valores_coherencia_extremadura_modificado")
