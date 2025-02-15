# Extrapolar as áreas de distribuição e gerar mapas de riqueza
# Rotina para estimar riqueza de espécies em áreas irregulares como domínios, biomas, países, cidades etc.
library(dplyr)
library(ggplot2)
library(devtools)
# install_github('ropensci/CoordinateCleaner', force = T)
library(CoordinateCleaner)
library(countrycode)
library(paleobioDB)
library(readr)
library(terra)
library(viridis)
library(tidyverse)
library(sf)
# install.packages("rgeos", repos="http://R-Forge.R-project.org")
# install.packages("rgdal", repos="http://R-Forge.R-project.org")
# devtools::install_github("azizka/speciesgeocodeR")
# install.packages('lwgeom')
# install.packages('tmap')

library(tmap)
library(lwgeom)  # Pacote necessário para st_make_valid()


library(speciesgeocodeR)
library(letsR)
library(dismo)
library(fossil)

# irregular bins: domínios, biomas, ecorregiões etc.
# MORRONE SHAPE FILE:
morrone <- vect('Americas_MORRONE_NEARTICO.shp') # objeto Spatial
crs(morrone) <- "+proj=longlat +datum=WGS84"
plot(morrone)

head(belos_bons_final)
bel_sf <- st_as_sf(x = belos_bons_final, coords = c('long', 'lat'))
bel_sf
st_crs(bel_sf) <- "+proj=longlat +datum=WGS84"
plot(bel_sf$geometry, add = T)

# devtools::install_github('ropensci/hoardr')
# devtools::install_github('wmgeolab/rgeoboundaries')

# library(rgeoboundaries)

head(belos_bons_final)
colnames(belos_bons_final) <- c('species', 'decimallongitude',
                                'decimallatitude')
morrone_sf <- as(morrone, 'Spatial')
plot(morrone_sf)

# Domínios
morrone_sf@data$Dominio

sp.class <- SpGeoCod(belos_bons_final, morrone_sf, areanames = "Dominio")
sp.class
summary(sp.class)
plot(sp.class)
plot(sp.class, type = "speciesrichness")
head(sp.class$spec_table) # tabela com as ocorrências
WriteOut(sp.class, type = "nexus") # produz um arquivo species_classification.nex
WriteOut(sp.class, type = 'BioGeoBEARS') # os nomes das espécies não apresentam '_',
# e o arquivo se chama BioGeoBEARS.txt

########################################################
# para termos um mapa:

# Transformar os pontos em um objeto espacial
belo_sf <- st_as_sf(belos_bons_final,
      coords = c('decimallongitude',
                 'decimallatitude'), crs = 4326) # Ajuste os nomes das colunas de coordenadas conforme necessário
belo_sf

morrone_sf <- st_read('Americas_MORRONE_NEARTICO.shp')
morrone_sf
pot_sf <- st_transform(belo_sf, st_crs(morrone_sf))
pot_sf

# Fazer a interseção para saber a qual bioma cada ocorrência pertence
pot_dominio <- st_join(pot_sf, morrone_sf)

# Verifica se há problemas de geometria
st_is_valid(morrone_sf) # existe

# Corrigir as geometrias inválidas
morrone_sf <- st_make_valid(morrone_sf)

# Verificar novamente
all(st_is_valid(morrone_sf))  # Deve retornar TRUE, mas não retornou


# você pode filtrar geometrias inválidas:
morrone_sf <- morrone_sf[st_is_valid(morrone_sf), ]

# Fazer a interseção para saber a qual bioma cada ocorrência pertence
pot_dominio <- st_join(pot_sf, morrone_sf)
pot_dominio
# Contar a diversidade por bioma (número de espécies únicas por bioma)
nspec_dominio <- pot_dominio %>%
  group_by(Dominio) %>%  # Ajuste para o nome correto da coluna de biomas
  summarise(nspec = n_distinct(species))

nspec_dominio

names(morrone_sf)

# Converter `morrone_sf` para um `data.frame` antes do join
morrone_data <- as.data.frame(morrone_sf) %>%
  dplyr::select(Dominio, geometry)
plot(morrone_data$geometry)

# Unir os dados ao shapefile dos dominios
sf_dominio <- left_join(morrone_data, nspec_dominio, by = "Dominio")
sf_dominio <- st_as_sf(sf_dominio)

# Visualizar no mapa
plot(sf_dominio["nspec"],
     main = "Diversidade de Espécies por Domínio",
     col = terrain.colors(10))

# Ativar modo de visualização interativa
tmap_mode("view")

# Criar mapa temático
tm_shape(sf_dominio) +
  tm_polygons("nspec", palette = "viridis", title = "Diversidade de Espécies") +
  tm_layout(legend.outside = TRUE)
