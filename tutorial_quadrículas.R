# Extrapolar as áreas de distribuição e gerar mapas de riqueza
# Tutorial para extrapolar áreas de distribuição baseadas em quadrículas
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
# 

rm(list = ls())

dat_bel <- read_tsv(file = 'Belostomatidae.txt')
dat_bel

dat_clB <- dat_bel %>% filter(taxonRank %in% "SPECIES" | is.na(taxonRank))

dat_bel <- dat_bel %>%
  dplyr::select(species, decimalLongitude, decimalLatitude)

head(dat_bel)

dat_bel <- as.data.frame(dat_bel)
dat_bel

###### Iremos verificar se as coordenadas são válidas
cl <- cc_val(dat_bel, lat = "decimalLatitude", lon = "decimalLongitude")
cl

fl <- cc_equ(cl, value = "flagged", lat = "decimalLatitude", lon = "decimalLongitude")
fl

cl <- cl[fl, ]
cl

# tirando os NAs 
belos <- na.omit(cl)
belos

# Excluindo duplicados
pontos_dups <- duplicated(belos)
pontos_dups

# Mas quantas coordenadas são duplicadas?
length(which(pontos_dups==TRUE)) # 4003

# E quantas não são?
length(which(pontos_dups==FALSE)) # 11537

# Uma forma de limpar os dados
# gerar um novo objeto e guardar todas as coordenadas duplicadas
dupl <- which(pontos_dups==TRUE)
dupl # obs.: o que foi gravado? as coordenadas?

# Criar um novo objeto sem os duplicados
nodupl <- belos[-dupl,]
dim(nodupl) # 11537 pontos
################################################################
################################################################

# para termos uma matriz com as espécies como nomes de linhas:
head(nodupl)

belostoD <- data.frame(id = 1:nrow(nodupl), species = nodupl$species)
belostoV <- vect(nodupl, geom = c("decimalLongitude", "decimalLatitude"),
                 crs = "+proj=longlat +datum=WGS84")
belostoV
plot(belostoV)

# neotropical SHAPE FILE:
neo <- vect('regneotropical.shp')
neo
crs(neo) <- "+proj=longlat +datum=WGS84"
neo
plot(neo)

# seleção dos pontos para uma região especificada
pts_belosto_bons <- belostoV[neo, ]
plot(pts_belosto_bons, cex = 0.8, col = "green")

# pegando as coordenadas desse vetor:
belos_bons <- geom(pts_belosto_bons)[, c('x', 'y')]

# extraindo os atributos
labels <- as.data.frame(pts_belosto_bons)[, ('species')]

# criando a matriz final combinada
belos_bons_final <- data.frame(species = labels, belos_bons)


# Plotando os pontos sobre os mapas
points(belos_bons_final[, c(2, 3)], col = 'red', pch = 19)
dim(belos_bons_final) # 862 pontos

# mínimos polígonos convexos
source('calcRange_quadricula_MPC.R')
colnames(belos_bons_final) <- c('species', 'long', 'lat')
rang1 <- CalcRange_quadricula_MPC(x = belos_bons_final, shape_file = neo, resol = 5)
rang1$pres_abs
cores <- viridis(n = 10, alpha = 0.5, option = 'B')
plot(neo)
plot(rang1$geometry, col = cores, add =T)

## shapefiles dos MSTs
nomes <- unique(labels)
nomes
# MPC_shapes_n <- list.files('out_MPC/', pattern = "^mst_.*\\.shp$", full.names = T)
resul <- table(belos_bons_final$species)
names(resul)
resul[1]

MPC_shapes_n <- list()
pasta <- 'MPC_'
for(k in 1:length(nomes)){
  if(resul[k] < 3){
    next
  }
  MPC_shapes_n[[k]] <- vect(paste0('out_MPC/', pasta, names(resul)[k], '/MPC_', names(resul)[k], '.shp'))
}
MPC_shapes_n <- MPC_shapes_n[-which(sapply(MPC_shapes_n, is.null))]
MPC_shapes_n
for(i in 1:length(MPC_shapes_n)){
  plot(MPC_shapes_n[[i]], add = T)
}


# buffers de área definida
source('calcRange_buffer_quadricula.R')
colnames(belos_bons_final) <- c('species', 'long', 'lat')
rang2 <- CalcRangeBuffer_q(x = belos_bons_final, shape_file = neo, resol = 5, buffer.width = 500000)
plot(neo)
plot(rang2$geometry, add = T, col = cores)

BUFF_shapes_n <- list()
pasta <- 'BUFF_'
for(k in 1:length(nomes)){
  if(resul[k] < 3){
    next
  }
  BUFF_shapes_n[[k]] <- vect(paste0('out_buffers/', pasta, names(resul)[k], '/BUFF_', names(resul)[k], '.shp'))
}
BUFF_shapes_n <- BUFF_shapes_n[-which(sapply(BUFF_shapes_n, is.null))]
BUFF_shapes_n
for(i in 1:length(BUFF_shapes_n)){
  plot(BUFF_shapes_n[[i]], add = T)
}

# MSTs
source('calcRange_MST_quadricula.R')
colnames(belos_bons_final) <- c('species', 'long', 'lat')
rang3 <- CalcRangeMST_q(x = belos_bons_final, shape_file = neo, resol = 5)
plot(neo)
plot(rang3$geometry, add = T, col = cores)

## shapefiles dos MSTs
MST_shapes_n <- list.files('out_MST/', pattern = "^mst_.*\\.shp$", full.names = T)
MST_shapes_n
MST_shapes <- list()
for(k in 1:length(MST_shapes_n)){
  MST_shapes[[k]] <- vect(MST_shapes_n[k])
}
MST_shapes

for(i in 1:length(MST_shapes)){
  plot(MST_shapes[[i]], add = T)
}

# nexus
source('range_nexus.R')
Range_nexus(rang1$pres_abs, meto = 'MPC') # o arquivo estará disponível no diretório /out_mat

# BioGeoBears:
source('range_BioGeoBears.R')
Range_BioGeoBears(rang1$pres_abs, meto = 'MPC') # o arquivo estará disponível no diretório /out_mat

