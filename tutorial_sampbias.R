###################################################
############ rotina para sample bias ##############
library(sampbias)
library(sf)

# convertendo para um vetor:
belo_sf <- st_as_sf(x = belos_bons_final,
                   coords = c("decimallongitude",
                              "decimallatitude"))
st_crs(belo_sf) <- "+proj=longlat +datum=WGS84"

# neotropical SHAPE FILE:
neo_sf <- st_read('regneotropical.shp')
neo_sf
st_crs(neo_sf) <- "+proj=longlat +datum=WGS84"
pot_sf <- st_transform(belo_sf, st_crs(neo_sf))
pot_sf

plot(neo_sf$geometry)

# fazendo a exclusão de outros pontos fora da região:
belo_sf_mod <- st_contains_properly(neo_sf$geometry, belo_sf) # subamostrar os pontos terrestres
belo_sf_mod

# DOWNLOAD ENVIRONMENTAL VARIABLES ####
library(devtools)
library(sdmpredictors)
library(terra)
# library(raster)


# we'll use functions of the 'sdmpredictors' package to access different online datasets
pred_datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
pred_datasets
names(pred_datasets)
pred_datasets[ , 1:4]  # these are the datasets currently available for download using the 'sdmpredictors' package; you can check their URLs for more info on their contents
pred_datasets[ , c("dataset_code", "citation")]  # remember to ALWAYS cite the actual data sources, not just the package you used for downloading!

pred_layers <- list_layers(datasets = pred_datasets)
unique(pred_layers$dataset_code)
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$name)  # example of terrestrial variables dataset
# unique(pred_layers[pred_layers$dataset_code == "Freshwater", ]$name)  # example of marine variables dataset

# let's choose one dataset (e.g. WorldClim) and one particular set of variables (e.g. altitude and the bioclimatic ones, which are in rows 1 to 20):
layers_choice1 <- unique(pred_layers[pred_layers$dataset_code == "WorldClim", c("name", "layer_code")])
layers_choice1
layers_choice1 <- layers_choice1[c(1, 13:20), ]
layers_choice1

# Freshwater:
layers_choice2 <- unique(pred_layers[pred_layers$dataset_code == "Freshwater", c("name", "layer_code")])
layers_choice2
layers_choice2 <- layers_choice2[c(127, 128), ]
layers_choice2

# define folder for downloading / fetching the variables' map layers:
options(sdmpredictors_datadir = "../outputs/sdmpredictors")
# load the layers to the current R session (downloading them if they aren't already in the folder defined above):
layers1 <- load_layers(layers_choice1$layer_code, rasterstack = FALSE)  # rasterstack=TRUE gives error when there are layers with different extent
layers1  # a list of raster maps
# see how many elements in 'layers':
length(layers1)

layers2 <- load_layers(layers_choice2$layer_code, rasterstack = FALSE)  # rasterstack=TRUE gives error when there are layers with different extent
layers2  # a list of raster maps


# convert each layer to 'SpatRaster' class (from package 'terra'), which is much faster to process
layers1 <- lapply(layers1, rast)
plot(layers1[[1]], main = names(layers1)[1])
plot(layers1[[5]], main = names(layers1)[5])

# layers2 <- lapply(layers2, rast)
# plot(layers2[[1]], main = names(layers2)[1])
# plot(layers2[[2]], main = names(layers2)[2])

# find out if your layers have different extents or resolutions:
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$cellsize_lonlat)  # in this case 0.08333333 - spatial resolution can then be coarsened as adequate for your species data and modelling region (see below)
sapply(layers1, ext)  # if you get different extents (which doesn't happen with WorldClim, but may happen with other datasets), you'll have to crop all layers to the minimum common extent before proceeding
# for example, if the first layer has the smallest extent:
#layers <- lapply(layers, crop, extent(layers[[1]]))
# unique(pred_layers[pred_layers$dataset_code == "Freshwater", ]$cellsize_lonlat)  # in this case 0.08333333 - spatial resolution can then be coarsened as adequate for your species data and modelling region (see below)
# sapply(layers2, ext)  # if you get different extents (which doesn't happen with WorldClim, but may happen with other datasets), you'll have to crop all layers to the minimum common extent before proceeding

# ATENÇÃO: AS LAYERS 2 APRESENTAM EXTENSÃO MENOR, MAS POSSUEM RESOLUÇÃO MAIS FINA!
layers1 <- lapply(layers1, crop, ext(layers2[[1]]))
sapply(layers1, ext)  # if you get different extents (which doesn't happen with WorldClim, but may happen with other datasets), you'll have to crop all layers to the minimum common extent before proceeding
plot(layers1[[1]])
# plot(layers2[[1]])

# # MORRONE SHAPE FILE:
# # morrone <- vect(paste0(diret1, '/morroneNew/NeotropicMap_Geo.shp'))
# morro_vect <- vect(paste0(diret1, '/morroneNew/NeotropicMap_Geo.shp'))
# # Asul_topo <- mask(crop(global_topo, morrone), morrone)
# # plot(Asul_topo)
datum <- "+proj=longlat +datum=WGS84"
# morro_vect <- project(morro_vect, datum)

camadas_cut1 <- rast(layers1)
plot(camadas_cut1)
# camadas_cut2 <- rast(layers2)

# camadas_cut2 <- mask(crop(camadas_cut2, morro_vect), morro_vect)
camadas_cut1 <- mask(crop(camadas_cut1, neo_sf), neo_sf)
plot(camadas_cut1[[1]])
# plot(camadas_cut2[[1]])

camadas_cut1 <- project(camadas_cut1, datum)
# camadas_cut2 <- project(camadas_cut2, datum)

camadas_cut1[[1]]@ptr$range_max # valor máximo

# pontos:
plot(neo_sf$geometry, border = 'tan')
plot(belo_sf[unlist(belo_sf_mod), ], col = 'darkblue', cex = 0.4, add = T)

# ------------- CV chuva -----------------
# plot(camadas_cut1[[5]])
CVLayerPrec <- camadas_cut1[[5]]
CVLayerPrec <- aggregate(CVLayerPrec, fact = 5)
plot(CVLayerPrec)

# threshold our CV raster into bands, and specify a more convenient
# 25%
thresh_25 <- classify(CVLayerPrec,
                      rcl = rbind(c(-Inf, camadas_cut1[[5]]@ptr$range_max*0.25, 1),
                                  c(camadas_cut1[[5]]@ptr$range_max*0.25, Inf,  2)))
plot(thresh_25)

# convert our thresholded raster to a set of sf polygons
CV_25 <- st_as_sf(as.points(thresh_25))
CV_25 <- st_cast(CV_25, "MULTIPOINT")
library(viridis)
plot(CV_25)
st_crs(CV_25) <- "+proj=longlat +datum=WGS84"
CV_25$WC_bio15[CV_25$WC_bio15==2] <- NA
plot(CV_25)

CV_v_25 <- vect(CV_25)
CV_v_25
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

######################################
# threshold our CV raster into bands, and specify a more convenient
# 50%
thresh_50 <- classify(CVLayerPrec,
                      rcl = rbind(c(-Inf, camadas_cut1[[5]]@ptr$range_max*0.50, 1),
                                  c(camadas_cut1[[5]]@ptr$range_max*0.50, Inf,  2)))
plot(thresh_50)

# convert our thresholded raster to a set of sf polygons
CV_50 <- st_as_sf(as.points(thresh_50))
CV_50 <- st_cast(CV_50, "MULTIPOINT")
library(viridis)
plot(CV_50)
st_crs(CV_50) <- "+proj=longlat +datum=WGS84"
CV_50$WC_bio15[CV_50$WC_bio15==2] <- NA
plot(CV_50)

CV_v_50 <- vect(CV_50)
CV_v_50
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

#############################
# threshold our CV raster into bands, and specify a more convenient
# 50%
thresh_75 <- classify(CVLayerPrec,
                      rcl = rbind(c(-Inf, camadas_cut1[[5]]@ptr$range_max*0.75, 1),
                                  c(camadas_cut1[[5]]@ptr$range_max*0.75, Inf,  2)))
plot(thresh_75)

# convert our thresholded raster to a set of sf polygons
CV_75 <- st_as_sf(as.points(thresh_75))
CV_75 <- st_cast(CV_75, "MULTIPOINT")
library(viridis)
plot(CV_75)
st_crs(CV_75) <- "+proj=longlat +datum=WGS84"
CV_75$WC_bio15[CV_75$WC_bio15==2] <- NA
plot(CV_75)

CV_v_75 <- vect(CV_75)
CV_v_75
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

# ---------------- pluviosidade anual -------------------
pluviAnual <- camadas_cut1[[2]]
plot(pluviAnual)
pluviAnual <- aggregate(pluviAnual, fact = 5)

#############################
# threshold our CV raster into bands, and specify a more convenient
# 25%
thresh_25 <- classify(pluviAnual,
                      rcl = rbind(c(-Inf, camadas_cut1[[2]]@ptr$range_max*0.25, 1),
                                  c(camadas_cut1[[2]]@ptr$range_max*0.25, Inf,  2)))
plot(thresh_25)

# convert our thresholded raster to a set of sf polygons
pluviAnual_25 <- st_as_sf(as.points(thresh_25))
pluviAnual_25 <- st_cast(pluviAnual_25, "MULTIPOINT")
library(viridis)
plot(pluviAnual_25)
st_crs(pluviAnual_25) <- "+proj=longlat +datum=WGS84"
pluviAnual_25$WC_bio12[pluviAnual_25$WC_bio12==2] <- NA
plot(pluviAnual_25)

pluviAnual_v_25 <- vect(pluviAnual_25)
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

# ---------------- altitude -------------------
altitude <- camadas_cut1[[1]]
plot(altitude)
altitude <- aggregate(altitude, fact = 5)

#############################
# threshold our CV raster into bands, and specify a more convenient
# 25%
thresh_25 <- classify(altitude,
                      rcl = rbind(c(-Inf, camadas_cut1[[1]]@ptr$range_max*0.25, 1),
                                  c(camadas_cut1[[1]]@ptr$range_max*0.25, Inf,  2)))
plot(thresh_25)

# convert our thresholded raster to a set of sf polygons
alt_25 <- st_as_sf(as.points(thresh_25))
alt_25 <- st_cast(alt_25, "MULTIPOINT")
library(viridis)
plot(alt_25)
st_crs(alt_25) <- "+proj=longlat +datum=WGS84"
alt_25$WC_alt[alt_25$WC_alt==2] <- NA
plot(alt_25)

alt_v_25 <- vect(alt_25)
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

#############################
# threshold our CV raster into bands, and specify a more convenient
# 50%
thresh_50 <- classify(altitude,
                      rcl = rbind(c(-Inf, camadas_cut1[[1]]@ptr$range_max*0.50, 1),
                                  c(camadas_cut1[[1]]@ptr$range_max*0.50, Inf,  2)))
plot(thresh_50)

# convert our thresholded raster to a set of sf polygons
alt_50 <- st_as_sf(as.points(thresh_50))
alt_50 <- st_cast(alt_50, "MULTIPOINT")
library(viridis)
plot(alt_50)
st_crs(alt_50) <- "+proj=longlat +datum=WGS84"
alt_50$WC_alt[alt_50$WC_alt==2] <- NA
plot(alt_50)

alt_v_50 <- vect(alt_50)
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)


#############################
# threshold our CV raster into bands, and specify a more convenient
# 75%
thresh_75 <- classify(altitude,
                      rcl = rbind(c(-Inf, camadas_cut1[[1]]@ptr$range_max*0.75, 1),
                                  c(camadas_cut1[[1]]@ptr$range_max*0.75, Inf,  2)))
plot(thresh_75)

# convert our thresholded raster to a set of sf polygons
alt_75 <- st_as_sf(as.points(thresh_75))
alt_75 <- st_cast(alt_75, "MULTIPOINT")
library(viridis)
plot(alt_75)
st_crs(alt_75) <- "+proj=longlat +datum=WGS84"
alt_75$WC_alt[alt_75$WC_alt==2] <- NA
plot(alt_75)

alt_v_75 <- vect(alt_75)
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

#############################
# threshold our CV raster into bands, and specify a more convenient
# 10%
thresh_10 <- classify(altitude,
                      rcl = rbind(c(-Inf, camadas_cut1[[1]]@ptr$range_max*0.10, 1),
                                  c(camadas_cut1[[1]]@ptr$range_max*0.10, Inf,  2)))
plot(thresh_10)

# convert our thresholded raster to a set of sf polygons
alt_10 <- st_as_sf(as.points(thresh_10))
alt_10 <- st_cast(alt_10, "MULTIPOINT")
library(viridis)
plot(alt_10)
st_crs(alt_10) <- "+proj=longlat +datum=WGS84"
alt_10$WC_alt[alt_10$WC_alt==2] <- NA
plot(alt_10)

alt_v_10 <- vect(alt_10)
# CV = as.data.frame(CV, xy = TRUE); sp::coordinates(CV) = ~x+y
# plot(CV_v)

# camadas gazeteers
getwd()
cit <- vect("sampbias/cities_gaz.shp")
# cit_sf <- st_read('sampbias/cities_gaz.shp')
# plot(cit_sf$geometry)
cit <- crop(cit, ext(neo_sf))
# cit_sf <- crop(cit_sf$geometry, ext(neo_sf$geometry))
# cit_sf <- st_contains_properly(cit_sf$geometry, neo_sf$geometry)
datum <- "+proj=longlat +datum=WGS84"
cit <- project(cit, datum)
plot(cit)

roads <- vect('sampbias/roads_gaz.shp')
roads
roads <- crop(roads, ext(neo_sf))
datum <- "+proj=longlat +datum=WGS84"
roads <- project(roads, datum)
plot(roads)

rivers <- vect("sampbias/rivers_gaz.shp")
rivers
rivers <- crop(rivers, ext(neo_sf))
datum <- "+proj=longlat +datum=WGS84"
rivers <- project(rivers, datum)
plot(rivers)

# coordenadas modificadas:
belo_sf_mod <- belo_sf[unlist(belo_sf_mod), ]

gazetteers <- list(cities = cit, roads = roads,
                   pluviCV = CV_v_25, rivers = rivers,
                   pluviYear = pluviAnual_v_25, alt_10 = alt_v_10,
                   alt_25 = alt_v_25)
gazetteers

gazetteers_1 <- list(cities = cit, roads = roads,
                     rivers = rivers)
gazetteers_1

# diret4 <- "D:/Dropbox/papers_books/mapas_paleoBot/josi_bot"
# usoTerra <- rast(paste0(diret4, '/human-modification/lulc-human-modification-terrestrial-systems_geographic.tif'))
# plot(usoTerra, col = viridis(100))
# res(usoTerra)

# prepare data
doms <- as.data.frame(st_coordinates(belo_sf_mod))
doms
colnames(doms) <- c("decimalLongitude", "decimalLatitude")
doms$species <- belo_sf_mod$species
doms

# make our raster
dmr <- crop(rast(resolution = 3), neo_sf)
dmr

# run the calculator, supplying our SpatVectors as a named list, and a buffer of
# zero as we do not need to worry about the raster edges. This may take 10 minutes to run
# options("download.file.method" = "libcurl")
# spatially project the modeled effects of our sampling biases
# devtools::install_github("ropensci/rnaturalearthhires")
# devtools::install_github("ropensci/rnaturalearth")
# devtools::install_github("ropensci/rnaturalearthdata")
# rnaturalearth::ne_download(scale = "medium", type = "urban_areas", category = "cultural", returnclass = "sf")

testbias1 <- calculate_bias(doms, gaz = gazetteers,
                           inp_raster = dmr, buffer = 0)

testbias2 <- calculate_bias(doms, gaz = gazetteers_1,
                           inp_raster = dmr, buffer = 0)

# library("rnaturalearth")
project_bias(testbias1)
map_bias(project_bias(testbias1), type = 'log_sampling_rate')
# map_bias(project_bias(testbias2))
summary(testbias2)
# plot(testbias2)

proj1 <- project_bias(testbias1)
plot(proj1)


# calculate the relative strengths of the factors at 0.1 to 10 km ranges. The
# code is just the internals of plot.sampbias()
means <- colMeans(testbias$bias_estimate)
means
dists <- seq(0.1, 10, length.out = 1000)
r_road  <- means["q"] * exp(-means["w_road"]  * dists)
r_rivers <- means["q"] * exp(-means["w_rivers"] * dists)
r_cities <- means["q"] * exp(-means["w_cities"] * dists)
r_pluriAnual <- means["q"] * exp(-means["w_pluviYear"] * dists)
r_alt_10 <- means["q"] * exp(-means["w_alt_10"] * dists)
r_alt_25 <- means["q"] * exp(-means["w_alt_25"] * dists)
r_pluviCV <- means['q'] * exp(-means['w_pluviCV'] * dists)

# plot summaries
cols <- c("#440154FF", "#FDE725FF", "red", 'green', 'grey', 'orange', 'pink')
par(mfrow = c(2, 1))
boxplot(testbias$bias_estimate[,c("w_cities", "w_roads", 
                                  "w_rivers", "w_pluviYear", "w_alt_10", "w_alt_25", "w_pluviCV")],
        ylab = "Posterior weight",
        col = cols, border = cols, names = c("Cities", "Roads", "Rivers", "PluviYear",
                                             "alt10", "alt25", "pluviCV"))

plot(dists, r_cities, type = "l", lwd = 3, col = cols[1], xlab = "Distance (km)",
     ylab = "Sampling rate")
lines(dists, r_road, lwd = 3, col = cols[2])
lines(dists, r_rivers, lwd = 3, col = cols[3])
lines(dists, r_pluriAnual, lwd = 3, col = cols[4])
lines(dists, r_alt_10, lwd = 3, col = cols[5])
lines(dists, r_alt_25, lwd = 3, col = cols[6])
lines(dists, r_pluviCV, lwd = 3, col = cols[7])
legend("topright", c("Cities", "Roads", "Rivers", "pluviA", "Alt10",
                     "Alt25", "pluviCV"), lty = 1,
       lwd = 3, bty = "n", col = cols)


