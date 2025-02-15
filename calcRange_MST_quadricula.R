# Range baseado num MST
CalcRangeMST_q <- function(x, shape_file, resol) {
  # projection
  wgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  warning("Assuming lat/long wgs84 coordinates")
  
  # check for geosphere package
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package 'geosphere' not found. Please install the package.", call. = FALSE)
  }
  
  # fix different input data types data.frame
  if (is.data.frame(x)) {
    names(x) <- tolower(names(x))
    dat <- x[, c("species", "long", "lat")]
  }
  ## spgeoOUt
  if (is.spgeoOUT(x)) {
    dat <- x$samples[, 1:3]
  }
  
  # check for species with less than 3 records
  filt <- table(dat$species)
  # filt
  sortout <- names(filt[filt <= 2])
  # sortout
  filt <- filt[filt > 2]
  
  dat.filt <- droplevels(subset(dat, dat$species %in% as.character(names(filt))))
  
  # check for species where all lat or long ar identical, or almost identical,
  # to prevent line polygons longitude
  test <- split(dat.filt, f = dat.filt$species)
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallongitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))
  
  # latitude
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallatitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))
  # 
  # # test for almost perfect fit
  # test2 <- sapply(test, function(k) {
  #   round(abs(cor(k[, "decimallongitude"], k[, "decimallatitude"])), 6)
  # })
  # sortout2 <- names(test2[test2 == 1])
  # sortout <- c(sortout, sortout2)
  # dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))
  # 
  # sortout <- sortout[!is.na(sortout)]
  # 
  # if (length(sortout) > 0) {
  #   warning("found species with < 3 occurrences:", paste("\n", sortout))
  # }
  
  # calculate buffers:
  # inp <- split(dat.filt, f = dat.filt$species)
  
  # # test for occurrences spanning > 180 degrees
  # test <- lapply(inp, function(k) {
  #   SpatialPoints(k[, 2:3])
  # })
  # test <- lapply(test, "extent")
  # test <- lapply(test, function(k) {
  #   (k@xmax + 180) - (k@xmin + 180)
  # })
  # test <- unlist(lapply(test, function(k) {
  #   k >= 180
  # }))
  # if (any(test)) {
  #   stop("data includes species spanning >180 degrees.")
  # }
  
  print('1) preparing the coordinates...')
  coordin <- matrix(as.matrix(dat.filt[, c(2, 3)]), nrow(dat.filt),
                    2, dimnames = list(dat.filt[,1], colnames(dat.filt)[c(2, 3)]))
  
  cols1 <- setNames(viridis(n = length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))
  
  # preparing species' msts...
  #######
  ### vetor de acúmulo de espécimes ###
  species <- unique(rownames(coordin))
  qde <- 0
  for(i in species){qde[i] <- length(which(rownames(coordin) == i))}
  qde <- cumsum(qde)
  tabela <- matrix(NA, nrow = dim(coordin), nc = 4)
  #######
  
  print('2) calculating mst...')
  tempo <- coordin
  
  # tempo <- as.data.frame(na.exclude(tempo))
  # colnames(tempo) <- c('species', colnames(coordin[, c(1, 2)]))
  
  Long <- tempo[,1]
  Lat <- tempo[, 2]
  names(Long) <- rownames(tempo)
  
  tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(tempo), 2, dimnames = list(rownames(tempo),
                                        colnames(tempo)[c(1,2)]))
  #######
  # preparando o raster com grid do shapefile escolhido
  grid <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                 crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  #########
  # producing a raster of the shapefile
  mask.raster <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(as(shape_file, 'Spatial'), mask.raster)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # mask.raster[is.na(mask.raster)] <- 0
  r <- merge(r, mask.raster)
  ncellras <- ncell(r)
  
  dir.create('out_MST/') # temporary folder
  
  ##############
  conta <- 0
  coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))),
                   dimnames = list(seq(1:ncellras), unique(rownames(tempo)))
                   [c(1:2)])
  lista_r <- list()
  
  for(j in unique(rownames(tempo))){
    conta <- conta + 1
    tempoo <- subset(tempo, rownames(tempo) == j)
    
    #### shapefile ###
    tempo.d <- as.data.frame(tempoo)
    tempo_shape <- lats2Shape(lats = tempo.d)
    write.shapefile(tempo_shape, paste0(c('out_MST/pointshape_', j), collapse = ''))
    resul1_shape <- vect(paste0(c('out_MST/pointshape_', j, '.shp'),
                                collapse = ''), crs = "+proj=longlat +datum=WGS84")
    # resul1_shape <- readShapeSpatial('tempshape1_out_MST.shp') # shapefile
    # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    terra::project(resul1_shape, "+proj=longlat +datum=WGS84")
    # resul1_shape <- vect(resul1_shape, crs = "+proj=longlat")
    
    ##### MST based on geographic distance #####
    # rownames(resul1_shape@coords) <- rownames(tempo.d)
    # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
    dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
    mst2 <- dino.mst(dista)
    rownames(mst2) <- rownames(tempo.d)
    colnames(mst2) <- rownames(tempo.d)
    lats <- cbind(resul1_shape$long,
                  resul1_shape$lat)
    rownames(lats) <- rownames(tempo.d)
    colnames(lats) <- c('longitude', 'latitude')
    mst_shape <- msn2Shape(msn = mst2, lats = lats)
    write.shapefile(mst_shape, paste0(c('out_MST/mst_', j), collapse = ''))
    
    print(paste0(c(conta + 2, ') calculating mst... Done'), collapse = ''))
    
    # preparing shapes
    pontos_linha <- shapefile(paste0(c('out_MST/mst_', j, '.shp'), collapse = ''),
                              warnPRJ = FALSE)
    proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
    # plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
    # plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)
    
    idx <- which(names(qde) == j)
    tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- tempo[which(rownames(tempo) == j), 1]
    tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- tempo[which(rownames(tempo) == j), 2]
    tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
    
    # back-transforming lines in points:
    suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))
    
    lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presenças
    proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84   
    # plot(neo)
    # plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
    #      alpha = 0.5)
    writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"), 
                                         collapse = ''), overwrite = TRUE)
    
    # teste <- tabelao[j,j] # posições dos táxons do nó
    teste <- conta
    # text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])
    
    # presence-absence matrix:
    # preparing the matrix
    for(i in 1:ncellras){
      if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
        coor.l[i, conta] <- 1
      } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
        coor.l[i, conta] <- 0
      }
    }
  }
  
  coor.l <- na.exclude(coor.l)
  # coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
  # rownames(coor.ll) <- rownames(coor.l)
  write.table(x = coor.l, file = 'out_MST/pres_abs_MST.txt', sep = '\t')
  
  ##################
  print(paste0(conta + 3, ') Loading the raster files and adding the presences...'))
  layers_MST <- rast(list.files('out_MST/', pattern = "\\.tif$", full.names = T))
  MST_stack <- stack(layers_MST)
  # MST_stack[is.na(MST_stack)] <- 0
  
  # Somar todas as camadas para obter riqueza de espécies por célula
  MST_soma <- sum(rast(MST_stack), na.rm = T)
  return(list(geometry = MST_soma, pres_abs = coor.l))
}