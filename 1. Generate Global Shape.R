## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

## Dependencies

setwd("/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0") # /Volumes/Laminaria Albacora Jellyfish

# ---------------------------------------------------------------------------------------

region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
crs(region.as.raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## --------------------
# ocean shape (polygon)

ocean <- shapefile(world.shape)
crs(ocean) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ocean <- as(ocean, "SpatialPolygons")

## Create clipping polygon

clipper <- as(extent(region.as.raster), "SpatialPolygons")
proj4string(clipper) <- CRS(proj4string(ocean))
ocean <- raster::intersect(ocean, clipper)

sections <- seq(extent(region.as.raster)[3],extent(region.as.raster)[4],length.out = 20)
cl.2 <- makeCluster(4)
registerDoParallel(cl.2)
ocean.surface.global <- foreach(s = 1:(length(sections)-1), .verbose=F, .packages=c("raster","gstat","ncdf4")) %dopar% {  
  section.i <- as(extent(xmin , xmax , sections[s] , sections[s+1]), "SpatialPolygons")
  proj4string(section.i) <- CRS(proj4string(ocean))
  ocean.r <- raster::intersect(ocean, section.i)
  section.i <- crop(region.as.raster, extent(xmin , xmax , sections[s] , sections[s+1]))
  ocean.r <- rasterize(ocean.r, section.i)
  ocean.r[!is.na(ocean.r)] <- -1
  return(ocean.r)
}
stopCluster(cl.2) ; rm(cl.2)

ocean.surface.global.t <- ocean.surface.global
ocean.surface.global.t$fun <- mean
ocean.surface.global.t$na.rm <- TRUE
ocean.surface.global.t <- do.call(mosaic, ocean.surface.global.t)

ocean.surface.global.t[is.na(ocean.surface.global.t)] <- 1
ocean.surface.global.t[ocean.surface.global.t == -1] <- NA
crs(ocean.surface.global.t) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# ----------------

plot(ocean.surface.global.t , col="black")
writeRaster(ocean.surface.global.t,filename="Data/ocean_raw",format="GTiff",overwrite=T)

# --------------------------------------------------------
# Remove isolated pixels

ocean.surface.global.t <- raster("Data/ocean_raw.tif")

rc <- raster::clump(ocean.surface.global.t)
plot(rc)

rcfreq <- raster::freq(rc)
rcfreq <- rcfreq[order(rcfreq[,2], decreasing = TRUE),]
head(rcfreq)

ocean.surface.global.t <- rc == 1
plot(ocean.surface.global.t)

writeRaster(ocean.surface.global.t,filename="Data/ocean_raw",format="GTiff",overwrite=T)

focal.matrix.1 <- matrix(c(0,1,0,1,1,1,0,1,0),nrow=3,ncol=3)
focal.matrix.2 <- matrix(1,nrow=3,ncol=3) # 1
blanks <- 1 ; previous.blanks <- 0

while( blanks != 0 & previous.blanks != blanks ) {
  
  previous.blanks <- blanks
  ocean.surface.global.t.f.1 <- focal(ocean.surface.global.t, w=focal.matrix.1, fun=sum, na.rm=TRUE , pad=FALSE)
  ocean.surface.global.t.f.2 <- focal(ocean.surface.global.t, w=focal.matrix.2, fun=sum, na.rm=TRUE , pad=FALSE)
  
  ocean.surface.global.t.f.values.1 <- which(getValues(ocean.surface.global.t.f.1) <= 2)
  ocean.surface.global.t.f.values.2 <- which(getValues(ocean.surface.global.t.f.2) <= 3)
  ocean.surface.global.t.f.values <- unique(ocean.surface.global.t.f.values.1, ocean.surface.global.t.f.values.2)
  blanks <- length(ocean.surface.global.t.f.values) 
  
  if( blanks > 0 ) { ocean.surface.global.t[ocean.surface.global.t.f.values] <- NA ; print(blanks) } 

}

rc <- raster::clump(ocean.surface.global.t)
rcfreq <- raster::freq(rc)
rcfreq <- rcfreq[order(rcfreq[,2], decreasing = TRUE),]
head(rcfreq)
ocean.surface.global.t <- rc == 1
writeRaster(ocean.surface.global.t,filename="Data/ocean_raw",format="GTiff",overwrite=T)

# --------------------------------------------------------

ocean.surface.global.t <- raster("Data/ocean_raw.tif")

clipper <- shapefile(super.clipper)
clipper <- rasterize(clipper, region.as.raster)
clipper[!is.na(clipper)] <- 1 ; clipper[is.na(clipper)] <- 2 ; clipper[clipper == 1] <- NA ; clipper[clipper == 2] <- 1
ocean.surface.global.t <- mask(ocean.surface.global.t,clipper)

if( exists("missing.islands") ) {
  
  missing.islands.r <- rasterize(as(shapefile(missing.islands),'SpatialPoints'), region.as.raster)
  missing.islands.r[!is.na(missing.islands.r)] <- 1
  
  ocean.surface.global.t[which(getValues(missing.islands.r) == 1)] <- NA
  
}

writeRaster(ocean.surface.global.t,filename="Data/ocean_raw",format="GTiff",overwrite=T)

# --------------------------------------------------------
# Last Clean

ocean.surface.global.t <- raster("Data/ocean_raw.tif")

focal.matrix.1 <- matrix(c(0,1,0,1,1,1,0,1,0),nrow=3,ncol=3)
focal.matrix.2 <- matrix(1,nrow=3,ncol=3) # 1
blanks <- 1 ; previous.blanks <- 0

while( blanks != 0 & previous.blanks != blanks ) {
  
  previous.blanks <- blanks
  ocean.surface.global.t.f.1 <- focal(ocean.surface.global.t, w=focal.matrix.1, fun=sum, na.rm=TRUE , pad=FALSE)
  ocean.surface.global.t.f.2 <- focal(ocean.surface.global.t, w=focal.matrix.2, fun=sum, na.rm=TRUE , pad=FALSE)
  
  ocean.surface.global.t.f.values.1 <- which(getValues(ocean.surface.global.t.f.1) <= 2)
  ocean.surface.global.t.f.values.2 <- which(getValues(ocean.surface.global.t.f.2) <= 3)
  ocean.surface.global.t.f.values <- unique(ocean.surface.global.t.f.values.1, ocean.surface.global.t.f.values.2)
  blanks <- length(ocean.surface.global.t.f.values) 
  
  if( blanks > 0 ) { ocean.surface.global.t[ocean.surface.global.t.f.values] <- NA ; print(blanks) } 
  
}

writeRaster(ocean.surface.global.t,filename="Data/ocean.tif",format="GTiff",overwrite=T)

# --------------------------------------------------------
# Coastline

#ocean.r.f <- Ocean.global.final
#region.as.raster <- Ocean.global.final
#region.as.raster[region.as.raster == 1] <- NA

ocean.r.f <- raster("Data/ocean.tif")

library(sdmpredictors)
lgm.shape <- load_layers("MS_bathy_21kya")
lgm.shape[!is.na(lgm.shape)] <- 1
ocean.r.f <- lgm.shape

# ---------------------------------
xmin <- extent(ocean.r.f)[1]
xmax <- extent(ocean.r.f)[2]
ymin <- extent(ocean.r.f)[3]
ymax <- extent(ocean.r.f)[4]
resolution <- res(ocean.r.f)
region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
region.as.raster <- raster(region.as.table)
extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
crs(region.as.raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# ---------------------------------

sections <- seq(extent(region.as.raster)[3],extent(region.as.raster)[4],length.out = 16)

cl.2 <- makeCluster(4)
registerDoParallel(cl.2)
Coastline.global <- foreach(s = 1:(length(sections)-1), .verbose=F, .packages=c("raster","gstat","ncdf4")) %dopar% {  
  
  section.i <- crop(region.as.raster, extent( xmin , xmax , sections[s] - 0.25 , sections[s+1] + 0.25))
  shape.ocean.r <- crop(ocean.r.f, section.i)    
  shape.r <- section.i
  
  ocean.cells.1 <- which(getValues(shape.ocean.r) == 1)
  ocean.cells.0 <- which(is.na(getValues(shape.ocean.r)))
  ocean.cells.adjacent.1 <- adjacent(shape.ocean.r, ocean.cells.1, directions=8)
  ocean.cells.adjacent.0 <- adjacent(shape.ocean.r, ocean.cells.0, directions=8)
  coast.line <- intersect(ocean.cells.adjacent.0,ocean.cells.adjacent.1)
  shape.r[coast.line] <- 1
  coast.line.r <- calc(stack(shape.ocean.r,shape.r), function(x) { ifelse( x[1] == 1 & x[2] == 1 , 1 , NA ) })
  coast.line.r <- mosaic(coast.line.r,region.as.raster,fun=mean,na.rm=TRUE)
  
  writeRaster(coast.line.r,file=paste0("Data/Temp/Temp.",s,".tif"),format="GTiff",overwrite=TRUE)

}
stopCluster(cl.2) ; rm(cl.2)

Coastline.global <- list.files(path="Data/Temp" , full.names=TRUE,pattern = "Temp" ) 

Coastline.global.final <- calc( stack(Coastline.global) , fun=mean , na.rm=TRUE)
crs(Coastline.global.final) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(Coastline.global.final,col="black")
n.cells <- as.data.frame(Coastline.global.final,xy=TRUE) ; n.cells <- n.cells[!is.na(n.cells[,3]),] ; nrow(n.cells)
file.remove(Coastline.global)

# Remove complex pixels 

clipper <- shapefile(super.clipper)
clipper <- rasterize(clipper, region.as.raster)
clipper[!is.na(clipper)] <- 1 ; clipper[is.na(clipper)] <- 2 ; clipper[clipper == 1] <- NA ; clipper[clipper == 2] <- 1
Coastline.global.final <- mask(Coastline.global.final,clipper)
n.cells <- as.data.frame(Coastline.global.final,xy=TRUE) ; n.cells <- n.cells[!is.na(n.cells[,3]),] ; nrow(n.cells)

writeRaster(Coastline.global.final,filename="Data/coast_line",format="GTiff",overwrite=T)
file.remove("data/ocean_raw.tif")

# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
## 0.05 BIO-ORACLE

# Version with all intertidal

bathy <- raster('/Volumes/Laminaria/Dropbox/Raw Data/Rasters/Bathymetry/Bathymetry _ depth Mean.tif')
rc <- bathy
rc[!is.na(rc)] <- 1
plot(rc)
Ocean.global.final <- rc

region.as.raster <- rc
ocean.r.f <- rc
ocean.r.f[is.na(ocean.r.f)] <- 1
ocean.r.f[ocean.r.f == 0] <- NA
sections <- seq(extent(region.as.raster)[3],extent(region.as.raster)[4],length.out = 32)

# ---------------------------------

# Version with open waters

bathy <- raster('/Volumes/Laminaria/Dropbox/Raw Data/Rasters/Bathymetry/Bathymetry _ depth Mean.tif')
xmin <- -180
xmax <- 180

## check for small closed seas

rc <- raster::clump(bathy)
rcfreq <- raster::freq(rc)
rcfreq <- rcfreq[order(rcfreq[,2], decreasing = TRUE),]

plot(bathy)
plot(rc)
plot(rc == 0)
rc <- rc == 0

region.as.raster <- rc
ocean.r.f <- rc
ocean.r.f[is.na(ocean.r.f)] <- 1
ocean.r.f[ocean.r.f == 0] <- NA
sections <- seq(extent(region.as.raster)[3],extent(region.as.raster)[4],length.out = 32)

Ocean.global.final <- ocean.r.f
Ocean.global.final[is.na(Ocean.global.final)] <- 2
Ocean.global.final[Ocean.global.final == 1] <- NA
Ocean.global.final[!is.na(Ocean.global.final)] <- 1
writeRaster(Ocean.global.final,filename="Data/ocean_oracle",format="GTiff",overwrite=T)

# ---------------------------------
# ---------------------------------

cl.2 <- makeCluster(4)
registerDoParallel(cl.2)
Coastline.global <- foreach(s = 1:(length(sections)-1), .verbose=F, .packages=c("raster","gstat","ncdf4")) %dopar% {  
  
  section.i <- crop(region.as.raster, extent( xmin , xmax , sections[s] - 0.25 , sections[s+1] + 0.25))
  shape.ocean.r <- crop(ocean.r.f, section.i)    
  shape.r <- section.i
  
  ocean.cells.1 <- which(getValues(shape.ocean.r) == 1)
  ocean.cells.0 <- which(is.na(getValues(shape.ocean.r)))
  ocean.cells.adjacent.1 <- adjacent(shape.ocean.r, ocean.cells.1, directions=8)
  ocean.cells.adjacent.0 <- adjacent(shape.ocean.r, ocean.cells.0, directions=8)
  coast.line <- intersect(ocean.cells.adjacent.0,ocean.cells.adjacent.1)
  shape.r[coast.line] <- 1
  coast.line.r <- calc(stack(shape.ocean.r,shape.r), function(x) { ifelse( x[1] == 1 & x[2] == 1 , 1 , NA ) })
  coast.line.r <- mosaic(coast.line.r,region.as.raster,fun=mean,na.rm=TRUE)
  
  writeRaster(coast.line.r,file=paste0("Data/Temp/T.",s,".tif"),format="GTiff",overwrite=TRUE)
  
}
stopCluster(cl.2) ; rm(cl.2)

Coastline.global <- list.files(path="Data/Temp" , full.names=TRUE,pattern = "T" ) 
Coastline.global.final <- stack(Coastline.global)
Coastline.global.final[Coastline.global.final == 0] <- NA

Coastline.global.final <- calc( Coastline.global.final , fun=mean , na.rm=TRUE)
crs(Coastline.global.final) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(Coastline.global.final,col="black")
file.remove(Coastline.global)

writeRaster(Coastline.global.final,filename="Data/coast_line_oracle",format="GTiff",overwrite=T)

clipper <- shapefile(super.clipper)
clipper <- rasterize(clipper, region.as.raster)
clipper[!is.na(clipper)] <- 1 ; clipper[is.na(clipper)] <- 2 ; clipper[clipper == 1] <- NA ; clipper[clipper == 2] <- 1
Coastline.global.final <- mask(Coastline.global.final,clipper)
n.cells <- as.data.frame(Coastline.global.final,xy=TRUE) ; n.cells <- n.cells[!is.na(n.cells[,3]),] ; nrow(n.cells)

writeRaster(Coastline.global.final,filename="Data/coast_line_oracle",format="GTiff",overwrite=T)

