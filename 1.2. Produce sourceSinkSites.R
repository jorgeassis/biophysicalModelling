## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2021
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("0. Config.R") 
source("Dependences.R")

## -----------------------------------------
## Produce Source Sink sites from Bathymetry (additionalSourceSinkRegions)

bathymetryRange <- c(0,40)
bathymetryRasterFile <- "../Data/Spatial/Emodnet Bathymetry Atlantic.tif"

additionalSourceSinkRegions

## ------------------------------------------------------------------------------------------------------------------------------
##
## 
## 
## ------------------------------------------------------------------------------------------------------------------------------

## Define region of interest

if( is.null(alternativeLandmass)) { worldmap <- getMap(resolution = "high") }
if( ! is.null(alternativeLandmass)) { worldmap <- shapefile(alternativeLandmass) }

worldmap <- gBuffer(worldmap, byid=TRUE, width=0)
worldmap <- crop(worldmap,extent(min.lon,max.lon,min.lat,max.lat))

worldmap <- st_as_sf(worldmap)
worldmap$ID <- 1:nrow(worldmap)

## --------------------------------------

## Produce hexagon lists

shapeVertex <- data.frame(Lon=c(min.lon,min.lon,max.lon,max.lon),
                          Lat=c(min.lat,max.lat,max.lat,min.lat))
shape <- spPolygons(as.matrix(shapeVertex))
crs(shape) <- dt.projection
shape <- st_as_sf(shape)
hexagons.address <- h3::polyfill(shape, sim.resolution)

# -------

polygons.all <- h3_to_geo_boundary_sf(hexagons.address)
hexagons.land <- st_intersects( worldmap,polygons.all) # works more or less...

polygons.all <- polygons.all[unique(unlist(hexagons.land)),]
hexagons.address.land <- hexagons.address[unique(unlist(hexagons.land))]
polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)
plot(polygons.land,col="red")

# polygons.all <- h3_to_geo_boundary_sf(hexagons.address)
# centroid.all <- st_centroid(polygons.all)
# centroid.all$ID <- 1:nrow(centroid.all)
# 
# hexagons.land <- st_join(centroid.all, worldmap, join = st_intersects)
# hexagons.land <- as.data.frame(hexagons.land)
# hexagons.address.land <- hexagons.address[which(!is.na(hexagons.land$ID.y))]
# 
# polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)
# plot(polygons.land)

hexagons.address.ocean <- unique( hexagons.address[ ! hexagons.address %in% hexagons.address.land ] )
plot(h3_to_geo_boundary_sf(hexagons.address.ocean),col="blue")

## --------------------------------------
## Get shore hexagons

shoreTest <- function(hexagon) {
  neighbors <- k_ring(hexagon, radius = 1)
  if( sum(hexagons.address.land %in% neighbors) > 0 & sum(hexagons.address.ocean %in% neighbors) > 0 ) { return(hexagon) } else { return(NULL) }
}

cl <- makeCluster(number.cores)
clusterExport(cl, c("hexagons.address.land","shoreTest","hexagons.address.ocean","k_ring"))
hexagons.address.shore <- unlist(unique(parLapply(cl, hexagons.address.land , function(x) { shoreTest(x) })))
stopCluster(cl)

hexagons.address.land <- hexagons.address.land[! hexagons.address.land %in% hexagons.address.shore]
hexagons.address.ocean <- hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.shore]

## Clean shore hexagons 

shoreTestCleaner <- function(hexagon) {
  neighbors <- k_ring(hexagon, radius = 1)
  if( sum(hexagons.address.land %in% neighbors) > 0 & sum(hexagons.address.ocean %in% neighbors) > 0 ) { return(NULL) } else { return(hexagon) }
}

cl <- makeCluster(number.cores)
clusterExport(cl, c("hexagons.address.shore","hexagons.address.land","shoreTestCleaner","hexagons.address.ocean","k_ring"))
hexagons.address.ocean.add <- unlist(unique(parLapply(cl, hexagons.address.shore , function(x) { shoreTestCleaner(x) })))
stopCluster(cl)

hexagons.address.ocean <- unique(c(hexagons.address.ocean,hexagons.address.ocean.add))

shoreTestCleaner <- function(hexagon) {
  neighbors <- k_ring(hexagon, radius = 1)
  if( sum(hexagons.address.land %in% neighbors) > 0 & sum(hexagons.address.ocean %in% neighbors) > 0 ) { return(hexagon) } else { return(NULL) }
}

cl <- makeCluster(number.cores)
clusterExport(cl, c("hexagons.address.shore","hexagons.address.land","shoreTestCleaner","hexagons.address.ocean","k_ring"))
hexagons.address.shore.add <- unlist(unique(parLapply(cl, hexagons.address.shore , function(x) { shoreTestCleaner(x) })))
stopCluster(cl)

hexagons.address.shore <- unique(c(hexagons.address.shore,hexagons.address.shore.add))
hexagons.address.land <- hexagons.address.land[! hexagons.address.land %in% hexagons.address.shore]
hexagons.address.ocean <- hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.shore]

polygons.shore <- h3_to_geo_boundary_sf(hexagons.address.shore)
polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)
polygons.ocean <- h3_to_geo_boundary_sf(hexagons.address.ocean)

## --------------------------------------
## Remove outer frame

centroid.land <- st_centroid(polygons.land)
buffer.remover <- which(st_coordinates(centroid.land)[,1] <= max.lon - buffer.val &
                          st_coordinates(centroid.land)[,1] >= min.lon + buffer.val &
                          st_coordinates(centroid.land)[,2] <= max.lat - buffer.val &
                          st_coordinates(centroid.land)[,2] >= min.lat + buffer.val )

hexagons.address.land <- hexagons.address.land[buffer.remover]
polygons.land <- polygons.land[buffer.remover,]

centroid.shore <- st_centroid(polygons.shore)
buffer.remover <- which(st_coordinates(centroid.shore)[,1] <= max.lon - buffer.val &
                          st_coordinates(centroid.shore)[,1] >= min.lon + buffer.val &
                          st_coordinates(centroid.shore)[,2] <= max.lat - buffer.val &
                          st_coordinates(centroid.shore)[,2] >= min.lat + buffer.val )

hexagons.address.shore <- hexagons.address.shore[buffer.remover]
polygons.shore <- polygons.shore[buffer.remover,]

if( sum(hexagons.address.shore %in% hexagons.address.ocean) > 0 | sum(hexagons.address.land %in% hexagons.address.ocean) > 0 | sum(hexagons.address.ocean %in% hexagons.address.land) > 0 ) { stop("Error :: 001") }

## ----------

ggplot() + geom_sf(data = polygons.land, fill=NA, colour="red", size=0.1)
ggplot() + geom_sf(data = polygons.shore, fill="Black", colour="Black", size=0.1)
ggplot() + geom_sf(data = polygons.ocean, fill=NA, colour="red", size=0.1)
