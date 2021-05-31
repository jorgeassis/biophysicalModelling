## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2021
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("0. Config.R") 
source("Dependences.R")

## -----------------------------------------
## Produce additional Source Sink sites from Bathymetry (additionalSourceSinkRegions)

bathymetryRange <- c(0,-40)
bathymetryRasterFile <- "../Data/Spatial/Emodnet Bathymetry Atlantic.tif"

## -------------

if( is.null(alternativeLandmass)) { worldmap <- getMap(resolution = "high") }
if( ! is.null(alternativeLandmass)) { worldmap <- shapefile(alternativeLandmass) }

worldmap <- gBuffer(worldmap, byid=TRUE, width=0)
worldmap <- crop(worldmap,extent(min.lon,max.lon,min.lat,max.lat))
worldmap <- st_as_sf(worldmap)
worldmap$ID <- 1:nrow(worldmap)

## -------------

bathymetryRaster <- raster(bathymetryRasterFile)
bathymetryRaster <- crop(bathymetryRaster,extent(min.lon,max.lon,min.lat,max.lat))
bathymetryRasterCells <- Which(bathymetryRaster <= bathymetryRange[1] & bathymetryRaster >= bathymetryRange[2], cells=TRUE)
additionalSourceSinkRegions <- xyFromCell(bathymetryRaster,bathymetryRasterCells)
additionalSourceSinkRegions <- as.data.frame(additionalSourceSinkRegions)
coordinates(additionalSourceSinkRegions)=~x+y
proj4string(additionalSourceSinkRegions)<- CRS("+proj=longlat +datum=WGS84") 
additionalSourceSinkRegions <- st_as_sf(additionalSourceSinkRegions)

additionalSourceSinkRegions.land <- st_intersects( worldmap[,"ID"],additionalSourceSinkRegions)
additionalSourceSinkRegions <- additionalSourceSinkRegions[-unique(unlist(additionalSourceSinkRegions.land)),]
additionalSourceSinkRegions <- sf:::as_Spatial(additionalSourceSinkRegions)
additionalSourceSinkRegions$ID <- 1:nrow(additionalSourceSinkRegions)

plot(additionalSourceSinkRegions)
raster::shapefile(additionalSourceSinkRegions, "../Data/Spatial/additionalSourceSink_0_40.shp",overwrite=TRUE)

## -----------------------------------------------------
## -----------------------------------------------------
