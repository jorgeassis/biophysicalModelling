## ------------------------------------------------------------------ ##
##                                                                    ##
##  PlankTonic v3.0                                                   ##
##  R Pipelines for biophysical modelling                             ##
##                                                                    ##
## --------------------------------------------                       ##
##                                                                    ##
##  Jorge Assis [ jorgemfa@gmail.com ]                                ##
##  biodiversityDS                                                    ##
##  biodiversityDataScience.com                                       ##
##                                                                    ##
## ------------------------------------------------------------------ ##

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "testLandMass" 
project.folder <- "../"
data.folder <- "/media/Jellyfish/Raw Data [Environment]/Copernicus/" # "../Data/" # 
results.folder <- paste0(project.folder,"Results/")
temp.folder <- paste0(project.folder,"tempFolder/")

# -----------------------------------

number.cores <- 20 # parallel::detectCores() / 2

# -----------------------------------
# Main region

sim.resolution <- 5 # >= 5 for HR # https://h3geo.org/docs/core-library/restable/
min.lon <- -180
max.lon <- 180
min.lat <- -90
max.lat <- 90
dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# -----------------------------------
# Source Sink sites

sourceSinkLocationType <- "centroid" # centroid or peripheral along H3 polygons
alternativeLandmass <- NULL # "Data/Shapefiles/Global Landmass.shp" # NULL or .shp file for specific landmass
removeLandmassSourceSinkSites <- FALSE ## landmass regions are unwanted Source Sink sites

maskSourceSinkSites <- NULL # "../Data/Spatial/innerSeas.shp" ## NULL or .shp file to mask Source Sink sites
maskSourceSinkSitesType <- "exclude" ## include or exclude Source Sink sites from sim

addSourceSinkRegions <- NULL # "../connectivityNEPacific/Data/sourceSinkSites.shp" ## NULL or .shp file for additional Points or Polygon regions
addSourceSinkRegionsForceCoast <- NULL ## added regions are strictly coastal distributed

addSourceSinkRegionsBathymetry <- NULL ## NULL or range
bathymetryRasterFile <- NULL ## NULL or .tif raster layer

# -----------------------------------
# Hycom config

rawDataDimensions <- 2
rawDataDepth <- c(0)

# -----------------------------------
# Traits

months.all <- 1:12
from.day <- 1
to.day <- 31
from.year <- 2000
to.year <- 2020

allow.back.to.origin <- FALSE                     # at t == t.start

n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE           # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 180                      # Days allowed to travel

# -----------------------------------
# Illustration (movie)

movie.year <- 2009
movie.sites.buffer <- 2 # Nearby cells to include, 0 for xy only

movie.sites.xy <- NULL # "../Data/Spatial/movie.shp" # NULL  "../Data/Spatial/movie.shp" # 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 ) , ncol=2 , byrow=TRUE) 

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------