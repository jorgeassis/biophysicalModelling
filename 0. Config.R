## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2021
## ------------------------------------------------------------------------------------------------------------------

## ------------------------------------
## Git

# library(credentials)
# ssh_keygen()
# set_github_pat()

## ------------------------------------
## ------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "Mangrove" 
project.folder <- "../"
data.folder <- paste0(project.folder,"Data/")
results.folder <- paste0(project.folder,"Results/")

# -----------------------------------

number.cores <- 32

# -----------------------------------
# Main region

min.lon <- -180 
max.lon <- 180
min.lat <- -50
max.lat <- 50
dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sim.resolution <- 4 # https://h3geo.org/docs/core-library/restable/

# -----------------------------------
# Source Sink sites

sourceSinkLocationType <- "centroid" # centroid or peripheral along H3 polygons

alternativeLandmass <- NULL ## .shp file for specific landmass
removeLandmassSourceSinkSites <- FALSE ## landmass regions are unwanted Source Sink sites

addSourceSinkRegions <- NULL ## .shp file for additional Points or Polygon regions

addSourceSinkRegionsBathymetry <- NULL ## NULL or range
bathymetryRasterFile <- NULL ## NULL or .tif raster layer

maskSourceSinkSites <- "../Data/Spatial/mangroveOccurrence.shp" ## .shp file to mask Source Sink sites
maskSourceSinkSitesType <- "include" ## include or exclude Source Sink sites from sim

# -----------------------------------
# Hycom config

rawDataDimensions <- 2
rawDataDepth <- c(0)

# -----------------------------------
# Traits

months.all <- 1:12 # c(9,10,11,12,1,2,3,4) 
from.day <- 1 ; to.day <- 31
from.year <- 2008 ; to.year <- 2016 #  2008:2017

allow.back.to.origin <- FALSE                     # at t == t.start

n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 180                       # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Ilustration (movie)

movie.year <- 2016
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- NULL # "../Data/Spatial/movie.shp" # 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 ) , ncol=2 , byrow=TRUE) 

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------