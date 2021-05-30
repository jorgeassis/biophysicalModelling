## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2021
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "01 Intertidal"
project.folder <- "/Users/jorgeassis/Google Drive/theMarineDataScientist/Projects/Portuguese EEZ Connectivity estimates/"
data.folder <- paste0(project.folder,"Data/")

# -----------------------------------

number.cores <- 20 # [!!]

# -----------------------------------
# Region

# Change 0
landmass.shp <- NULL # NULL for default landmass // .shp file for specific landmass
alternativeLandmassPolygon <- NULL

# Change 1
additionalLandmassPolygon <- 
  additional.landmass.shp <- NULL # NULL for none // .shp file for addiitonal landmass regions

# Change 1.2
additional.source.sink.shp <- "../Data/sourceSinkSites.shp" # "../Data/sourceSinkPolygons_0.shp" # "../Data/Shapefiles/rockyHabitats" # NULL
additional.source.sink.shp.force.shore <- TRUE # If additional.source.sink.shp are new regions
additionalSourceSinkSites <- #POLY? XY? ??
  
  # Change 2
  source.sink.loc.type <- "centroid" # centroid peripheral along H3 polygons
locationSourceSinkSites <- "centroid" # centroid peripheral along H3 polygons

# Change 3
unwanted.release.coastline <- TRUE # If TRUE, landmass regions are unwanted sourceSink sites
removeLandmassSourceSinkSites <- TRUE # If TRUE, landmass regions are unwanted sourceSink sites

# Change 4
unwanted.release.sites.shp <- NULL 
maskSourceSinkPolygon <- NULL # .shp file to mask sourceSink sites
maskSourceSinkType <- "include" # include exclude

maskSourceSinkBathymetry <- NULL # range of bathymetry
bathymetryRasterFile <- NULL # .tif bathymetry file


dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -105 
max.lon <- 25
min.lat <- -40
max.lat <- 41

buffer <- TRUE
buffer.val <- 0.1

sim.resolution <- 6 # https://h3geo.org/docs/core-library/restable/

# -----------------------------------
# Traits

months.all <- 1:12 # c(9,10,11,12,1,2,3,4) 
from.day <- 1 ; to.day <- 31
from.year <- 2008 ; to.year <- 2017 #  2008:2017

depth.range <- c(0)

allow.back.to.origin <- FALSE                     # at t == t.start

n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 60                       # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Hycom config

final.dimensions <- 2

# -----------------------------------
# Ilustration (movie)

movie.year <- 2017
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- NULL # "../Data/Movie.shp" 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------