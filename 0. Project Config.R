## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

# setwd("/media/Bathyscaphe/Transport Simulation in the Atlantic Halodule/Git")

## ------------------------------------
## Files and folders

project.name <- "Halodule"
project.folder <- "/media/Bathyscaphe/Transport Simulation in the Atlantic Halodule/"
data.folder <- paste0(project.folder,"Data/")

landmass.shp <- NULL # "../Data/mainLandAzores.shp" 
bathymetry.tif <- NULL

additional.landmass.shp <- NULL # "../Data/Dispersal simulations/Shapefiles/additionalSites.shp" 
additional.source.sink.shp <- "../Data/sourceSinkSites.shp" # "../Data/sourceSinkPolygons_0.shp" # "../Data/Shapefiles/rockyHabitats" # NULL
additional.source.sink.shp.force.shore <- TRUE # If additional.source.sink.shp are new regions

source.sink.loc.type <- "centroid" # centroid peripheral

unwanted.release.coastline <- TRUE
unwanted.release.sites.shp <- NULL # "Data/Dispersal simulations/Shapefiles/unwantedSites.shp" # NULL

# -----------------------------------

number.cores <- 20 # [!!]

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -105 
max.lon <- 25
min.lat <- -40
max.lat <- 41

buffer <- TRUE
buffer.val <- 0.1

sim.resolution <- 5 # https://github.com/uber/h3/blob/master/docs/core-library/restable.md

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