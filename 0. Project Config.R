## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "Perna"
project.folder <- "/media/Bathyscaphe/Transport Simulation in South Africa Perna perna/"

coastline.shp <- "Data/Shapefiles/Global Coastline.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp"

additional.sourcesink.shp <- NULL
landmass.shp.2 <- NULL

unwanted.release.coastline <- FALSE
unwanted.release.sites.shp <- "Data/Shapefiles/Unwanted.shp" # NULL

bathymetry.tif <- NULL

## ------------------------------------

source("Dependences.R")

## ------------------------------------

number.cores <- 10
parallel.computational.sections <- 10
parallel.computational.buffer <- 2 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- 16
max.lon <- 34
min.lat <- -38
max.lat <- -28 # 
source.sink.dist <- 1 # km

# -----------------------------------
# Traits

months.all <- 1:12 # Spawning 5:10 (30 days off)
from.day <- 1 ; to.day <- 31
from.year <- 2008 ; to.year <- 2012
depth.range <- c(0)

kill.by.raft <- TRUE                              # Will eliminate particles that got to another cell - first raft event # May need a new particle every day
n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 30                        # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Hycom config

buffer <- TRUE 
buffer.val <- 0.05 
final.dimensions <- 2

# -----------------------------------
# Ilustration (movie)

movie.year <- -999
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- "Data/Shapefiles/Movie.shp" # NULL
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 

# --------------------------------------------------------------
# --------------------------------------------------------------

# source("1. Get Data.R")

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
