## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

# Updated
## ------------------------------------
## Files and folders

project.name <- "Test"
project.folder <- "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/TestScript"

coastline.shp <- "Data/Shapefiles/Global Coastline.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp"
bathymetry.tif <- NULL
unwanted.release.sites.shp <- NULL

## ------------------------------------

number.cores <- 6      
parallel.computational.sections <- 6
parallel.computational.buffer <- 2 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -20
max.lon <- 10
min.lat <- 20
max.lat <- 50
source.sink.dist <- 50 # km

# -----------------------------------
# Traits

months.all <- 3:6
from.day <- 1 ; to.day <- 31
from.year <- 2003 ; to.year <- 2003
depth.range <- c(0)

kill.by.raft <- TRUE                              # Will eliminate particles that got to another cell - first raft event # May need a new particle every day
n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- TRUE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 90                       # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Hycom config

buffer <- TRUE 
buffer.val <- 0.05 
final.dimensions <- 2

# -----------------------------------
# Ilustration (movie)

movie.year <- 2003
movie.sites.xy <- matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 ) , ncol=2 , byrow=TRUE) 
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

# --------------------------------------------------------------
# --------------------------------------------------------------

# source("1. Get Data.R")

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------