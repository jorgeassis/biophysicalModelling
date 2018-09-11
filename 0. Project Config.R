## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

# Updated
## ------------------------------------
## Files and folders

project.name <- "NAtlantic"
project.folder <- "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/Plankton"

coastline.shp <- "Data/Shapefiles/Global Coastline.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp"
bathymetry.tif <- NULL
unwanted.release.sites.shp <- "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/Plankton/Data/Shapefiles/unwanted.shp"

## ------------------------------------

number.cores <- 12      
parallel.computational.sections <- 12
parallel.computational.buffer <- 2 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -25
max.lon <- 32
min.lat <- 22
max.lat <- 81.5
source.sink.dist <- 10 # km

# -----------------------------------
# Traits

months.all <- 5:11 # Spawning 5:10 (30 days off)
from.day <- 1 ; to.day <- 31
from.year <- 2003 ; to.year <- 2012
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

movie.year <- 2010
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/Plankton/Data/Shapefiles/movie.shp" 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 


# --------------------------------------------------------------
# --------------------------------------------------------------

source("Dependences.R")

# --------------------------------------------------------------

# source("1. Get Data.R")

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------