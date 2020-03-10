## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "WAMPA"
project.folder <- "/media/Jellyfish/Dropbox/Manuscripts/West Africa MPA connectivity/"

coastline.shp <- "Data/Shapefiles/Global Coastline.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp"
bathymetry.tif <- NULL

additional.landmass.shp <- "Data/Spatial/surf_marine_2_new130520192.shp" # NULL

unwanted.release.coastline <- TRUE
unwanted.release.sites.shp <- NULL # "Data/Shapefiles/Unwanted.shp"

## ------------------------------------

source("Dependences.R")

## ------------------------------------

number.cores <- 40
parallel.computational.sections <- 10
parallel.computational.buffer <- 2 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -30
max.lon <- -7.5
min.lat <- 4.75
max.lat <- 24.5
source.sink.dist <- 1 # km

# -----------------------------------
# Traits

months.all <- 1:12 # Spawning 5:10 (30 days off)
from.day <- 1 ; to.day <- 31
from.year <- 2003 ; to.year <- 2012
depth.range <- c(0)

kill.by.raft <- TRUE                              # Will eliminate particles that got to another cell - first raft event # May need a new particle every day
n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 200                       # Days allowed to travel
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

movie.sites.xy <- "Data/Spatial/movie.shp" 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 

# --------------------------------------------------------------
# --------------------------------------------------------------

# source("1. Get Data.R")

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
