## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "CentralityMPA"
project.folder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Network analyses reveal weak and strong links in the European network of Marine Protected Areas/"

coastline.shp <- "Data/Shapefiles/Global Coastline.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp"
bathymetry.tif <- NULL

additional.landmass.shp <- "Data/notake_merged_Med_P.shp" 

unwanted.release.coastline <- TRUE
unwanted.release.sites.shp <- NULL # "Data/Shapefiles/Unwanted.shp"

## ------------------------------------

source("Dependences.R")

## ------------------------------------

number.cores <- 8
parallel.computational.sections <- 8
parallel.computational.buffer <- 0.5 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -9
max.lon <- 38
min.lat <- 29
max.lat <- 47
source.sink.dist <- 1 # km

# -----------------------------------
# Traits

months.all <- 1:12 # c(9,10,11,12,1,2,3,4) # Spawning 5:10 (30 days off)
from.day <- 1 ; to.day <- 31
from.year <- 2008 ; to.year <- 2017

depth.range <- c(0)

kill.by.raft <- TRUE                              # Will eliminate particles that got to another cell - first raft event # May need a new particle every day
n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 125                       # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Hycom config

buffer <- TRUE 
buffer.val <- 0.05 
final.dimensions <- 2

# -----------------------------------
# Ilustration (movie)

movie.year <- 2017
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- "Data/Movie.shp" 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 

# --------------------------------------------------------------
# --------------------------------------------------------------

# source("1. Get Data.R")

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
