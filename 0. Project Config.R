## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2020
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

project.name <- "NAtlantic"
project.folder <- "/media/Bathyscaphe/Transport Simulations Explain Genetic Differention of North Atlantic Marine Forests/"

currentsVel.data.folder <- paste0(project.folder,"Data/")

# [decap] coastline.shp <- "Data/Shapefiles/Global Coastline.shp" # "../Data/mainLandAzores.shp"
landmass.shp <- "Data/Shapefiles/Global Landmass.shp" # "../Data/shoreLineAzores.shp"

bathymetry.tif <- NULL

additional.landmass.shp <- "../Data/Dispersal simulations/Shapefiles/additionalSites.shp" 
additional.source.sink.shp <- NULL 

source.sink.generation.type <- "centroid" # peripherical [decap] centroid

unwanted.release.coastline <- FALSE
unwanted.release.sites.shp <- "Data/Dispersal simulations/Shapefiles/unwantedSites.shp" # NULL

## ------------------------------------

number.cores <- 1
# [decap] parallel.computational.sections <- 1
# [decap] parallel.computational.buffer <- 0.5 # degrees

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -10 # -43.5
max.lon <- 10 # 38
min.lat <- 35 # 21
max.lat <- 50 # 78

buffer <- TRUE
buffer.val <- 0.1

sim.resolution <- 6 # https://github.com/uber/h3/blob/master/docs/core-library/restable.md

# -----------------------------------
# Traits

months.all <- 1:12 # 1:12 c(9,10,11,12,1,2,3,4) # Spawning 5:10 (30 days off)
from.day <- 1 ; to.day <- 31
from.year <- 2013 ; to.year <- 2013

depth.range <- c(0)

allow.retention <- TRUE
allow.back.to.origin <- FALSE

n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE            # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 30                       # Days allowed to travel
behaviour <- FALSE                                # Only settle after period

# -----------------------------------
# Hycom config

final.dimensions <- 2

# -----------------------------------
# Ilustration (movie)

movie.year <- 2013
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only

movie.sites.xy <- "Data/Dispersal simulations/Shapefiles/movie.shp" 
# matrix( c(  -8.892305, 37.956704 , -9.225347 , 38.411873 , -9.489235 , 38.708553 , -5 , 50 , - 2 , 45 , -5 , 44 , 5 , 40 , -10 , 30 ) , ncol=2 , byrow=TRUE) 

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------