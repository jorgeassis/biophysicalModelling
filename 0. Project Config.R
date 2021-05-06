## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

# setwd("/media/Bathyscaphe/Transport Simulation in the Atlantic Halodule/Git")

## ------------------------------------
## Files and folders

project.name <- "Selvagens"
project.folder <- "/media/Bathyscaphe/Transport Simulations in Selvagens/"
data.folder <- paste0(project.folder,"Data/")

landmass.shp <- "Data/Shapefiles/Global Landmass.shp" # NULL "../Data/mainLandAzores.shp" "Data/Shapefiles/Global Landmass.shp"
bathymetry.tif <- NULL

# Revise to be polygons or points from source sink 

additional.landmass.shp <- NULL # "../Data/Dispersal simulations/Shapefiles/additionalSites.shp" 
additional.source.sink.shp <- NULL # "../Data/sourceSinkSites.shp" # "../Data/sourceSinkPolygons_0.shp" # "../Data/Shapefiles/rockyHabitats" # 
additional.source.sink.shp.force.shore <- NULL # TRUE # If additional.source.sink.shp are new regions

source.sink.loc.type <- "centroid" # centroid peripheral

unwanted.release.coastline <- FALSE
unwanted.release.sites.shp <- NULL # "Data/Dispersal simulations/Shapefiles/unwantedSites.shp" # NULL

# -----------------------------------

number.cores <- 40

# -----------------------------------
# Region

dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
min.lon <- -35.5
max.lon <- 12
min.lat <- 10
max.lat <- 51.5

buffer <- TRUE
buffer.val <- 0.1

sim.resolution <- 5 # https://h3geo.org/docs/

# -----------------------------------
# Traits

months.all <- 1:12 # c(9,10,11,12,1,2,3,4) 
from.day <- 1 ; to.day <- 31
from.year <- 2013 ; to.year <- 2017

depth.range <- c(0)

allow.back.to.origin <- FALSE                     # at t == t.start

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

movie.year <- 2017
movie.sites.buffer <- 0 # Nearby cells to include, 0 for xy only
movie.sites.xy <-  matrix( c(  
  -15.8646,30.156, -15.706,28.1755,-8.746,41.049,-16.2856,33.088,-31.276,39.354, -25.0687,36.89,-5.825,43.7255,-1.6839,46.341,-4.76741201706258,48.0318,-9.21962136117718,43.1807847122748,-9.60139625249854,38.7215600343441,-9.26082016150725,32.6135115371994,-16.6436324483549,32.465324224514,-13.7432419070933,29.1061425199271,-17.9619917852919,28.6443637594147,-26.9707805876953,38.678689734983,-9.04109360535107,36.8726319197159,-3.67976563514028,43.5740615601251,-3.81160156883399,35.3460902665874,-1.83406256342836,36.7670941667282,0.319257686902193,38.678689734983,8.14152308606219,40.9069099535121,6.42765594804399,37.1533536252814,-6.09675775285826,35.8107285375731,-24.8174603373647,37.7116791418826 ) , ncol=2 , byrow=TRUE) 
# "../Data/Movie.shp" 

## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------