## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------


## ------------------------------------

getwd()

project.name <- "Test"
folder.data <- ""

## ------------------------------------

land.poly <- "/Volumes/Jellyfish/Dropbox/Raw Data/Shapefiles/World Present/LandMass Polygon HR.shp"
missing.islands.poly <- "Data/Shapefiles/missing_islands.shp"
unwanted.release.sites.poly <- "Data/Shapefiles/land_clipper_2.shp"

## ------------------------------------

xmin <- -30
xmax <- 45
ymin <- 20
ymax <- 85
resolution <- 0.01

## ------------------------------------
## Hycom

months.all <- 3:4
from.day <- 1 ; to.day <- 31
from.year <- 2003 ; to.year <- 2003
depth.range <- c(0)
buffer <- TRUE 
buffer.val <- 0.05 
final.dimensions <- 2

# -----------------------------------

source("Dependences.R")

source("2. Get Data.R")

# Test dimensions of data (currents) as function
