## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

getwd()

## ------------------------------------

xmin <- -30
xmax <- 45
ymin <- 20
ymax <- 85
resolution <- 0.01

world.shape <- "Data/Shapefiles/gshhs/world_hd.shp" # gshhs/world_hd.shp naturalearthdata/world_md.shp
missing.islands <- "Data/Shapefiles/missing_islands.shp" # Add new regions
super.clipper <- "Data/Shapefiles/land_clipper_2.shp"

# -----------------------------------
