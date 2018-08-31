
# Points over Polygon

library(spatialEco)
new_shape <- point.in.poly(pnts, ind_adm)

library(sp)
pt.in.poly <- sp::over(ind_adm, pnts, fn = NULL)

You could also use the st_intersection function from the sf package:
  
  Load the library

library(sf)

Create a simple feature geometry (polygon) from your polygon

ind_adm <- st_as_sf(ind_adm)

Create a simple feature geometry (point) from your points of interest

(24047 is the EPSG code for India)

pnts <- st_as_sf(pnts) %>% st_set_crs(., 24047)

Keep only the points inside the polygon

kept_points <- st_intersection(ind_adm, pnts)

