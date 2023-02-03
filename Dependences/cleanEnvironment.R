## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

options(warn=-1)

rm("polygons.land.sp","particles.reference.bm","polygons.plot.ocean","polygons.all","hexagons.land","polygons.land","polygons.plot.land", "centroid.land","hexagons.address.1","hexagons.address.2","hexagons.sourcesink.shp","centroid.shore","polygons.sourceSink","polygons.plot.sourcesink")

closeAllConnections(); gc(reset=TRUE)
print(head(list.memory()))

options(warn=0)
