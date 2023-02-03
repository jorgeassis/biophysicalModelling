## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##

options(warn=-1)

# Open and crop Landmass

if( is.null(alternativeLandmass)) { worldmap <- getMap(resolution = "high") }
if( ! is.null(alternativeLandmass)) { worldmap <- readOGR(dsn = substr(alternativeLandmass,1,unlist(gregexpr("/",alternativeLandmass))[2]), layer = substr(alternativeLandmass,unlist(gregexpr("/",alternativeLandmass))[2]+1,nchar(alternativeLandmass)-4)) }

cat("Landmass from external spatialPolygons:",as.character(!is.null(alternativeLandmass)),"\n")

worldmap <- gBuffer(worldmap, byid=TRUE, width=0)
worldmap <- crop(worldmap,extent(min.lon,max.lon,min.lat,max.lat))
worldmap <- st_as_sf(worldmap[,1])
st_crs(worldmap) <- 4326
worldmap <- st_union(worldmap)
worldmap <- st_wrap_dateline(worldmap, options = "WRAPDATELINE=YES", quiet = TRUE)

## --------------------------------------
# Get hexagons addresses

shapeVertex <- data.frame(Lon=c(min.lon,min.lon,max.lon,max.lon),Lat=c(min.lat,max.lat,max.lat,min.lat))
shapeVertex <- st_as_sf(shapeVertex, coords = c("Lon", "Lat"), crs = 4326)
shapeVertex <- sf::st_make_grid(shapeVertex, n=10)

hexagons.address <- h3jsr::polyfill(shapeVertex, sim.resolution) 
hexagons.address <- unique(unlist(hexagons.address))
names(hexagons.address) <- NULL

hexagons.address.poly <- h3_to_geo_boundary_sf(hexagons.address)
st_crs(hexagons.address.poly) <- 4326

intersectionWhich <- st_intersects(hexagons.address.poly, worldmap)
intersectionWhich <- which(lengths(intersectionWhich) > 0)
intersectionWhich.area <- st_area(hexagons.address.poly[intersectionWhich,])

intersection <- st_intersection(hexagons.address.poly, worldmap)
intersection <- st_area(intersection)

hexagons.address.land <- hexagons.address[intersectionWhich[which((as.numeric(intersection) / as.numeric(intersectionWhich.area)) > 0.5)]]
hexagons.address.ocean <- hexagons.address[ ! hexagons.address %in% hexagons.address.land ]

## --------------------------------------
## Get sourceSink sites

shoreTest <- function(hexagon) {
  neighbors <- h3_k_ring(hexagon, 1)
  if( sum( neighbors %in% hexagons.address.land ) > 0 & sum( neighbors %in% hexagons.address.ocean ) > 0 ) { return(hexagon) } else { return(NULL) }
}

cl <- makeCluster(number.cores-2)
clusterExport(cl, c("hexagons.address.land","shoreTest","hexagons.address.ocean","h3_k_ring"))
hexagons.address.sourcesink <- unlist(unique(parLapply(cl, hexagons.address.land , function(x) { shoreTest(x) })))
stopCluster(cl)
closeAllConnections()
gc(reset=TRUE) 

hexagons.address.land <- hexagons.address.land[! hexagons.address.land %in% hexagons.address.sourcesink]
hexagons.address.ocean <- hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.sourcesink]

## --------------------------------------
## Remove inshore sourceSink sites

clusterPolygons <- function(polys, distance){
  dist_matrix = st_distance(polys, by_element = FALSE)
  class(dist_matrix) = NULL
  connected = dist_matrix == distance
  g = igraph::graph_from_adjacency_matrix(connected)
  return(components(g)$membership)
}

hexagons.address.sourcesink.poly <- h3_to_geo_boundary_sf(hexagons.address.sourcesink)
hexagons.address.sourcesink.poly$cluster = clusterPolygons(hexagons.address.sourcesink.poly, 0)
  
convexHull <- data.frame(id=unique(hexagons.address.sourcesink.poly$cluster),count=sapply(unique(hexagons.address.sourcesink.poly$cluster), function(x) { sum( hexagons.address.sourcesink.poly$cluster == x)} ))
convexHullTest <- convexHull[convexHull$count > 4 & convexHull$count < 500,"id"]

hexagons.address.sourcesink.remove <- character(0)
for( i in convexHullTest ){
  convexHullTest.i <- hexagons.address.sourcesink.poly[hexagons.address.sourcesink.poly$cluster == i,1]
  convexHullTest.i.sum <- nrow(hexagons.address.sourcesink.poly[hexagons.address.sourcesink.poly$cluster == i,1])
  convexHullTest.i <- st_remove_holes(st_cast(st_union(convexHullTest.i)))
  convexHullTest.i <- unlist(h3jsr::polyfill(convexHullTest.i, sim.resolution) )
  if( sum(convexHullTest.i %in% hexagons.address.ocean ) > sum(convexHullTest.i %in% hexagons.address.land ) )  { hexagons.address.sourcesink.remove <- unique(c(hexagons.address.sourcesink.remove,convexHullTest.i)) }
}

hexagons.address.sourcesink <- unique(hexagons.address.sourcesink[!hexagons.address.sourcesink %in% hexagons.address.sourcesink.remove])
hexagons.address.land <- unique(c(hexagons.address.sourcesink.remove,hexagons.address.land))
hexagons.address.ocean <- unique(hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.sourcesink.remove])

# plot(h3_to_geo_boundary_sf(hexagons.address.ocean),col="blue")
# plot(h3_to_geo_boundary_sf(hexagons.address.sourcesink),col="blue")
# plot(h3_to_geo_boundary_sf(hexagons.address.land),col="blue")

## --------------------------------------
## Remove outer frame

if( min.lon != -180 & max.lon != 180 & min.lat != -90 & max.lat != 90 ) {
  
  polygons.sourceSink <- h3_to_geo_boundary_sf(hexagons.address.sourcesink)
  polygons.land <- h3_to_geo_boundary_sf(hexagons.address.land)
  polygons.ocean <- h3_to_geo_boundary_sf(hexagons.address.ocean)
  
  centroid.land <- st_centroid(polygons.land)
  centroid.shore <- st_centroid(polygons.sourceSink)
  centroid.ocean <- st_centroid(polygons.ocean)
  
  buffer.val <- 0.5
  centroid.min.lon <- min(c(extent(polygons.sourceSink)[1],extent(polygons.land)[1],extent(polygons.ocean)[1]))
  centroid.max.lon <- max(c(extent(polygons.sourceSink)[2],extent(polygons.land)[2],extent(polygons.ocean)[2])) 
  centroid.min.lat <- min(c(extent(polygons.sourceSink)[3],extent(polygons.land)[3],extent(polygons.ocean)[3])) 
  centroid.max.lat <- max(c(extent(polygons.sourceSink)[4],extent(polygons.land)[4],extent(polygons.ocean)[4])) 
    
  buffer.remover <- which(st_coordinates(centroid.land)[,1] <= centroid.max.lon - buffer.val &
                          st_coordinates(centroid.land)[,1] >= centroid.min.lon + buffer.val &
                          st_coordinates(centroid.land)[,2] <= centroid.max.lat - buffer.val &
                          st_coordinates(centroid.land)[,2] >= centroid.min.lat + buffer.val )
  hexagons.address.land <- hexagons.address.land[buffer.remover]
  
  buffer.remover <- which(st_coordinates(centroid.shore)[,1] <= centroid.max.lon - buffer.val &
                          st_coordinates(centroid.shore)[,1] >= centroid.min.lon + buffer.val &
                          st_coordinates(centroid.shore)[,2] <= centroid.max.lat - buffer.val &
                          st_coordinates(centroid.shore)[,2] >= centroid.min.lat + buffer.val )
  hexagons.address.sourcesink <- hexagons.address.sourcesink[buffer.remover]
  
  buffer.remover <- which(st_coordinates(centroid.ocean)[,1] <= centroid.max.lon - buffer.val &
                          st_coordinates(centroid.ocean)[,1] >= centroid.min.lon + buffer.val &
                          st_coordinates(centroid.ocean)[,2] <= centroid.max.lat - buffer.val &
                          st_coordinates(centroid.ocean)[,2] >= centroid.min.lat + buffer.val )
  hexagons.address.ocean <- hexagons.address.ocean[buffer.remover]
  
  hexagons.address <- unique(c(hexagons.address.land,hexagons.address.sourcesink,hexagons.address.ocean))
  
}

## --------------------------------------

if( sum(hexagons.address.sourcesink %in% hexagons.address.ocean) > 0 | sum(hexagons.address.land %in% hexagons.address.ocean) > 0 | sum(hexagons.address.ocean %in% hexagons.address.land) > 0 ) { stop("Error :: 001") }
if( sum( length(hexagons.address.sourcesink) , length(hexagons.address.ocean) , length(hexagons.address.land) ) != length(hexagons.address) ) { stop("Error :: 003") }

## --------------------------------------
## --------------------------------------
## Produce source sink sites

cat("Type of source sink sites:",sourceSinkLocationType,"\n")

if(sourceSinkLocationType == "peripheral") {
  cl <- makeCluster(number.cores)
  clusterExport(cl, c("hexagons.address.sourcesink","h3_to_geo_boundary"))
  source.sink.hexagons <- parApply(cl,data.frame(hexagons.address.sourcesink),1, function(x) { data.frame(address=x,h3_to_geo_boundary(x[[1]])[[1]][,c("lng","lat")]) })
  stopCluster(cl) 
  source.sink.hexagons <- do.call(rbind, source.sink.hexagons)
}

if(sourceSinkLocationType == "centroid") {
  polygons.sourceSink <- h3_to_geo_boundary_sf(hexagons.address.sourcesink)
  source.sink.hexagons <- st_centroid(polygons.sourceSink)
  source.sink.hexagons <- data.frame(lng=st_coordinates(source.sink.hexagons)[,1],lat=st_coordinates(source.sink.hexagons)[,2],address=hexagons.address.sourcesink)
}
source.sink.hexagons <- source.sink.hexagons[which(!duplicated(source.sink.hexagons[,c("lng","lat")])),]

# plot(polygons.sourceSink[1,])
# points( source.sink.hexagons[,c("lng","lat")])

## -------------------------
## Test if over land

source.sink.hexagons.t <- source.sink.hexagons
coordinates( source.sink.hexagons.t) = ~lng+lat
crs(source.sink.hexagons.t) <- dt.projection
polygons.land.sp <- sf:::as_Spatial(h3_to_geo_boundary_sf(hexagons.address.land))
crs(polygons.land.sp) <- dt.projection
source.sink.hexagons.over.land <- over(source.sink.hexagons.t, polygons.land.sp)

if( ! 0 %in% dim(source.sink.hexagons.over.land) > 0 ) { stop("Error :: 994")}

## --------------------------------------

cat("Remove landmass source sink sites from model:",as.character(removeLandmassSourceSinkSites),"\n")

if( removeLandmassSourceSinkSites ){ 
  source.sink.xy <- data.frame(cells.id=source.sink.hexagons$address,x=source.sink.hexagons$lng,y=source.sink.hexagons$lat,source=0,stringsAsFactors = FALSE) 
}

if(! removeLandmassSourceSinkSites ){ 
  source.sink.xy <- data.frame(cells.id=source.sink.hexagons$address,x=source.sink.hexagons$lng,y=source.sink.hexagons$lat,source=1,stringsAsFactors = FALSE) 
}

if(sum(duplicated(source.sink.xy[,2:3])) > 0 ) { stop("Error :: PT001")}

cat("Hexagons listing and addresses: TRUE","\n")

## --------------------------------------
## Add additional source sink sites from bathymetric range

cat("Additional source sink sites from bathymetric range:",ifelse(!is.null(addSourceSinkRegionsBathymetry),"TRUE","FALSE"),"\n")

if( ! is.null(addSourceSinkRegionsBathymetry)) {
  
  stop("Review :: 001")
  
  bathymetryRaster <- raster(bathymetryRasterFile)
  bathymetryRaster <- crop(bathymetryRaster,extent(min.lon,max.lon,min.lat,max.lat))
  
  # bathymetryRaster <- aggregate(bathymetryRaster,6)
  
  bathymetryRasterCells <- Which(bathymetryRaster < addSourceSinkRegionsBathymetry[1] & bathymetryRaster >= addSourceSinkRegionsBathymetry[2], cells=TRUE)
  additionalSourceSinkRegions <- xyFromCell(bathymetryRaster,bathymetryRasterCells)
  additionalSourceSinkRegions <- as.data.frame(additionalSourceSinkRegions)
  
  cl <- makeCluster(number.cores)
  clusterExport(cl, c("geo_to_h3","additionalSourceSinkRegions","sim.resolution"))
  hexagons.address.additional <- unlist(unique(parLapply(cl, 1:nrow(additionalSourceSinkRegions) , function(x) { geo_to_h3(c(additionalSourceSinkRegions[x,2],additionalSourceSinkRegions[x,1]), sim.resolution ) })))
  stopCluster(cl)
  hexagons.address.additional <- unique(hexagons.address.additional)
  hexagons.address.additional <- hexagons.address.additional[hexagons.address.additional %in% hexagons.address.ocean]
  
  hexagons.address.ocean <- hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.additional]
  
  hexagons.address.sourcesink <- c(hexagons.address.sourcesink,hexagons.address.additional)
  
  # Add new land sites
  
  bathymetryRasterCells <- Which(bathymetryRaster >= addSourceSinkRegionsBathymetry[1], cells=TRUE)
  additionalSourceSinkRegions <- xyFromCell(bathymetryRaster,bathymetryRasterCells)
  additionalSourceSinkRegions <- as.data.frame(additionalSourceSinkRegions)
  
  cl <- makeCluster(number.cores)
  clusterExport(cl, c("h3","additionalSourceSinkRegions","sim.resolution"))
  hexagons.address.land.additional <- unlist(unique(parLapply(cl, 1:nrow(additionalSourceSinkRegions) , function(x) { geo_to_h3(c(additionalSourceSinkRegions[x,2],additionalSourceSinkRegions[x,1]), sim.resolution ) })))
  stopCluster(cl)
  hexagons.address.land.additional <- unique(hexagons.address.land.additional)
  hexagons.address.land.additional <- hexagons.address.land.additional[hexagons.address.land.additional %in% hexagons.address.ocean]
  hexagons.address.land <- c(hexagons.address.land,hexagons.address.land.additional)
  
  hexagons.address.ocean <- hexagons.address.ocean[! hexagons.address.ocean %in% hexagons.address.land.additional]
  
  ## ----------
  
  if(sourceSinkLocationType == "peripheral") {
    source.sink.hexagons.additional <- lapply(hexagons.address.additional ,function(x) { data.frame(address=x,h3_to_geo_boundary(x)[[1]][,c("lng","lat")]) })
    source.sink.hexagons.additional <- do.call(rbind, source.sink.hexagons.additional)
    source.sink.hexagons.additional <- source.sink.hexagons.additional[!duplicated(source.sink.hexagons.additional[,c("lng","lat")]),]
  }
  
  if(sourceSinkLocationType == "centroid") {
    source.sink.hexagons.additional <- h3_to_polygon(hexagons.address.additional)
    source.sink.hexagons.additional <- st_centroid(source.sink.hexagons.additional)
    source.sink.hexagons.additional <- data.frame(lng=st_coordinates(source.sink.hexagons.additional)[,1],lat=st_coordinates(source.sink.hexagons.additional)[,2],address=hexagons.address.additional)
  }
  
  additional.pts <- data.frame(cells.id=source.sink.hexagons.additional$address,x=source.sink.hexagons.additional$lng,y=source.sink.hexagons.additional$lat,source=1,stringsAsFactors = FALSE)
  source.sink.xy <- source.sink.xy[ ! source.sink.xy$cells.id %in% additional.pts$cells.id,]
  source.sink.xy <- rbind(source.sink.xy,additional.pts)
  
  ## ----------
}

## --------------------------------------
## Add additional source sink sites from shp

cat("Additional source sink sites from spatialPolygons:",ifelse(!is.null(addSourceSinkRegions),"TRUE","FALSE"),"\n")

if( ! is.null(addSourceSinkRegions) ) {
  
  stop("Review :: 001")
  
  additional.pts <- data.frame()
  
  for(i in 1:length(addSourceSinkRegions)){
    
    additional.shp <- shapefile(addSourceSinkRegions[i])
    additional.shp$ID <- 1:nrow(additional.shp)
    crs(additional.shp) <- dt.projection
    additional.shp <- st_as_sf(additional.shp)
    
    # If polygons
    
    if( as.character(st_geometry_type(additional.shp)[1]) != "POINT") {
      hexagons.address.additional <- unique(sapply(1:nrow(additional.shp) , function(x) { h3_to_parent(  polyfill(additional.shp[x,"ID"], sim.resolution + 1), sim.resolution  ) })) # ?? h3_to_parent :: sim.resolution + 1 
    }
    
    # If points
    
    if( as.character(st_geometry_type(additional.shp)[1]) == "POINT") {
      additional.shp.coords <- st_coordinates(additional.shp)
      hexagons.address.additional <- unique(sapply(1:nrow(additional.shp.coords) , function(x) { geo_to_h3(c(additional.shp.coords[x,2],additional.shp.coords[x,1]), sim.resolution  ) }))
    }
    
    if( TRUE %in% (hexagons.address.additional %in% hexagons.address.land) ) { 
      
      repositioning <- which(hexagons.address.additional %in% hexagons.address.land)
      for(repo in repositioning) {
        repo.dist <- h3_distance(hexagons.address.additional[repo], hexagons.address.sourcesink)
        hexagons.address.additional[repo] <- hexagons.address.sourcesink[which(repo.dist > 0)][which.min(repo.dist[which(repo.dist > 0)])]
      }
    }
    
    if( addSourceSinkRegionsForceCoast & FALSE %in% (hexagons.address.additional %in% hexagons.address.sourcesink) ) { 
      
      repositioning <- which(! hexagons.address.additional %in% hexagons.address.sourcesink)
      for(repo in repositioning) {
        repo.dist <- h3_distance(hexagons.address.additional[repo], hexagons.address.sourcesink)
        hexagons.address.additional[repo] <- hexagons.address.sourcesink[which(repo.dist > 0)][which.min(repo.dist[which(repo.dist > 0)])]
      }
    }
    
    if( TRUE %in% (hexagons.address.additional %in% hexagons.address.land) ) { stop("Error :: 138") }
    
    # plot(h3_to_polygon(hexagons.address.additional))
    
    if(sourceSinkLocationType == "peripheral") {
      source.sink.hexagons.additional <- lapply(hexagons.address.additional ,function(x) { data.frame(address=x,h3_to_geo_boundary(x)[[1]][,c("lng","lat")]) })
      source.sink.hexagons.additional <- do.call(rbind, source.sink.hexagons.additional)
      source.sink.hexagons.additional <- source.sink.hexagons.additional[!duplicated(source.sink.hexagons.additional[,c("lng","lat")]),]
    }
    
    if(sourceSinkLocationType == "centroid") {
      source.sink.hexagons.additional <- h3_to_polygon(hexagons.address.additional)
      source.sink.hexagons.additional <- st_centroid(source.sink.hexagons.additional)
      source.sink.hexagons.additional <- data.frame(lng=st_coordinates(source.sink.hexagons.additional)[,1],lat=st_coordinates(source.sink.hexagons.additional)[,2],address=hexagons.address.additional)
    }
    
    additional.pts.i <- data.frame(cells.id=source.sink.hexagons.additional$address,x=source.sink.hexagons.additional$lng,y=source.sink.hexagons.additional$lat,source=1,stringsAsFactors = FALSE)
    additional.pts <- rbind(additional.pts,additional.pts.i)
    
    ## ----------
    
  }
  
  additional.pts <- additional.pts[which(!duplicated(additional.pts[,c("x","y")])),]
  
  source.sink.xy <- rbind(source.sink.xy[ which(! source.sink.xy$cells.id %in% additional.pts$cells.id),],additional.pts)
  
}

## -----------------------------------------------------
## Polygon to Remove or Add release sites

cat("Mask source sink sites with spatialPolygons:",ifelse(!is.null(maskSourceSinkSites),"TRUE","FALSE"),"\n")

if( ! is.null(maskSourceSinkSites) ) {
  
  stop("Review :: 001")
  
  source.sink.xy.t <- source.sink.xy
  coordinates(source.sink.xy.t) <- c("x","y")
  crs(source.sink.xy.t) <- crs(worldmap)
  
  maskSourceSinkSites <- shapefile(maskSourceSinkSites)
  maskSourceSinkSites <- as(maskSourceSinkSites,"SpatialPolygons")
  crs(maskSourceSinkSites) <- crs(worldmap)
  
  points.over.polygon <- as.vector(which( ! is.na( sp::over( source.sink.xy.t , maskSourceSinkSites , fn = NULL) )) )
  
  if( length(points.over.polygon) > 0 ) { 
    
    if(maskSourceSinkSitesType == "include") { source.sink.xy[(1:nrow(source.sink.xy))[!(1:nrow(source.sink.xy)) %in% points.over.polygon],"source"] <- 0 }
    
    if(maskSourceSinkSitesType == "exclude") { source.sink.xy[points.over.polygon,"source"] <- 0 }
    
  }
  
}

## ---------------

if(sum(duplicated(source.sink.xy[,c("x","y")])) > 0) { stop("Error :: 007") }
if(sum(duplicated(source.sink.xy[,"cells.id"])) > 0) { stop("Error :: 008") }

## -----------------------------------------------------

hexagons.address.sourcesink <- unique(source.sink.xy[,"cells.id"])
hexagons.address.ocean <- unique(hexagons.address.ocean[ ! hexagons.address.ocean %in% hexagons.address.sourcesink ])
hexagons.address.land <- unique(hexagons.address.land[ ! hexagons.address.land %in% hexagons.address.sourcesink ])

if( length(hexagons.address.ocean) == 0 ) { hexagons.address.ocean <- geo_to_h3(c(max.lat,min.lon), sim.resolution  ) }

## -----------------------------------------------------
## Give new ids

unique.cells.id <- data.frame(old.id=unique(source.sink.xy$cells.id),new.id=1:length(unique(source.sink.xy$cells.id)))
source.sink.xy[,"cells.id"] <- sapply(source.sink.xy[,"cells.id"],function(x) { unique.cells.id[which(unique.cells.id$old.id == x),"new.id"]   })

# ggplot() + geom_sf(data=worldmap) + geom_point(data = source.sink.xy[source.sink.xy$source == 0,c("x","y")], aes(x = x, y = y), size = 1, shape = 1, fill = "darkred")
# ggplot() + geom_sf(data=worldmap) + geom_point(data = source.sink.xy[source.sink.xy$source == 1,c("x","y")], aes(x = x, y = y), size = 1, shape = 1, fill = "darkred")

# head(source.sink.xy)
# tail(source.sink.xy)

## ------------------

if( sum(hexagons.address.sourcesink %in% hexagons.address.ocean) > 0 | sum(hexagons.address.land %in% hexagons.address.ocean) > 0 | sum(hexagons.address.ocean %in% hexagons.address.land) > 0 ) { stop("Error :: 001") }

## ------------------

if( length(hexagons.address.ocean) <= 250000 ) { 
  polygons.plot.ocean <- h3_to_polygon(hexagons.address.ocean)
  polygons.plot.ocean <- st_wrap_dateline(polygons.plot.ocean, options = "WRAPDATELINE=YES", quiet = TRUE)
}

polygons.plot.land <- h3_to_polygon(hexagons.address.land)
polygons.plot.land <- st_wrap_dateline(polygons.plot.land, options = "WRAPDATELINE=YES", quiet = TRUE)

## ------------------

polygons.plot.sourcesink <- h3_to_polygon(hexagons.address.sourcesink)
polygons.plot.sourcesink <- st_wrap_dateline(polygons.plot.sourcesink, options = "WRAPDATELINE=YES", quiet = TRUE)

if( exists("polygons.plot.ocean")) {
  plot1 <- ggplot() + 
    geom_sf(data = polygons.plot.ocean, fill="#6EADEC", colour="#6EADEC", size=0.1) + 
    geom_sf(data = polygons.plot.land, fill="#E1BF6F", colour="#E1BF6F", size=0.1) + 
    geom_sf(data = polygons.plot.sourcesink, fill="Red", colour="#FFFFFF", size=0.1) + theme_bw()
  
}

if(! exists("polygons.plot.ocean")) {
  plot1 <- ggplot() + 
    geom_sf(data = polygons.plot.land, fill="#E1BF6F", colour="#E1BF6F", size=0.1) + 
    geom_sf(data = polygons.plot.sourcesink, fill="Red", colour="Red", size=0.1) + theme_bw()
}

cat("Export region of interest and hexagon list to shapefile: TRUE","\n")

hexagons.sourcesink.shp <- as(polygons.plot.sourcesink, "Spatial")
hexagons.sourcesink.shp$ID <- source.sink.xy[source.sink.xy$source == 1,"cells.id"]
hexagons.sourcesink.shp$IDSEQ <- 1:nrow(hexagons.sourcesink.shp)
hexagons.sourcesink.shp$HEX <- hexagons.address.sourcesink
shapefile(hexagons.sourcesink.shp, paste0(results.folder,"/sourceSinkSites.shp"), overwrite=TRUE)

## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------

cat("Export region of interest and hexagon list to RData: TRUE","\n")

initial.coords <- source.sink.xy[source.sink.xy$source == 1 , c("x","y") ]
source.cells.id <- source.sink.xy[source.sink.xy$source == 1,1]
save(source.sink.xy, file = paste0(results.folder,"/","sourceSinkSites.RData"))

options(warn=0)