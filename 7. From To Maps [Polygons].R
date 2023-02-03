## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

if( ! exists("pipeLiner") ) {
  
  closeAllConnections()
  rm(list=(ls()[ls()!="v"]))
  gc(reset=TRUE)
  source("0. Config.R")
  source("Dependences/mainFunctions.R")
  
  list.dirs(path = paste0("../Results"), recursive = FALSE)
  n.season <- "" # c("YearRound","SeasonSummer","SeasonWinter")
  spawn.p <- months.all
  pld.period <- c <- 30
  
}

## --------------------------

if( is.null(alternativeLandmass) ) { worldMap <- ne_countries(scale = "medium", returnclass = "sp") }
if( ! is.null(alternativeLandmass) ) { worldMap <- shapefile(alternativeLandmass) }

load(paste0(results.folder,"/","sourceSinkSites.RData"))
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
# source.sink.xy <- source.sink.xy[source.sink.xy$Source == 1,]

hexagons.sourcesink.shp <- shapefile(paste0(results.folder,"/sourceSinkSites.shp"))
if( "SOURCE" %in% names(hexagons.sourcesink.shp) ) { hexagons.sourcesink.shp <- hexagons.sourcesink.shp[hexagons.sourcesink.shp$SOURCE==1,"ID"] } 
names(hexagons.sourcesink.shp)
plot(hexagons.sourcesink.shp)

# Open connectivity

load(file=paste0(results.folder,"/particlePairedConnectivityAveraged",n.season,"TableNames.RData"))

Connectivity.desc <- paste0(results.folder,"/particlePairedConnectivityAveraged",n.season,"Table.desc")
Connectivity <- attach.big.matrix(Connectivity.desc)
Connectivity <- as.data.table(Connectivity[])
colnames(Connectivity) <- particles.connectivity.names

Connectivity <- Connectivity[Connectivity$Max.Time <= pld.period,]
Connectivity
max(Connectivity$Max.Time)

Connectivity <- Connectivity[Connectivity$Pair.from %in% source.sink.xy$Pair,]
Connectivity <- Connectivity[Connectivity$Pair.to %in% source.sink.xy$Pair,]

## ---------------------------------

worldMap <- crop(worldMap,extent(c(min.lon-5,max.lon+5,min.lat-5,max.lat+5))) 

## ---------------------------------

additional.source.sink <- shapefile(maskSourceSinkSites)
additional.source.sink$ID <- additional.source.sink$id
additional.source.sink <- additional.source.sink[,"ID"]
coordRef <- crs(additional.source.sink)

closest.source.sink.sites <- coordinates(additional.source.sink)
closest.source.sink.sites <- spDists(as.matrix(source.sink.xy[,c("Lon","Lat")]),closest.source.sink.sites, longlat = TRUE)
closest.source.sink.sites <- apply(closest.source.sink.sites,1,which.min)
source.sink.xy$Poly <- closest.source.sink.sites
source.sink.xy[source.sink.xy$Source == 0, "Poly"] <- 0

Connectivity$Poly.from <- sapply( Connectivity$Pair.from,function(x) { as.numeric(source.sink.xy[which(source.sink.xy[,"Pair"] == x),"Poly"]) })
Connectivity$Poly.to <- sapply( Connectivity$Pair.to,function(x) { as.numeric(source.sink.xy[which(source.sink.xy[,"Pair"] == x),"Poly"]) })

# Retention degree

load(paste0(results.folder,"/modelParameters.RData"))
n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell

# 1 Mau; 2 Tej;  2 Cad; 3 4 Fr

plot(worldMap)
points(source.sink.xy[source.sink.xy$Poly == 4,2:3])

polygonNames <- c("Mau","Por","Fra","Spa")
resulConnectivityProb <- matrix(NA,ncol=length(polygonNames),nrow=length(polygonNames))
resulConnectivityTime <- matrix(NA,ncol=length(polygonNames),nrow=length(polygonNames))
colnames(resulConnectivityTime) <- polygonNames
colnames(resulConnectivityProb) <- polygonNames
rownames(resulConnectivityTime) <- polygonNames
rownames(resulConnectivityProb) <- polygonNames
  
for( i in 1:length(polygonNames)) {
  for( j in 1:length(polygonNames)) {
    resulConnectivityProb[i,j] <- mean(unlist(Connectivity[ Poly.from == i & Poly.to == j,"Mean.Probability"]))
    resulConnectivityTime[i,j] <- mean(unlist(Connectivity[ Poly.from == i & Poly.to == j,"Mean.Time"]))
  }
}

write.csv(resulConnectivityProb,paste0(results.folder,"/resulConnectivityPolygonsProb.csv"), row.names = FALSE)
write.csv(resulConnectivityTime,paste0(results.folder,"/resulConnectivityPolygonsTime.csv"), row.names = FALSE)

## ------------------------------------------------------------------------------------------------------------

additional.source.sink <- gBuffer(additional.source.sink, byid=TRUE, width=0)

## -----------------

plot(worldMap , col="Gray",border="Gray")
plot(additional.source.sink , col="Black",border="Black", add=TRUE)

## --------------------------------------------------------------------------------

# GGPLOT with connections

theme_map <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank(),
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.position="bottom", 
    legend.box = "horizontal",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,30,-10),
    legend.key.height= unit(0.25, 'cm'),
    legend.key.width= unit(0.75, 'cm') )

# extent(worldMap) + c(20,-20,20,-20)

source.sink.xy <- as.data.table(source.sink.xy)

coords.0 <- source.sink.xy[,2:3]
coords.1 <- source.sink.xy[Pair %in% unique(unlist(Connectivity[ Poly.from == 1 ,"Pair.to"])) , .(Lon,Lat)]
coords.2 <- source.sink.xy[Pair %in% unique(unlist(Connectivity[ Poly.from == 2 ,"Pair.to"])) , .(Lon,Lat)]
coords.3 <- source.sink.xy[Pair %in% unique(unlist(Connectivity[ Poly.from == 3 ,"Pair.to"])) , .(Lon,Lat)]
coords.4 <- source.sink.xy[Pair %in% unique(unlist(Connectivity[ Poly.from == 4 ,"Pair.to"])) , .(Lon,Lat)]

coords.0.h3 <- h3_to_geo_sf(apply(coords.0,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.1.h3 <- h3_to_geo_sf(apply(coords.1,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.2.h3 <- h3_to_geo_sf(apply(coords.2,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.3.h3 <- h3_to_geo_sf(apply(coords.3,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.4.h3 <- h3_to_geo_sf(apply(coords.4,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))

# ---

coords.1.s <- source.sink.xy[Source == 1 & Poly == 1 , .(Lon,Lat) ] 
coords.2.s <- source.sink.xy[Source == 1 & Poly == 2 , .(Lon,Lat) ] 
coords.3.s <- source.sink.xy[Source == 1 & Poly == 3 , .(Lon,Lat) ] 
coords.4.s <- source.sink.xy[Source == 1 & Poly == 4 , .(Lon,Lat) ] 

coords.1.s.h3 <- h3_to_geo_sf(apply(coords.1.s,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.2.s.h3 <- h3_to_geo_sf(apply(coords.2.s,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.3.s.h3 <- h3_to_geo_sf(apply(coords.3.s,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))
coords.4.s.h3 <- h3_to_geo_sf(apply(coords.4.s,1,function(coords) { geo_to_h3( c(unlist(coords)[2],unlist(coords)[1]), res=sim.resolution) }))

colors <- c("#6FBBE8","#42AE57","#CB911C","#9914B3")

mapRegion <- ggplot() +
  geom_polygon(data = worldMap , fill = "#DADADA", colour = "#ffffff" , size=0.2 ,  aes(long, lat, group = group)) +

  geom_sf(data = coords.0.h3 , fill = "#A5A5A5", colour = "#A5A5A5" , size=0.1) +
  
  geom_sf(data = coords.1.h3 , fill = "#6FBBE8", colour = "#6FBBE8" , size=1) +
  geom_sf(data = coords.2.h3 , fill = "#42AE57", colour = "#42AE57" , size=1) +
  geom_sf(data = coords.3.h3 , fill = "#CB911C", colour = "#CB911C" , size=1) +
  geom_sf(data = coords.4.h3 , fill = "#9914B3", colour = "#9914B3" , size=1) +

  geom_sf(data = coords.1.s.h3 , fill = "#6FBBE8", colour = "Black" , size=0.1) +
  geom_sf(data = coords.2.s.h3 , fill = "#42AE57", colour = "Black" , size=0.1) +
  geom_sf(data = coords.3.s.h3 , fill = "#CB911C", colour = "Black" , size=0.1) +
  geom_sf(data = coords.4.s.h3 , fill = "#9914B3", colour = "Black" , size=0.1) +
  
  geom_polygon(data = additional.source.sink_df , size=0.5 , color="Black" , aes(long, lat, group = group, fill = factor(id) )) +
  scale_fill_manual(guide = guide_legend(title="Source / Sink regions of oceanographic connectivity", direction = "horizontal", title.position = "top", title.hjust = 0.5) , values =colors, labels=c("MAU","PT","FR","SP"), aesthetics = c( "fill")) +
  
  coord_sf(crs = "+init=epsg:3035") +
  #coord_map('lambert', lat0=extent(worldMap)[3], lat1=extent(worldMap)[4], xlim=c(extent(worldMap)[1], extent(worldMap)[2]), ylim=c(extent(worldMap)[3], extent(worldMap)[4])) + 
  theme_map

pdf(file=paste0(results.folder,"/resulConnectivityPolygons.pdf"), width=12)
plot(mapRegion)
dev.off()

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------