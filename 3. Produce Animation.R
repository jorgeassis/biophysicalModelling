## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Config.R")
source("Dependences/mainFunctions.R")
 
## ------------------------------------------------------------------------------------------------------------------------------
## 
## # Video with Particle Flow
## # Video with assignments
## 
## ------------------------------------------------------------------------------------------------------------------------------

# Video with Particle Flow

load(paste0(results.folder,"/","particles.video.location.RData"))
head(particles.video.location.dt)

mainTitle <- paste0("Potential connectivity [Year ",movie.year,"]")
simulation.name <- paste0(particle.max.duration," days propagule duration")

load(paste0(results.folder,"/","sourceSinkSites.RData"))
load(paste0(results.folder,"/","modelParameters.RData"))

load(paste0(results.folder,"/particleReferenceTableNames.RData"))
particles.reference.bm.desc <- paste0(results.folder,"/particleReferenceTable.desc")
particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
particles.reference.bm <- as.data.frame(particles.reference.bm[])
colnames(particles.reference.bm) <- particles.reference.names
head(particles.reference.bm)

sim.extent <-unique(as.numeric(unlist(strsplit(global.simulation.parameters$extent, split=","))))
movie.year <- global.simulation.parameters$movie.year
months <- unique(as.numeric(unlist(strsplit(global.simulation.parameters$sim.months , split=","))))
movie.sites.id <- unique(as.numeric(unlist(strsplit(global.simulation.parameters$movie.sites.id , split=","))))
n.hours.per.day <- global.simulation.parameters$n.hours.per.day
n.particles.per.cell <- global.simulation.parameters$n.particles.per.cell

# --------------------------------------

show.polygon.region.interest <- TRUE
sim.every.hours <- 24 / n.hours.per.day

min.lon <- sim.extent[1] ; max.lon <- sim.extent[2] ; min.lat <- sim.extent[3] ; max.lat <- sim.extent[4]
ratio <- abs(sim.extent[1]) +  abs(sim.extent[2]) : abs(sim.extent[4]) - abs(sim.extent[3])

# ----------------------

movie.sites.ids <- particles.reference.bm[particles.reference.bm$start.cell %in% movie.sites.id, "id" ]

# ----------------------

particles.video.location.dt[,"start.cell"] <- particles.reference.bm[match(particles.video.location.dt$particle.id,particles.reference.bm$id),"start.cell"]

# ---------------------------------------------------------------------------------------------------------

videoBuffer <- 1

shapeVertex <- data.frame(Lon=c(min.lon,min.lon,max.lon,max.lon),
                          Lat=c(min.lat,max.lat,max.lat,min.lat))
shape <- spPolygons(as.matrix(shapeVertex))
crs(shape) <- dt.projection

if( is.null(alternativeLandmass) ) { land.polygon <- getMap(resolution = "high") }
if( ! is.null(alternativeLandmass) ) { land.polygon <- shapefile(alternativeLandmass) }

land.polygon <- gBuffer(land.polygon, byid=TRUE, width=0)
crs(land.polygon) <- dt.projection 
land.polygon <- crop(land.polygon,shape)

land.polygon@bbox <- as.matrix(extent(c(min.lon + ifelse(min.lon < 0 , videoBuffer * (-1) , videoBuffer), max.lon + ifelse(max.lon < 0 , videoBuffer * (-1) , videoBuffer),  min.lat + ifelse(min.lat < 0 , videoBuffer * (-1) , videoBuffer),max.lat + ifelse(max.lat < 0 , videoBuffer * (-1) , videoBuffer) )))
plot(land.polygon, col="grey")

# ------------------

if( ! "Video" %in% list.files(paste0(results.folder,"/")) ) { dir.create(file.path(paste0(results.folder,"/Video"))) }

# ------------------------------------------------------

t.steps <- max(particles.video.location.dt$t.step.movie)
change.day.vect <- rep(1:unique(n.hours.per.day),length.out=t.steps+unique(sim.every.hours))[1:t.steps]
change.day <- rep(FALSE,length(change.day.vect))
change.day[change.day.vect == 1] <- TRUE

# ---------------------------------

cells.colors <- movie.sites.id
cells.colors <- data.frame(cell=cells.colors,color=distinctColors(length(cells.colors)), stringsAsFactors = FALSE)

# ---------------------------------
# Aggregate colors

particles.video.location.dt[,"color"] <- cells.colors[match(particles.video.location.dt$start.cell,cells.colors$cell),"color"]
 
# dist <- spDists( as.matrix( coords.cell[,2:3] ) , as.matrix(coords.cell[,2:3] ) )
# dist <- as.dist(dist)
# mds.coor <- cmdscale(dist)
# plot(mds.coor)
# plot(hclust(dist(1-dist), method="single"))
# 
# k <- length(unique(cells.colors$cell))
# hc <- hclust(dist, "single")
# 
# cells.colors.i <- distinctColors( k )
# cells.colors <- data.frame(cell=coords.cell$cells.id,color=sapply(cutree(hc, k = k),function(x) cells.colors.i[x]))

# ---------------------------------

days.months <- data.frame(
  day = unlist( sapply(sort(months),function(x) { 1:numberOfDays( as.Date(  paste(unique(movie.year),"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) } )) ,
  month = unlist(sapply(sort(months),function(x) { rep(x[1],( numberOfDays( as.Date(  paste(unique(movie.year),"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) ) ) } ))
)

# ---------------------------------

show.additional.landmass.shp <- FALSE

if(show.additional.landmass.shp) {
  
  additional.landmass.shp <- shapefile(paste0("../",alternativeLandmass))
  
}

theme_map <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
    panel.grid.major = element_line(color = "#979797", size = 0.05),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#ffffff", color = NA), 
    panel.background = element_rect(fill = "#ffffff", color = NA), 
    legend.background = element_rect(fill = "#ffffff", color = NA),
    panel.border = element_blank()
  )

boxes <- data.frame(maxlat = max.lat+2,minlat = min.lat-2,maxlong = max.lon+2,minlong = min.lon-2, id="1")
boxes <- transform(boxes, laby=(min.lat + min.lat )/2, labx=(max.lon+min.lon )/2)

boxesSim <- data.frame(maxlat = max.lat,minlat = min.lat,maxlong = max.lon,minlong = min.lon, id="1")
boxesSim <- transform(boxesSim, laby=(min.lat + min.lat )/2, labx=(max.lon+min.lon )/2)

map <- ggplot() +
  geom_polygon(data = land.polygon , fill = "#A6A6A6", colour = "#ffffff" , size=0.15 ,  aes(long, lat, group = group)) +
  geom_rect(data=boxes, aes(xmin=min.lon-2 , xmax=max.lon+2, ymin=min.lat-2, ymax=max.lat+2 ), color="transparent", fill="transparent") +
  geom_rect(data=boxesSim, aes(xmin=min.lon , xmax=max.lon, ymin=min.lat, ymax=max.lat ), color="Black", fill="transparent", linetype = "dashed",size=0.12) +
  coord_equal() + theme_map # coord_map(projection = "mercator")
map

print.day <- 0

for( t in 1:t.steps) {
  
  cat('\014')
  cat('\n')
  cat('\n',paste0(rep("-",100),collapse = ""))
  cat('\n Generating animation files step #',t,'out of',t.steps)
  cat('\n',paste0(rep("-",100),collapse = ""))

  if( change.day[t] ) { print.day <- print.day + 1 }
  print.date.sim <- format(as.Date(  paste(movie.year,"-",days.months[print.day,2],"-",days.months[print.day,1],sep="")  , "%Y-%m-%d"), "%d %b %Y")

  points.plot <- particles.video.location.dt[particles.video.location.dt$t.step.movie == t,]
  points.plot <- data.frame(Lon = points.plot$pos.lon , Lat = points.plot$pos.lat , color=points.plot$color )
  points.plot <- points.plot[complete.cases(points.plot),]

  map.t <- map + geom_point(data = points.plot ,  aes(x = Lon, y = Lat) , shape = 21, colour = points.plot$color, fill = points.plot$color, size = 1.5, stroke = 0.35, alpha = 0.9) +
    annotate(geom="text", x=min.lon, y=max.lat + 5, label=mainTitle,size=5.5,family="Helvetica", color = "#000000",hjust = 0) +
    annotate(geom="text", x=min.lon, y=max.lat + 3, label=paste0("Simulation: ",print.date.sim),size=4.5,family="Helvetica", color = "#000000",hjust = 0)
  
  png(file=paste0(results.folder,"/Video/Video map_",t,".png"), width=1280, height=720, bg = "#f5f5f2")
  print(map.t)
  dev.off()
  
}
 
png_files <- sprintf(paste0(results.folder,"/Video/Video map_%d.png"), 1:t.steps)
av_encode_video(png_files, paste0(results.folder,"/animation.mp4"), framerate = 24)
unlink( paste0(results.folder,"/Video/") , recursive=TRUE)

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------