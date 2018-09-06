## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

source("Dependences.R")

## ------------------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------------------

project.folder <- "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/TestScript"

## ------------------------------------

# Produce Animation

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",project.name,"SimulationResults.sql"))

project.name <- dbReadTable(sql, "Parameters")$project.name
sim.extent <- as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$extent, split=",")))

particles.lat <- dbReadTable(sql, "MovieLat")
particles.lon <- dbReadTable(sql, "MovieLon") # dbReadTable(sql, "MovieAlt")

movie.year <- dbReadTable(sql, "Parameters")$movie.year
months <- as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$sim.months , split=",")))
particles.to.sql.id <- as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$particles.to.sql.id , split=","))) 
movie.sites.id <- as.numeric(unlist(strsplit(dbReadTable(sql, "Parameters")$movie.sites.id , split=",")))

n.hours.per.day <- dbReadTable(sql, "Parameters")$n.hours.per.day
n.particles.per.cell <- dbReadTable(sql, "Parameters")$n.particles.per.cell
source.sink.xy <- dbReadTable(sql, "ReleaseSites")

dbDisconnect(sql)

# --------------------------------------

show.polygon.region.interest <- TRUE
sim.every.hours <- 24 / n.hours.per.day

min.lon <- sim.extent[1] ; max.lon <- sim.extent[2] ; min.lat <- sim.extent[3] ; max.lat <- sim.extent[4]
ratio <- abs(sim.extent[1]) +  abs(sim.extent[2]) : abs(sim.extent[4]) - abs(sim.extent[3])

source.sink.id <- 1:nrow(source.sink.xy)
particles.reference <- data.table( id = 1:(n.particles.per.cell * nrow(source.sink.xy) ) )
particles.reference[ , start.cell := as.numeric( sapply( source.sink.id ,function(x) { rep(x,n.particles.per.cell) })) ]

# ---------------------------------------------------------------------------------------------------------

land.polygon <- shapefile("Data/Shapefiles/Global Landmass.shp")
crs(land.polygon) <- "+proj=longlat +ellps=WGS84"
land.polygon <- crop(land.polygon, extent(sim.extent) )
plot(land.polygon, col="grey")

legend.x <- 3
legend.y <- 21
  
# ------------------

if( ! "Video" %in% list.files(paste0(project.folder,"/Results")) ) { dir.create(file.path(paste0(project.folder,"/Results/Video"))) }

# ------------------------------------------------------

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  while (format(date, format="%m") == m) { date <- date + 1 }
  return(as.integer(format(date - 1, format="%d")))
}

# --------------------------------------

days.months <- data.frame(
  day = unlist( sapply(months,function(x) { 1:numberOfDays( as.Date(  paste(movie.year,"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) } )) ,
  month = unlist(sapply(months,function(x) { rep(x[1],( numberOfDays( as.Date(  paste(movie.year,"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) ) ) } ))
)

# ------------------------------------------------------------------------------------------------------

t.steps <- nrow(particles.lat)
change.day.vect <- rep(1:n.hours.per.day,length.out=t.steps+sim.every.hours)[1:t.steps]
change.day <- rep(FALSE,length(change.day.vect))
change.day[change.day.vect == 1] <- TRUE

# ---------------------------------

colors <- c("#EB0505","#AB4343","#8C43AB" ,"#438CAB","#8CAC43" ,"#EB8005","#F0F009","#18DC2C","#189BDC","Black","#F587DB","#CB19A2","#725F2C","#725F2C","#D07FA7","#6683E2","#EB0505","#AB4343","#8C43AB" ,"#438CAB","#8CAC43","#EB8005","#F0F009","#18DC2C","#189BDC","Black")
cells.colors <- unique(particles.reference[ id %in% particles.to.sql.id , start.cell ])
cells.colors <- data.frame(cell=cells.colors,color=colors[1:length(cells.colors)], stringsAsFactors = FALSE)

# ---------------------------------

polygon.region.interest.xx <-  c( min(particles.lon[particles.lon != 0],na.rm=T) , min(particles.lon[particles.lon != 0],na.rm=T) , max(particles.lon[particles.lon != 0],na.rm=T) , max(particles.lon[particles.lon != 0],na.rm=T) )
polygon.region.interest.yy <-  c( min(particles.lat[particles.lat != 0],na.rm=T) , max(particles.lat[particles.lat != 0],na.rm=T) , max(particles.lat[particles.lat != 0],na.rm=T) , min(particles.lat[particles.lat != 0],na.rm=T) )

# ---------------------------------

show.polygon.region.interest <- FALSE
print.day <- 0

png(file=paste0(project.folder,"/Results/Video/Video map_%02d.png"), width=1280, height=720)
par( mar=c(0,0,0,0) , bg="#ffffff")

for( t in 1:t.steps) {
  
  if( change.day[t] ) { print.day <- print.day + 1 }
  
  print.date.sim <- format(as.Date(  paste(movie.year,"-",days.months[print.day,2],"-",days.months[print.day,1],sep="")  , "%Y-%m-%d"), "%d %b %Y")

  moving.ids <- particles.to.sql.id[which(particles.lon[t,] !=0)]
  moving.cell.ids <- particles.reference[ id %in% moving.ids , start.cell]
  moving.colors <- as.character( sapply(moving.cell.ids, function(x) { cells.colors[ cells.colors$cell %in% x , "color"] } ) )
  
  moving.lons <- particles.lon[t,] ; moving.lons <- moving.lons[moving.lons != 0]
  moving.lats <- particles.lat[t,] ; moving.lats <- moving.lats[moving.lats != 0]
  
  points.plot <- data.frame(Lon = moving.lons , Lat = moving.lats , color=moving.colors )
  points.plot <- points.plot[complete.cases(points.plot),]
  
  print(  plot(land.polygon , col="grey" , border="grey") )
  print(  text(legend.x,legend.y, paste0("t: ",print.date.sim), cex = 1, pos=4 , col="black") )
  print(  points(points.plot[,1], points.plot[,2], pch=16 , col=as.character(points.plot[,3]),cex=0.9) )
  
  if(show.polygon.region.interest) {
    
    polygon(polygon.region.interest.xx,polygon.region.interest.yy, lty =2)
    
  }
  
}

dev.off()

# ------------------------------------------------------

paste0(project.folder,"/Results/Video/Video map_%02d.png")

system( 'ffmpeg -s 1280x720 -i "/Volumes/Laminaria/Dropbox/Manuscripts/Transport simulations explain genetic differention of North Atlantic marine forests/TestScript/Results/Video/Video map_%02d.png" -vcodec libx264 -r 32 -pix_fmt yuv420p output.mp4 -y' )
file.remove( list.files(paste0(project.folder,"/Results/Video"),pattern="png",full.names=TRUE) )

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------