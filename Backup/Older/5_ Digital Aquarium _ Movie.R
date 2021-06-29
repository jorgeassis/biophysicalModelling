## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2016 Versio 2.0
## Simulation of Dispersal using ocean fields
##
## Assis, et al. 2015
##
## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
##

gist.directory <- "/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/"
raw.data.dir <- "/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/Data/"
results.directory <- "/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/Results/"
results.file <- "Test"
source(paste0(gist.directory,"Dependences.R"))

# --------------------------------------

show.polygon.region.interest <- TRUE
results.folder <- "SQLite"
file.name <- "Test"

min.lon <- -13.569048 ; max.lon <- -6.713579 ; min.lat <- 36.631013 ; max.lat <- 40.210395

parcticles.to.sql.years <- 2003
sim.every.hours <- 2

# --------------------------------------

new.extent <- extent( min.lon , max.lon , min.lat , max.lat ) # extent( -35.0 , 2.5 , 20.0 , 54.0 ) # NULL
ratio <- abs(new.extent[1]) +  abs(new.extent[2]) : abs(new.extent[4]) - abs(new.extent[3])

# ---------------------------------------------------------------------------------------------------------

raw.data.currents.files <- list.files(paste0(raw.data.dir,"/currents/"),full.names = TRUE,pattern="nc")
nc <- nc_open( raw.data.currents.files[1] , verbose=FALSE )
nc.date <- as.Date(ncvar_get( nc, "Time"), origin = "1970-01-01")
nc_close( nc )
months <- min(as.numeric(substr(nc.date, 6, 7))):max(as.numeric(substr(nc.date, 6, 7)))

land.polygon <- shapefile("/Volumes/Jellyfish/Dropbox/Raw Data/Shapefiles/World Present/GSHHS/GSHHS_f_L1.shp")
crs(land.polygon) <- "+proj=longlat +ellps=WGS84"
land.polygon <- crop(land.polygon, new.extent )
plot(land.polygon, col="grey")

# ------------------

if( ! "Video" %in% list.files(results.directory) ) { dir.create(file.path(paste0(results.directory,"/Video"))) }

# ------------------------------------------------------

numberOfDays <- function(date) {
  
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}

# --------------------------------------

days.months <- data.frame(
  
  day = unlist( sapply(months,function(x) { 1:numberOfDays( as.Date(  paste(parcticles.to.sql.years,"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) } )) ,
  month = unlist(sapply(months,function(x) { rep(x[1],( numberOfDays( as.Date(  paste(parcticles.to.sql.years,"-",x[1],"-","01",sep="")  , "%Y-%m-%d") ) ) ) } ))
  
)

# ------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/",results.folder,"/",paste0(file.name,".reference.particles.sql")))
particles.lat <- dbGetQuery(sql, "SELECT * FROM lat") ; dim(particles.lat)
particles.lon <- dbGetQuery(sql, "SELECT * FROM lon") ; dim(particles.lon)
particles.reference <- dbGetQuery(sql, "SELECT * FROM Particles_reference_table") ; dim(particles.reference)
dbDisconnect(sql)

sql <- dbConnect(RSQLite::SQLite(), paste0(results.directory,"/",results.folder,"/",paste0(file.name,".reference.particles.sql")))
particles.reference.source <- dbGetQuery(sql, "SELECT * from Particles_reference_table" ) ; dim(particles.reference.source)
dbDisconnect(sql)

particles.reference.source <- particles.reference.source[particles.reference.source$id %in% as.vector(unlist(particles.reference)), ]
particles.reference.source <- particles.reference.source[particles.reference.source$Year == parcticles.to.sql.years, ]

# ------------------------------------------------------

particles.time.x <- matrix(NA,ncol=max(particles.reference.source$t.finish),nrow=nrow(particles.lat))
particles.time.y <- matrix(NA,ncol=max(particles.reference.source$t.finish),nrow=nrow(particles.lat))

for(i in 1:nrow(particles.lat)) {
  
  id.i <- as.vector(unlist(particles.reference))[i]
  path.i <- particles.lon[i,][!is.na(particles.lon[i,])]
  particles.time.x[ i , particles.reference.source$t.start[particles.reference.source$id == id.i]:(particles.reference.source$t.start[particles.reference.source$id == id.i]+length(path.i)-1) ] <- path.i
  path.i <- particles.lat[i,][!is.na(particles.lat[i,])]
  particles.time.y[ i , particles.reference.source$t.start[particles.reference.source$id == id.i]:(particles.reference.source$t.start[particles.reference.source$id == id.i]+length(path.i)-1) ] <- path.i
  
}

# ------------------------------------------------------

t.steps.all <- ncol(particles.time.x)
change.day.vect <- rep(1:(24/sim.every.hours),length.out=t.steps.all+sim.every.hours)[1:t.steps.all]
change.day <- rep(FALSE,length(change.day.vect))
change.day[change.day.vect == 1] <- TRUE

particles.color <- c("#EB0505","#AB4343","#8C43AB","#438CAB","#8CAC43" ,"#EB8005","#F0F009","#18DC2C","#189BDC","Black","#F587DB","#CB19A2","#725F2C","#725F2C","#D07FA7","#6683E2","#EB0505","#AB4343","#8C43AB" ,"#438CAB","#8CAC43","#EB8005","#F0F009","#18DC2C","#189BDC","Black")

particles.color <- particles.color[1:length(unique(particles.reference.source$cell))]
length(particles.color) == length(unique(particles.reference.source$cell))

particles.color <- data.frame(cell=unique(particles.reference.source$cell),color=particles.color[1:length(unique(particles.reference.source$cell))])

polygon.region.interest.xx <-  c( min(particles.time.x,na.rm=T) , min(particles.time.x,na.rm=T) , max(particles.time.x,na.rm=T) , max(particles.time.x,na.rm=T) )
polygon.region.interest.yy <-  c( min(particles.time.y,na.rm=T) , max(particles.time.y,na.rm=T) , max(particles.time.y,na.rm=T) , min(particles.time.y,na.rm=T) )

global.extent <- extent( min(min(polygon.region.interest.xx),extent(land.polygon)[1]) , max(max(polygon.region.interest.xx),extent(land.polygon)[2]) , min(min(polygon.region.interest.yy),extent(land.polygon)[3]) , max(max(polygon.region.interest.yy),extent(land.polygon)[4])  )
global.extent <- as(global.extent, 'SpatialPolygons')  

png(file=paste0(results.directory,"Video/Video map_%02d.png"), width=1280, height=720)
par( mar=c(0,0,0,0) , bg="#ffffff")

print.day <- 0

for( t in 1:t.steps.all) {
  
  if( change.day[t] ) { print.day <- print.day + 1 }

  print.date.sim <- format(as.Date(  paste(parcticles.to.sql.years,"-",days.months[print.day,2],"-",days.months[print.day,1],sep="")  , "%Y-%m-%d"), "%d %b %Y")
  particles.color.t <- sapply(  as.vector(unlist(particles.reference)) , function(x) {  as.character(particles.color[particles.color[,1] == particles.reference.source[particles.reference.source$id == x,]$cell ,2]) } )

  points.plot <- data.frame(Lon = particles.time.x[,t] , Lat = particles.time.y[,t],color=particles.color.t )
  points.plot <- points.plot[complete.cases(points.plot),]
  
  print(  plot(global.extent , col="white" , border="white") )
  print(  plot(land.polygon , col="grey" , border="grey" , add=TRUE) )
  print(  text(140,18, paste0("t: ",print.date.sim), cex = 1.2, pos=4 , col="black") )
  print( points(points.plot[,1:2], pch=16 , col=as.character(points.plot[,3]),cex=0.9) )
  
  if(show.polygon.region.interest) {
     
    polygon(polygon.region.interest.xx,polygon.region.interest.yy, lty =2)
        
  }
}

dev.off()

# ------------------------------------------------------

results.directory
system( 'ffmpeg -s 1280x720 -i "/Volumes/Laminaria/Dropbox/Manuscripts/Genetic diversity drivers of Sargassum thunbergii/Results/Video/Video map_%02d.png" -vcodec libx264 -r 32 -pix_fmt yuv420p output.mp4 -y' )
file.remove( list.files(paste0(results.directory,"/Video"),pattern="png",full.names=TRUE) )

# ------------------------------------------------------
