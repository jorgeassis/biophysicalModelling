## ------------------------------------------------------------------------
## Get data from Hycom do Netcdf
## Extracting region, depth and time
## ----------------------------------------------------------------------------------------------------------------

## Important

## 1. If needed download an extra month to allow driffting for the period
## 2. Expand ocean region (to allow full navigation), and use exlusion region for no new particles

## ----------------------------------------------------------------------------------------------------------------

## Dependences

setwd("/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/") # Albacora Jellyfish
source("Dependences.R")

dump.folder <- "/Volumes/Jellyfish/Dropbox/Gist/One Aquarium V2.0/"
dump.folder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Phylogeographic patterns in the North Atlantic and Adjacent Seas/Dispersal simulations/Data"

## ----------------

ocean.extent <- raster("Data/_ Atlantic/ocean.tif")
plot(ocean.extent , col="black")
min.lon <- extent(ocean.extent)[1] ; max.lon <- extent(ocean.extent)[2] ; min.lat <- extent(ocean.extent)[3] ; max.lat <- extent(ocean.extent)[4]

## ----------------

# Test
min.lon <- -13.569048  ; max.lon <- -6.713579  ; min.lat <- 36.631013 ; max.lat <- 40.210395
# Asia
min.lon <- 111.5 ; max.lon <- 145 ; min.lat <- 18 ; max.lat <- 48.5
# Atlantic
min.lon <- -25 ; max.lon <- 32 ; min.lat <- 22 ; max.lat <- 81.5

## ----------------

months.all <- 3:4
from.day <- 1 ; to.day <- 31
from.year <- 2003 ; to.year <- 2003
depth.range <- c(0)
buffer <- TRUE
buffer.val <- 0.05
results.file <- "Test" # Sargassum Atlantic
final.dimensions <- 2
variable <- "currents"

## ----------------------------------------------------------------------------------------------------------------

files <- "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"

nc <- NULL
while(is.null(nc)) {
  sink("/dev/null") 
  try( nc <- nc_open( files , verbose=FALSE ) , TRUE)
  sink()
}

time.corrected <- as.Date(ncvar_get( nc, "time")/24, origin = "2000-01-01") 

years <- as.numeric(substr(time.corrected,1,4))
months <- as.numeric(substr(time.corrected,6,7))
days <- as.numeric(substr(time.corrected,9,10))
from.year <- min((from.year:to.year)[from.year:to.year %in% years])
to.year <- max((from.year:to.year)[from.year:to.year %in% years])

# -----------------------------------------------------------------

if( buffer == TRUE ) {  
  min.lon <- min.lon - buffer.val
  max.lon <- max.lon + buffer.val
  min.lat <- min.lat - buffer.val
  max.lat <- max.lat + buffer.val
}

depths <- ncvar_get( nc, "depth")
depth <-  which(depths %in% depth.range)
Latitude <- ncvar_get( nc, "lat")
Longitude <- ncvar_get( nc, "lon")

array.region.lon <- which(Longitude >= min.lon & Longitude <= max.lon)
array.region.lat <- which(Latitude >= min.lat & Latitude <= max.lat)
Longitude.array <- Longitude[array.region.lon]
Latitude.array <- Latitude[array.region.lat]

i.min <- min(array.region.lon)
j.min <- min(array.region.lat)
i.max <- max(array.region.lon)
j.max <- max(array.region.lat)

# -----------------------------------------------------------------

if( variable == "currents" ) {
          
          for (y in from.year:to.year){
            
            time.window <- time.corrected[ which( years == y & months %in% months.all & days >= from.day & days <= to.day )]          
            time.file.array <- which(time.corrected %in% time.window)
            
            time.stamp <- vector()
            
            for( k in 1:length(time.file.array)) {
              
              t <- time.file.array[k]
              u <- array(data = NA, dim = c( i.max-i.min+1 , j.max-j.min+1 , length(depth) ))
              v <- array(data = NA, dim = c( i.max-i.min+1 , j.max-j.min+1 , length(depth) ))

              for( d in 1:length(depth) ) {
                
                values.to.place.u <- NULL
                values.to.place.v <- NULL
                
                while(is.null(values.to.place.u)) {
                  
                  try( values.to.place.u <- ncvar_get( nc, "water_u", start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1)))
                  
                }
                
                while(is.null(values.to.place.v)) {
                  
                  try( values.to.place.v <- ncvar_get( nc, "water_v", start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1)))
                  
                }
                
                u[,,d] <- values.to.place.u
                v[,,d] <- values.to.place.v
                
              }
              
              if(dim(u)[3] == 1) { u <- u[,,1] }
              if(dim(v)[3] == 1) { v <- v[,,1] }
              
              # For 2d data
              
              if( k == 1 & is.na(dim(u)[3]) )  {  
                
                data.u <- array(data = NA, dim = c(dim(u)[1],dim(u)[2],length(time.window)) )
                data.v <- array(data = NA, dim = c(dim(v)[1],dim(v)[2],length(time.window)) ) 
                
                }  
              
              
              if( is.na(dim(u)[3]) )  {  
                
                data.u[,,k] <- u
                data.v[,,k] <- v 
                
                }     
              
              # For 3d data
              
              if(  k == 1 & !is.na(dim(u)[3]) )  { 
                
                      data.u <- array(data = NA, dim = c(dim(u)[1],dim(u)[2],dim(u)[3],length(time.window)), dimnames = c("Lon","Lat","Depth","Time"))
                      data.v <- array(data = NA, dim = c(dim(v)[1],dim(v)[2],dim(u)[3],length(time.window)), dimnames = c("Lon","Lat","Depth","Time")) 
                      
                      }
              
              if( !is.na(dim(u)[3]) )  {  data.u[,,,k] <- u
                                          data.v[,,,k] <- v }
            }
            
            #--------------------------------
            # Export to netcdf
            
            dimX <- ncdim_def( "X", "unit", 1:dim(data.u)[1] )
            dimY <- ncdim_def( "Y", "unit", 1:dim(data.u)[2] )
            dimT <- ncdim_def( "Time", "days", as.numeric(time.window) )
            dimD <- ncdim_def( "Depths", "meters", depth )
            
            mv <- -999 # missing value to use
            
            var1d <- ncvar_def( "Longitude", "Degrees", dimX, mv)
            var2d <- ncvar_def( "Latitude", "Degrees", dimY, mv)
            var3d <- ncvar_def( "Date", "days", dimT, mv)
            var4d <- ncvar_def( "Depth", "meters", dimD, mv)
            
            # 2D
            
            if(final.dimensions == 2) {
              var4d <- ncvar_def( "UComponent", "m", list(dimX,dimY,dimT), mv, prec="double")
              var5d <- ncvar_def( "VComponent", "m", list(dimX,dimY,dimT), mv,  prec="double")
              nc.file <- nc_create( paste(dump.folder,"/currents_2d_",results.file,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d) )
              ncvar_put( nc.file, var1d, Longitude.array )  
              ncvar_put( nc.file, var2d, Latitude.array )   
              ncvar_put( nc.file, var3d, as.numeric(time.window) )
              ncvar_put( nc.file, var4d, data.u )  
              ncvar_put( nc.file, var5d, data.v )
              nc_close(nc.file)
            }
            
            if(final.dimensions == 3) {
              var5d <- ncvar_def( "UComponent", "m", list(dimX,dimY,dimD,dimT), mv, prec="double")
              var6d <- ncvar_def( "VComponent", "m", list(dimX,dimY,dimD,dimT), mv,  prec="double")
              nc.file <- nc_create( paste(dump.folder,"/currents_3d_",results.file,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d,var6d) )
              ncvar_put( nc.file, var1d, Longitude.array )  
              ncvar_put( nc.file, var2d, Latitude.array )   
              ncvar_put( nc.file, var3d, as.numeric(time.window) )
              ncvar_put( nc.file, var4d, as.numeric(depths[depth]) )
              ncvar_put( nc.file, var5d, data.u )  
              ncvar_put( nc.file, var6d, data.v )
              nc_close(nc.file)
            }
            
          }
          
}

# -----------------------------------------------------------------

if(variable != "currents") {
          
          for (y in from.year:to.year){
            
            time.window <- time.corrected[ which( years == y & months %in% months.all & days >= from.day & days <= to.day )]          
            time.file.array <- which(time.corrected %in% time.window)
            
            time.stamp <- vector()
            
            for( k in 1:length(time.file.array)) {
              
              t <- time.file.array[k]
              
              var.data <- array(data = NA, dim = c( i.max-i.min+1 , j.max-j.min+1 , length(depth) ))

              for( d in 1:length(depth) ) {
                
                values.to.place <- NULL
                
                while(is.null(values.to.place)) {
                  
                  try( values.to.place <- ncvar_get( nc, variable, start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1))
                       , TRUE)
                  
                }

                var.data[,,d] <- values.to.place
                values.to.place <- NULL
                
              }
              
              # For 2d data
              
              if( k == 1 & is.na(dim(var.data)[3]) )  {  data.u <- array(data = NA, dim = c(dim(u)[1],dim(u)[2],length(time.window)), dimnames = c("Lon","Lat","Time"))
              data.v <- array(data = NA, dim = c(dim(v)[1],dim(v)[2],length(time.window)), dimnames = c("Lon","Lat","Time")) }  
              
              
              if( is.na(dim(var.data)[3]) )  {  data.variable[,,k] <- var.data }     
              
              # For 3d data
              
              if( k == 1 & !is.na(dim(var.data)[3]) )  {  data.variable <- array(data = NA, dim = c(dim(var.data)[1],dim(var.data)[2],dim(var.data)[3],length(time.window)), dimnames = c("Lon","Lat","Depth","Time"))  }

              if( !is.na(dim(var.data)[3]) )  { data.variable[,,,k] <- var.data  }
              
            }
            
            #--------------------------------
            # Export to netcdf
            
            dimX <- ncdim_def( "X", "unit", 1:dim(data.variable)[1] )
            dimY <- ncdim_def( "Y", "unit", 1:dim(data.variable)[2] )
            dimT <- ncdim_def( "Time", "days", as.numeric(time.window) )
            dimD <- ncdim_def( "Depths", "meters", depth )
            
            mv <- -999 # missing value to use
            
            var1d <- ncvar_def( "Longitude", "Degrees", dimX, mv)
            var2d <- ncvar_def( "Latitude", "Degrees", dimY, mv)
            var3d <- ncvar_def( "Date", "days", dimT, mv)
            var4d <- ncvar_def( "Depth", "meters", dimD, mv)
            
            # 2D
            
            if(final.dimensions == 2) {
              var4d <- ncvar_def( variable, "unit", list(dimX,dimY,dimT), mv, prec="double")
              nc.file <- nc_create( paste("Data/",variable,"_2d_",results.file,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d) )
              ncvar_put( nc.file, var1d, Longitude.array )  
              ncvar_put( nc.file, var2d, Latitude.array )   
              ncvar_put( nc.file, var3d, as.numeric(time.window) )
              ncvar_put( nc.file, var4d, data.variable )  
              nc_close(nc.file)
            }
            
            if(final.dimensions == 3) {
              var5d <- ncvar_def( variable, "unit", list(dimX,dimY,dimD,dimT), mv, prec="double")
              nc.file <- nc_create( paste(dump.folder,"/",variable,"_3d_",results.file,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d) )
              ncvar_put( nc.file, var1d, Longitude.array )  
              ncvar_put( nc.file, var2d, Latitude.array )   
              ncvar_put( nc.file, var3d, as.numeric(time.window) )
              ncvar_put( nc.file, var4d, as.numeric(depths[depth]) )
              ncvar_put( nc.file, var5d, data.variable )  
              nc_close(nc.file)
            }
            
          }
          
}

## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------