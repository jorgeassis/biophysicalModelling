## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------



## CHECK FOR WRONG DATA > 100




source("0. Project Config.R")
source("Dependences.R")

## -----------------------------------------

files <- "http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1"

nc <- NULL
increment <- 0
while(is.null(nc)) {
  sink("/dev/null") 
  try( nc <- nc_open( files , verbose=FALSE ) , TRUE)
  sink()
  increment <- increment + 1
  if(increment == 1000) { stop("1000 attempts failed!")}
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

      values.to.place.u[values.to.place.u == -30000 ] <- NA
      
      while(is.null(values.to.place.v)) {
        
        try( values.to.place.v <- ncvar_get( nc, "water_v", start=c(i.min,j.min,depth[d],t), count=c((i.max-i.min+1),(j.max-j.min+1),1,1)))
        
      }
      
      values.to.place.v[values.to.place.v == -30000 ] <- NA
      
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
    nc.file <- nc_create( paste(project.folder,"/Data","/currents_2d_",project.name,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d) )
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
    nc.file <- nc_create( paste(project.folder,"/Data","/currents_3d_",project.name,"_",y,".nc",sep=""), list(var1d,var2d,var3d,var4d,var5d,var6d) )
    ncvar_put( nc.file, var1d, Longitude.array )  
    ncvar_put( nc.file, var2d, Latitude.array )   
    ncvar_put( nc.file, var3d, as.numeric(time.window) )
    ncvar_put( nc.file, var4d, as.numeric(depths[depth]) )
    ncvar_put( nc.file, var5d, data.u )  
    ncvar_put( nc.file, var6d, data.v )
    nc_close(nc.file)
  }
  
}
          
# -----------------------------------------------------------------
# -----------------------------------------------------------------
