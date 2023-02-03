## Prepare video (animation) points

plot2 <- NULL

if( ! is.null(movie.sites.xy) ) {
  
  if(class(movie.sites.xy) == "character") { movie.sites.xy <- as.data.frame(shapefile(movie.sites.xy))[,2:3] }
  movie.sites.xy <- movie.sites.xy[complete.cases(movie.sites.xy),]
  movie.sites.id <- unique(sort( as.vector(get.knnx( source.sink.xy[ source.sink.xy[,"source"] == 1,c("x","y") ] , movie.sites.xy , k = 1 + movie.sites.buffer , algorithm="kd_tree" )$nn.index) ))
  movie.sites.id <- source.sink.xy[ source.sink.xy[,"source"] == 1, "cells.id" ][movie.sites.id]
  movie.sites.xy <- source.sink.xy[ source.sink.xy[,"cells.id"] %in% movie.sites.id ,c("x","y") ]
  
  plot2 <- plot1 + geom_point(data = movie.sites.xy, aes(x = x, y = y), size = 2.5, shape = 16, color = "darkred")
  
  particles.reference.bm <- attach.big.matrix(particles.reference.bm.desc)
  particles.reference <- as.data.table(particles.reference.bm[,c(1,2,9)])
  colnames(particles.reference) <- c("id","start.cell","state")
  setkey(particles.reference, id )
  
  particles.video.id <- particles.reference[ start.cell %in% movie.sites.id , id ]
  # particles.video.location.dt <- data.table()
  
  rm(particles.reference.bm)
  rm(particles.reference)

  cat("Animation defined:", as.character(! is.null(movie.sites.xy)),"\n")
  
}

if(! exists("movie.sites.id")) { movie.sites.id <- NULL}
if(! exists("particles.video.id")) { particles.video.id <- NULL}
