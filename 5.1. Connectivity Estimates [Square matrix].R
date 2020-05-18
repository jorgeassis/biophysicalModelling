## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("0. Project Config.R")

## --------------------------------------------------------------------------------------------------------------
##
##
## 
## --------------------------------------------------------------------------------------------------------------

## Display Square matrices

distance.probability <- read.big.matrix(paste0(project.folder,"/Results/Connectivity.Distance.bm"))
distance.probability <- as.data.frame(distance.probability[,])
colnames(distance.probability) <- c("Pair.from","Pair.to","Probability","SD.Probability","Max.Probability","Mean.Time","SD.Time","Time.max","Mean.events","SD.events","Max.events","Distance")

source.sink.xy <- read.big.matrix(paste0(project.folder,"/Results/source.sink.bm"))
source.sink.xy <- as.data.frame(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

distance.probability.m <- acast(distance.probability, Pair.from ~ Pair.to , value.var = "Probability",fill=0)
distance.probability.m <- distance.probability.m[,as.numeric(rownames(distance.probability.m))]

## --------------------

sorting.sites <- which.max(source.sink.xy[,3])
sites.to.sort <- source.sink.xy[-sorting.sites,2:3]

for (i in 2:nrow(source.sink.xy)) {
  
  related.site <- sorting.sites[i-1]
  
  colsest <- sort(spDists(as.matrix(source.sink.xy[2:3]), as.matrix(source.sink.xy[related.site,2:3]), longlat = TRUE),index.return=TRUE)$ix
  sorting.sites <- c(sorting.sites , colsest[!colsest %in% sorting.sites][1] )
  
}

source.sink.xy <- source.sink.xy[sorting.sites,]

# sorting.sites <- sorting.sites[-which(sorting.sites == 707)]
# sorting.sites[sorting.sites > 707 ] <- sorting.sites[sorting.sites > 707 ] - 1

distance.probability.m <- distance.probability.m[ sorting.sites , rev( sorting.sites ) ]
  
## --------------------

file.sampling.sites<- shapefile("../locationsFinal.shp")
sampling.sites <- as.data.frame(file.sampling.sites)[,c("coords.x1","coords.x2")]
colnames(sampling.sites) <- c("Lon","Lat")
sampling.sites.names <- as.character(file.sampling.sites$Site.Code)

# file.sampling.sites<- "/Volumes/Jellyfish/Dropbox/Manuscripts/_ Under Revision/Rejection of the Abundant Centre Hypothesis in marine mussels/Data/Sampling_sites.txt"
# sampling.sites <- read.table(file.sampling.sites,sep=",",header=T)
# sampling.sites.names <- as.character(sampling.sites$Site)
# sampling.sites.names <- c(paste0("P",1:9),"M3")

sampling.sites <- sapply(1:nrow(sampling.sites), function(x) { which.min(spDistsN1(as.matrix(source.sink.xy[,c("Lon","Lat")]),as.matrix(sampling.sites[x,c("Lon","Lat")],ncol=2),longlat = TRUE)) } )

## --------------------

plot(raster(-log(distance.probability.m)))

distance.probability.m <- distance.probability.m / max(distance.probability.m)
distance.probability.m.reclass <- distance.probability.m

distance.probability.m.reclass[distance.probability.m > 0 & distance.probability.m <= 0.1] <- 1
distance.probability.m.reclass[distance.probability.m > 0.1 & distance.probability.m <= 0.25] <- 2
distance.probability.m.reclass[distance.probability.m > 0.25 & distance.probability.m <= 0.5] <- 3
distance.probability.m.reclass[distance.probability.m > 0.5 & distance.probability.m <= 0.75] <- 4
distance.probability.m.reclass[distance.probability.m > 0.75 & distance.probability.m <= 1] <- 5

distance.probability.m.reclass <- raster(distance.probability.m.reclass)
extent(distance.probability.m.reclass) <- c(-nrow(distance.probability.m.reclass),-1,-nrow(distance.probability.m.reclass),-1)

plot(distance.probability.m.reclass,col=c("white","#E8C25A","#E3A577","#D56932","#B41B1C","#8A0000"),axes=FALSE, box=FALSE)
box(which = "plot", lwd = 0.5)

for(p in 1:length(sampling.sites.names)) {
  site <- -sampling.sites[p]
  lines(c(-738,site),c(site,site),col="gray" , lty = 5)
  lines(c(site,site),c(-738,site),col="gray" , lty = 5)
  points(site,site,pch=16,col="black",bg="black")
  
  text(site, y = site, labels = sampling.sites.names[p],pos=2)

}
     
