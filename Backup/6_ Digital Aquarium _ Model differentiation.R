## ------------------------------------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------------------------------------
##
## Digital Aquarium 2016 Versio 3.0
## Simulation of Dispersal using ocean fields
##

## ------------------------------------------------------------------------------------------------------------------
## Main Configuration

gist.directory <- "/Volumes/Laminaria/Dropbox/Gist/One Aquarium V2.0" # Albacora Jellyfish
results.folder <- "/Volumes/Laminaria/Dropbox/Manuscripts/Phylogeographic patterns in the North Atlantic and Adjacent Seas/Dispersal simulations/Results/"

results.files <- "Bifurcaria.bifurcata.fst.60day"
distance.geographic.file <- "distance.Bifurcaria.bifurcata.txt"
distance.currents.file <- "conn.ss.Bifurcaria.bifurcata.60.days.txt"
distance.differentition.file <- "fst.Bifurcaria.bifurcata.txt"

## ------------------------


distance.geographic.file <- paste0(results.folder,distance.geographic.file)
distance.currents.file <- paste0(results.folder,distance.currents.file)
distance.differentition.file <- paste0(results.folder,distance.differentition.file)

## ------------------------------------------------------------------------------------------------------------------

distance.geographic <- read.table(distance.geographic.file,dec=",",sep=";" , stringsAsFactors=FALSE ) 
distance.currents <- read.table(distance.currents.file,dec=".",sep=";" , stringsAsFactors=FALSE, header=FALSE)
distance.differentition <- read.table(distance.differentition.file,dec=",",sep=";",header=TRUE )[,-1] 

distance.currents[distance.currents == 0] <- NA

matrix.cor <- data.frame()

for(i in 1:(nrow(distance.geographic)-1)) {
  for(j in (i+1):nrow(distance.geographic)) {
    if(i!=j) { matrix.cor <- rbind(matrix.cor , data.frame( Differentition=as.numeric(as.character(distance.differentition[i,j])),
                                                            Distance = as.numeric(as.character(distance.geographic[i,j])),
                                                            Ocean.currents.min = min( c( as.numeric( as.character(distance.currents[i,j])) , as.numeric( as.character(distance.currents[j,i])) ) , na.rm=T  ) ,
                                                            Ocean.currents.mean = mean( c( as.numeric( as.character(distance.currents[i,j])) , as.numeric( as.character(distance.currents[j,i])) ) , na.rm=T ) ,
                                                            Ocean.currents.max = max( c( as.numeric( as.character(distance.currents[i,j])) , as.numeric( as.character(distance.currents[j,i])) ) , na.rm=T )
                                                            ) )
      
    } } }

matrix.cor[,1] <- matrix.cor[,1]  / (1-matrix.cor[,1] ) # CHECK THIS AFTERWARDS
matrix.cor <- matrix.cor[ complete.cases(matrix.cor) , ]

matrix.cor[,3] <- log(matrix.cor[,3])
matrix.cor[,4] <- log(matrix.cor[,4])
matrix.cor[,5] <- log(matrix.cor[,5])

## ------------------------------------------------------------------------------------------------------------------

# SP
# matrix.cor[matrix.cor$Differentition>1.5,3] <- matrix.cor[matrix.cor$Differentition>1.5,3] -25
# matrix.cor[matrix.cor$Differentition>1.5,4] <- matrix.cor[matrix.cor$Differentition>1.5,4] -25
# matrix.cor[matrix.cor$Differentition>1.5,5] <- matrix.cor[matrix.cor$Differentition>1.5,5] -25

# FV
# matrix.cor[matrix.cor$Distance<40,2] <- matrix.cor[matrix.cor$Distance<40,2] + 50
# matrix.cor[matrix.cor$Differentition>0.75,3] <- matrix.cor[matrix.cor$Differentition>0.75,3] -50
# matrix.cor[matrix.cor$Differentition>0.75,4] <- matrix.cor[matrix.cor$Differentition>0.75,4] -50
# matrix.cor[matrix.cor$Differentition>0.75,5] <- matrix.cor[matrix.cor$Differentition>0.75,5] -50

# FC 68
# matrix.cor[matrix.cor$Differentition>7.5,3] <- matrix.cor[matrix.cor$Differentition>7.5,3] -75
# matrix.cor[matrix.cor$Differentition>7.5,4] <- matrix.cor[matrix.cor$Differentition>7.5,4] -75
# matrix.cor[matrix.cor$Differentition>7.5,5] <- matrix.cor[matrix.cor$Differentition>7.5,5] -75

# CC 94
# matrix.cor[matrix.cor$Differentition>1,3] <- matrix.cor[matrix.cor$Differentition>1,3] -27
# matrix.cor[matrix.cor$Differentition>1,4] <- matrix.cor[matrix.cor$Differentition>1,4] -27
# matrix.cor[matrix.cor$Differentition>1,5] <- matrix.cor[matrix.cor$Differentition>1,5] -27

# PP cox
# matrix.cor[matrix.cor$Differentition>10,3] <- matrix.cor[matrix.cor$Differentition>10,3] -47
# matrix.cor[matrix.cor$Differentition>10,4] <- matrix.cor[matrix.cor$Differentition>10,4] -47
# matrix.cor[matrix.cor$Differentition>10,5] <- matrix.cor[matrix.cor$Differentition>10,5] -47

# AN
# matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.min > -60,3] <- matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.min > -60,3] -47
# matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.mean > -60,4] <- matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.mean > -60,4] -47
# matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.max > -60,5] <- matrix.cor[matrix.cor$Differentition>0.3 & matrix.cor$Ocean.currents.max > -60,5] -47

# BB
# matrix.cor[matrix.cor$Differentition>2.5 ,3] <- matrix.cor[matrix.cor$Differentition>2.5 ,3] -37
# matrix.cor[matrix.cor$Differentition>2.5 ,4] <- matrix.cor[matrix.cor$Differentition>2.5 ,4] -37
# matrix.cor[matrix.cor$Differentition>2.5 ,5] <- matrix.cor[matrix.cor$Differentition>2.5 ,5] -37

## ------------------------------------------------------------------------------------------------------------------

# 8 8

pdf( file=paste0(results.folder,"Images/",results.files,".ibd.pdf") , width = 8, height = 8 )
plot(matrix.cor$Distance,matrix.cor$Differentition,axes=FALSE,pch = 16,bg="grey",col="#A5A5A5",xlab="Geographic distance (km)",ylab="Genetic differentiation (Fst / 1 - Fst)")
box() ; axis(2,las=2) ; axis(1,las=0)
data.x.axis <- sort(matrix.cor$Distance[!is.na(matrix.cor$Distance)],decreasing=FALSE)
new.data <- data.frame(Distance= data.x.axis )
lines( data.x.axis , predict.lm( lm(Differentition ~ Distance, data=matrix.cor) ,new.data  ) , lty=2 , col="Black")
dev.off()

cor.ibd <- cor(matrix.cor$Distance,matrix.cor$Differentition , use = "complete.obs",method="pearson")
fit.ibd <- lm(Differentition ~ Distance, data=matrix.cor)

## -----------------

pdf( file=paste0(results.folder,"Images/",results.files,".currents.min.pdf") , width = 8, height = 8 )
plot(matrix.cor$Ocean.currents.min,matrix.cor$Differentition , axes=FALSE,pch = 16,bg="grey",col="#A5A5A5",xlab="Ocean connectivity (min)",ylab="Genetic differentiation (Fst / 1 - Fst)")
box() ; axis(2,las=2) ; axis(1,las=0)
data.x.axis <- sort(matrix.cor$Ocean.currents.min[!is.na(matrix.cor$Ocean.currents.min)],decreasing=FALSE)
new.data <- data.frame(Ocean.currents.min= data.x.axis )
lines( data.x.axis , predict.lm( lm(Differentition ~ Ocean.currents.min, data=matrix.cor) ,new.data  ) , lty=2 , col="Black")
dev.off()

cor.min <- cor(matrix.cor$Ocean.currents.min,matrix.cor$Differentition , use = "complete.obs",method="pearson")
fit.min <- lm(Differentition ~ Ocean.currents.min, data=matrix.cor)

## -----------------

pdf( file=paste0(results.folder,"Images/",results.files,".currents.mean.pdf") , width = 8, height = 8 )
plot(matrix.cor$Ocean.currents.mean,matrix.cor$Differentition,axes=FALSE,pch = 21,bg="grey",col="#A5A5A5",xlab="Ocean connectivity (mean)",ylab="Genetic differentiation (Fst / 1 - Fst)")
box() ; axis(2,las=2) ; axis(1,las=0)
data.x.axis <- sort(matrix.cor$Ocean.currents.mean[!is.na(matrix.cor$Ocean.currents.mean)],decreasing=FALSE)
new.data <- data.frame(Ocean.currents.mean= data.x.axis )
lines( data.x.axis , predict.lm( lm(Differentition ~ Ocean.currents.mean, data=matrix.cor) ,new.data  ) , lty=2 , col="black")
dev.off()

cor.mean <- cor(matrix.cor$Ocean.currents.mean,matrix.cor$Differentition , use = "complete.obs",method="pearson")
fit.mean <- lm(Differentition ~ Ocean.currents.mean, data=matrix.cor)

## -----------------

pdf( file=paste0(results.folder,"Images/",results.files,".currents.max.pdf") , width = 8, height = 8 )
plot(matrix.cor$Ocean.currents.max,matrix.cor$Differentition,axes=FALSE,pch = 21,bg="grey",col="#A5A5A5",xlab="Ocean connectivity (max)",ylab="Genetic differentiation (Fst / 1 - Fst)")
box() ; axis(2,las=2) ; axis(1,las=0)
data.x.axis <- sort(matrix.cor$Ocean.currents.max[!is.na(matrix.cor$Ocean.currents.max)],decreasing=FALSE)
new.data <- data.frame(Ocean.currents.max= data.x.axis )
lines( data.x.axis , predict.lm( lm(Differentition ~ Ocean.currents.max, data=matrix.cor) ,new.data  ) , lty=2 , col="black")
dev.off()

cor.max <- cor(matrix.cor$Ocean.currents.max,matrix.cor$Differentition , use = "complete.obs",method="pearson")
fit.max <- lm(Differentition ~ Ocean.currents.min, data=matrix.cor)

## -----------------

data.table( AIC(fit.ibd) , summary(fit.ibd)$adj.r.squared , cor.ibd ,
            AIC(fit.min) , summary(fit.min)$adj.r.squared , cor.min ,
            AIC(fit.mean) , summary(fit.mean)$adj.r.squared , cor.mean ,
            AIC(fit.max) , summary(fit.max)$adj.r.squared , cor.max )

## ----------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------