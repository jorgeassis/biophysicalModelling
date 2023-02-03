## ------------------------------------------------------------------------------------------------------------------
## 
## 
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)
closeAllConnections()

library(plyr)
library(gridExtra)
library(grid)
library(AICcmodavg)
library(data.table)
library(rnaturalearth)
library(lme4)

theme_map2 <- theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),  
                    axis.title.x = element_blank(),  axis.title.y = element_blank(), panel.grid.major = element_blank() ,
                    text = element_text(size=13) ,
                    legend.title = element_blank() ,
                    legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                    legend.background = element_rect(fill="white", size=0.2, linetype="solid",  colour ="#979797"))

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=13) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)) ,
                   legend.title = element_blank() ,
                   legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                   legend.background = element_rect(fill="white", size=0.2, linetype="solid",  colour ="#979797"))

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------------------------

source.sink.xy <- loadRData(paste0("../Results/","sourceSinkSites.RData"))
source.sink.xy <- data.table(source.sink.xy[,])
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )

worldMap <- ne_countries(scale = 10, returnclass = "sp")
resultsDF <- read.csv(file="../Results/Genetics/mainData.csv")

regressionDF <- data.frame()
overallModelDF <- data.frame()
discardedEntries <- data.frame()

for( study in sort(unique(resultsDF$Study)) ) {
  
  resultsDF.i <- resultsDF[which(resultsDF$Study == study) , ]
  
  for( species in sort(unique(resultsDF.i$Species)) ) {
    
    resultsDF.ii <- resultsDF.i[which(resultsDF.i$Species == species) , ]
    
    for( marker in sort(unique(resultsDF.ii$Marker)) ) {
      
      for( index in sort(unique(resultsDF.ii$Index)) ) {
        
        resultsDF.iii <- resultsDF.ii[which(resultsDF.ii$Marker == marker & resultsDF.ii$Index == index) , ]
        
        pd <- 0
        
        harmonicIndex <- read.csv(paste0("../Results/connectivityExport/Sim000Days/harmonicIndex.csv"))
        betweennessIndex <- read.csv(paste0("../Results/connectivityExport/Sim000Days/betweennessIndex.csv"))
        
        resultsDF.iii <- resultsDF.iii[which(resultsDF.iii$propaguleDuration == pd ), ]
        
        if( nrow(resultsDF.iii) < 8) { discardedEntries <- rbind(discardedEntries,data.frame(study=study,species=species,marker=marker)); next }
        if( sum(resultsDF.iii$Differentiation != 0) == 0 ) { next }
        
        resultsDF.iiii <- resultsDF.iii
        
        betweenness <- apply(resultsDF.iiii[,c("from","to")],1,function(x) { mean( c(ifelse(length(betweennessIndex[betweennessIndex$Pair == x[1],2])>0,betweennessIndex[betweennessIndex$Pair == x[1],2],NA) , ifelse(length(betweennessIndex[betweennessIndex$Pair == x[2],2])>0,betweennessIndex[betweennessIndex$Pair == x[2],2],NA)) , na.rm=T )  })
        harmonic <- apply(resultsDF.iiii[,c("from","to")],1,function(x) { mean( c(ifelse(length(harmonicIndex[harmonicIndex$Pair == x[1],2])>0,harmonicIndex[harmonicIndex$Pair == x[1],2],NA) , ifelse(length(harmonicIndex[harmonicIndex$Pair == x[2],2])>0,harmonicIndex[harmonicIndex$Pair == x[2],2],NA)) , na.rm=T )  })
        
        resultsDF.iiii <- data.frame(Differentiation = resultsDF.iiii$Differentiation,
                                     connectivity = resultsDF.iiii$Connectivity.mean,
                                     directionality = (resultsDF.iiii$Connectivity.min - resultsDF.iiii$Connectivity.max) / resultsDF.iiii$Connectivity.min,
                                     betweenness = betweenness,
                                     harmonic = harmonic,
                                     site = resultsDF.iiii$from
        )
        
        resultsDF.iiii <- resultsDF.iiii[complete.cases(resultsDF.iiii),]
        
        # Overall matrix 
        overallModelDF <- rbind(overallModelDF,data.frame(Species=species, Marker=marker, Study=study, resultsDF.iiii))
        
        model.null <- lmer(Differentiation ~ 1  + (1|site), resultsDF.iiii, REML=F)
        model <- lmer(Differentiation ~ connectivity + directionality + (1|site), resultsDF.iiii, REML=F) # random effect on population 
        model.harmonic <- lmer(Differentiation ~ connectivity + directionality + harmonic + (1|site), resultsDF.iiii, REML=F) # random effect on population 
        model.betweenness <- lmer(Differentiation ~ connectivity + directionality + betweenness + (1|site), resultsDF.iiii, REML=F) # random effect on population 
        
        model.null <- lm(observed ~ predicted , data=data.frame(observed=resultsDF.iiii$Differentiation,predicted=predict(model.null)) )
        model <- lm(observed ~ predicted , data=data.frame(observed=resultsDF.iiii$Differentiation,predicted=predict(model)) )
        model.harmonic <- lm(observed ~ predicted , data=data.frame(observed=resultsDF.iiii$Differentiation,predicted=predict(model.harmonic)) )
        model.betweenness <- lm(observed ~ predicted , data=data.frame(observed=resultsDF.iiii$Differentiation,predicted=predict(model.betweenness)) )
        
        model.null.aic <- AIC(model.null)
        model.null.r2 <- summary(model.null)$adj.r.squared
        tryCatch( model.null.pVal <- lmp(model.null) , error=function(e) { model.null.pVal <<- 1 } )
        model.null.corr <- cor(predict(model.null),resultsDF.iiii$Differentiation , use = "complete.obs",method="pearson")
        
        model.aic <- AIC(model)
        model.r2 <- summary(model)$adj.r.squared
        model.pVal <- lmp(model)
        model.corr <- cor(predict(model),resultsDF.iiii$Differentiation , use = "complete.obs",method="pearson")
        
        model.harmonic.aic <- AIC(model.harmonic)
        model.harmonic.r2 <- summary(model.harmonic)$adj.r.squared
        model.harmonic.pVal <- lmp(model.harmonic)
        model.harmonic.corr <- cor(predict(model.harmonic),resultsDF.iiii$Differentiation , use = "complete.obs",method="pearson")
        
        model.betweenness.aic <- AIC(model.betweenness)
        model.betweenness.r2 <- summary(model.betweenness)$adj.r.squared
        model.betweenness.pVal <- lmp(model.betweenness)
        model.betweenness.corr <- cor(predict(model.betweenness),resultsDF.iiii$Differentiation , use = "complete.obs",method="pearson")
        
        best.model <- c("model.null.aic","model.aic","model.harmonic.aic","model.betweenness.aic")[which.min(c(model.null.aic,model.aic,model.harmonic.aic,model.betweenness.aic))]
        
        regressionDF <- rbind(regressionDF,
                              data.frame(Study=study,
                                         Species=species,
                                         Marker=marker,
                                         Index=index,
                                         n=length(unique(c(resultsDF.iii$from,resultsDF.iii$to))),
                                         nPairs=nrow(resultsDF.iiii),
                                         pd=pd,
                                         
                                         best.model=best.model,
                                         
                                         model.null.aic = model.null.aic,
                                         model.null.r2 = model.null.r2,
                                         model.null.pVal = model.null.pVal,
                                         model.null.corr = model.null.corr,
                                         
                                         model.aic = model.aic,
                                         model.r2 = model.r2,
                                         model.pVal = model.pVal,
                                         model.corr = model.corr,
                                         
                                         model.harmonic.aic = model.harmonic.aic,
                                         model.harmonic.r2 = model.harmonic.r2,
                                         model.harmonic.pVal = model.harmonic.pVal,
                                         model.harmonic.corr = model.harmonic.corr,
                                         
                                         model.betweenness.aic = model.betweenness.aic,
                                         model.betweenness.r2 = model.betweenness.r2,
                                         model.betweenness.pVal = model.betweenness.pVal,
                                         model.betweenness.corr = model.betweenness.corr,
                                         
                                         averageDifferentiation=mean(resultsDF.iiii$Differentiation)
                                         
                              ))
        
        resultsDFTest <- data.frame(observed=resultsDF.iiii$Differentiation,predictedOT=predict(model.betweenness))
        averageDifferentiation <- mean(resultsDFTest$observed)
        range.y <- range(resultsDFTest$observed)
        range.x <- range(resultsDFTest$predictedOT)
        
        plot1 <- ggplot(resultsDFTest, aes(x=predictedOT, y=observed, group = 1)) +
          geom_point(size=1.75,color="#5B5B5B", alpha = 0.5) +
          scale_y_continuous(limits = range.y) +
          geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=FALSE,size=0.35, linetype = "dashed") +
          ylab(paste0("Observed genetic differentiation")) + xlab("Predicted genetic differentiation [Oceanographic transport]") + mainTheme +
          theme_minimal(base_size = 14) + mainTheme + 
          annotate("label", alpha = 0.65, label.padding=unit(0.75, "lines"), x = range.x[1], y = range.y[2], hjust=0,vjust=1 , 
                   label = paste0("AICc: ",format(round(model.betweenness.aic, 3), nsmall = 3), "\nAdjusted R2: ",format(round(model.betweenness.r2, 3), nsmall = 3), "\np-Value: ",format.pval(model.betweenness.pVal, digits=3, eps = 0.001), "\nPearson's Corr.: ",format(round(model.betweenness.corr, 3), nsmall = 3)))
        
        resultsDFTest <- data.frame(observed=resultsDF.iiii$Differentiation,predictedOT=predict(model.harmonic))
        averageDifferentiation <- mean(resultsDFTest$observed)
        range.y <- range(resultsDFTest$observed)
        range.x <- range(resultsDFTest$predictedOT)
        
        plot2 <- ggplot(resultsDFTest, aes(x=predictedOT, y=observed, group = 1)) +
          geom_point(size=1.75,color="#5B5B5B", alpha = 0.5) +
          scale_y_continuous(limits = range.y) +
          geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=FALSE,size=0.35, linetype = "dashed") +
          ylab(paste0("Observed genetic differentiation")) + xlab("Predicted genetic differentiation [Oceanographic transport]") + mainTheme +
          theme_minimal(base_size = 14) + mainTheme + 
          annotate("label", alpha = 0.65, label.padding=unit(0.75, "lines"), x = range.x[1], y = range.y[2], hjust=0,vjust=1 , 
                   label = paste0("AICc: ",format(round(model.harmonic.aic, 3), nsmall = 3), "\nAdjusted R2: ",format(round(model.harmonic.r2, 3), nsmall = 3), "\np-Value: ",format.pval(model.harmonic.pVal, digits=3, eps = 0.001), "\nPearson's Corr.: ",format(round(model.harmonic.corr, 3), nsmall = 3)))
        
        if( ! dir.exists( paste0("../Results/Genetics/Models/") ) ) { dir.create(file.path( paste0("../Results/Genetics/Models/") ), showWarnings = FALSE, recursive=TRUE) } 
        
        # -----------------------------------
        
        # Map with populations and links
        
        source.sink.xy.i <- source.sink.xy[source.sink.xy$Pair %in% unique(c(resultsDF.iii$from,resultsDF.iii$to)),]
        
        x.min <- min(source.sink.xy.i$Lon) 
        x.max <- max(source.sink.xy.i$Lon) 
        y.min <- min(source.sink.xy.i$Lat) 
        y.max <- max(source.sink.xy.i$Lat) 
        
        x.min <- x.min - min(c(10, round( abs(diff(c(x.min , x.max))) + 2 )  ))
        x.max <- x.max + min(c(10, round( abs(diff(c(x.min , x.max))) + 2 )  ))
        y.min <- y.min - min(c(10, round( abs(diff(c(y.max , y.min))) + 2 )  ))
        y.max <- y.max + min(c(10, round( abs(diff(c(y.max , y.min))) + 2 )  ))
        
        if(abs(diff(c(x.max , x.min))) > abs(diff(c(y.max , y.min))) ) { 
          y.max <- ( y.min + (abs(diff(c(y.max , y.min))) / 2) ) + ( abs(diff(c(x.max , x.min))) / 2 )
          y.min <- ( y.min + (abs(diff(c(y.max , y.min))) / 2) ) - ( abs(diff(c(x.max , x.min))) / 2 )
        }
        
        if(abs(diff(c(x.max , x.min))) < abs(diff(c(y.max , y.min))) ) { 
          
          x.max <- ( x.min + (abs(diff(c(x.max , x.min))) / 2) ) + ( abs(diff(c(y.max , y.min))) / 2 )
          x.min <- ( x.min + (abs(diff(c(x.max , x.min))) / 2) ) - ( abs(diff(c(y.max , y.min))) / 2 )
          
        }
        
        regionMap <- crop(worldMap,extent(x.min,x.max,y.min,y.max))
        
        resultsDF.iii$fromCoordX <-  sapply(resultsDF.iii$from, function(x) { source.sink.xy[source.sink.xy$Pair == x , "Lon"]  } )
        resultsDF.iii$fromCoordY <-  sapply(resultsDF.iii$from, function(x) { source.sink.xy[source.sink.xy$Pair == x , "Lat"]  } )
        resultsDF.iii$toCoordX <-  sapply(resultsDF.iii$to, function(x) { source.sink.xy[source.sink.xy$Pair == x , "Lon"]  } )
        resultsDF.iii$toCoordY <-  sapply(resultsDF.iii$to, function(x) { source.sink.xy[source.sink.xy$Pair == x , "Lat"]  } )
        
        lineConnections <- list()
        lineStrenght <- numeric(0)
        
        for( l.i in sort(resultsDF.iii$Connectivity.mean, index.return=T)$ix ){
          lineStrenght <- c(lineStrenght,(resultsDF.iii[l.i,"Connectivity.mean"] ) * (-1) )
          pointFrom <- c(as.numeric(as.character(resultsDF.iii[l.i,"fromCoordX"])),as.numeric(as.character(resultsDF.iii[l.i,"fromCoordY"]))) 
          pointTo <- c(as.numeric(as.character(resultsDF.iii[l.i,"toCoordX"])),as.numeric(as.character(resultsDF.iii[l.i,"toCoordY"]))) 
          routes_sl <- gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 100, addStartEnd = TRUE, sp = TRUE, breakAtDateLine=TRUE)
          lineConnections = c(lineConnections,sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = l.i), match.ID = F))
        }
        
        lineConnectionsSp <- do.call(rbind, lineConnections)
        lineStrenght <- (lineStrenght - min(lineStrenght)) / max( lineStrenght - min(lineStrenght) )
        
        for( l.s in 1:length(lineStrenght) ) {
          
          if( lineStrenght[l.s] >= 0.8 ) { lineStrenght[l.s] <- 1 }
          if( lineStrenght[l.s] >= 0.6 & lineStrenght[l.s] < 0.8 ) { lineStrenght[l.s] <- 2 }
          if( lineStrenght[l.s] >= 0.4 & lineStrenght[l.s] < 0.6 ) { lineStrenght[l.s] <- 3 }
          if( lineStrenght[l.s] >= 0.2 & lineStrenght[l.s] < 0.4 ) { lineStrenght[l.s] <- 4 }
          if( lineStrenght[l.s] < 0.2 ) { lineStrenght[l.s] <- 5 }
          
        }
        
        myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
        
        plot3 <- ggplot() + 
          geom_polygon(data = regionMap, aes(x = long, y = lat, group = group), fill="#CDCDCD", colour = "#9E9E9E" , size=0.25 ) +
          geom_sf(data = st_as_sf(lineConnectionsSp) , size= lineStrenght / 10 , colour = myColors[lineStrenght]) +
          geom_point(data = source.sink.xy.i, aes(x = Lon, y = Lat), colour = "#000000",size=2.5) +
          theme_minimal(base_size = 14) + theme_map2 
        
        pdf(file=paste0("../Results/Genetics/Models/Map_",gsub(" / "," ",species)," _ ",marker," _ ",index," _ ",study,"_betweenness.pdf"),width=17,height=6)
        print(grid.arrange(plot3,plot1, ncol=2, top = textGrob(paste0(species, " (", marker,")"),hjust = 0,gp=gpar(fontsize=17,font=1))))
        dev.off()
        
        pdf(file=paste0("../Results/Genetics/Models/Map_",gsub(" / "," ",species)," _ ",marker," _ ",index," _ ",study,"_harmonic.pdf"),width=17,height=6)
        print(grid.arrange(plot3,plot2, ncol=2, top = textGrob(paste0(species, " (", marker,")"),hjust = 0,gp=gpar(fontsize=17,font=1))))
        dev.off()
        
      }
    }
  }
}

## ---------------

nrow(regressionDF)
mean(regressionDF$model.harmonic.r2)
sum(regressionDF$model.harmonic.pVal > 0.05)
discardedEntries
write.csv(regressionDF, file="../Results/Genetics/mainResults.csv")

## -----------------------------------
## -----------------------------------
## Global GLMM

head(overallModelDF)

overallModelDF$Species <- as.factor(overallModelDF$Species)
overallModelDF$Marker <- as.factor(overallModelDF$Marker)
overallModelDF$Study <- as.factor(overallModelDF$Study)
overallModelDF$site <- as.factor(overallModelDF$site)

overallModelDF <- overallModelDF[overallModelDF$Differentiation != 0,]

model.null <- lmer(Differentiation ~ 1  + (1|site) + (1|Study) + (1|Species) + (1|Marker), overallModelDF, REML=F)
model <- lmer(Differentiation ~ connectivity + directionality + (1|site) + (1|Study) + (1|Species) + (1|Marker), overallModelDF, REML=F) # random effect on population 
model.harmonic <- lmer(Differentiation ~ connectivity + directionality + harmonic + (1|site) + (1|Study) + (1|Species) + (1|Marker), overallModelDF, REML=F) # random effect on population 
model.betweenness <- lmer(Differentiation ~ connectivity + directionality + betweenness + (1|site) + (1|Study) + (1|Species) + (1|Marker), overallModelDF, REML=F) # random effect on population 

model.null <- lm(observed ~ predicted , data=data.frame(observed=overallModelDF$Differentiation,predicted=predict(model.null)) )
model <- lm(observed ~ predicted , data=data.frame(observed=overallModelDF$Differentiation,predicted=predict(model)) )
model.harmonic <- lm(observed ~ predicted , data=data.frame(observed=overallModelDF$Differentiation,predicted=predict(model.harmonic)) )
model.betweenness <- lm(observed ~ predicted , data=data.frame(observed=overallModelDF$Differentiation,predicted=predict(model.betweenness)) )

model.null.aic <- AIC(model.null)
model.null.r2 <- summary(model.null)$adj.r.squared
tryCatch( model.null.pVal <- lmp(model.null) , error=function(e) { model.null.pVal <<- 1 } )
model.null.corr <- cor(predict(model.null),overallModelDF$Differentiation , use = "complete.obs",method="pearson")

model.aic <- AIC(model)
model.r2 <- summary(model)$adj.r.squared
model.pVal <- lmp(model)
model.corr <- cor(predict(model),overallModelDF$Differentiation , use = "complete.obs",method="pearson")

model.harmonic.aic <- AIC(model.harmonic)
model.harmonic.r2 <- summary(model.harmonic)$adj.r.squared
model.harmonic.pVal <- lmp(model.harmonic)
model.harmonic.corr <- cor(predict(model.harmonic),overallModelDF$Differentiation , use = "complete.obs",method="pearson")

model.betweenness.aic <- AIC(model.betweenness)
model.betweenness.r2 <- summary(model.betweenness)$adj.r.squared
model.betweenness.pVal <- lmp(model.betweenness)
model.betweenness.corr <- cor(predict(model.betweenness),overallModelDF$Differentiation , use = "complete.obs",method="pearson")

resultsDFTest <- data.frame(observed=overallModelDF$Differentiation,predictedOT=predict(model.harmonic))
averageDifferentiation <- mean(resultsDFTest$observed)
range.y <- range(resultsDFTest$observed)
range.x <- range(resultsDFTest$predictedOT)

plot1 <- ggplot(resultsDFTest, aes(x=predictedOT, y=observed, group = 1)) +
  geom_point(size=1.75,color="#5B5B5B", alpha = 0.5) +
  scale_y_continuous(limits = range.y) +
  geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=FALSE,size=0.35, linetype = "dashed") +
  ylab(paste0("Observed genetic differentiation")) + xlab("Predicted genetic differentiation [Oceanographic transport]") + mainTheme +
  theme_minimal(base_size = 14) + mainTheme + 
  annotate("label", alpha = 0.65, label.padding=unit(0.75, "lines"), x = range.x[1], y = range.y[2], hjust=0,vjust=1 , 
           label = paste0("AICc: ",format(round(model.harmonic.aic, 3), nsmall = 3), "\nAdjusted R2: ",format(round(model.harmonic.r2, 3), nsmall = 3), "\np-Value: ",format.pval(model.harmonic.pVal, digits=3, eps = 0.001), "\nPearson's Corr.: ",format(round(model.harmonic.corr, 3), nsmall = 3)))

plot1
