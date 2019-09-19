## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------
##
## ------------------------------------------------------------------------------------------------------------------

source("0. Project Config.R")

sql.project.name <- "Atlantic"

number.cores <- 5

## ------------------------------------------------------------------------------------------------------------------

sql <- dbConnect(RSQLite::SQLite(), paste0(sql.directory,"/",sql.project.name,"SimulationResults.sql"))
Connectivity <- data.table(dbReadTable(sql, "Connectivity"))
source.sink.xy <- dbReadTable(sql, "SourceSinkSites")
dbDisconnect(sql)

Connectivity <- Connectivity[ , j=list(mean(Probability, na.rm = TRUE) , sd(Probability, na.rm = TRUE) , max(Probability, na.rm = TRUE) , mean(Time.mean, na.rm = TRUE) , sd(Time.mean, na.rm = TRUE) , max(Time.mean, na.rm = TRUE) , mean(Number.events, na.rm = TRUE) , sd(Number.events, na.rm = TRUE) , max(Number.events, na.rm = TRUE) ) , by = list(Pair.from,Pair.to)]
colnames(Connectivity) <- c("Pair.from" , "Pair.to" , "Mean.Probability" , "SD.Probability" , "Max.Probability" , "Mean.Time" , "SD.Time" , "Max.Time" , "Mean.events" , "SD.events" , "Max.events" )
Connectivity[is.na(Connectivity)] <- 0

Connectivity ; gc()

## -------------------

source.sink.xy <- source.sink.xy[source.sink.xy$cells.id %in% unique(c(Connectivity$Pair.from,Connectivity$Pair.to)),]

## -------------------

results <- data.frame(day=numeric(30),clusters=numeric(30),modul=numeric(30),signif=numeric(30))

for (n.days in 1:30) {
        
      print(n.days)
  
      comb <- Connectivity
      comb[ which(comb[,"Max.Time"] > n.days) , "Mean.Probability" ] <- 0
      comb <- comb[,c("Pair.from","Pair.to","Mean.Probability")]
      comb <- comb[comb$Pair.from != comb$Pair.to,]
      comb <- as.data.frame( comb[ sort( as.vector(unlist(comb[,"Mean.Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )
      
      net.function <<- prod
      graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
      
      E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,1/-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
      graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
      graph.obj <- as.undirected(graph.obj, mode = "collapse", edge.attr.comb = "min") # min / mean / max
      graph.obj <- simplify(graph.obj, remove.multiple = TRUE, remove.loops = TRUE)
      
      clustering.method <- "walktrap.community" # Uni: fastgreedy.community** walktrap.community leading.eigenvector.community Bi: walktrap.community edge.betweenness.community(slow)

      membership.graph <- get(clustering.method)(graph.obj)$membership

      modul.t <- modularity(graph.obj,membership.graph) 
      n.clusters.t <- length(unique(membership.graph))

      n.cells <- length(membership.graph)
      unique.members <- sort(unique( membership.graph ))
      # 
      # permutat <- 4999
      # mods <- sapply(1:permutat, function(i){
      #   modularity( graph.obj , sample( unique.members , n.cells , replace = TRUE) )
      # })
      
      # signif.t <- signif(sum( mods > modularity(graph,get(clustering.method)(graph)$membership) ) / permutat,digits=4)
    
      results[n.days,1] <- n.days
      results[n.days,2] <- n.clusters.t
      results[n.days,3] <- modul.t
      results[n.days,4] <- 0
      
}


pdf( file=paste0(project.folder,"Results/Days vs Clusters.pdf") , width = 9, height = 9 )

par(mfrow=c(1,1),mar = c(1, 1, 1, 1))

x <- results$day
y <- results$clusters/1000

par(mar = c(5, 5.5, 3, 3))
plot(bezierCurve(x,y,100)$x,bezierCurve(x,y,100)$y,col="#000000",type="l", lwd=1,axes=FALSE, lty=2,ylab="Clusters (x1000)",xlab="Dispersal period (day)")
points(x,y,pch=19, col="#A2A2A2", bg="#A2A2A2")

axis(2,las=2,col="White",col.ticks="Black")
axis(1,las=0,col="White",col.ticks="Black")
box()

dev.off()


## -------------------

bezierCurve <- function(x, y, n)
{
  outx <- NULL
  outy <- NULL
  
  i <- 1
  for (t in seq(0, 1, length.out=n))
  {
    b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    
    i <- i+1
  }
  
  return (list(x=outx, y=outy))
}

bez <- function(x, y, t)
{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  for (i in 0:n)
  {
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
  }
  
  return (list(x=outx, y=outy))
}
## -----------------------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------
