
# Points over Polygon

library(spatialEco)
new_shape <- point.in.poly(pnts, ind_adm)

library(sp)
pt.in.poly <- sp::over(ind_adm, pnts, fn = NULL)

# Points along line

require(sp)

x <- c(18.25721, 18.25763,18.25808,18.25846,18.25864,18.25886,18.25892,18.25913,18.25940,18.25962,
       18.25976,18.25997,18.26021,18.26048,18.26061,18.26086,18.26107,18.26128,18.26154,18.26219,
       18.26276,18.26350,18.26445,18.26510,18.26584,18.26668,18.26704,18.26807,18.26850,18.26944,
       18.27020,18.27080,18.27111,18.27134,18.27168,18.27191,18.27217,18.27254,18.27309,18.27345,
       18.27368,18.27389,18.27398,18.27400,18.27392,18.27383,18.27370) 

y <- c(44.69540,44.69539,44.69544,44.69552,44.69563,44.69586,44.69608,44.69644,44.69672,44.69687
       ,44.69701,44.69718,44.69737,44.69763,44.69771,44.69778,44.69781,44.69781,44.69782,44.69776
       ,44.69772,44.69778,44.69794,44.69805,44.69814,44.69822,44.69824,44.69826,44.69821,44.69805
       ,44.69775,44.69737,44.69728,44.69717,44.69701,44.69687,44.69671,44.69649,44.69616,44.69598
       ,44.69578,44.69560,44.69539,44.69513,44.69490,44.69476,44.69453)

river<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
proj4string(river) <- CRS("+init=epsg:4326")


norm_vec <- function(x) sqrt(sum(x^2))
new_point <- function(p0, p1, di) { # Finds point in distance di from point p0 in direction of point p1
    v = p1 - p0
    u = v / norm_vec(v)
    return (p0 + u * di)
}

find <- function(river, M) {

  result = river[1,,drop=FALSE] 
  # for all subsequent points p1, p2 in this data.frame norm_vec(p2 - p1) = M at all times
  equidistantPoints = river[1,,drop=FALSE] 
  river = tail(river, n = -1)
  accDist = 0


  while (length(river) > 0) {
    point = river[1,]
    lastPoint = result[1,]

    dist = norm_vec(point - lastPoint)    

    if ( accDist + dist > M ) {
      np = new_point(lastPoint, point, M - accDist)
      equidistantPoints = rbind(np, equidistantPoints) # add np to equidistantPoints
      result = rbind(np, result) # add np to result
      accDist = 0 # reset accDist
    } else {
      #move point from river to result  
      river = tail(river, n = -1)
      result = rbind(point, result)    
      #update accDist  
      accDist = accDist + dist
    }
  }
  allPoints = result[NROW(result):1,] # reverse result
  return(list(newPoints = equidistantPoints, allPoints = allPoints))
}

r = cbind(x,y) 
result = find(r, 0.003)
plot(result$allPoints, type="l", col="red", asp = 1)
points(result$allPoints, col="red")
points(result$newPoints, col="cyan")
