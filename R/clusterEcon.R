#' Cluster county-level economic data
library(rgeoda)
library(sf)
library(spatialEco)
library(terra)
library(cluster)
library(factoextra)

setwd("C:/evans/GITS/ecoecon")
  data.dir <- file.path(getwd(), "data")

#*********************************
# Read tabular data
econ <- read.csv(file.path(data.dir, "Oregon_CSD_Table.csv"))
  ch.idx <- which(lapply(1:ncol(econ), function(i) is.character(econ[,i])) == TRUE)[-c(1:2)]  
    if(length(ch.idx) > 0) econ <- econ[-ch.idx]
names(econ)[1:2] <- c("AFFGEOID", "NAME")

#*********************************
# Read census subdivision polygons
subdiv <- st_read(file.path(data.dir, "OR_subdivisions.shp"))
  subdiv <- st_transform(subdiv, st_crs("+proj=longlat +datum=WGS84 +no_defs"))
    subdiv <- merge(subdiv, econ, by="AFFGEOID")
  d <- st_drop_geometry(subdiv[,15:42])

#*********************************
# Read ecoregion data 
er <- rast(file.path(data.dir, "ecoregions.tif"))
er.dat <- read.csv(file.path(data.dir, "ecoregion_data.csv"))

#*********************************
# Screen collinearity
( cl <- collinear(d, p=0.75) )
  if(length(cl) > 0) 
    d <- d[,-which(names(d) %in% cl)]

( cl.test <- rfUtilities::multi.collinear(d, perm = TRUE, leave.out = TRUE, n = 999) )
    cl <- cl.test[cl.test$frequency > 0,]$variables

#*********************************
# PCA clustering
res.pca <- prcomp(d, scale = TRUE)
  fviz_eig(res.pca)

fviz_pca_var(res.pca, col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)  
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969")

# Coordinates, cos2 and contributions of 
#   individual observations
ind.coord <- res.pca$x	
center <- res.pca$center
scale <- res.pca$scale
getdistance <- function(ind_row, center, scale){
  return(sum(((ind_row-center)/scale)^2))
  }
d2 <- apply(d, 1, getdistance, center, scale)
cos2 <- function(ind.coord, d2){return(ind.coord^2/d2)}
ind.cos2 <- apply(ind.coord, 2, cos2, d2)

				
contrib <- function(ind.coord, comp.sdev, n.ind){
  100*(1/n.ind)*ind.coord^2/comp.sdev^2
}
ind.contrib <- t(apply(ind.coord, 1, contrib, 
                 res.pca$sdev, nrow(ind.coord)))

#head(ind.coord[,1:4])
#head(ind.cos2[,1:4])
#head(ind.contrib[,1:4])
				
subdiv$pca <- apply(ind.contrib, MARGIN=1, which.max)
  plot(subdiv["pca"])

pca <- rasterize(vect(subdiv), rast(er), field="pca",  background=NA, 
                touches=FALSE, update=FALSE, sum=FALSE, cover=FALSE)
writeRaster(pca, file.path(data.dir, "results", "pca_cluster.tif"), 
            overwrite=TRUE)

#*********************************
# Clustering
( opt <- optimal.k(daisy(d), nk = 20) )
n = 10

# agglomerative hierarchical clustering
ag <- agnes(daisy(d), diss = TRUE, metric = "manhattan", 
            method="complete", stand = TRUE)
subdiv$agnes <- cutree(as.hclust(ag), k = n) * 10000
  plot(subdiv["agnes"])

# Spatially structured clustering (2nd order contingency) 
queen_w <- queen_weights(subdiv, order=2)
subdiv$azp <- azp_greedy(n, queen_w, d)$Clusters * 10000
  plot(subdiv["azp"])

rc <- rasterize(vect(subdiv), rast(er), field="agnes",  background=NA, 
                touches=FALSE, update=FALSE, sum=FALSE, cover=FALSE)
writeRaster(rc, file.path(data.dir, "results", "econclust.tif"), 
            overwrite=TRUE)

#*********************************
# Merge ecoregions and clusters
r <- c(er, rc)
  r <- as(stack(r), "SpatialPixelsDataFrame")
  ( cts <- table(r@data[,1], r@data[,2]) )

ecoecon <- rc + er
writeRaster(ecoecon, file.path(data.dir, "results", "econecon.tif"), 
            overwrite=TRUE)



