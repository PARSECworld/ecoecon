#' Cluster county-level economic data
library(rgeoda)
library(sf)
library(spatialEco)
library(raster)

setwd("C:/evans/GITS/ecoecon")
  data.dir <- file.path(getwd(), "data")

#*********************************
# Read tabular data
econ <- read.csv(file.path(data.dir, "Oregon_CSD_Table.csv"))
  econ$delta.gini <- econ[,"Gini_Index_2019"]  - econ[,"Gini_Index_2010"]                                
  econ$delta.pop <- econ[,"Total_Population_2019"] - econ[,"Total_Population_2010"]
  ch.idx <- which(lapply(1:ncol(econ), function(i) is.character(econ[,i])) == TRUE)[-c(1:2)]  
    if(length(ch.idx) > 0) econ <- econ[-ch.idx]
names(econ)[1:2] <- c("AFFGEOID", "NAME")

( cl <- collinear(econ[,-c(1:2)], p=0.75) )
  if(length(cl) > 0) 
    econ <- econ[,-which(names(econ) %in% cl)]

#*********************************
# Read spatial data and merge
subdiv <- st_read(file.path(data.dir, "OR_subdivisions.shp"))
  subdiv <- st_transform(subdiv, st_crs("+proj=longlat +datum=WGS84 +no_defs"))
    subdiv <- merge(subdiv, econ, by="AFFGEOID")

#*********************************
# Clustering
d <- st_drop_geometry(subdiv[,15:22])
( opt <- optimal.k(daisy(d), nk = 20) )
n = 10

# agglomerative hierarchical clustering
ag <- agnes(daisy(d), diss = TRUE, metric = "manhattan", 
            method="complete", stand = TRUE)
subdiv$agnes <- cutree(as.hclust(ag), k = n)
  plot(subdiv["agnes"])

# Spatially structured clustering (2nd order contingency) 
queen_w <- queen_weights(subdiv, order=2)
subdiv$azp <- azp_greedy(n, queen_w, d)$Clusters
  plot(subdiv["azp"])

#*********************************
# Read ecoregion data 
er <- rast(file.path(data.dir, "ecoregions.tif"))
er.dat <- read.csv(file.path(data.dir, "ecoregion_data.csv"))
rc <- rasterize(vect(subdiv), rast(er), field="agnes",  background=NA, 
                touches=FALSE, update=FALSE, sum=FALSE, cover=FALSE)
