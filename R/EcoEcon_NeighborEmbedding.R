#' kNN neighbor embedding of county-level economic data
suppressMessages(invisible(lapply(c("sf", "terra", "spatialEco", 
                 "umap", "ggplot2", "rgl", "e1071", "viridis", 
				 "factoextra", "cluster", "gridExtra"), 
				 require, character.only=TRUE)))

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
  econ <- st_drop_geometry(subdiv[,c(1,15:49)])

econ.names <- names(econ)[-c(1,34:36)]
  post <- sort(econ.names[grep("2019", econ.names)]) 
  pre <- sort(econ.names[which(!econ.names %in% post)])  
  new <- sort(unlist(lapply(strsplit(pre, "_"), function(x) paste(x[-length(x)], 
              collapse = "_") ) ) )
pecon <- data.frame(parm=new, pre=pre, post=post)
  delta <- list()
    for(i in 1:nrow(pecon)) {
      y1 <- pecon$pre[i]
  	  y2 <- pecon$post[i] 
      delta[[i]] <- 
  	    as.numeric(econ[,y2]) -
        as.numeric(econ[,y1])
    }
  names(delta) <- pecon$parm
econ.delta <- as.data.frame(do.call(cbind, delta))
  names(econ.delta) <- paste(names(econ.delta), "_delta")

econ <- data.frame(AFFGEOID = econ[,1], econ[,post], 
                        econ.delta, econ[,c(34:36)])
 
#*********************************
#   Index and remove NA's
#   Screen collinearity
na.idx <- which(is.na(econ), arr.ind = TRUE)
  cat(paste(names(econ)[unique(na.idx[,2])], "\n"),
      "columns have NA values", "\n")

if(nrow(na.idx) > 0){
  cat("removing rows:", unique(na.idx[,1]), "\n")
  econ <- econ[-unique(na.idx[,1]),]   
}

#( cl <- collinear(na.omit(econ[,-1]), p=0.75) )
#  if(length(cl) > 0) 
#    econ <- econ[,-which(names(econ) %in% cl[-c(6,7)])]

#*********************************
# Read ecoregional data 
er <- rast(file.path(data.dir, "ecoregions.tif"))
  er.dat <- read.csv(file.path(data.dir, "ecoregion_data.csv"))
    er.dat <- er.dat[which(er.dat$Value %in% values(er)),]

#*********************************
# Aggregate ecoregional data to political boundaries
#   not what we want but, for visualization 
( classes <- sort(terra::unique(er)[,1]) )
eco <- extract(er, vect(subdiv))
  d <- as.data.frame(do.call(rbind,tapply(eco[,2], eco$ID, function(x) { 
                     prop.table(table(factor(x, levels=classes)))})))
  names(d) <- paste0("ER", names(d)) 
  
er.maj  <- apply(d, 1, function(x) names(d)[which.max(x)])  
subdiv$Value <- as.numeric(substring(er.maj, 3, 6))
subdiv$er.prop <- apply(d, 1, function(x) x[which.max(x)] ) 

ecoecon <- merge(st_drop_geometry(subdiv), er.dat, by="Value") 
# apply(ecoecon[,64:72], 2, function(x) length(unique(x)) )

#***********************************************************
# Create UMAP Manifold then project to  
#   new statistical space
d = ncol(econ)
k = 12

custom.settings <- umap.defaults 
  custom.settings$metric = c("euclidean", "manhattan", "cosine", "pearson")[2] 
  custom.settings$n_components = d
  custom.settings$random_state = 42
  custom.settings$n_neighbors = k  

eumap <- umap(econ[,-1], config=custom.settings)

#***********************************************************
# Cluster data, assign clusters to subdiv

# Search for optimal k
asw <- numeric(10)
  for (k in 2:10) {
    asw[k] <- cluster::pam(scale(eumap[["layout"]]), k)$silinfo$avg.width
  }
k.opt <- which.max(asw)

p <- cmeans(scale(eumap[["layout"]]), centers = k.opt, iter.max = 100, 
            dist = "euclidean", m = 2)
  pval <- p$membership			
	p <- data.frame(cluster=p$cluster, pval)

# New table with clusters, pvalue, majority realm, fraction of realm,
#   first two UMAP embeddings and, covariates used in UMAP embeddings
lab = ecoecon[,c("Realm_World_Ecosystem", "LF_ClassName")[1]]
  lab = gsub("Nearctic ", "", lab)
pidx <- c(1,2)

xy <- data.frame(subdiv$AFFGEOID, st_coordinates(st_centroid(subdiv))[,1:2])
  names(xy) <- c("AFFGEOID", "long", "lat") 
econ <- data.frame(AFFGEOID = econ[,1], cluster=p$cluster, pvalue=apply(pval, 1, max),
                   realm = lab[-unique(na.idx[,1])],
				   er_fraction = ecoecon$er.prop[-unique(na.idx[,1])],
				   eumap[["layout"]][,pidx],
                   econ[,-1])
  econ <- merge(xy[-unique(na.idx[,1]),], econ, by="AFFGEOID")

write.csv(econ, file.path(getwd(), "results", "SubdivEconResults.csv"),
          row.names = FALSE)

subdiv$cluster <- insert.values(econ$cluster, NA, sort(unique(na.idx[,1])))
  subdiv$cluster[which(is.na(subdiv$cluster))] <- 0 
 
#***********************************************************
# Focal ecological regime assignment 
sdiv <- rasterize(vect(subdiv), er, field="cluster")

fclass <- function(x) {
  xp <- prop.table(table(factor(x)))
  as.numeric(names(xp)[which.max(xp)])
 }
er.class <- focal(er, matrix(1, 5, 5), fclass)

fprop <- function(x) {
  xp <- prop.table(table(factor(x)))
  xp[which.max(xp)]
 }
er.prop <- focal(er, matrix(1, 5, 5), fprop)

# erdiv <- c(sdiv, er.class, er.prop)
#   cls <- as.points(erdiv)
#   cls$ecoecon <- (cls$USGS_ELU_Global - 5000) * cls$cluster
# ecoecon <- rasterize(cls, er.class, "ecoecon")

ecoecon <- (er.class - 5000) * sdiv
  ecoecon <- c(ecoecon, sdiv, er.class, er.prop)
    names(ecoecon) <- c("ecoecon", "subdiv", "realm", "prealm")
    writeRaster(ecoecon, file.path(getwd(), "results",  
                "ecoecon.tif"), overwrite = TRUE)

#***********************************************************
# Calculate per-group PCA

group <- sort(unique(econ$cluster))
  contributions <- list()

pdf(file.path(getwd(), "results", "GroupPCA.pdf"), height=10, width=10)
  for(i in 1:length(group)) {
    d <- econ[econ$cluster == group[i],][c(1,10:29)]
      res.pca <- prcomp(d[,-1], scale = TRUE)
    plt <- fviz_eig(res.pca, main=paste0("Cluster ", i, " PCA variances"))
	  print(plt)
    plt <- fviz_pca_var(res.pca, col.var = "contrib",	
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE) + 
				 ggtitle(paste0("Cluster ", i, " PCA contribution"))
      print(plt)	  
    plt <- fviz_pca_biplot(res.pca, repel = TRUE,
                    col.var = "#2E9FDF", col.ind = "#696969") +
					ggtitle(paste0("Cluster ", i, " PCA biplot"))
      print(plt)
    ind.coord <- res.pca$x	
    center <- res.pca$center
    scale <- res.pca$scale
    getdistance <- function(ind_row, center, scale){
      return(sum(((ind_row-center)/scale)^2))
      }
    d2 <- apply(d[,-1], 1, getdistance, center, scale)
    cos2 <- function(ind.coord, d2){return(ind.coord^2/d2)}
    ind.cos2 <- apply(ind.coord, 2, cos2, d2)	
    contrib <- function(ind.coord, comp.sdev, n.ind){
      100*(1/n.ind)*ind.coord^2/comp.sdev^2
    }
    ind.contrib <- t(apply(ind.coord, 1, contrib, 
                     res.pca$sdev, nrow(ind.coord)))
    contributions[[i]] <- data.frame(AFFGEOID=d[,1], cluster=i, ind.contrib) 
  }
dev.off()
  
#***********************************************************
# PLot results

#### Plot by ecological realm
# 2D plot of first two embeddings color coded by realm, 
#   point size represents fraction of realm within subdiv
ggplot(econ, aes(x = X1, y = X2, color = realm)) + 
  geom_point(aes(size = er_fraction)) +
    scale_size_continuous(range = c(0.5, 6)) +
	  guides(size = "none", color=guide_legend(ncol=1)) +
	    scale_color_viridis(discrete = TRUE, option = "A") +
          theme_bw()

#### Plot by [X,Y]Z and volume of 4th dimension  
# 3D plot of first four embeddings where dim 1-3 are [X,Y]Z
#   and point color/size are volume of 4th embedding
rescale <- function(x) ( x - min(x) ) / ( max(x) - min(x) ) 
x <- rescale(eumap[["layout"]][,1])
y <- rescale(eumap[["layout"]][,2])
z <- rescale(eumap[["layout"]][,3])
v <- rescale(eumap[["layout"]][,4])
rbPal <- colorRampPalette(rev(c("red", "yellow", "green", "darkorchid4")))
vclr <- rbPal(30)[as.numeric(cut(v, breaks = 30))]
plot3d(x, y, z, type="s", size=v*2,5, col=vclr,
      xlab="", ylab="", zlab="") 

#### Plot Fuzzy C-means cluster
# 2D plot of first two embeddings color coded by cluster  
#   point size represents probability of cluster assignment 
ggplot(econ, aes(x = X1, y = X2, color = factor(cluster))) + 
  geom_point(aes(size = pvalue)) +
    scale_colour_brewer(palette = "Spectral") +
      scale_size_continuous(range = c(1, 6)) +  
	    guides(size = "none", color=guide_legend(ncol=1))

# plot clusters by subdivision	 
plot(subdiv["cluster"], pal=RColorBrewer::brewer.pal(k.opt+1,"Spectral"))
    
#### 3D Plot by [X,Y]Z and cluster
# 3D plot of first three embeddings where point colors
#   represent cluster and size pvalue of membership 
rescale <- function(x) ( x - min(x) ) / ( max(x) - min(x) ) 
x <- rescale(eumap[["layout"]][,1])
y <- rescale(eumap[["layout"]][,2])
z <- rescale(eumap[["layout"]][,3])
v <- econ$pvalue 
cls <- factor(econ$cluster)
  levels(cls) <- rev(RColorBrewer::brewer.pal(k.opt,"Spectral"))
    cls <- as.character(cls)	
plot3d(x, y, z, type="s", size=v*2,5, col=cls,
      xlab="", ylab="", zlab="") 
 
#***********************************************************
# Exploratory plots of raw parameter distributions
pcol = 10:44 # column indices for model parameters

pdf(file.path(getwd(), "results", "ParameterClusteredDistributions.pdf"), height=8, width=11)
  lapply(split(pcol, ceiling(seq_along(pcol)/4)), function(x) {
    plist <- vector(mode = "list", length = 4)
    i=0
    for(j in x) {
      i=i+1
      plist[[i]] <- 
        ggplot(data = econ, aes_string(x = names(econ)[j], fill = "cluster")) + 
          geom_histogram(aes(y =..density..), bins=8, alpha=0.7, position="dodge") + 
            geom_density(aes_string(x = names(econ)[j]), inherit.aes = FALSE) +
              scale_fill_brewer(palette = "BuPu") +			  
                geom_rug() +
                  labs(x = names(econ)[j]) +
                    theme_minimal()  
     }
    g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], 
                 nrow=2)  
      print(g) })    
dev.off()

pdf(file.path(getwd(), "results", "ParameterStackedDistributions.pdf"), height=8, width=11)
  lapply(split(pcol, ceiling(seq_along(pcol)/4)), function(x) {
    plist <- vector(mode = "list", length = 4)
    i=0
    for(j in x) {
      i=i+1
      plist[[i]] <- 
        ggplot(data = econ, aes_string(x = names(econ)[j], fill = "cluster")) + 
          geom_density(alpha=0.3) +
            scale_fill_brewer(palette = "BuPu") +			  
              geom_rug() +
                labs(x = names(econ)[j]) +
                  theme_minimal() 
     }
    g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], 
                 nrow=2)  
      print(g) })    
dev.off()

#***********************************************************
# Cluster aggregation statistics
clust.med <- aggregate(econ[,7:ncol(econ)], by=list(econ[,"cluster"]), median)
clust.mad <- aggregate(econ[,7:ncol(econ)], by=list(econ[,"cluster"]), mad)

# plot distribution of Gini_Index._delta
i = Gini_Index._delta
ggplot(data = clust.med, mapping = aes_string(x = i)) + 
  geom_histogram(aes(y=..density..), bins=k.opt, fill="blue", 
                 color="white", alpha=0.7) + 
  geom_density() +
  geom_rug() +
  labs(x=i) +
  theme_minimal()  
 