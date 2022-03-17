#' kNN neighbor embedding of county-level economic data
suppressMessages(invisible(lapply(c("sf", "terra", "spatialEco", 
                 "umap", "ggplot2", "rgl", "e1071", "viridis", 
				 "factoextra", "cluster", "gridExtra"), 
				 require, character.only=TRUE)))

setwd("C:/evans/GITS/ecoecon")
  data.dir <- file.path(getwd(), "data")

set.seed(42)

mdl = c("all", "deltas")[2]

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
  names(econ.delta) <- paste0(names(econ.delta), "_delta")

if(mdl == "all") {
  econ <- data.frame(AFFGEOID = econ[,1], econ[,post], 
                     econ.delta, econ[,c(34:36)])
	cidx = 0 
} else if(mdl == "deltas") { 
  econ <- data.frame(AFFGEOID = econ[,1], econ.delta, 
                     econ[,c(34:36)]) 
	cidx = 0 					 
}

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
d = ncol(econ)-1
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
	ecoecon[[4]] <- mask(ecoecon[[4]], ecoecon[[1]])
    writeRaster(ecoecon, file.path(getwd(), "results",  
                "ecoecon.tif"), overwrite = TRUE)

#***********************************************************
# Plot results

pdf(file.path(getwd(), "results", "UMAP_Cluster_resutls.pdf"),
    height=8.5, width=11) 

#### Plot by ecological realm
# 2D plot of first two embeddings color coded by realm, 
#   point size represents fraction of realm within subdiv
p <- ggplot(econ, aes(x = X1, y = X2, color = realm)) + 
  geom_point(aes(size = er_fraction)) +
    scale_size_continuous(range = c(0.5, 6)) +
	  guides(size = "none", color=guide_legend(ncol=1)) +
	    scale_color_viridis(discrete = TRUE, option = "A") +
          theme_bw()
print(p)

#### Plot Fuzzy C-means cluster
# 2D plot of first two embeddings color coded by cluster  
#   point size represents probability of cluster assignment 
p <- ggplot(econ, aes(x = X1, y = X2, color = factor(cluster))) + 
  geom_point(aes(size = pvalue)) +
    scale_colour_brewer(palette = "Spectral") +
      scale_size_continuous(range = c(1, 6)) +  
	    guides(size = "none", color=guide_legend(ncol=1))
print(p) 

# plot clusters by subdivision
plot(st_geometry(subdiv))	 
  plot(subdiv["cluster"][-which(subdiv$cluster == 0),], 
       pal=RColorBrewer::brewer.pal(k.opt,"Spectral"), 
	   add=TRUE)

plot(ecoecon)

dev.off() 
 
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
    
# rgl::rgl.snapshot(file.path(getwd(), "results", "4D_embeddings.png"))	
	
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
 
# rgl::rgl.snapshot(file.path(getwd(), "results", "3D_embeddings_cluster_p.png"))	

#***********************************************************
# Calculate per-group PCA
#   loadings / eigenvectors are in the eigload object
#   feature contributions are in the contribution object
pca.out = file.path(getwd(), "results/pca") 

getdistance <- function(ind_row, center, scale){
  return(sum(((ind_row-center)/scale)^2))
}
cos2 <- function(ind.coord, d2){
  return(ind.coord^2/d2)
}
contrib <- function(ind.coord, comp.sdev, n.ind){
  100*(1/n.ind)*ind.coord^2/comp.sdev^2
}
	
group <- sort(unique(econ$cluster))
  contributions <- list()
  eigload <- list()

pdf(file.path(pca.out, "GroupPCA.pdf"), height=10, width=10)
  for(i in 1:length(group)) {
    d <- econ[econ$cluster == group[i],][c(1,10:ncol(econ))]
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

    d2 <- apply(d[,-1], 1, getdistance, center, scale)
    ind.cos2 <- apply(ind.coord, 2, cos2, d2)	
    ind.contrib <- t(apply(ind.coord, 1, contrib, 
                     res.pca$sdev, nrow(ind.coord)))
    contributions[[i]] <- data.frame(AFFGEOID=d[,1], cluster=i, ind.contrib) 
	eigload[[i]] <- data.frame(cluster=group[i], round(res.pca$rotation, 6))
  }
dev.off()

rn <- rownames(eigload[[1]])
  eigload <- do.call(plyr::rbind.fill, eigload) 
    eigload <- data.frame(parm=rn, eigload) 
      write.csv(eigload, file.path(pca.out, "eigenvectors.csv"), 
                row.names = FALSE)

contributions <- do.call(plyr::rbind.fill, contributions)   
  write.csv(contributions, file.path(pca.out, "FeatureContribution.csv"), 
            row.names = FALSE)
  
#***********************************************************
# Exploratory plots of raw parameter distributions
pcol = 10:ncol(econ) # column indices for model parameters
econ$cluster <- factor(econ$cluster)

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
              scale_color_viridis(discrete = TRUE, option = "A") +			  
                geom_rug() +
                  labs(x = names(econ)[j]) +
                    theme_minimal()  
     }
	if(length(x) == 1) { 
      g <- grid.arrange(plist[[1]], nrow=1)  
	} else if(length(x) == 2) {	
      g <- grid.arrange(plist[[1]], plist[[2]], nrow=1)  
	} else if(length(x) == 3) {	
      g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], nrow=2)  
	} else if(length(x) == 4) {
      g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], nrow=2)  
    }				 				 
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
            scale_color_viridis(discrete = TRUE, option = "A")  +			  
              geom_rug() +
                labs(x = names(econ)[j]) +
                  theme_minimal() 
     }
	if(length(x) == 1) { 
      g <- grid.arrange(plist[[1]], nrow=1)  
	} else if(length(x) == 2) {	
      g <- grid.arrange(plist[[1]], plist[[2]], nrow=1)  
	} else if(length(x) == 3) {	
      g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], nrow=2)  
	} else if(length(x) == 4) {
      g <- grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], nrow=2)  
    }				 				 
    print(g) })    
dev.off()

#***********************************************************
# Cluster aggregation statistics
clust.med <- aggregate(econ[,pcol], by=list(econ[,"cluster"]), median)
clust.mad <- aggregate(econ[,pcol], by=list(econ[,"cluster"]), mad)

clr <- RColorBrewer::brewer.pal(nrow(clust.med), "PRGn")

pdf(file.path(getwd(), "results", "AggregatedMedian.pdf"), 
    height=10, width=10)
  par(mfrow=c(2,2))
    hist(clust.med[,2], breaks=nrow(clust.med), freq=FALSE, col=clr,
         main=paste0("Median ", names(clust.med)[2]),
    	 xlab=names(clust.med)[2], ylab="pdf")
    	 lines(density(clust.med[,2]))
         legend("topright", legend=paste0("clust", 1:nrow(clust.med)), 
    	       fill=clr, cex=0.75)
  for(j in 3:ncol(clust.med)){
    hist(clust.med[,j], breaks=nrow(clust.med), freq=FALSE, col=clr,
        main=paste0("Median ", names(clust.med)[j]),
    	  xlab=names(clust.med)[j], ylab="pdf")
    	  lines(density(clust.med[,j]))
  }
dev.off()  

pdf(file.path(getwd(), "results", "AggregatedMAD.pdf"), 
    height=10, width=10)
  par(mfrow=c(2,2))
    hist(clust.mad[,2], breaks=nrow(clust.mad), freq=FALSE, col=clr,
         main=paste0("MAD ", names(clust.mad)[2]),
    	 xlab=names(clust.mad)[2], ylab="pdf")
    	 lines(density(clust.mad[,2]))
         legend("topright", legend=paste0("clust", 1:nrow(clust.mad)), 
    	       fill=clr, cex=0.75)
  for(j in 3:ncol(clust.mad)){
    hist(clust.mad[,j], breaks=nrow(clust.mad), freq=FALSE, col=clr,
        main=paste0("MAD ", names(clust.mad)[j]),
    	  xlab=names(clust.mad)[j], ylab="pdf")
    	  lines(density(clust.mad[,j]))
  }
dev.off()  
