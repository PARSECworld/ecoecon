#' ESDA of county-level economic data
invisible(lapply( c("rgeoda", "sf", "sp", "spatialEco", "spdep", 
          "RColorBrewer", "ggplot2", "dplyr", "EcoGenetics",
		  "rfUtilities"), library, character.only = TRUE))  

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

#*********************************
# Identify collinearity and multicollinearity 
# in data
d <- sf::st_drop_geometry(subdiv[,15:(ncol(subdiv)-1)])
( cl <- collinear(d, p=0.75) )
  
( cl.test <- rfUtilities::multi.collinear(d, perm = TRUE, leave.out = TRUE, n = 999) )
 
#*********************************
# Explore neighbor matrix (Wij)
subdiv.sp <- as(subdiv, "Spatial")
or.nb <- poly2nb(pl = subdiv.sp, queen = TRUE)
  summary(or.nb)
  
plot(subdiv.sp, border = "gray", main="1st order neighbors (Wij)")
  plot(or.nb, coordinates(subdiv.sp),
       col = "red", add = TRUE)

#***************************
# Smoothing using Spatial Moving Average
# (Weighted mean using Wij)
Wij <- spdep::nb2listw(spdep::poly2nb(pl = as(subdiv, "Spatial"), queen = TRUE))
subdiv <- dplyr::transmute(subdiv, Density_z = Gini_Index_2019 - mean(Gini_Index_2019), 
                           SMA_z = lag.listw(Wij, Density_z)) 
  subdiv$gini_sma <- lag.listw(x = Wij, subdiv.sp$Gini_Index_2019)

# Scatterplot of GINI 2019 and its spatial moving average
ggplot(data = subdiv, aes(x = Density_z, y = SMA_z)) +
  geom_point(color = "gray") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  coord_equal()

#*****************************************************
# Univariate LISA Local Autocorrelation
subdiv <- st_read(file.path(data.dir, "OR_subdivisions.shp"))
  subdiv <- st_transform(subdiv, st_crs("+proj=longlat +datum=WGS84 +no_defs"))
    subdiv <- merge(subdiv, econ, by="AFFGEOID")

subdiv.nb <- queen_weights(subdiv)

lisa <- local_moran(subdiv.nb, subdiv["Gini_Index_2019"])
  lisa.df <- data.frame(clust=0:6, label=lisa$labels, 
                        col=lisa$colors)  
    cat("False detection rate: ", lisa_fdr(lisa, 0.05), "\n")  
    lms <- lisa_values(gda_lisa = lisa)
    pvals <- lisa_pvalues(lisa)

subdiv$clust <- lisa_clusters(lisa, cutoff = 0.05) 
subdiv <- left_join(subdiv, lisa.df, by ="clust")
  lcols <- st_drop_geometry(subdiv[,"col"])[,1] 

#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) inset=c(-0.2,0)
plot(st_geometry(subdiv), col=lcols, main="GINI 2019 SMA LISA")
  legend("bottomright", legend=lisa.def$label, bg="white", 
         fill=lisa.def$col, horiz = FALSE, cex=0.75)

#*****************************************************
# Bivariate (Lapse Rate) LISA Local Autocorrelation
subdiv <- st_read(file.path(data.dir, "OR_subdivisions.shp"))
  subdiv <- st_transform(subdiv, st_crs("+proj=longlat +datum=WGS84 +no_defs"))
    subdiv <- merge(subdiv, econ, by="AFFGEOID")

subdiv.nb <- queen_weights(subdiv)
qsa <- local_bimoran(subdiv.nb, subdiv[c("Gini_Index_2010", "Gini_Index_2019")])
    lisa.df <- data.frame(clust=0:6, label=qsa$labels, 
                        col=qsa$colors)

subdiv$clust <- lisa_clusters(qsa, cutoff = 0.05) 
subdiv <- left_join(subdiv, lisa.df, by ="clust")
  lcols <- st_drop_geometry(subdiv[,"col"])[,1] 

plot(st_geometry(subdiv), col=lcols, main="GINI 2010 & 2019 Bivariate LISA")
  legend("bottomright", legend=lisa.def$label, bg="white", 
         fill=lisa.df$col, horiz = FALSE, cex=0.75)

#*****************************************************
# Moran's-I univariate and correlation correlogram 
or.nb <- poly2nb(pl = as(subdiv, "Spatial"), 
                 queen = TRUE)
  summary(or.nb)
  
( gini2019_cor <- sp.correlogram(or.nb, st_drop_geometry(subdiv[,"Gini_Index_2019"])[,1], 
                                 order=8, method="corr", zero.policy=TRUE) )
  plot(gini2019_cor, main="Lagged Correlation of GINI 2018")

( gini2019_acor <- sp.correlogram(or.nb, st_drop_geometry(subdiv[,"Gini_Index_2019"])[,1], 
                                 order=8, method="I", zero.policy=TRUE) )
  plot(gini2019_acor, main="Lagged Spatial Autocorrelation of GINI 2018")


