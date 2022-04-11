# Percent of bachelors degree or higher 
#   "pct_bachelors_degree_higher._delta"
# total retired
#   "RetireTotal"
# poverty status 
#   "Percent_Population_Poverty_Status_Determined._delta"   
# total population
#   "Total_Population._delta"
# GINI
# 
# median_income 
#
#
setwd("C:/evans/GITS/ecoecon")
  data.dir <- "C:/evans/GITS/ecoecon/data"
  out.dir <-  "C:/evans/GITS/ecoecon/results"
load(file.path(out.dir, "UMAP_all.RData"))   
econ.dat <- econ

suppressMessages(invisible(lapply(c("sf", "terra", "spatialEco", 
                 "rgeoda", "ggplot2", "viridis"), 
				 require, character.only=TRUE)))

#*******************************************
# Data prep
ecoecon <- rast(file.path(out.dir, "ecoecon.tif"))
  ecoecon <- ecoecon[[1]]
    ecoecon[ecoecon == 0] <- NA

# convert to vector polygons  
  ecoecon_poly <- st_as_sf(as.polygons(ecoecon, na.rm = TRUE))
    idx <- which(ecoecon_poly$ecoecon == 0)
      if(length(idx) > 0) ecoecon_poly <- ecoecon_poly[-idx,]
	#plot(st_geometry(ecoecon_poly))
  #st_write(ecoecon_poly, file.path(out.dir, "ecoecon_polys.shp"))
	
deltas <- data.frame(AFFGEOID=econ[,1], econ[,grep("delta", names(econ))])
  econ <- merge(subdiv[,1:11], deltas, by="AFFGEOID") 

gini <- rasterize(vect(econ), ecoecon, field="Gini_Index._delta") 
med_income <- rasterize(vect(econ), ecoecon, field="Median_Income_Households._delta") 
pct_ba <- rasterize(vect(econ), ecoecon, field="pct_bachelors_degree_higher._delta") 
pct_poverty <- rasterize(vect(econ), ecoecon, 
                 field="Percent_Population_Poverty_Status_Determined._delta") 

ecoecon_poly <- st_cast(st_read(file.path(out.dir, "ecoecon_polys.shp")), "POLYGON")
  ecoecon_poly$gini <- as.numeric(extract(gini, vect(ecoecon_poly), "median")[,2])
  ecoecon_poly$med_income <- as.numeric(extract(med_income, vect(ecoecon_poly), "median")[,2])
  ecoecon_poly$pct_ba <- as.numeric(extract(pct_ba, vect(ecoecon_poly), "median")[,2])
  ecoecon_poly$pct_poverty <- as.numeric(extract(pct_poverty, vect(ecoecon_poly), "median")[,2])

#*******************************************
# Wij matrix
wij <- queen_weights(ecoecon_poly)

#*******************************************
# Calculate univaraite LISA 
pdf(file.path(out.dir, "ecoecon_LISA.pdf"), height=10, width=10)

# Univaraite GINI delta
lisa <- local_moran(wij, ecoecon_poly['gini'])
  lisa_colors <- lisa_colors(lisa)
  lisa_labels <- lisa_labels(lisa)
  lisa_clusters <- lisa_clusters(lisa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Gini delta, Univaraite Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Univaraite Median Income delta
lisa <- local_moran(wij, ecoecon_poly['med_income'])
  lisa_colors <- lisa_colors(lisa)
  lisa_labels <- lisa_labels(lisa)
  lisa_clusters <- lisa_clusters(lisa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Median Income delta, Univaraite Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Univaraite percent bachelors degree delta
lisa <- local_moran(wij, ecoecon_poly['pct_ba'])
  lisa_colors <- lisa_colors(lisa)
  lisa_labels <- lisa_labels(lisa)
  lisa_clusters <- lisa_clusters(lisa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "percent bachelors degree delta, Univaraite Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Univaraite percent poverty delta
lisa <- local_moran(wij, ecoecon_poly['pct_poverty'])
  lisa_colors <- lisa_colors(lisa)
  lisa_labels <- lisa_labels(lisa)
  lisa_clusters <- lisa_clusters(lisa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Percent poverty delta, Univaraite Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

dev.off()


#  # Bivaraite GINI and Median Income deltas
#  qsa <- local_bimoran(wij, ecoecon_poly[c('gini', 'med_income')])
#    lisa_colors <- lisa_colors(qsa)
#    lisa_labels <- lisa_labels(qsa)
#    lisa_clusters <- lisa_clusters(qsa)
#    plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
#         function(x){return(lisa_colors[[x+1]])}), 
#         border = "#333333", lwd=0.2)
#    title(main = "Median Income and GINI Deltas, Bivaraite Local Moran")
#    legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
#           border = "#eeeeee", ncol = 4, cex=0.65)
#  
