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
subdiv.class <- read.csv(file.path(out.dir, "SubdivEconResults.csv"))
  econ <- merge(subdiv[,1:11], subdiv.class, by="AFFGEOID") 

ecoecon.all <- rast(file.path(out.dir, "ecoecon.tif"))
  ecoecon <- ecoecon.all[[1]]
    ecoecon[ecoecon == 0] <- NA

# convert to vector polygons  
ecoecon_poly <- st_as_sf(as.polygons(ecoecon, na.rm = TRUE))   
    idx <- which(ecoecon_poly$ecoecon == 0)
      if(length(idx) > 0) ecoecon_poly <- ecoecon_poly[-idx,]
	#plot(st_geometry(ecoecon_poly))
  #st_write(ecoecon_poly, file.path(out.dir, "ecoecon_polys.shp"))

ecoecon_poly <- st_cast(st_read(file.path(out.dir, "ecoecon_polys.shp")), "POLYGON")
  ecoecon_poly$cluster <- extract(ecoecon.all[[2]], vect(ecoecon_poly), "max")[,2]  
 
#*******************************************
# variables
# "Gini_Index_2010" and "Gini_Index_2019"  
# "Median_Income_Households_2012" and "Median_Income_Households_2019"  
# "Total_Population_2010" and "Total_Population_2019"   
# "pct_bachelors_degree_higher_2012" and "pct_bachelors_degree_higher_2019"  
# "Population_Poverty_Status_Determined_2012" and "Population_Poverty_Status_Determined_2019" 
# "RetireTotal"                                      
# "RetirementIncome"    

gini2010 <- rasterize(vect(subdiv), ecoecon, field="Gini_Index_2010")
gini2019 <- rasterize(vect(subdiv), ecoecon, field="Gini_Index_2019")
  ecoecon_poly$gini2010 <- as.numeric(extract(gini2010, vect(ecoecon_poly), 
                                      "median", na.rm=TRUE)[,2])
  ecoecon_poly$gini2019 <- as.numeric(extract(gini2019, vect(ecoecon_poly), 
                                      "median", na.rm=TRUE)[,2])

pct_ba2012 <- rasterize(vect(subdiv), ecoecon, field="pct_bachelors_degree_higher_2012") 
pct_ba2019 <- rasterize(vect(subdiv), ecoecon, field="pct_bachelors_degree_higher_2019") 
  ecoecon_poly$pct_ba2012 <- as.numeric(extract(pct_ba2012, vect(ecoecon_poly), 
                                        "median", na.rm=TRUE)[,2])
  ecoecon_poly$pct_ba2019 <- as.numeric(extract(pct_ba2019, vect(ecoecon_poly), 
                                        "median", na.rm=TRUE)[,2])

med_income2012 <- rasterize(vect(subdiv), ecoecon, field="Median_Income_Households_2012") 
med_income2019 <- rasterize(vect(subdiv), ecoecon, field="Median_Income_Households_2019") 
  ecoecon_poly$med_income2012 <- as.numeric(extract(med_income2012, vect(ecoecon_poly), 
                                            "median", na.rm=TRUE)[,2])
  ecoecon_poly$med_income2019 <- as.numeric(extract(med_income2019, vect(ecoecon_poly), 
                                            "median", na.rm=TRUE)[,2])

poverty2012 <- rasterize(vect(subdiv), ecoecon, field="Percent_Population_Poverty_Status_Determined_2012") 
poverty2019 <- rasterize(vect(subdiv), ecoecon, field="Percent_Population_Poverty_Status_Determined_2019") 
  ecoecon_poly$poverty2012 <- as.numeric(extract(poverty2012, vect(ecoecon_poly), 
                                         "median", na.rm=TRUE)[,2])
  ecoecon_poly$poverty2019 <- as.numeric(extract(poverty2019, vect(ecoecon_poly), 
                                         "median", na.rm=TRUE)[,2])

#*******************************************
# Wij matrix
wij <- queen_weights(ecoecon_poly)

#*******************************************
# Calculate biivaraite LISA 
pdf(file.path(out.dir, "ecoecon_bivariate_LISA.pdf"), height=10, width=10)

# Bivaraite GINI 2010 & 2019
qsa <- local_bimoran(wij, ecoecon_poly[c('gini2010', 'gini2019')])
  lisa_colors <- lisa_colors(qsa)
  lisa_labels <- lisa_labels(qsa)
  lisa_clusters <- lisa_clusters(qsa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "GINI 2010 & 2019 Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Bivaraite Median Income 2012 & 2019
qsa <- local_bimoran(wij, ecoecon_poly[c('med_income2012', 'med_income2019')])
  lisa_colors <- lisa_colors(qsa)
  lisa_labels <- lisa_labels(qsa)
  lisa_clusters <- lisa_clusters(qsa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Median Income 2012 & 2019 Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Bivaraite Percent Bachelor’s degree 2012 & 2019
qsa <- local_bimoran(wij, ecoecon_poly[c('pct_ba2012', 'pct_ba2019')])
  lisa_colors <- lisa_colors(qsa)
  lisa_labels <- lisa_labels(qsa)
  lisa_clusters <- lisa_clusters(qsa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Bachelor’s degree 2012 & 2019 Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

# Bivaraite Poverty 2012 & 2019
qsa <- local_bimoran(wij, ecoecon_poly[c('poverty2012', 'poverty2019')])
  lisa_colors <- lisa_colors(qsa)
  lisa_labels <- lisa_labels(qsa)
  lisa_clusters <- lisa_clusters(qsa)
  plot(st_geometry(ecoecon_poly), col=sapply(lisa_clusters, 
       function(x){return(lisa_colors[[x+1]])}), 
       border = "#333333", lwd=0.2)
  title(main = "Percent Poverty 2012 & 2019 Local Moran")
  legend('bottomleft', legend = lisa_labels, fill = lisa_colors, 
         border = "#eeeeee", ncol = 4, cex=0.65)

dev.off()

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
  title(main = "percent Bachelor’s degree delta, Univaraite Local Moran")
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


#### Deltas
# gini <- rasterize(vect(econ), ecoecon, field="Gini_Index._delta") 
# med_income <- rasterize(vect(econ), ecoecon, field="Median_Income_Households._delta") 
# pct_ba <- rasterize(vect(econ), ecoecon, field="pct_bachelors_degree_higher._delta") 
# pct_poverty <- rasterize(vect(econ), ecoecon, 
#                  field="Percent_Population_Poverty_Status_Determined._delta") 
#   ecoecon_poly$gini <- as.numeric(extract(gini, vect(ecoecon_poly), "median")[,2])
#   ecoecon_poly$med_income <- as.numeric(extract(med_income, vect(ecoecon_poly), "median")[,2])
#   ecoecon_poly$pct_ba <- as.numeric(extract(pct_ba, vect(ecoecon_poly), "median")[,2])
#   ecoecon_poly$pct_poverty <- as.numeric(extract(pct_poverty, vect(ecoecon_poly), "median")[,2])


