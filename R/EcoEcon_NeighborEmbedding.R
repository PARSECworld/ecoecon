#' Cluster county-level economic data
library(sf)
library(terra)
library(spatialEco)
library(RColorBrewer)
library(umap)
library(rsvd) 
library(plotly) 
library(ggplot2)

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
  econ <- st_drop_geometry(subdiv[,15:49])

econ.names <- names(econ)[-c(33:35)]
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
econ.delta <- do.call(cbind, delta)

econ <- data.frame(econ.delta, econ[,c(33:35)])
 
#*********************************
#   Index and remove NA's
#   Screen collinearity
na.idx <- which(is.na(econ), arr.ind = TRUE)
  cat(paste(names(econ)[unique(na.idx[,2])], "\n"),
      "have NA values", "\n")

if(nrow(na.idx) > 0){
  cat("removing rows:", unique(na.idx[,1]), "\n")
  econ <- econ[-unique(na.idx[,1]),]   
}

( cl <- collinear(na.omit(econ), p=0.75) )
  if(length(cl) > 0) 
    econ <- econ[,-which(names(econ) %in% cl[-c(6,7)])]

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

eumap <- umap(econ, config=custom.settings)

pidx <- c(1,2)
lab = ecoecon[,c("Realm_World_Ecosystem", "LF_ClassName")[1]]
  lab = gsub("Nearctic ", "", lab)


edat <- data.frame(eumap[["layout"]][,pidx], realm = lab[-unique(na.idx[,1])],
                   er_fraction = ecoecon$er.prop[-unique(na.idx[,1])]) 

ggplot(edat, aes(x = X1, y = X2, color = realm)) + 
  geom_point(aes(size = er_fraction)) +
    scale_size_continuous(range = c(0.5, 6)) +
	  guides(size = "none", color=guide_legend(ncol=1)) +
	    scale_color_viridis(discrete = TRUE, option = "A") +
          theme_bw()
	  

# # Create interactive plot
#   fig <- plot_ly(edat, x = ~X1, y = ~X2, split = ~label,
#                  type = 'scatter', mode = 'markers', colors=clr) %>%  
#   layout(  
#     plot_bgcolor = "#e5ecf6",
#     legend = list(title = list(text = "Ecoregions")),
#     xaxis = list(title = "0"),  
#     yaxis = list(title = "1")) 
# 
# fig
