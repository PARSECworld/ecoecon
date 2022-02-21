# ecoecon development 0.0-2

ecoecon R code for development of a Ecological-Economic spatial metric  

# Available code and functions in ecoecon 0.0-2 are:

		clusterEcon - Spatial clustering of economic data
        EcoEcon_NeighborEmbedding - model for neighbor embedding clustering method
		
# Available data in ecoecon 0.0-2 are:
	
Model development for Oregon (in data directory)

		ecoregion_data.csv - Tabular ecoregions data (Value is key field) 
		ecoregions.tif - USGS 2015 World Ecological Land Units (OR subset)
		OR_subdivisions.shp - Shapefile of Oregon census subdivisions 
		Oregon_CSD_Table.csv - Tabular census data

Preliminary results (in data/results directory)

		ecoecon.tif - Raster stack (ecoecon, subdivision, realm class, realm prop) 
                of combined cluster and ecoregion realms, encoding 
                is ((realm – 5000) * subdivision cluster 1:8) 	
		SubdivEconResults.csv - Economic variables and clustering results
		ParameterClusteredDistributions.pdf - plots of each variables histogram 
                colored by UMAP cluster with overlay PDF
		ParameterStackedDistributions.pdf - Overlaid probability density functions for 
                each variable and cluster.		
		
**Bugs**: Users are encouraged to report bugs here. Go to [issues](https://github.com/PARSECworld/ecoecon/issues) in the menu above, and press new issue to start a new bug report, documentation correction or feature request. You can direct questions to <jeffrey_evans@tnc.org>.

**This will eventually be an R package but, is currently not able to compile so please download manually.**
