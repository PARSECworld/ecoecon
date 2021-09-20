# ecoecon development 0.0-1

ecoecon R package for development of a Ecological-Economic spatial metric  

# Available functions in ecoecon 0.0-1 are:

		clusterEcon - Spatial clustering of economic data

# Available data in ecoecon 0.0-1 are:
	
Model development for Oregon (in data directory)

		ecoregion_data.csv - Tabular ecoregions data (Value is key field) 
		ecoregions.tif - USGS 2015 World Ecological Land Units (OR subset)
		OR_subdivisions.shp - Shapefile of Oregon census subdivisions 
		Oregon_CSD_Table.csv - Tabular census data

Preliminary results (in data/results directory)

		econclust.tif - Raster of clustered census data
		ecoecon.tif - Raster of combined cluster and ecoregion data, encoding 
						is cluster (1:n * 10000) + ecoregion value		

**Bugs**: Users are encouraged to report bugs here. Go to [issues](https://github.com/PARSECworld/ecoecon/issues) in the menu above, and press new issue to start a new bug report, documentation correction or feature request. You can direct questions to <jeffrey_evans@tnc.org>.

**or, for the development version, run the following (requires the remotes package):**
`remotes::install_github("PARSECworld/ecoecon")`
