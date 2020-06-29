#Coat-color-phenology-analysis

Species like hares, ptarmigans, Arctic fox and other turn white for winter. With less and less snow in warming climate, these species are mismatched and noticable to predators. Populations in which at least some animals retain summer colors have better chances to survive, especially if evolution works fast enough to propagate these 'summer color' genes. In this work, we outline hot spots for evolution to do its work. Protection of these spots is what needs to be done to help winter-white species to survive in a quickly changing environment. 

For details, see our paper: 'Winter Coat Color Polymorphisms Identify Global Hotspots for Evolutionary Rescue from Climate Change', Mills, Bragina, Kumar, Zimova, Lafferty, Feltner, Davis, Hacklander, Alves, Good, Melo-Ferreira, Dietz, Abramov, Lopatina, Fay.

The following files are included:

1)  "Data_final_2018.xlsx" 
This is the data used to build models predicting the probability of winter white for each species across their range. It is needed to run the R script. Column descriptions are as follows: 

OBJECTID_1 - unique object identifier
spsCode - five letter species code of specimen 
scient_name - scientific name of specimen
comName - common name of speciemn
hisSpsCode - historical species code of specimen
country - country where specimen was collected
locality - more specific information on location where specimen was collected
source - source of specimen 
observation - type of observation used to obtain specimen
observer - observer from current study who recored information on specimen
verbDate - date specimen was collected
year - year specimen was collected
colorSymbol - winter coat color of specimen; 1 = winter white; 2 = winter brown
femaleNum - number of females in the observation
maleNum - number of males in the observation
unkNum - number of individuals of unknown sex in the observation
decLat - latitude in decimal degrees of location where specimen was collected 
decLong- longitude in decimal degrees of location where specimen was collected 
numberObs - number of individuals in the observation
textExcerp - excerpt from source of observation if source was a written account 
snow - snow cover duration based on the Global SnowPack product of the German Aerospace Center (www.DLR.de/eoc) for all hydrological years (Sept. 1 – Aug. 31) between 2000/2001 and 2014/2015
bio_2 - mean diurnal temperature range (average monthly max temp - min temp) derived from Bioclim data set [http://www.worldclim.org]
bio_3 - isothermality in diurnal temperature (bio_2/temperature annual range) derived from Bioclim data set [http://www.worldclim.org]
bio_15 - precipitation seasonality (coefficient of variation) derived from Bioclim data set [http://www.worldclim.org]
alt - altitude in meters
dist_to_coast - distance to coast in meters
humanfp - human footprint indexed from Global Human Influence Index, v2 [1995 – 2004]: http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-indexgeographic/data-download


2)  "rasterstack4cov.tif"
This is the a raster stack of the three covariates (bio_2, bio_3 & snow) used in making the winter white probability maps. It is needed to run the R script. It can be downloaded from a repository at https://datadryad.org/stash/dataset/doi:10.5061/dryad.8m0p1.


3)  "global_phenology_final_script.R"
This is the R script associated with the analyses from this publication. Note in order to make species maps you need to download and read into R their range maps from the IUCN red list website: http://www.iucnredlist.org/  
