# ForestDrought
This R code runs the core analyses for the paper: 

Crockett ETH, Qingfeng G, Atkins JE, Sun G, Potter KM, Coztanza J,
Ollinger S, Woodall C, McNulty  S, Trettin C, Holgerson, J, and Xiao J.
Influences of structural and species diversity on forest resistance to drought. Ecology Letters.

(c) Erin Crockett, 2025
erin.crockett@unbc.ca

Files:
- Run_Analyses.R - provides code to run the core analyses of the paper based on the ForestDroughtData.csv dataset
- Functions_for_Analyses.R - provides helper functions necessary to run the main code
- Preprocessing.R - provides code to combine and clean the datasets listed below
- ForestDroughtData.csv - provides example data useful to run the analyses. These are not the real data; these are simulated data with similar properties to the real data, and thus may provide slightly different results, because as noted in the paper the forest plot coordinates and derivative products are protected for privacy reasons. Folks interested in the exact plot coordinates may contact the US Forest Service to apply for special permissions and go through the data access process.
- SpeciesList_Conifer.csv - provides information on which species found within the FIA plot data are conifers (used during Preprocessing.R script)

The code in this repository depends on several R packages that are listed in the code. It also draws from Jarret Byrnes' excellent work for lavaan spatial corrections.

We used the following publicly available data:  
- FIA plot data: https://apps.fs.usda.gov/fia/datamart/datamart.html (note this provides fuzzed coordinates rather than the actual coordinates used in the data analyses)
- Landsat-based maps of net primary productivity from Google Earth Engine (Robinson et al., 2018, doi:10.1002/rse2.74)
- SoilGrids database 1000m2 product: https://files.isric.org/soilgrids/latest/data_aggregated/
- EPA Ecoregions: https://www.epa.gov/eco-research/ecoregions-north-america 
- Canopy cover: https://data.fs.usda.gov/geodata/rastergateway/treecanopycover
- Climate variables from the Daymet database: https://doi.org/10.3334/ORNLDAAC/1852
