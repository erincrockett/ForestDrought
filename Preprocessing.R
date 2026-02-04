# This R code outlines the preprocessing steps to combine and clean data sources
# to generate the input dataset used for the paper: 

# Crockett ETH, Qingfeng G, Atkins JE, Sun G, Potter KM, Coztanza J,
# Ollinger S, Woodall C, McNulty  S, Trettin C, Holgerson, J, and Xiao J.
# Influences of structural and species diversity on forest resistance to drought.
# Ecology Letters.

# The R Code "Run_Analyses.R" provides the code to run the main statistical analyses.
# The .csv file provides example code to run these analyses. The code in this file will
# a dataset similar to the .csv provided (this code retains extra columns).

# (c) Erin Crockett, 2025
# erin.crockett@unbc.ca

#R version 4.3.1
library(raster)      #v3.6-23
library(rgdal)       #v1.6-7
library(reshape2)    #v1.4.4
library(spatialEco)  #v1.3-7
library(SPEI)        #v1.8.1
library(dplyr)       #v1.1.2 
library(rFIA)        #v1.0.0
library(vegan)       #v2.6-4
library(MASS)        #v7.3-60
library(foreign)     #v0.8-84
library(ggplot2)     #v3.4.2
library(viridis)     #v0.6.3



###[1] FOREST PLOT DATA ########################################################

###(1.1) Input Parameters ------------------------------------------------------
states <- c("MO","KS","NE","SD","ND","NM","AZ","NV",
            "CO","UT","WY","MT","ID","WA","OR","CA",
            "TN","TX","KY","VA","OK",
            "WV","MD","DE","NJ","PA","NY",
            "CT","RI","MA","VT","NH","ME",
            "OH","MI","IN","IL","WI","MN","IA",
            "SC","GA","FL","AL","MS","AR","LA","NC")          
states <- sort(states)

#Set the Base Directory
#>User to change their own
base.directory <- "F:/DataFiles"

#Conifer Species List based mainly from: https://www.srs.fs.usda.gov/pubs/misc/ag_654/volume_1/vol1_table_of_contents.htm
# and also from https://herbaria.plants.ox.ac.uk/bol/conifers/Explore
setwd(base.directory)
fia.species <- read.csv("SpeciesList_Conifer.csv")
fia.species$SpeciesName <- paste( fia.species$Genus, fia.species$Species)

#Conifer or Not - change to simply to 1/0  (other numbers based on which data source list came from)
fia.species[ is.na(fia.species$Conifer), "Conifer" ] <- 0
fia.species[ fia.species$Conifer > 1 , "Conifer"]    <- 1

height.NA.cutoffF <- 0.10  #What proportion of NAs to allow in height data before discard that plot (0.10 means 10%)

#Note how many FIA plots there are in each state
plots.by.state <- data.frame( State=states, nPlot=NA, nPlotOrg=NA)

#Source the structural properties/diversity functions
setwd("F:/DataFiles")
source("Functions_for_Analyses.R")


###(1.2) DOWNLOAD FIA DATA for given states -------------------------------------
#for the PLOT, CONDITION, TREE, and SUBPLOT tables
#note: slow download, large files, several GB required. Need to create folder for each state
options(timeout=3600)
for( ii in 1:length( states ) ){
   print(ii)
   setwd( paste( base.directory, states[ii] , sep="/" ) )
   getFIA(states = states[ii], dir = paste(base.directory, states[ii], sep="/"),
          tables= c("PLOT","COND","TREE","SUBPLOT") )
}


#Add ID for Each Unique Plot
fia.ec.plots <- as.data.frame( matrix(NA, nrow=0, ncol=2) )
colnames(fia.ec.plots) <- c("CN","EC_ID")
for( ii in 1:length( states ) ){
   setwd( paste0( "F:/FIA_Public_Data/", states[ii] ) )
   fia.plot <- read.csv( paste0( states[ii], "_PLOT.csv") )
   fia.plot$EC_ID <- paste( states[ii] , fia.plot$UNITCD , 
                            fia.plot$COUNTYCD ,fia.plot$PLOT , sep="_")
   fia.ec.plots <- rbind( fia.ec.plots, fia.plot[ ,c("CN","EC_ID")] )
}
setwd("F:/DataFiles")
save(fia.ec.plots, file="EC_addIDs.RData")  


###(1.3) CALCULATE METRICS -----------------------------------------------------

## Slow to run  (tree file is large to load)
for( ii in 1:length( states ) ){
   
   #Read in the Data
   setwd( paste( "F:/FIA_Public_Data", states[ii] , sep="/" ) )
   fia.plot <- read.csv( paste0( states[ii], "_PLOT.csv") )
   fia.cond <- read.csv( paste0( states[ii], "_COND.csv") )
   fia.tree <- read.csv( paste0( states[ii], "_TREE.csv") )
   
   ### STEP 1 - DATA CLEANING ##################################################
   ##[S1.1] PLOT TABLE -----------------------------------------------------------
   plots.by.state[ii,"nPlotOrg"] <- nrow(fia.plot)
   #Check Plot Stats Code (if sampled, then = 1)
   fia.plot <- fia.plot[ fia.plot$PLOT_STATUS_CD == 1, ]
   #Check for same plot types (plot code)  [DESIGNCD]
   #Note: Most are code "1"  (but WA,OR,CA - which have many 501,502, ...)
   #Some OK codes not listed here since not contained in the data
   fia.plot <- fia.plot[ fia.plot$DESIGNCD %in% c( 1, 311,312,313,314,315,316,328,
                                                   220,221, 501,502,503,504,505,506 ),  ]
   #Check for Quality Status/type of measurement  ==1 or ==7
   fia.plot <- fia.plot[ fia.plot$QA_STATUS %in% c(1,7),  ]
   #Check Sample Method  1==all physically done in the field (Good)
   fia.plot <- fia.plot[ fia.plot$SAMP_METHOD_CD == 1, ]
   #Extract Useful Columns
   fia.plot <- fia.plot[ ,c("CN","MEASYEAR","LAT","LON","ELEV",
                            "ECOSUBCD","DESIGNCD", "QA_STATUS" )]
   
   ##[S1.2] CONDITION TABLE ------------------------------------------------------
   #Extract based on the plot table above
   fia.cond <- fia.cond[ fia.cond$PLT_CN %in% fia.plot$CN, ]
   ## Determine the % of each plot on forested (1) land and other land conditions
   #Aggregate by PLT_CN and Condition number
   cond.agg <- aggregate( CONDPROP_UNADJ ~ PLT_CN + COND_STATUS_CD, data=fia.cond, FUN=sum  )
   #Switch from long to wide format
   cond.agg.wide <- stats::reshape( cond.agg, idvar="PLT_CN", timevar="COND_STATUS_CD", direction="wide" )
   #If not all conditions present, need to add them to the matrix
   if( ncol(cond.agg.wide) != 6 ){
      all.col.names <- c( "PLT_CN", paste0("CONDPROP_UNADJ.", 1:5) )
      cols.to.add <- setdiff( all.col.names, colnames(cond.agg.wide) )  #find which not present
      cond.agg.wide[ ,cols.to.add] <- NA                #add the columns
      cond.agg.wide <- cond.agg.wide[ ,all.col.names]   #reorder to match 
   }
   #Change colnames and switch NA to 0
   colnames(cond.agg.wide)[2:6] <- paste( "COND", 1:5, sep="_")
   cond.agg.wide[ is.na(cond.agg.wide) ] <- 0
   #If Field Stand Age Code is  NA,0,998,999 (not measured), use STDAGE column
   for(aa in 1:nrow(fia.cond) ){
      #Find if its one of these four values
      if( fia.cond[aa, "FLDAGE"] %in% c(NA,999,998,0) ){
         #Look to see if a non NA value for STDAGE
         #If so, use that value ,#If not, mark as NA
         if( !is.na(fia.cond[aa, "STDAGE"]) ){
            fia.cond[aa, "FLDAGE"] <- fia.cond[aa, "STDAGE"]
         }else{
            fia.cond[aa, "FLDAGE"] <- NA
         }
      }
   }
   #Extract Useful columns
   condition.columns <- c("PLT_CN","CONDPROP_UNADJ", "STDORGCD",
                          "FLDAGE","STDAGE", "LIVE_CANOPY_CVR_PCT", 
                          "DSTRBCD1", "DSTRBYR1", 
                          "SLOPE","ASPECT","BALIVE", "CARBON_DOWN_DEAD", 
                          "CARBON_LITTER", "CARBON_SOIL_ORG", "CARBON_STANDING_DEAD", 
                          "CARBON_UNDERSTORY_AG", "CARBON_UNDERSTORY_BG" )
   fia.cond <- fia.cond[  ,condition.columns]
   #Shorten column names
   colnames(fia.cond) <- gsub( "CARBON","C", colnames(fia.cond) )
   colnames(fia.cond) <- gsub( "STANDING","STND", colnames(fia.cond) )
   colnames(fia.cond) <- gsub( "UNDERSTORY","UNDSTY", colnames(fia.cond) )
   
   ##[S1.3] TREE TABLE -----------------------------------------------------------
   #Extract Out Individual Trees in Plots based on selected above for Plot and Condition Tables
   fia.tree <- fia.tree[ fia.tree$PLT_CN %in% fia.cond$PLT_CN, ]
   #Extract trees of at least 5" in size  (These largeish trees are sampled in four subplots)
   fia.tree <- fia.tree[ fia.tree$DIA >= 5.0, ]
   #Extract trees that are still living
   fia.tree <- fia.tree[ fia.tree$STATUSCD == 1, ]
   #Remove any trees of unknown/unlisted species
   fia.tree <- fia.tree[ !is.na( fia.tree$SPCD) , ]
   #For Western Areas -- Extract Just those within the subplot, rather than the macroplot
   #ie within 24ft of the subplot centre
   if( states[ii] %in% c("CA","OR","WA") ){
      fia.tree <- fia.tree[ fia.tree$DIST <= 24,  ]
   }
   #Extract Useful Columns
   tree.columns <- c("PLT_CN","PLOT","SUBP","TREE","SPCD",
                     "DIA","HT","ACTUALHT","HTCD", "DIACHECK",
                     "CARBON_AG","CARBON_BG")
   fia.tree <- fia.tree[ ,tree.columns]
   #Add a row for abundance - used in Height NAs and others later on
   fia.tree$Abundance <- 1   
   
   ##(S1.3.2) Determine Height NA values and any Plots to Remove -----------------
   #Source my function
   prop.NAsF <- height.na.search( fia.tree )
   #Only keep plots with no/few NA values
   prop.NAsF <- prop.NAsF[ prop.NAsF$Height.nas < height.NA.cutoffF,  ]
   fia.tree <- fia.tree[ fia.tree$PLT_CN %in% prop.NAsF$PLT_CN, ]
   ##(S1.3.3) Add 4 Letter Species Codes & Conifer -------------------------------
   ## Add Species Code (use 4 letter code rather than number since will add to column names)
   #Replace name on Tree df to match that of the Species df
   colnames(fia.tree)[ which( colnames(fia.tree) == "SPCD" ) ] <- "FIACode"
   species.df <- fia.species[ ,c("FIACode","SPECIESCode","Conifer", "SpeciesName", "Genus", "Species")]
   fia.tree <- merge( x=fia.tree, y=species.df, by="FIACode", all.x=TRUE)
   ##(S1.3.4) Extract only plots that have trees in all four subplots -----------------
   plot.ids <- fia.plot$CN
   subplot.df <- data.frame( PLT_CN = plot.ids, nSubplots=NA)
   #Find number subplots in each plot
   for( tt in 1:length(plot.ids) ){
      tree.tt <- fia.tree[ fia.tree$PLT_CN == plot.ids[tt], ]
      subplot.df[tt,2] <- length( unique( tree.tt$SUBP ) )
   }
   #Extract these as vector, and then from the tree dataframe
   four.subplot.ids <- subplot.df[ subplot.df$nSubplots == 4, "PLT_CN"]
   fia.tree <- fia.tree[ fia.tree$PLT_CN %in% four.subplot.ids, ]
   
   ##[S1.4] MATCH TABLES - keeping matching plot ids -----------------------------
   #Extract only those plots that meet all the criteria above - and extract for each table
   matching.plot.ids <- intersect( fia.plot$CN, intersect( fia.tree$PLT_CN, fia.cond$PLT_CN) )
   fia.tree <- fia.tree[ fia.tree$PLT_CN %in% matching.plot.ids, ]
   fia.cond <- fia.cond[ fia.cond$PLT_CN %in% matching.plot.ids, ]
   fia.plot <- fia.plot[ fia.plot$CN %in% matching.plot.ids, ]   
   
   ### STEP 2 - CALC METRICS ###################################################
   ##[S2.1] Structural Properties & Diversity ----------------------------------
   ## Number of Stems in the Plot
   nstems.df <- stats::aggregate(Abundance ~ PLT_CN, fia.tree, FUN=sum )  
   colnames(nstems.df)[2] <- "nStems"
   ## Proportion of Conifer Trees in the Plot
   #Note: if species Code is listed as NA -- then conifer.df will not have a row for that PLT_CN
   conifer.df <- stats::aggregate(Conifer ~ PLT_CN, fia.tree, FUN=mean )  
   conifer.df$ForestType <- NA
   for(rr in 1:nrow(conifer.df) ){
      if( conifer.df[rr,"Conifer"] >= 0.70 ){
         conifer.df[rr,"ForestType"] <- "Conifer"
      }else if( conifer.df[rr,"Conifer"] <= 0.30 ){
         conifer.df[rr,"ForestType"] <- "Broadleaf"
      }else{
         conifer.df[rr,"ForestType"] <- "Mixed"
      }
   }
   ## Stand Density Index
   sdi.df <- data.frame( PLT_CN = matching.plot.ids,
                         SDI.cm = NA,
                         SDI.usa = NA )
   for(tt in 1:length(matching.plot.ids) ){
      #Extract EC_ID
      trees.tt <- fia.tree[ which(fia.tree$PLT_CN %in% matching.plot.ids[tt]), ]
      n.trees.tt <- nrow(trees.tt)
      #Metric Units
      plot.area.ha <- 0.0672  #ha per plot
      nStem.per.ha <- n.trees.tt / plot.area.ha
      trees.tt$DIA_cm <- 2.54 * trees.tt$DIA 
      quad.sum.cm.tt <- sum( trees.tt$DIA_cm^2 )
      quadratic.mean.cm.tt <- sqrt( quad.sum.cm.tt / n.trees.tt )
      #USA/Imperial Units
      area.sqft.4 <- 4 * (pi * 24^2) #with 4 subplots
      plot.area.acre <- area.sqft.4 / 43560    #acre = 43560 sq ft
      nStem.perAcre <- n.trees.tt / plot.area.acre
      quad.sum.usa.tt <- sum( trees.tt$DIA^2 )
      quadratic.mean.usa.tt <- sqrt( quad.sum.usa.tt / n.trees.tt )
      ## Calculate SDI (stand density index)
      b.exp <- -1.605  #Exponent Value
      sdi.df[tt, "SDI.cm"]  <- nStem.per.ha  * (quadratic.mean.cm.tt / 25.4)^b.exp
      sdi.df[tt, "SDI.usa"] <- nStem.perAcre * (quadratic.mean.usa.tt / 10)^b.exp
   }

   ## Basal Area Properties  (area = sq inches)
   basal.df <- calc.basal.metrics( fia.tree )
   ## DBH Size Classes (in inches)
   dbh.classes.df.2  <- calc.size.classes( fia.tree, size.typeF="DBH", class.widthF=2 ) 
   all.dbh.data <- merge( dbh.classes.df.2, basal.df, by="PLT_CN" )

   ## Height Size Classes & Propertiees (in feet)
   height.df <- calc.height.metrics( fia.tree )
   height.classes.df  <- calc.size.classes( fia.tree, size.typeF="Height", class.widthF=5  )
   all.height.data <- merge( height.classes.df, height.df, by="PLT_CN" )

   ## Biomass
   biomass.df <- calc.biomass( fia.tree )
   
   ## Bind datasets together
   all.equal( nstems.df$PLT_CN, conifer.df$PLT_CN )
   all.equal( nstems.df$PLT_CN, sdi.df$PLT_CN )
   all.equal( nstems.df$PLT_CN, all.dbh.data$PLT_CN )
   all.equal( nstems.df$PLT_CN, all.height.data$PLT_CN )
   all.equal( nstems.df$PLT_CN, biomass.df$PLT_CN )
   nrow(nstems.df) == nrow(sdi.df )
   structure.df <- merge( nstems.df,    sdi.df,          by="PLT_CN" )
   structure.df <- merge( structure.df, all.dbh.data,    by="PLT_CN" )
   structure.df <- merge( structure.df, all.height.data, by="PLT_CN" )
   structure.df <- merge( structure.df, biomass.df,      by="PLT_CN" )  
   structure.df <- merge( structure.df, conifer.df,      by="PLT_CN" , all.x=TRUE) 
   structure.df <- merge( structure.df, cond.agg.wide,   by="PLT_CN" , all.x=TRUE, all.y=FALSE) 

   ##[S2.2] Species Diversity Metrics ---------------------------------------------------
   tree.comm <- fia.tree[ ,c("PLT_CN","SPECIESCode","Abundance")] #Extract only columns needed
   commMat <- commMatFunction( tree.comm )
   diversity.df <- data.frame( PLT_CN = as.numeric(rownames(commMat) ) )
   commMat.PA <- commMat
   commMat.PA[commMat.PA > 0] <- 1 #create presence-absence
   diversity.df$q0 <- rowSums(commMat.PA)					                     #species richness
   diversity.df$SpShannon <- vegan::diversity(x=commMat, index="shannon")
   diversity.df$Sp_Even <- diversity.df$SpShannon / log( diversity.df$q0 )
   
   ##[S2.3] Covariates ----------------------------------------------------------
   #Indicate which columns to calculate weighted average
   #(if numeric variable and have multiple conditions in the same plot)
   #Note: for stand origin code -- get partial values if some plots planted and others are natural
   w.avg.cols <- c( "FLDAGE", "STDORGCD", "LIVE_CANOPY_CVR_PCT",  "SLOPE","ASPECT",
                    "C_DOWN_DEAD", "C_LITTER", "C_SOIL_ORG", "C_STND_DEAD", 
                    "C_UNDSTY_AG", "C_UNDSTY_BG")
   #Note: not using "STDAGE" since have "FLDAGE"
   #Create (nearly) empty df to add rows to
   condition.df <- data.frame( PLT_CN = structure.df$PLT_CN )
   #Add new column to the growing dataframe
   for( cc in 1:length(w.avg.cols) ){
      add.col <- get.weighted.avg( fia.cond, w.avg.cols[cc] )
      condition.df <- merge( condition.df, add.col, by="PLT_CN", all.x=TRUE)
   }
   colnames(fia.plot)[1] <- "PLT_CN"
   colnames(fia.plot)[3:4] <- c("LATx","LONx")  #indicate that fuzzed coordinates
   
   ##[S2.4] MERGE & SAVE Files --------------------------------------------------
   structure.div.df <- merge( structure.df, diversity.df, by= "PLT_CN")
   covariates.df <- merge( fia.plot, condition.df, by="PLT_CN")
   setwd( base.directory )
   write.csv( structure.div.df, paste0( states[ii], "_Structure_Diversity.csv"), row.names=FALSE )
   write.csv( covariates.df,    paste0( states[ii], "_Covariates.csv"), row.names=FALSE  )
}#Close for ii loop

#Make into one file
setwd( base.directory )
structure.df <- read.csv( paste0( states[1], "_Structure_Diversity.csv") )
covariates.df <- read.csv(  paste0( states[1], "_Covariates.csv")  )
for( ii in 2:length(states) ){
   str.ii <- read.csv( paste0( states[ii], "_Structure_Diversity.csv") )
   cov.ii <- read.csv(  paste0( states[ii], "_Covariates.csv")  )
   structure.df  <- rbind(structure.df, str.ii)
   covariates.df <- rbind(covariates.df, cov.ii)
}
setwd("F:/DataFiles")
write.csv( structure.df,  file="Structure.csv",  row.names=FALSE )
write.csv( covariates.df, file="Covariates.csv", row.names=FALSE  )


## STEP 3 - Spatial Files -------------------------------------------------------
#Projection information for different datasets
nad83.prj <-    "+proj=longlat +datum=NAD83 +no_defs"
wgs84.prj <-    "+proj=longlat +datum=WGS84 +no_defs"
lamb.con.prj <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
soil.prj <-     "+proj=igh +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
albers.prj <-   "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

## Read in the real FIA Coordinates 
 #Note: file not provided to protect the privacy of landowners
setwd("F:/DataFiles")
fia.realC.df <- read.csv("ActualCoordinates.csv")   
colnames(fia.realC.df)[1] <- "PLT_CN"
fia.realC.df <- fia.realC.df[ complete.cases(fia.realC.df$ACTUAL_LAT), ]

##Bind together into df 
fia.full <- merge( fia.realC.df, fia.ec.plots[ ,c("PLT_CN","EC_ID")],
                   by="PLT_CN", all.x=TRUE, all.y=FALSE)
fia.full <- fia.full[ complete.cases( fia.full$EC_ID), ]
fia.df <- fia.full[ order( fia.full$MEASYEAR, decreasing=TRUE), ]
fia.df <- fia.df[ !duplicated( fia.df$EC_ID), ]
fia.df$State <- substr( fia.df$EC_ID, start=1, stop=2 )

##Make into Spatial Pts df and Project
   #FIA coords originally in NAD 83
fia.sp <- SpatialPointsDataFrame( coords= fia.df[ ,c("ACTUAL_LON","ACTUAL_LAT")], data=fia.df, 
                                  proj4string=CRS(nad83.prj) )
#Projection to WGS as per the Productivity dataset
fia.prj.wgs84 <- spTransform(fia.sp, CRS(wgs84.prj) )
#Projection to Lambert Conformal Conic as per the Climate Dataset
fia.prj.lambcon <- spTransform(fia.sp, CRS(lamb.con.prj) )

#Function to Extract Raster Values from Files in the Folder 
rasterStack <- function( setwdF ){
   setwd(setwdF)
   my.filesF <- list.files() 
   my.rasterF <- raster::stack(my.filesF)
   return(my.rasterF)
}


###[2] PRODUCTIVITY ###########################################################
## Download Productivity 
#Landsat-based maps of net primary productivity from Google Earth Engine (Robinson et al., 2018, doi:10.1002/rse2.74),  (Accessed 2022 May & Jun)
#> User needs to do separately - download to folder: "F:/Productivity"

#Function - to get Productivity values at each location
get.Productivity.3x3 <- function( raster.fileF, sptl.ptsF, dir.prodF="F:/Productivity" ){
   #ras.folderF: name folder to get stack of raster files
   #sptl.ptsF:   spatial points dataset with coordinates of points that want to extract
   #dir.prodF:   file folder for productivity files
   setwd( dir.prodF )
   rasterF <- raster( raster.fileF )
   #Add year of Productivity
   yearF <- substr( raster.fileF, start=nchar(raster.fileF)-7, stop= nchar(raster.fileF)-4 )
   #Extract Focal Cell Numbers
   ras.valsF <- raster::extract( rasterF, sptl.ptsF, cellnumbers=TRUE)
   ras.valsF <- as.data.frame(ras.valsF)
   colnames(ras.valsF)[2] <- paste0("Prod_", yearF )
   ras.valsF$EC_ID <- sptl.ptsF$EC_ID
   ras.valsF$PLT_CN <- sptl.ptsF$PLT_CN
   # Extract Adjacent Cell Numbers
   ras.adjF <- raster::adjacent( rasterF, cells= unlist(ras.valsF$cells), 
                                 pairs=TRUE, directions=8 )
   colnames(ras.adjF)[1] <- "cells"
   # Add the original (focal) cell to the bottom of this df
   add.rows <- matrix( c(ras.valsF$cells, ras.valsF$cells), ncol=2, byrow=FALSE )
   ras.adjF <- rbind( ras.adjF, add.rows)
   #Merge Cell Values with the EC_ID code
   ras.adj.dfF <- merge( ras.adjF, ras.valsF, by="cells", all.x=TRUE)
   #Extract Productivity Values for all cells
   ras.adj.dfF$ProdNearby <- rasterF[ unlist(ras.adj.dfF$to) ]
   #For year 2020, NA values listed as 0
   if( grepl( "2020", raster.fileF) ){
      ras.adj.dfF[ which(ras.adj.dfF$ProdNearby  == 0 ),"ProdNearby" ] <- NA
   }
   #Aggregate/Summarise by Group (EC_ID)
   sumry.dataF <- stats::aggregate( ras.adj.dfF$ProdNearby, by=list(EC_ID = ras.adj.dfF$EC_ID), FUN=mean , na.rm=TRUE)
   colnames(sumry.dataF)[2] <- paste0("Prod_3x3_", yearF)
   #Number of NA values
   na.dataF <- stats::aggregate( ras.adj.dfF$ProdNearby, by=list(EC_ID = ras.adj.dfF$EC_ID), FUN=function(x){ sum( is.na(x)) } )
   colnames(na.dataF)[2] <- paste0("nNA_", yearF)
   sumry.dataF <- merge( sumry.dataF, na.dataF, by="EC_ID")
   sumry.dataF <- merge( ras.valsF[ ,c("PLT_CN", paste0("Prod_", yearF ), "EC_ID")], sumry.dataF, by="EC_ID")
   return( sumry.dataF )
}

## Extract productivity Values for each year 
setwd("F:/Productivity")
raster.files <- list.files()
raster.files <- raster.files[ grep("tif", raster.files)]       #Select tif files
raster.files <- raster.files[ -grep(".ovr", raster.files)]     #Remove .ovr (pyramid files)
raster.files <- raster.files[ -grep("xml", raster.files)]      #Remove xml file
#Create Template  
ptm <- proc.time()
prod.df <- get.Productivity.3x3( raster.files[1] , fia.prj.wgs84 )
proc.time() - ptm
#Add to template during loop  (takes long time)
for( rr in 2:length(raster.files)){
   prod.rr <- get.Productivity.3x3( raster.files[rr], fia.prj.wgs84 )
   prod.rr$PLT_CN <- NULL
   prod.df <- merge( prod.df, prod.rr, by="EC_ID")
   print( paste( "Completed rr", rr))
}
setwd("F:/DataFiles")
save( prod.df, file= "Productivity_3x3.RData") 


####[3] SOIL ###################################################################
#Download Soil data, 1km version
#https://files.isric.org/soilgrids/latest/data_aggregated/  (Accessed 2022 May & Jun)
#> User needs to do separately - download bulk density (bdod) and nitrogen to folder: "F:/Soils"

#Function - to get Soil values at each location  
get.SoilData <- function( sptl.ptsF, soil.propertyF, dir.soilF=soil.dir ){
   #sptl.ptsF:  spatial points dataset with coordinates of points that want to extract
   #soil.propertyF: which soil aspect to consider:  eg "bdod"
   #dir.prodF: file folder for productivity files
   #Load file
   soil.rastersF <- rasterStack( paste( dir.soilF, soil.propertyF, sep="/" ) )
   if( proj4string(soil.rastersF) != proj4string(sptl.ptsF) ){
      sptl.ptsF <- spTransform(sptl.ptsF, CRS( proj4string(soil.rastersF) ) )
   }
   #Extract Soil values at each point
   soil.valsF <- raster::extract( soil.rastersF, sptl.ptsF )
   soil.valsF <- as.data.frame(soil.valsF)
   #Rename columns   #Note: list.files() sorted colnames from 0-6, so not in descending order
   depth.valuesF <- c("0_5","100_200","15_30","30_60","5_15","60_100")
   colnames(soil.valsF) <- paste( soil.propertyF, depth.valuesF, sep="_" )
   #Calculate Average Values (not weighted by depth distances)
   soil.valsF[ ,paste0(soil.propertyF, "_0_30_avg") ] <- rowMeans(soil.valsF[ ,c(1,3,5)])
   soil.valsF[ ,paste0(soil.propertyF, "_0_60_avg") ] <- rowMeans(soil.valsF[ ,c(1,3,4,5)])
   soil.valsF[ ,paste0(soil.propertyF, "_0_200_avg") ] <- rowMeans(soil.valsF)
   #Add ID Columns    
   soil.valsF$EC_ID  <- sptl.ptsF@data$EC_ID
   rownames(soil.valsF) <- sptl.ptsF@data$PLT_CN
   return(soil.valsF)
}
#Project- make same as soils data
fia.prj.soils <- spTransform(fia.sp, CRS(soil.prj) )

soil.dir <- "F:/SoilGrids"
setwd(soil.dir)
soil.bdod <- get.SoilData( fia.prj.soils, "bdod" , soil.dir )
soil.nitrogen <- get.SoilData( fia.prj.soils, "nitrogen" , soil.dir )
soil.df <- merge(soil.bdod, soil.nitrogen, by="EC_ID")
setwd("F:/DataFiles")
save(soil.df, file="Soil_Values.RData")


####[4] EPA ECOREGIONS ########################################################
#Download EPA Ecoregions file: https://www.epa.gov/eco-research/ecoregions-north-america  (Accessed 2021 Nov)
#> User needs to do separately and save to folder F:/EPA_Ecoregions
#Note: initially downloaded Level 3 file, then subsequently reduced detail to level 1 ecoregions

setwd("F:/EPA_Ecoregions")
ecoReg <- readOGR( dsn=".", "NA_CEC_Eco_Level3")
#Define Projection as Lambert 
ecoReg.prj <- spTransform(ecoReg, CRS(lamb.prj) )
#Replace NA with EPA in the column names
colnames(ecoReg.prj@data) <- gsub( "NA","EPA", colnames(ecoReg.prj@data) )
#Extract just the Code Numbers (so its not so large a file)
ecoReg.prj <- ecoReg.prj[ ,c("EPA_L3CODE","EPA_L2CODE")]
#Overlay points in Polygons
if( proj4string(ecoReg.prj) != proj4string(fia.prj.lambcon) ){
   print( "check EPA projection" )
}
fia.epa <- spatialEco::point.in.poly( fia.prj.lambcon, ecoReg.prj)
#Remove Coords And Make Non-Spatial
fia.epa <- fia.epa@data
#Save
setwd("F:/DataFiles")
save(fia.epa, file="FIA_EPA.RData")


####[5] CANOPY COVER ############################################################
#Need to download canopy cover data 
#https://data.fs.usda.gov/geodata/rastergateway/treecanopycover  (Accessed 2022 Aug)
#> User needs to do separately and save to folder F:/NLCD_CanopyCover

get.canopycover <- function( sptl.ptsF,  canopy.rasF ){
   #Load the raster file
   canopy.rasF <- raster( canopy.dirF )
   #Check Projection Info 
   if( proj4string(canopy.rasF) != proj4string(sptl.ptsF) ){
      print( "check projection - pts" )
      break
   }
   #Extract Canopy Cover values at each point
   canopy.valsF <- raster::extract( canopy.rasF, sptl.ptsF )
   canopy.valsF <- as.data.frame(canopy.valsF)
   #Add ID Columns    
   canopy.valsF$EC_ID  <- sptl.ptsF@data$EC_ID
   rownames(canopy.valsF) <- sptl.ptsF@data$PLT_CN
   return(canopy.valsF)
}

#Project- make same as canopy cover data
fia.prj.albers <- spTransform(fia.sp, CRS(albers.prj) )
setwd("F:/CanopyCover")
canopy.cover <- get.canopycover( fia.prj.albers, "nlcd_2011_treecanopy_2019_08_31.img")
colnames(canopy.cover)[1] <- "CanopyCov2011"
#Save to Folder:
setwd("F:/DataFiles")
save(canopy.cover, file="Canopy_Cover.RData")


####[6] CLIMATE ################################################################
#Download Monthly Climate Data from Daymet
#https://doi.org/10.3334/ORNLDAAC/1852 (version 4)  (Accessed 2022 May)
#> User needs to do this separately
#Save to folders 'Prcp' 'Tmax' 'Tmin' within folder "F:/Climate_Monthly"

#Function - Extract climate values for each point, and then transform from wide to long data
get.ClimateValues <- function( varF, sptl.ptsF , dir.climF="F:/Climate_Monthly"){
   #varF: "Prcp" "Tmin" or "Tmax:  names of climate variables on the folder
   #sptl.ptsF:  spatial points dataset with coordinates of points that want to extract
   #dir.climF: base directory for climate files
   #Load Climate Data
   clim.rasterF <- rasterStack( paste( dir.climF, varF, sep="/") )
   ## Get Climate Values at Each of the Points
   if( proj4string(clim.rasterF) != proj4string(sptl.ptsF) ){
       sptl.ptsF <- spTransform(sptl.ptsF, CRS( proj4string(clim.rasterF) ) )
       print("TRANSFORM")
   }
   #Extract Climate values at each point
   clim.valsF <- raster::extract( clim.rasterF, sptl.ptsF )
   #Dataframe with Shorter column names
   clim.valsF <- as.data.frame(clim.valsF)
   if( varF == "Prcp"){ colnames(clim.valsF) <- gsub( pattern="daymet_v4_prcp_monttl_na", "prcp", colnames(clim.valsF) )  }
   if( varF == "Tmax"){ colnames(clim.valsF) <- gsub( pattern="daymet_v4_tmax_monavg_na", "tmax", colnames(clim.valsF) )  }
   if( varF == "Tmin"){ colnames(clim.valsF) <- gsub( pattern="daymet_v4_tmin_monavg_na", "tmin", colnames(clim.valsF) )  }
   #Add ID Column
   clim.valsF$PLT_CN <- sptl.ptsF@data$PLT_CN
   clim.valsF$EC_ID <- sptl.ptsF@data$EC_ID
   return( clim.valsF )
}
## Extract Climate Data Each Point
clim.Type.all <- c("Prcp","Tmax","Tmin")
for(cc in 1:length(clim.Type.all) ){
   clim.Type <- clim.Type.all[cc]
   clim.df <- get.ClimateValues(clim.Type, fia.prj.lambcon) 
   #Save initial dataset
   setwd( paste0( "F:/DataFiles") )
   save( clim.df, file= paste0( clim.Type, ".RData" ) )
}  

## Calculate Average 20 year values
dir.create("F:/DataFiles/Avg20yr")
#Function - get average over past 20 years
clim.20yr <- function( clim.dataF ){
   clim.20yr.avgF <- data.frame( PLT_CN=clim.dataF$PLT_CN, EC_ID=clim.dataF$EC_ID )
   for( ff in 240:(ncol(clim.dataF)-2) ){  #start at 240=20yr*12mo, last 2 columns are IDs
      col.ff <- rowMeans( clim.dataF[  , (ff-240+1):ff ] )   
      new.row.name <- colnames( clim.dataF)[ff]
      clim.20yr.avgF[ , new.row.name] <- col.ff
   }
   return(clim.20yr.avgF)
}

setwd("F:/DataFiles")
load( "Prcp.RData" )
clim.prcp <- clim.20yr( clim.df )
rm(clim.df) 

load( "Tmax.RData" )
clim.tmax <- clim.20yr( clim.df )
rm(clim.df) 

load( "Tmin.RData" )
clim.tmin <- clim.20yr( clim.df )
rm(clim.df) 

setwd("F:/DataFiles")
save( clim.prcp , file= "Prcp_Avg.RData" )
save( clim.tmax , file= "Tmax_Avg.RData" )
save( clim.tmin , file= "Tmin_Avg.RData" )


####[7] CALCULATE SPEI #########################################################

## Function - Change Wide to Long Format and add Year and Month Columns --------
wide.to.longF <- function( dfF, varF ){
   dfF$EC_ID <- NULL
   long.climF <- reshape2::melt( dfF, id.vars="PLT_CN")
   var.name <- strsplit(as.character(long.climF$variable), split="_" )
   long.climF$Year <- NA
   long.climF$Month <- NA
   for(jj in 1:length(var.name)){
      name.jj <- var.name[[jj]]
      name.jj <- strsplit(name.jj, split="\\.")
      long.climF[jj,"Year"]  <- as.numeric( name.jj[[2]][1] )
      long.climF[jj,"Month"] <- as.numeric( name.jj[[3]][2] )
   }
   colnames(long.climF)[3] <- varF
   long.climF$variable <- NULL
   #Create ID to merge between Prcp/Tmax/Tmin datasets
   long.climF$MergeID <- paste( long.climF$PLT_CN, long.climF$Year, long.climF$Month, sep="_")
   return(long.climF)
}

## Function - SPEI drought index at various time steps 
calc.Drought <- function( site.dfF, latF, monthly.stepsF=c(3,6,12)  ){
   #site.dfF: df which has precip and temperature values for each time step of interest
   #latF: latitude of the site
   #monthly.stepsF:  which time steps to calculate for
   #Note: #PRCP - monthly precipitation total (mm)
          #TMED - monthly mean temperature (C)

   #Calc Potential Evapotranspiration
   PET_tho <- thornthwaite( Tave = site.dfF$Tavg, lat=latF ,  verbose=F)
   site.dfF$PET <- as.numeric( PET_tho[ ,1] )  #convert time series column into reg
   #Calc Water Balance 
   site.dfF$BAL <- site.dfF$Prcp - site.dfF$PET
   #Create dataframe to add do depending on how many intervals are desired
   sp.valuesF <- data.frame( Year = site.dfF$Year, Month = site.dfF$Month )
   #Calc SPEI for different Intervals (time step of data and here in months)
   for(ff in 1:length(monthly.stepsF)){
      spei.ff <- spei( site.dfF[ ,"BAL"], scale= monthly.stepsF[ff] , verbose=F )
      new.col.name <- paste0("SPEI_", monthly.stepsF[ff])
      sp.valuesF[ , new.col.name] <- as.numeric( spei.ff$fitted )
   }
   return( sp.valuesF)
}

#Read in Climate Data
setwd("F:/DataFiles")
load( "Prcp.RData" ) 
clim.prcp <- clim.df
load( "Tmax.RData" )
clim.tmax <- clim.df
load( "Tmin.RData" )
clim.tmin <- clim.df
rm(clim.df)

#Create empty list
drought.list <- list()

ptm <- proc.time() 
for(ii in 1:nrow(fia.df) ){
   #Extract latitude for the evapotranspiration metrics
   latitude.ii <- fia.df[ii,"ACTUAL_LAT"]
   #Extract mini matrix of climate values for that site ID
   prcp.ii <- clim.prcp[ which(clim.prcp$PLT_CN == fia.df[ii,"PLT_CN"]), ]
   tmax.ii <- clim.tmax[ which(clim.tmax$PLT_CN == fia.df[ii,"PLT_CN"]), ]
   tmin.ii <- clim.tmin[ which(clim.tmin$PLT_CN == fia.df[ii,"PLT_CN"]), ]
   prcp.ii.long <- wide.to.longF( prcp.ii, "Prcp" )
   tmax.ii.long <- wide.to.longF( tmax.ii, "Tmax" )
   tmin.ii.long <- wide.to.longF( tmin.ii, "Tmin" )
   climate.ii <- merge( prcp.ii.long, tmax.ii.long[ ,c("Tmax","MergeID")], by="MergeID" )
   climate.ii <- merge( climate.ii,   tmin.ii.long[ ,c("Tmin","MergeID")], by="MergeID" )
   #Create Average Temperature based on Tmax and Tmin
   climate.ii$Tavg <- (climate.ii$Tmax + climate.ii$Tmin) / 2
   #Order the data from beginning to end -- since the merge orders by ID
   climate.ii <- climate.ii[ order( climate.ii$Year, climate.ii$Month), ]
   if( nrow(climate.ii) != 504 ){  #how many months of climate data
      print("Check Climate.ii")
      break
   }
   #Calculate SPEI Drought index
   drought.list[[ii]] <- calc.Drought( climate.ii, latitude.ii, monthly.stepsF= c(3,6,9,12,15,18,21,24) )  #c(3,6,12) ) #c(1,2,3,5,6,8,12)  )
   rm( latitude.ii, prcp.ii, prcp.ii.long, tmax.ii, tmax.ii.long, tmin.ii, tmin.ii.long, climate.ii  )
   
   if( ii %in% seq( 100, 50000, 100) ){
      print( paste( "ii", ii, " -- of", nrow(fia.df)) )
   }
}#close ii
proc.time() - ptm  
   
#Add Plot IDs as the names & #Save File
names(drought.list) <- fia.df$PLT_CN
setwd( "F:/DataFiles" )
save(drought.list, file="SPEI.RData" ) 
   

####[8] SPEI-PRODUCTIVITY ######################################################

#Set Month to Extract SPEI data for Events
#> User to Repeat this section #8 for each month from May-Oct
month.num <- 8      # 5,6,7,8,9,10   (May-Oct)
months.df <- data.frame( MonthCh=c("May","June","July","Aug","Sep","Oct"),
                         MonthNum=5:10)
setwd("F:/DataFiles")
load( "Productivity_3x3.RData")   #object: prod.df

## Read in FIA Plot Information
setwd("F:/DataFiles")
str.df <- read.csv("Structure.csv")

#Combine - to just include those sites that met FIA criteria
fia.ec.plots <- merge( fia.ec.plots, str.df[ ,c("PLT_CN","DBH_Div_q0_2")], by="PLT_CN" )

## Function - Find Climate Events, that were also 'normal' before and after the designated period
search.for.events <- function( speiF , xYrsF = 1 , col.speiF= "SPEI_12" ,
                                 fia.yrsSampledF = fia.state.ii$MEASYEAR , yrs.beforeF=2 , yrs.afterF=0 , 
                                 fia.rangeF=3 , normal.boundsF=1  ){
   #speiF: df: with the spei values for that site for the month of inquiry (jun,jul,aug,sep)
   #xYrsF: length of drought event (consecutive yrs with SPEI value < threshold)
   #col.speiF: ch: which timestep of spei to use
   #fia.yrsSampledF: which years the FIA plot was sampled
   #yrs.beforeF: how many years before event to look for normal SPEI conditions
   #yrs.afterF: how many years before event to look for normal SPEI conditions
   #fia.rangeF: number of years before the event that can consider FIA measurements
   #normal.boundsF: what values define normal conditions (eg -1 and +1 spei values)

   #Create Event Vector to Fill - 1s indicate climate events (start with all 0s)
   event.vecF <- rep( 0, nrow(speiF) )
   max.fia.yearF <- 2020 - yrs.afterF - xYrsF        #eg if sampled in 2016, could have event 2017, and have values for 2018,2019,2020
   min.fia.yearF <- 1990 + yrs.beforeF - fia.rangeF  #eg productivity start 1990, so min FIA year is 1988 (if 5 years before event is OK)
   fia.yrsSampledF <- fia.yrsSampledF[ which( fia.yrsSampledF <= max.fia.yearF) ]
   fia.yrsSampledF <- fia.yrsSampledF[ which( fia.yrsSampledF >= min.fia.yearF) ]    #min year 1988 since start of the data
   if( length(fia.yrsSampledF) > 0 ){
      for( ff in 1:length(fia.yrsSampledF) ){
         #Extract SPEI values for the possible date range
         fia.yrF <- which( speiF[ ,"Year"] == fia.yrsSampledF[ff] )
         event.yr.rangeF <- (fia.yrF+1):(fia.yrF+fia.rangeF)
         #Make sure that falls within the possible years of SPEI data (if fia period long andn after period is short)
         event.yr.rangeF <- intersect( event.yr.rangeF, 1:(nrow(speiF)-yrs.afterF) )  #Cant go beyond range of the data
         #Extract the SPEI values of the event year
         spei.event.valsF <- speiF[ event.yr.rangeF ,col.speiF] 
    
         norm.vecF <- NULL 
         spei.event.valsF <- NULL
         spei.catF <- NULL
          for( yy in 1:length(event.yr.rangeF) ){
            year.yyF <- event.yr.rangeF[yy]
            vals.yrs.yyF <- speiF[ year.yyF:(year.yyF+xYrsF-1) , col.speiF]
            #If all years considered drought years  OR  All are normal range years        
            if( all( vals.yrs.yyF <= -1 ) || all( vals.yrs.yyF < 1 & vals.yrs.yyF > -1)  ){
               #Check normal before
               vals.bef.yyF <- speiF[ (year.yyF-yrs.beforeF) : (year.yyF-1) , col.speiF]
               if( any( is.na(vals.bef.yyF)) ){
                  norm.vecF[yy] <- "No"
                  spei.event.valsF[yy] <- 0
               }else{
                   if(  all( (vals.bef.yyF > -normal.boundsF) , (vals.bef.yyF < normal.boundsF ) ) ){
                     norm.vecF[yy] <- "Yes"
                     if( all(vals.yrs.yyF < -2 ) ){
                        spei.catF[yy] <- -2
                     }else if( all(vals.yrs.yyF < -1.75 ) ){
                        spei.catF[yy] <- -1.75
                     }else if( all(vals.yrs.yyF < -1.5 ) ){
                        spei.catF[yy] <- -1.5
                     }else if( all(vals.yrs.yyF < -1.25 ) ){
                        spei.catF[yy] <- -1.25
                     }else if( all(vals.yrs.yyF < -1.0 ) ){
                        spei.catF[yy] <- -1
                     }else{
                        spei.catF[yy] <- -0.01   #for mixed (normal range) plots
                     }
                     #Add actual values
                     spei.event.valsF[yy] <- paste( vals.yrs.yyF, collapse="," )
                   }else{
                           norm.vecF[yy] <- "No"
                           spei.event.valsF[yy] <- 0    #Add 0 to indicate that not interested in this one
                   }
               }
            }else{ #close if all <-1
               norm.vecF[yy] <- "No"
               spei.event.valsF[yy] <- 0
            }
         }#close for yy
         #If at least one period had normal conditions before/after
         if( any(norm.vecF == "Yes") ){
            max.yyF <-  which.min( spei.catF )
            #Add Event Year SPEI Category
            event.yrF <- event.yr.rangeF[ max.yyF ]
            event.vecF[event.yrF] <- min( spei.catF , na.rm=TRUE )           
         }
      }#close for ff
   }#close if have fia years sampled in the right time frame
   return(event.vecF)
}

#Create Empty Dataframe to fill during the loops
cnames <- c("Year","Month","Prod","Prod_DT", "Recent_PLT_CN","EC_ID","State","YrSampled")
results.df <- as.data.frame( matrix( NA, nrow=0, ncol= length(cnames)) )
colnames(results.df) <- cnames

r2.names <- c("EC_ID","State", "MeanProd","nNA","n0", paste0( "spei_", seq(3,24,3) ) )
r2.df <- as.data.frame( matrix( NA, nrow=0, ncol= length(r2.names)) )
colnames(r2.df) <- r2.names
r2.df.dt <- r2.df

states <- c("MO","KS","NE","SD","ND","NM","AZ","NV",
            "CO","UT","WY","MT","ID","WA","OR","CA",
            "TN","TX","KY","VA","OK",
            "WV","MD","DE","NJ","PA","NY",
            "CT","RI","MA","VT","NH","ME",
            "OH","MI","IN","IL","WI","MN","IA",
            "SC","GA","FL","AL","MS","AR","LA","NC")     

setwd( "F:/DataFiles" )
load("SPEI.RData" )  #object: drought.list

#Change Drought.list names
my.names <- names(drought.list)
n.times <- 0
for( dd in 1:length(my.names)){
   id <- fia.ec.plots[ fia.ec.plots$PLT_CN == my.names[dd], "EC_ID"]
   if( length(id) > 0 ){
      if( length( unique(id)) > 1 ){
         print("Check dd -- multiple entries")
         break
      }else{
         id <- id[1]
      }
      names(drought.list)[[dd]] <- id
   }else{
      n.times <- c(n.times, 1)
      names(drought.list)[[dd]] <- "none"
   }
   if( dd %in% seq(200,15000,200)){
      print( paste( dd, "of", length(my.names)))
   }
}

for( ss in 1:length(states)){
   #Extract Data for the state
   state.ss <- states[ss]
   fia.state <- fia.ec.plots[ fia.ec.plots$State == state.ss, ]
   fia.state <- fia.state[ fia.state$EC_ID %in% prod.df$EC_ID , ]
   prod.ss <- prod.df[ prod.df$EC_ID %in% fia.state$EC_ID, ]
   #Extract most recent year
   fia.ss <- fia.state[ order(fia.state$MEASYEAR), ]
   fia.ss <- fia.ss[ !duplicated( fia.ss$EC_ID), ]
   #Change Productivity df to have each year as column
   rownames(prod.ss) <- prod.ss$EC_ID
   prod.ss$EC_ID <- NULL
   prod.ss <- prod.ss[ ,sort( colnames(prod.ss))]

   for( ii in 1:nrow(prod.ss)){
      ec.id.ii <- rownames(prod.ss)[ii]
      #Productivity - Normalized
      prod.ii <- unlist( prod.ss[ii, ] )
      if( !is.na( tail(prod.ii,1)  ) ){
         if( tail(prod.ii,1) == 0 ){ prod.ii[ length(prod.ii) ] <- NA }  #Replace 0 with NA value (different since mosaiced this one)
      }
      prod_nNA.ii <-  prod.ii[ 1:31 ] 
      prod_center.ii <-  prod.ii[ 33:63 ]   
      prod_3x3.ii <-  prod.ii[ 64:94 ] 
      #Add NA values to Prod3x3 for two conditions
      center.na <- which( is.na(prod_center.ii) )  #if focal cell is NA
      if( length( center.na) >  0){ prod_3x3.ii[ center.na] <- NA }
      prod.nas <- which( prod_nNA.ii >= 3  )  #if at least 3 NAs
      if( length( prod.nas) >  0){ prod_3x3.ii[ prod.nas] <- NA }
      #Extract PLT_CN - and use to extract right dataframe from the SPEI list
      plt.cn.ii <- fia.ss[ fia.ss$EC_ID == ec.id.ii, "PLT_CN"]
      plt.cn.ii <- as.character(plt.cn.ii)
      #Note: some sites do not show up if did not meet the conditions
      #for including the plot in structural diversity calculations (eg not every subplot had trees) 
      if( ec.id.ii %in% names(drought.list) ){
         #Make sure Productivity has least 15 records  (a number of NA and 0s)
         if( length( which(prod_3x3.ii > 0)) > 15 ){  
            spei.ii <- drought.list[[ ec.id.ii ]] 
            #SPEI - Extract Correct Months & Years
            spei.ii <- spei.ii[ spei.ii$Month == month.num , ]
            spei.ii <- spei.ii[ spei.ii$Year >= 1990, ]  #Productivity data start at 1990
            spei.ii <- spei.ii[ spei.ii$Year <= 2020, ]  #Productivity ends at 2020
            #Watch for any "Inf" values in SPEI  >> replace with NA
            for( cc in 3:10){
               spei.ii[ ,cc][ is.infinite(spei.ii[ ,cc]) ] <- NA
            }
            #Get years Sampled
            fia.state.ii <- fia.state[which(fia.state$EC_ID == ec.id.ii), ]    
            fia.state.ii <- fia.state.ii[ order(fia.state.ii$MEASYEAR), ]
            #Search for Climate Events
            spei.ii$spei_3mo_3be_3af_3fia  <- search.for.events( spei.ii , col.speiF= "SPEI_3" ,  yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_6mo_3be_3af_3fia  <- search.for.events( spei.ii , col.speiF= "SPEI_6" ,  yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_9mo_3be_3af_3fia  <- search.for.events( spei.ii , col.speiF= "SPEI_9" ,  yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_12mo_3be_3af_3fia <- search.for.events( spei.ii , col.speiF= "SPEI_12" , yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_15mo_3be_3af_3fia <- search.for.events( spei.ii , col.speiF= "SPEI_15" , yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_18mo_3be_3af_3fia <- search.for.events( spei.ii , col.speiF= "SPEI_18" , yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_21mo_3be_3af_3fia <- search.for.events( spei.ii , col.speiF= "SPEI_21" , yrs.beforeF=3, fia.rangeF=3 )
            spei.ii$spei_24mo_3be_3af_3fia <- search.for.events( spei.ii , col.speiF= "SPEI_24" , yrs.beforeF=3, fia.rangeF=3 )
            #Add Productivity
            spei.ii$ProdCenter <- prod_center.ii  
            spei.ii$Prod_nNA <- prod_nNA.ii
            spei.ii$Prod_3x3 <- prod_3x3.ii 
            ## Add Detrended Productivity - Linear --- 3x3 productivity
            #Linear regression
            lm.ii <- lm( prod_3x3.ii ~ c(1990:2020) )
            #Check if significant -- if not, use regular values
            lm.ii.sum <- summary(lm.ii)
            if( lm.ii.sum$coefficients[2,4] < 0.05 ){  #Pvalue
               prod.ii.detrend <- residuals(lm.ii)
               #Add in NAs if needed (NA in productivity dataset)
               if( length(prod.ii.detrend) != nrow(spei.ii) ){
                  #Find out which years have NA values
                  prod.years <- paste0( "Prod_3x3_", c(1990:2020))
                  na.cols <- setdiff( prod.years, names(prod.ii.detrend) )
                  #Add NA values to the end of the vector
                  na.vals <- rep(NA, length(na.cols)) 
                  names(na.vals) <- na.cols
                  prod.ii.detrend <- c(prod.ii.detrend, na.vals)
                  #Sort into the right order (if have NAs in the middle)
                  prod.ii.detrend <- prod.ii.detrend[ sort( names(prod.ii.detrend)) ]
                  print( paste( "NA issue", ii, state.ss ) )
               }
               spei.ii$Prod_3x3_DT <- prod.ii.detrend
            }else{
               spei.ii$Prod_3x3_DT <- prod_3x3.ii
            }#close if significant
             #Add Identifiers
            spei.ii$Recent_PLT_CN <- as.numeric( plt.cn.ii )
            spei.ii$EC_ID <- ec.id.ii
            spei.ii$State <- state.ss
            ### Add Which Years Site was Sampled
            fia.state.ii <- fia.state[ fia.state$EC_ID == ec.id.ii,  ]
            colnames(fia.state.ii)[ which(colnames(fia.state.ii)=="MEASYEAR")] <- "Year"
            fia.state.ii$YrSampled <- 1
            spei.ii <- merge( spei.ii, fia.state.ii[ ,c("Year","YrSampled","PLT_CN")], by="Year", all.x=TRUE )
            spei.ii$YrSampled[ is.na( spei.ii$YrSampled) ] <- 0
            #Bind to growing dataframe
            results.df <- rbind( results.df, spei.ii )
            ### Examine Productivity - SPEI relationships across different timesteps
            get.R2 <- function( spei.colF,  prod.colF="Prod_3x3", spei.dfF=spei.ii ){
               lm.spF <- lm( spei.ii[ ,prod.colF] ~ spei.ii[, spei.colF] )
               return( summary(lm.spF)$r.squared )
            }
            nNA <- length( which( is.na(prod_3x3.ii)))
            n0 <- length( which( prod_3x3.ii == 0 ))
            
            r2.vals <- c( ec.id.ii, state.ss, mean(prod.ii, na.rm=TRUE) , nNA, n0, 
                          get.R2("SPEI_3") , get.R2("SPEI_6"),  get.R2("SPEI_9"),  get.R2("SPEI_12"),   
                          get.R2("SPEI_15"), get.R2("SPEI_18"), get.R2("SPEI_21"), get.R2("SPEI_24")    )
            r2.vals <- matrix( r2.vals, nrow=1)
            colnames(r2.vals) <- r2.names
            
            r2.vals.dt <- c( ec.id.ii, state.ss, mean(prod.ii, na.rm=TRUE) , nNA, n0, 
                             get.R2("SPEI_3",  "Prod_3x3_DT"), get.R2("SPEI_6",  "Prod_3x3_DT"), 
                             get.R2("SPEI_9",  "Prod_3x3_DT"), get.R2("SPEI_12", "Prod_3x3_DT"),   
                             get.R2("SPEI_15", "Prod_3x3_DT"), get.R2("SPEI_18", "Prod_3x3_DT"), 
                             get.R2("SPEI_21", "Prod_3x3_DT"), get.R2("SPEI_24", "Prod_3x3_DT")    )
            r2.vals.dt <- matrix( r2.vals.dt, nrow=1)
            colnames(r2.vals.dt) <- r2.names
            
            r2.df <- rbind( r2.df,   r2.vals )
            r2.df.dt <- rbind( r2.df.dt,  r2.vals.dt ) 
         }#close if at least 15 values
      }#close if found in SPEI names
      rm( ec.id.ii, prod.ii, plt.cn.ii, spei.ii, r2.vals, r2.vals.dt )
   }#close ii
   print( paste( state.ss, ss) )
}#close ss

#Change to numeric cols
for(cc in 3:13){
   r2.df[ ,cc] <- as.numeric( r2.df[ ,cc] )
   r2.df.dt[ ,cc] <- as.numeric( r2.df.dt[ ,cc] )
}

## Save Intermediary Files
mon.ch <- months.df[ which(months.df$MonthNum == month.num), "MonthCh"]
setwd("F:/DataFiles")
save( results.df, file= paste0("SPEI_events_", month.num, ".RData") )
save(r2.df, file= paste0("r2_vals_", month.num, ".RData") )
save(r2.df.dt, file= paste0("r2_detrended_vals_", month.num, ".RData") )


#Function: Create vector with productivity values before, after, during, for a SPEI columns
create.ProdSPEI.vec <- function( colF="spei_9mo_3be_3af_3fia" , spei.colF="SPEI_9" ,
                                              xYrsF=1 , nYrs.beforeF=3 , nYrs.afterF=0 , 
                                              yearColumnF= "YrSampled" , 
                                              iiF=nn , results.dfF=results.df , 
                                              df.newF=event.df , column.namesF=new.names ){
   if( results.dfF[ iiF , colF] != 0 ){
      ## Productivity Values
      #Average Productivity Before
      normal.beforeF <- mean( results.dfF[ (iiF-nYrs.beforeF):(iiF-1) , "Prod_3x3"] )
      #Average Productivity After
      normal.afterF <- mean( results.dfF[ (iiF+xYrsF):(iiF+nYrs.afterF+xYrsF-1) , "Prod_3x3"] )
      #Combine into productivity vector
      prod.vecF <- c(  results.dfF[iiF+xYrsF-1,"Prod_3x3"] ,  normal.beforeF, normal.afterF,
                       results.dfF[ iiF+xYrsF,"Prod_3x3"], results.dfF[ iiF+xYrsF+1,"Prod_3x3"], results.dfF[ iiF+xYrsF+1,"Prod_3x3"] )
      ## Check for FIA Sampling Before/After Drought
      fia.3yr.beforeF <- sum( results.dfF[ (iiF-3):(iiF-1) , yearColumnF] )     #Years Before   
      fia.5yr.beforeF <- sum( results.dfF[ (iiF-5):(iiF-1) , yearColumnF] )
      fia.3yr.afterF  <- sum( results.dfF[ (iiF+xYrsF):(iiF+xYrsF+2) , yearColumnF] )     #Years After
      fia.5yr.afterF  <- sum( results.dfF[ (iiF+xYrsF):(iiF+xYrsF+4) , yearColumnF] )
      fia.7yr.afterF  <- sum( results.dfF[ (iiF+xYrsF):(iiF+xYrsF+6) , yearColumnF] )
      fia.yr.eventF  <- sum( results.dfF[ iiF , yearColumnF] )                  #Year of Event
      if( is.na(fia.5yr.afterF ) ){ fia.5yr.afterF <- 0 }
      if( is.na(fia.7yr.afterF ) ){ fia.7yr.afterF <- 0 }
      if( fia.5yr.beforeF > 0 ){ 
         before.plt.cn <- results.dfF[ (iiF-5):(iiF-1) , "PLT_CN"] 
         before.plt.cn <- before.plt.cn[ which( !is.na( before.plt.cn) ) ]
               before.plt.cn <- before.plt.cn[ length(before.plt.cn)] #Take the most recent (last) one if there are multiple
      }else{
         before.plt.cn <- NA
      }
      if( fia.5yr.afterF > 0 ){ 
         after.plt.cn <- results.dfF[ (iiF+1):(iiF+5) , "PLT_CN"] 
         after.plt.cn <- after.plt.cn[ which( !is.na( after.plt.cn) ) ]
         after.plt.cn <- after.plt.cn[1] #Take the first one if there are multiple
      }else{
         after.plt.cn <- NA
      }
      if( fia.yr.eventF > 0 ){ 
         event.plt.cn <- results.dfF[ iiF , "PLT_CN"] 
         event.plt.cn <- event.plt.cn[ which( !is.na( event.plt.cn) ) ]
      }else{
         event.plt.cn <- NA
      }
      ## Create Values to Add to Growing datafame
      new.vecF <- c( prod.vecF  ,
                     #0,0, 0,0, 0,0 ,
                     fia.3yr.beforeF, fia.5yr.beforeF, fia.3yr.afterF, fia.5yr.afterF, fia.7yr.afterF, 
                     before.plt.cn, after.plt.cn, event.plt.cn )
      new.vecF <- as.data.frame( matrix(new.vecF, nrow=1 ) )
      new.vecF <- cbind( results.dfF[iiF, c("PLT_CN","EC_ID","Year", colF)],  
                         results.df[ iiF+xYrsF-1,  c(spei.colF) ] ,      #add +Xyrs-1 to get SPEI value of the 2nd (or 3rd) year of drought
                         colF, new.vecF )
      names(new.vecF) <- column.namesF
      df.newF <- rbind( df.newF, new.vecF)
   }
   return(df.newF)
}


##Load Results Df 
mon.ch <- months.df[ which(months.df$MonthNum == month.num), "MonthCh"]
setwd("F:/DataFiles") 
load( paste0("SPEI_events_", month.num, ".RData")  )  #object - results.df

## Create Productivity df with values before, during, and after the event
new.names <- c("Recent_PLT_CN", "EC_ID", "Year", "SPEIcat", "SPEIvalue", "Event", "Prod_3x3", "Prod_3x3Before","Prod_3x3After",
               "Prod_3x3_1a","Prod_3x3_2a","Prod_3x3_3a",
               "FIA_3yr_before","FIA_5yr_before","FIA_3yr_after","FIA_5yr_after", "FIA_7yr_after",
               "Before_PLT_CN","After_PLT_CN","Event_PLT_CN")
event.df <- as.data.frame( matrix( NA, nrow=0, ncol=length(new.names)))
#Params
event.time <- 1   #Event duration (eg 1 yr)
fia.numyears <- 3 #FIA measurements in X years before climate event
yrs.bef <- 3
yrs.aft <- 3

#Loop through all rows to extract values for the climate events
for( nn in 1:row(results.df) ){
   #For SPEI columns
   for(dd in seq(3,24,3)){
      colname.thisrun <- paste0("spei_",dd,"mo_",yrs.bef,"be_",yrs.aft,"af_",fia.numyears,"fia" )
      spei.col <- paste0( "SPEI_", dd)
      event.df <- create.ProdSPEI.vec( colname.thisrun,  spei.colF=spei.col, 
                        xYrsF=event.time, nYrs.beforeF=yrs.bef , nYrs.afterF=yrs.aft, df.newF=event.df ) 
   }#close dd
   if( nn %in% seq(10000, 100000000, 10000) ){
      print( paste( "nn", nn, "of", nrow(results.df) , "  -  ", Sys.time()) )
   }
}
## Save Values
setwd( "F:/DataFiles")
save(event.df,  file=paste0("SPEI_EventDf_", month.num, ".RData") )  


## Create Fig S3 - Heatmap Figure

months <- c("May","June","July","Aug","Sep","Oct")
month.num <- c(5,6,7,8,9,10)

## Create Template from First Month (May)
setwd("F:/DataFiles")
load( paste0( "r2_detrended_vals_", month.num[1] , ".RData") )  #r2.df.dt
#Calculate Mean & Median R2 Values
mat.avg <- apply( r2.df.dt[ ,3:13], MARGIN=2, FUN=mean)
mat.med <- apply( r2.df.dt[ ,3:13], MARGIN=2, FUN=median)
df.avg <- as.data.frame( matrix( mat.avg, nrow=1))
colnames(df.avg) <- names(mat.avg)
df.med <- as.data.frame( matrix( mat.med, nrow=1))
colnames(df.med) <- names(mat.med)
df.avg$Month <- months[1]
df.med$Month <- months[1]
df.avg$MonthNum <- month.num[1]
df.med$MonthNum <- month.num[1]
## Write for All
for(mm in 2:length(months)){
   setwd("F:/DataFiles")
   load( paste0( "r2_detrended_vals_", month.num[mm] , ".RData") )  #r2.df.dt
   #Calculate Mean & Median R2 Values
   mat.avg.mm <- apply( r2.df.dt[ ,3:13], MARGIN=2, FUN=mean)
   mat.med.mm <- apply( r2.df.dt[ ,3:13], MARGIN=2, FUN=median)
   df.avg.mm <- as.data.frame( matrix( mat.avg.mm, nrow=1))
   colnames(df.avg.mm) <- names(mat.avg.mm)
   df.med.mm <- as.data.frame( matrix( mat.med.mm, nrow=1))
   colnames(df.med.mm) <- names(mat.med.mm)
   df.avg.mm$Month <- months[mm]
   df.med.mm$Month <- months[mm]
   df.avg.mm$MonthNum <- month.num[mm]
   df.med.mm$MonthNum <- month.num[mm]
   #Bind to Previous
   df.avg <- rbind( df.avg, df.avg.mm)
   df.med <- rbind( df.med, df.med.mm)
}
create.heatmap <- function( df.R2F , monthsF=months ){
   #Create Long Data (from Wide) 
   R2F <- unlist(df.R2F[ 4:11])
   long.dfF <- as.data.frame( R2F )
   long.dfF$SPEI <- rep( seq(3,24,3)  , each=length(monthsF) )
   long.dfF$Month <- rep( monthsF  , length( seq(3,24,3) ) )
   long.dfF$SPEI <- factor( long.dfF$SPEI, levels=seq(3,24,3))
   long.dfF$Month <- factor( long.dfF$Month, levels=monthsF)
   long.dfF$R2 <- round( long.dfF$R2 , digits=3)
   #Plot
   my.plotF <- ggplot( long.dfF, aes(SPEI, Month, fill=R2)) +
      geom_tile( ) + 
      coord_fixed() +
      geom_text( aes(label=R2) , color="black", size=2 ) +
      scale_fill_gradientn( colors= mako(20 ))
   return(my.plotF)
}
setwd("F:/DataFiles")
png( filename="FigS3_Heatmap_ProdSPEI.png", height=14, width=14, units="cm", res=300)
   create.heatmap( df.med )
dev.off()



####[9] PREP STATS DATA ########################################################

###(9.1) Read in and Merge Data -----------------------------------------------
#Run this section of code for different drought thresholds:  (-1, -1.25, -1.5, -1.75, -2)
#and different SPEI months and timesteps  (not set up as a loop)
drought.thresh <- -2.00
drought.thresh.chr <- "-2.00"
#SPEI info
month.ch <- "Aug"   #spei month chr
month.num <- 8      #spei month num
spei.timestep <- 9  #spei timestep (# of months)
n.yr.drought <- 1   #how many years for the drought event
n.yr.before <- 3    #how many years of "normal" spei years beforehand
n.yr.fia.msr <- 3   #how many years before drought fia sampling must have occurred
##

## Drought Events df by month/timestep
setwd("F:/DataFiles")
load( paste0("SPEI_EventDf_", month.num, ".RData") )  #dim(event.df)
drought.event.col <- paste0( "spei_", spei.timestep, "mo_", n.yr.before,"be_","3af_", n.yr.fia.msr, "fia" )
spei.df.basic <- event.df[ which(event.df$Event == drought.event.col), ]

str.df <- read.csv("Structure.csv")
cov.df <- read.csv("Covariates.csv")
load("Prcp_Avg.RData")
load("Tmax_Avg.RData")
load("Tmin_Avg.RData")
load("FIA_EPA.RData")       #object: fia.epa
load("Soil_Values.RData")   #object: soil.df
load("Canopy_Cover.RData")

#Merge - covariate, structure, epa, soils, canopy cover
mod.data <- merge( cov.df,  str.df,        by="PLT_CN" )
mod.data <- merge( mod.data, fia.ec.plots[ ,c("PLT_CN","EC_ID")] , by="PLT_CN", all.x=T, all.y=F)
mod.data <- merge( mod.data, fia.epa,      by="EC_ID", all.x=T )
mod.data <- merge( mod.data, soil.df,      by="EC_ID")
mod.data <- merge( mod.data, canopy.cover, by="EC_ID")
colnames(mod.data)[2] <- "Before_PLT_CN"
mod.data$EC_ID <- NULL  #Remove before the merge (so dont have duplicates)
spei.df <- merge( spei.df.basic, mod.data, by="Before_PLT_CN" , all.x=TRUE, all.y=FALSE )

##Add EPA L1 Code
spei.df$EPA_L1CODE <- NA
for(ii in 1:nrow(spei.df)){
   temp.ii <- strsplit( as.character( spei.df[ii,"EPA_L2CODE"] ) , "\\." )
   spei.df[ii,"EPA_L1CODE"] <- temp.ii[[1]][1]
}
##Extract appropriate climate for each site  
spei.df$PrcpAvg <- NA
spei.df$TmaxAvg <- NA
spei.df$TminAvg <- NA
for(ii in 1:nrow(spei.df)){
   id.ii <- spei.df[ii,"EC_ID"]   
   event.yr.ii <- spei.df[ii,"Year"]   
   #Climate Data starts 1980, 20yr avg starts 2000
   if( event.yr.ii >= 2000 ){
      col.num.ii <- grep(event.yr.ii, colnames(clim.prcp) )[8]  #Using August = 8
      prcp.ii <- clim.prcp[ which( clim.prcp$EC_ID == id.ii) , col.num.ii]
      spei.df[ii,"PrcpAvg"] <- prcp.ii
      tmax.ii <- clim.tmax[ which( clim.tmax$EC_ID == id.ii) , col.num.ii]
      spei.df[ii,"TmaxAvg"] <- tmax.ii
      tmin.ii <- clim.tmin[ which( clim.tmin$EC_ID == id.ii) , col.num.ii]
      spei.df[ii,"TminAvg"] <- tmin.ii
   }
}
spei.df$Tavg <- ( spei.df$TmaxAvg + spei.df$TminAvg) / 2

###(9.2) Filtering & Cleaning ---------------------------------------------------
#Just Drought Years - at least at/beyond the SPEI threshold
spei.df <- spei.df[ which(spei.df$SPEIvalue <= drought.thresh), ]
#Productivity - Remove NA values 
spei.df <- spei.df[ complete.cases(spei.df$Prod_3x3), ]
spei.df <- spei.df[ complete.cases(spei.df$Prod_3x3Before), ]
#Canopy Cover - must be at least 10% (2011 year)
spei.df <- spei.df[ which(spei.df$CanopyCov2011 >= 10), ]
#Productivity - must be at least 500 "Before" and 100 "during event"
spei.df <- spei.df[ which(spei.df$Prod_3x3Before >= 500), ]
spei.df <- spei.df[ which(spei.df$Prod_3x3 >= 100), ]
#A few FLDAGE have 0 - remove these values
spei.df[ which(spei.df$FLDAGE == 0),    "FLDAGE"] <- NA
spei.df <- spei.df[ complete.cases(spei.df$FLDAGE), ]
#Extract most extreme event for plots that have multiple events
spei.df <- spei.df[ order(spei.df$SPEIvalue), ]
spei.df <- spei.df[ -which(duplicated(spei.df$EC_ID) ) , ]

#For FigS10 - keep multiple drought years
#These 3 lines commented out since these are supplementary, not the main analyses
#Instead of running the previous two lines, run the next three lines and then continue with the remainder of the code
  #id.dups.table <- table(spei.df$EC_ID)
  #id.dups <- names(id.dups.table[id.dups.table>1])
  #spei.df <- spei.df[ which(spei.df$EC_ID %in% id.dups), ] 

##Transformations & Calcs
spei.df$Resistance <- log( spei.df$Prod_3x3 / spei.df$Prod_3x3Before )
spei.df$SpRich_log <- log(spei.df$q0)
spei.df$FLDAGE_log <- log(spei.df$FLDAGE)
spei.df$nitrogen_0_200_avg_log <- log( spei.df$nitrogen_0_200_avg )
spei.df$SDI_log <- log(spei.df$SDI.cm)
spei.df$Biomass_log <- log(spei.df$Biomass_AG) 
spei.df$BasalArea_log <- log(spei.df$BasalArea)
spei.df$nStems_log <- log(spei.df$nStems)

##Scale the Variables
scale.vars <- c( "SpRich_log" , "q0", "SpShannon" , "Sp_Even" , "Resistance",
                 "FLDAGE_log" , "BasalArea_log", "Biomass_log", 
                 "nStems_log" , "Height_mean" , "SDI_log",
                 "PrcpAvg" , "Tavg", "bdod_0_200_avg", "nitrogen_0_200_avg_log" )
for(ii in 1:length(scale.vars)){
   new.vals.ii <- base::scale( spei.df[ , scale.vars[ii] ] )
   new.vals.ii <- as.numeric(new.vals.ii)                
   spei.df[ ,paste0(scale.vars[ii],"_z")] <- new.vals.ii
}

##Structural Diversity - as One Variable
pca.str.div <- princomp( spei.df[  ,c("Height_Div_q0_5","DBH_Div_q0_2")] )
summary(pca.str.div)
spei.df$StrRich <- pca.str.div$scores[ ,1]  
spei.df$StrRich_z <- base::scale( spei.df$StrRich, center=FALSE, scale=TRUE)

pca.str.div <- princomp( spei.df[  ,c("Height_Shannon_5","DBH_Shannon_2")] )
summary(pca.str.div)
spei.df$StrShannon <- pca.str.div$scores[ ,1]  
spei.df$StrShannon_z <- base::scale( spei.df$StrShannon, center=FALSE, scale=TRUE)

if( any( is.na(spei.df$Sp_Even) ) ){
   spei.df[ is.na(spei.df$Sp_Even) , "Sp_Even"] <- 0
}
if( any( is.na(spei.df$Height_Evenness_5) ) ){
   spei.df[ is.na(spei.df$Height_Evenness_5) , "Height_Evenness_5"] <- 0
}
if( any( is.na(spei.df$DBH_Evenness_2) ) ){
   spei.df[ is.na(spei.df$DBH_Evenness_2) , "DBH_Evenness_2"] <- 0
}
pca.str.div <- princomp( spei.df[  ,c("Height_Evenness_5","DBH_Evenness_2")] )
summary(pca.str.div)
spei.df$Str_Even <- pca.str.div$scores[ ,1]  
spei.df$Str_Even_z <- base::scale( spei.df$Str_Even, center=FALSE, scale=TRUE)

##Create Factor for SPEI category & L1Code
spei.df$SPEI_factor <- factor( as.character(spei.df$SPEIcat), levels=c("-1", "-1.25", "-1.5", "-1.75","-2"))
spei.df$EPA_L1CODE <- as.factor( spei.df$EPA_L1CODE)

##Save Dataset  
  #Repeat section 9 for other months & timesteps if desired & rename files
setwd( "F:/DataFiles" )
write.csv( spei.df, file="ForestDroughtData.csv")
write.dbf( spei.df[ ,c("Before_PLT_CN","LATx","LONx","Year","EPA_L1CODE")], file="ForestDroughtGIS_v3.dbf" )  #Format for ArcGIS Visualizations

## END
