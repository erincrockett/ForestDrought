# This R code runs the core analyses for the paper: 

# Crockett ETH, Qingfeng G, Atkins JE, Sun G, Potter KM, Coztanza J,
# Ollinger S, Woodall C, McNulty  S, Trettin C, Holgerson, J, and Xiao J.
# Influences of structural and species diversity on forest resistance to drought.
# Ecology Letters.

# (c) Erin Crockett, 2025
# erin.crockett@unbc.ca

#R version 4.3.1
library(lavaan)   #v0.6-15
library(sf)       #v1.0-14
library(sp)       #v2.0-0
library(ape)      #v5.7-1

rm( list=ls() )


### (1) Set & Create File Directories ###################################################

#Set Your Own Base Folder Directory From Which to Read in the Function Files and Data
base.folder <- "..."

setwd( base.folder )
data.df.full <- read.csv( "ForestDroughtData.csv" )
source( "Functions_for_Analyses.R" )

## Create Folders to Save Files
base.folder.stats <- paste0(base.folder, "/StatsMods")
if( !dir.exists(base.folder.stats) ){ dir.create(base.folder.stats) }
base.folder.figs <- paste0(base.folder, "/Figs")
if( !dir.exists(base.folder.figs) ){ dir.create(base.folder.figs) }
dd.fol <- c("T-2.00", "T-1.00" )
mm.fol <- c("m11_basic_z", "m12_basic_z",  "m13_basic_z", 
            "m15_basic_z", "m17_basic_z", "m19_basic_z")
ss.fol <- c("SEM_reg", "SEM_SptlCoeffs", "SEM_EPA_reg", "SEM_EPA_SptlCoeffs",
            "SEM_ForestType_reg","SEM_ForestType_SptlCoeffs",
            "SEM_MultiThresh", "SEM_MultiThresh_SptlCoeffs")
for(dd in 1:length(dd.fol) ){
   dir.x1 <- paste0( base.folder.stats, "/", dd.fol[dd] )
   if( !dir.exists(dir.x1) ){ dir.create(dir.x1) }
   for(mm in 1:length(mm.fol) ){
      dir.x2 <- paste0( dir.x1, "/", mm.fol[mm] )
      if( !dir.exists(dir.x2) ){ dir.create(dir.x2) }
      for(ss in 1:length(ss.fol) ){
         dir.x3 <- paste0( dir.x2, "/", ss.fol[ss] )
         if( !dir.exists(dir.x3) ){ dir.create(dir.x3) }
      }#ss
   }#mm
}#dd


### (2) Run The Stats Analyses ####################################################

#Drought Threshold
full.drought.name <- "T-2.00"
data.df <- data.df.full[ which(data.df.full$SPEI_factor == "-2"), ]
table(data.df$EPA_L1CODE)

#N Pts Per Group (eg EPA region)
n.per.group <- 20

## Make into Spatial Pts df and Project
nad83.prj <- "+proj=longlat +datum=NAD83 +no_defs"
lamb.prj  <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
mod.sp <- SpatialPointsDataFrame( coords= data.df[ ,c("LONx","LATx")], data=data.df, 
                                  proj4string=CRS(nad83.prj) )
mod.prj <- spTransform(mod.sp, CRS(lamb.prj))

my.y.var <- "Resistance_z"

#Base Mod formulas - from which to add Str and Sp variables
m11.formula <- c(
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z" ,
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z" ,
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z + nitrogen_0_200_avg_log_z" )

m12.formula <- c(
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z" ,
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z" ,
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + bdod_0_200_avg_z + nitrogen_0_200_avg_log_z" ,
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z" ,   
   "FLDAGE_log_z + Tavg_z + PrcpAvg_z + nitrogen_0_200_avg_log_z" )


### [m11] Model a: StrRich + SpRich ---------------------------------------------------------------
m11.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m11_basic_z")
run.SEM( model.nameF="m11.basic" , model.formulaF=m11.formula , 
          str.varF = "StrRich_z"  , sp.varF = "SpRich_log_z"  ,
          base.folderF=m11.folder , groupF="EPA_L1CODE"  )
write.mod.formula( paste0(m11.folder, "/SEM_reg") , newNameF="m11.MultiT")
run.SEM( model.nameF="m11.basic" , model.formulaF=m11.formula , 
          str.varF = "StrRich_z"  , sp.varF = "SpRich_log_z"  ,
          base.folderF=m11.folder ,  groupF="ForestType"  )


### [m15] Model c: StrRich + StrEven ---------------------------------------------------------------
m15.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m15_basic_z")
run.SEM( model.nameF="m15.basic" , model.formulaF=m11.formula , 
          str.varF = "StrRich_z"  , sp.varF = "Str_Even_z"  ,
          base.folderF=m15.folder , groupF="EPA_L1CODE" )
write.mod.formula( paste0(m15.folder, "/SEM_reg") , newNameF="m15.MultiT")
run.SEM( model.nameF="m15.basic" , model.formulaF=m11.formula , 
          str.varF = "StrRich_z"  , sp.varF = "Str_Even_z"  ,
          base.folderF=m15.folder ,  groupF="ForestType" )


### [m17] Model d: SpRich + SpEven ---------------------------------------------------------------
m17.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m17_basic_z")
run.SEM( model.nameF="m17.basic" , model.formulaF=m11.formula , 
          str.varF = "SpRich_log_z"  , sp.varF = "Sp_Even_z"  ,
          base.folderF=m17.folder ,  groupF="EPA_L1CODE" )
write.mod.formula( paste0(m17.folder, "/SEM_reg") , newNameF="m17.MultiT")
run.SEM( model.nameF="m17.basic" , model.formulaF=m11.formula , 
          str.varF = "SpRich_log_z"  , sp.varF = "Sp_Even_z"  ,
          base.folderF=m17.folder ,  groupF="ForestType"  )


### [m19] Model e: StrEven + SpEven ---------------------------------------------------------------
m19.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m19_basic_z")
run.SEM( model.nameF="m19.basic" , model.formulaF=m11.formula , 
          str.varF = "Str_Even_z"  , sp.varF = "Sp_Even_z"  ,
          base.folderF=m19.folder ,  groupF="EPA_L1CODE"  )
write.mod.formula( paste0(m19.folder, "/SEM_reg") , newNameF="m19.MultiT")
run.SEM( model.nameF="m19.basic" , model.formulaF=m11.formula , 
          str.varF = "Str_Even_z"  , sp.varF = "Sp_Even_z"  ,
          base.folderF=m19.folder ,   groupF="ForestType" )


###[m12] Model f: StrRich + SpRich + StrEven + SpEven -----------------------------------------
m12.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m12_basic_z")
run.SEM.wEven( model.nameF="m12.RichEven" , model.formulaF=m12.formula , 
               str.varF = "StrRich_z"  , sp.varF = "SpRich_log_z"  ,
               str.evenF = "Str_Even_z" , sp.evenF = "Sp_Even_z" ,
               restrictPathsF="ForceCovary" ,
               base.folderF=m12.folder ,  groupF="EPA_L1CODE" )
write.mod.formula( paste0(m12.folder, "/SEM_reg") , newNameF="m12.MultiT")
run.SEM.wEven( model.nameF="m12.RichEven" , model.formulaF=m12.formula , 
                str.varF = "StrRich_z"  , sp.varF = "SpRich_log_z"  ,
               str.evenF = "Str_Even_z" , sp.evenF = "Sp_Even_z" ,
               restrictPathsF="ForceCovary" ,
               base.folderF=m12.folder , groupF="ForestType" )


###[m13] Model a2: StrShannon + SpShannon --------------------------------------------
m13.folder <- paste0( base.folder.stats,"/", full.drought.name,  "/m13_basic_z")
run.SEM( model.nameF="m13.basic" , model.formulaF=m11.formula , 
         str.varF = "StrShannon_z"  , sp.varF = "SpShannon_z"  ,
         base.folderF=m13.folder , groupF="EPA_L1CODE" )
write.mod.formula( paste0(m13.folder, "/SEM_reg") , newNameF="m13.MultiT")
run.SEM( model.nameF="m13.basic" , model.formulaF=m11.formula , 
         str.varF = "StrShannon_z"  , sp.varF = "SpShannon_z"  ,
         base.folderF=m13.folder ,  groupF="ForestType" )


##[m12 Multi SPEI] With Multiple SPEI Coefficients ----------------------------------------
setwd(base.folder.stats)
m12.MultiT <- readLines("m12.MultiT.txt")
d.name <- "T-1.00"

## Make into Spatial Pts df and Project
nad83.prj <- "+proj=longlat +datum=NAD83 +no_defs"
lamb.prj  <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
mod.sp.full <- SpatialPointsDataFrame( coords= data.df.full[ ,c("LONx","LATx")], data=data.df.full, 
                                  proj4string=CRS(nad83.prj) )
mod.prj.full <- spTransform(mod.sp.full, CRS(lamb.prj))

m12extra.folder <- paste0( base.folder.stats,"/", d.name,  "/m12_basic_z")
SEM.MultiThresh( model.nameF="m12.MultiT" , model.formulaF=m12.MultiT , 
               str.varF = "StrRich_z"  , sp.varF = "SpRich_log_z"  ,
               base.folderF=m12extra.folder ,  dataF=data.df.full ,
               groupF="SPEI_factor", sptl.ptsF=mod.prj.full )


##[m12 MultiReg] Multiple Regression ----------------------------------------------------
m12.multi.folder <- paste0(m12.folder, "/MultiReg")
if( !dir.exists(m12.multi.folder) ){ dir.create(m12.multi.folder) }

explan.vars <- c("StrRich_z", "SpRich_log_z" , "Str_Even_z", "Sp_Even_z",
                 "FLDAGE_log_z", "Tavg_z" , "PrcpAvg_z", 
                 "bdod_0_200_avg_z" , "nitrogen_0_200_avg_log_z" )
data.multi <- data.df[ , c("Resistance_z", explan.vars) ]

my.lm <- lm( Resistance_z ~ . , data=data.multi )
my.lm.sumry <- summary(my.lm)              
my.lm.sel <- step(my.lm, direction = "backward")
my.lm.sel.sumry <- summary(my.lm.sel)

setwd(m12.multi.folder)
write.csv(my.lm.sumry$coefficients,     file= "MultiReg_coeffs.csv" )
write.csv(my.lm.sel.sumry$coefficients, file= "MultiReg_coeffs_BWsel.csv" )


### (3) Create Figures ###############################################################

m.folder <- "m12_basic_z" 
m.name <-    "m12" 
str.metric <- "StrRich_z" 
sp.metric  <- "SpRich_log_z" 
str.even.metric <- "Str_Even_z"
sp.even.metric <- "Sp_Even_z"

#Figure Options
title.size <- 1.5
axis.size <- 1.15
mtext.size <- 0.9
label.size <- 1.0
side.spacing <- 0.1

##(3.1) Fig1 - Path Analysis Diagram 
#Examine coefficients within the final file of the following folder extract values. Diagram created in powerpoint.
#setwd( paste( base.folder.stats, full.drought.name , m.folder, "SEM_SptlCoeffs", sep="/") )                        


##(3.2) Fig2 - Coeffs by EPA Ecoregion ------------------------------------------------
epa.df <- data.frame(EPA_L1CODE = c(7,6,10,9,5,8) ,
                     EPA_Name = c("Marine West Coast","NW Forested Mtns",
                                  "N. American Deserts", "Great Plains",
                                  "Northern Forests" , "Eastern Temperate" ))
#Read in Sptl Files
epa.folder.sptl <- paste( base.folder.stats, full.drought.name , m.folder, "SEM_EPA_SptlCoeffs", sep="/")                         
setwd(epa.folder.sptl)
all.files <- list.files()

x.int <- strsplit(tail(all.files,1), "\\.")[[1]]
mod.num <- strsplit(x.int[3], "_")[[1]][1]
mod.files <- all.files[ grep(mod.num, all.files) ]
#Add EPA Zone to df and combine
get.epa.zone <- function( epaF ){
   epa.intF <- strsplit( epaF, "_" )
   epa.numF <- epa.intF[[1]][2]
   epa.numF <- substr( epa.numF , start=4, stop=5)
   return( as.numeric( epa.numF) )
}
mod.results <- read.csv( mod.files[1] )
mod.results$EPA_L1CODE <- get.epa.zone( mod.files[1] )
for(ii in 2:nrow(epa.df) ){
   mod.ii <- read.csv( mod.files[ii] )
   mod.ii$EPA_L1CODE <- get.epa.zone(mod.files[ii])
   mod.results <- rbind( mod.results, mod.ii)
}

#Get text names for figure labels
mod.results <- merge( mod.results, epa.df , by="EPA_L1CODE")
#Get Data for Coeffs
mod.strDiv <- mod.results[ which(mod.results$rhs== str.metric) , ]
mod.strDiv$EPA_L1CODE <- factor( mod.strDiv$EPA_L1CODE , levels=c(epa.df$EPA_L1CODE) )
mod.strDiv <- mod.strDiv[ order(mod.strDiv$EPA_L1CODE), ]
mod.spRich <- mod.results[ which(mod.results$rhs== sp.metric ) , ]
mod.spRich <- mod.spRich[ which(mod.spRich$yvar== "Resistance_z" ) , ]
mod.spRich$EPA_L1CODE <- factor( mod.spRich$EPA_L1CODE , levels=c(epa.df$EPA_L1CODE) )
mod.spRich <- mod.spRich[ order(mod.spRich$EPA_L1CODE), ]
mod.spEven <- mod.results[ which(mod.results$rhs== sp.even.metric ) , ]
mod.spEven <- mod.spEven[ which(mod.spEven$yvar== "Resistance_z" ) , ]
mod.spEven$EPA_L1CODE <- factor( mod.spEven$EPA_L1CODE , levels=c(epa.df$EPA_L1CODE) )
mod.spEven <- mod.spEven[ order(mod.spEven$EPA_L1CODE), ]

## Create Figure
#Panels a,b,c - created in ArcGIS
 #uses the values just created above from mod.strDiv, mod.spRich, and mod.spEven
 #to colour the EPA level 1 ecoregion polygons for the USA
#This code create the figures for panels c,d,e
setwd(base.folder.figs)
png( "Fig2_EPA.png" , width=8, height=3.5 , units="in" , res=900)

   par(mfrow=c(1,3))
   par( oma=c(0,4,0,0))
   par( mar=c(10,3,2,2))
   threshs <- 1:nrow(mod.strDiv) 

   ## Structural Richness
   str.main <- mod.strDiv[ , "Estimate"]
   str.up95 <- str.main + 1.96 * (mod.strDiv[ , "Std.err"] )  #upper 95% confidence limits
   str.lo95 <- str.main - 1.96 * (mod.strDiv[ , "Std.err"] )
   y.max <- ceiling( max(str.up95) *20)/20
   y.min <- floor( min(str.lo95) *20)/20
   plot( x=threshs, y=str.main, pch=15, cex=2, col="blue4", cex.lab=axis.size , 
         ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   mtext("Path Coefficient", side = 2, line = 3 , cex = mtext.size)
   axis(1, at=threshs, labels=mod.strDiv$EPA_Name , las=2)
   axis(2, at=seq(-0.8, 0.4, 0.2) ) #, labels=c("-0.4", "-0.2", "0.0", "0.2") )
   abline(h=0, lty=5)
   title( "(D) Structural Richness" , line=0.5, adj=0.02 , cex.main=title.size )
   box()
   arrows(x0=threshs, y0=unlist(str.lo95) , y1=unlist(str.up95),
          length=0, lwd=3 , col="blue4")

   ## Sp Richness
   sp.main <- mod.spRich[ , "Estimate"]
   sp.up95 <- sp.main + 1.96 * (mod.spRich[ , "Std.err"] )  #upper 95% confidience limits
   sp.lo95 <- sp.main - 1.96 * (mod.spRich[ , "Std.err"] )
   y.max <- ceiling( max(sp.up95) *20)/20
   y.min <- floor( min(sp.lo95) *20)/20
   plot( x=threshs, y=sp.main, pch=15, cex=2, col="blue4",cex.lab=axis.size , 
      ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   axis(1, at=threshs, labels=mod.spRich$EPA_Name , las=2)
   axis(2, at=seq(-1.4, 0.2, 0.2) ,
        labels=c(-1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2) ) #, 
   abline(h=0, lty=5)
   title( "(E) Species Richness" , line=0.5, adj=0.02 , cex.main=title.size )
   box()
   arrows(x0=threshs, y0=unlist(sp.lo95) , y1=unlist(sp.up95),
          length=0, lwd=3 , col="blue4")

   ## Sp Evenness
   sp.even.main <- mod.spEven[ , "Estimate"]
   speven.up95 <- sp.even.main + 1.96 * (mod.spEven[ , "Std.err"] )  #upper 95% confidience limits
   speven.lo95 <- sp.even.main - 1.96 * (mod.spEven[ , "Std.err"] )
   y.max <- ceiling( max(speven.up95) *20)/20
   y.min <- floor( min(speven.lo95) *20)/20
   plot( x=threshs, y=sp.even.main, pch=15, cex=2, col="blue4",cex.lab=axis.size , 
         ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   axis(1, at=threshs, labels=mod.spEven$EPA_Name , las=2)
   ylab.seq <- seq(-0.4, 0.4, 0.2)
   axis(2, at= ylab.seq ) # seq(-1.2, 0.2, 0.2) ,
   abline(h=0, lty=5)
   title( "(F) Species Evenness" , line=0.5, adj=0.02 , cex.main=title.size )
   box()
   arrows(x0=threshs, y0=unlist(speven.lo95) , y1=unlist(speven.up95),
          length=0, lwd=3 , col="blue4")
dev.off()


##(3.3) Fig3 - Coeffs by Forest Type ------------------------------------------------
forest.folder.sptl <- paste( base.folder.stats, full.drought.name , m.folder, "SEM_ForestType_SptlCoeffs", sep="/")                         
setwd(forest.folder.sptl)
all.files <- list.files()
mod.files <- tail( all.files , 3)  #last 3 for three forest types

mod.results <- read.csv(mod.files[1] )
mod.results$ForestType <- "Broadleaf"
for (mm in 2:length(mod.files) ) {
   file.mm <- read.csv(mod.files[mm] )
   if(mm==2){ file.mm$ForestType <- "Conifer" }
   if(mm==3){ file.mm$ForestType <- "Mixed" }
   mod.results <- rbind(mod.results, file.mm)
}
mod.results$ForestType <- factor(mod.results$ForestType, 
                                 levels=c("Broadleaf","Mixed","Conifer") )
mod.results <- mod.results[ order(mod.results$ForestType), ]

mod.results <- mod.results[ which(mod.results$yvar== "Resistance_z") , ]
mod.strDiv <- mod.results[ which(mod.results$rhs== str.metric) , ]
mod.spDiv <- mod.results[ which(mod.results$rhs== sp.metric) , ]
mod.spEven <- mod.results[ which(mod.results$rhs== sp.even.metric ) , ]

#Create Figure
setwd(base.folder.figs)
png( "Fig3_ForestType.png" , width=8, height=3.5 , units="in" , res=900)

   par(mfrow=c(1,3))
   par( oma = c(0, 4, 0, 0))
   par( mar=c(7,2,3,2))

   ##(A) Structural Richness
   str.main <- mod.strDiv[ , "Estimate"]
   str.up95 <- str.main + 1.96 * (mod.strDiv[ , "Std.err"] )  #upper 95% confidience limits
   str.lo95 <- str.main - 1.96 * (mod.strDiv[ , "Std.err"] )
   ##(B) Sp Richness
   sp.main <- mod.spDiv[ , "Estimate"]
   sp.up95 <- sp.main + 1.96 * (mod.spDiv[ , "Std.err"] )  #upper 95% confidience limits
   sp.lo95 <- sp.main - 1.96 * (mod.spDiv[ , "Std.err"] )
   ##(C) Sp Evenness
   sp.even.main <- mod.spEven[ , "Estimate"]
   speven.up95 <- sp.even.main + 1.96 * (mod.spEven[ , "Std.err"] )  #upper 95% confidience limits
   speven.lo95 <- sp.even.main - 1.96 * (mod.spEven[ , "Std.err"] )
   #Round to nearest 0.05 -- set same axes for both
   y.max <- ceiling( max(str.up95 , sp.up95, speven.up95) *20)/20
   y.min <- floor( min(str.lo95 , sp.lo95, speven.lo95) *20)/20
   threshs <- 1:nrow(mod.strDiv) 
   y.ax.labs <- c(-0.6, -0.4, -0.2, 0, 0.2, 0.4)

   ##Str Rich
   plot( x=threshs, y=str.main, pch=15, cex=2, col="blue4", cex.lab=axis.size , 
      xlim=c(0.25, 3.75),  ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   mtext("Path Coefficient", side = 2, line = 3 , cex = mtext.size)
   axis(1, at=threshs, labels=mod.strDiv$ForestType , las=2 , cex.axis = axis.size )
   axis(2, at=y.ax.labs , cex.axis = label.size ) #, 
   abline(h=0, lty=5)
   title( "(A) Structural Richness" , line=0.5, adj=0.02 , cex.main=title.size)
   box()
   arrows(x0=threshs, y0=unlist(str.lo95) , y1=unlist(str.up95),
          length=0, lwd=3 , col="blue4")
   #(B) Sp Richness
   plot( x=threshs, y=sp.main, pch=15, cex=2, col="blue4",cex.lab=axis.size , 
      xlim=c(0.25, 3.75), ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   axis(1, at=threshs, labels=mod.spDiv$ForestType , las=2 , cex.axis = axis.size )
   axis(2, at=y.ax.labs , cex.axis = label.size  )
   abline(h=0, lty=5)
   title( "(B) Species Richness" , line=0.5, adj=0.02 , cex.main=title.size)
   box()
   arrows(x0=threshs, y0=unlist(sp.lo95) , y1=unlist(sp.up95),
          length=0, lwd=3 , col="blue4")
   #(C) Sp Evenness
   plot( x=threshs, y=sp.even.main, pch=15, cex=2, col="blue4", cex.lab=axis.size , 
         xlim=c(0.25, 3.75), ylim=c(y.min ,  y.max), axes=FALSE, xlab="" , ylab="Path Coefficient")
   axis(1, at=threshs, labels=mod.spEven$ForestType , las=2 , cex.axis = axis.size )
   axis(2, at=y.ax.labs , cex.axis = label.size )
   abline(h=0, lty=5)
   title( "(C) Species Evenness" , line=0.5, adj=0.02 , cex.main=title.size)
   box()
   arrows(x0=threshs, y0=unlist(speven.lo95) , y1=unlist(speven.up95),
          length=0, lwd=3 , col="blue4")
dev.off()


##(3.4) Fig4 - Coeffs by Spei Thresholds -----------------------------------------------------------
spei.folder.sptl <- paste( base.folder.stats, "T-1.00" , m.folder, "SEM_MultiThresh_SptlCoeffs", sep="/")                         
setwd(spei.folder.sptl)
sptl.files <- list.files()
sptl.files <- sptl.files[ grep("SptlCoeffs", sptl.files)]

get.DThresh <- function(xF){
   xF <- strsplit(xF, "-")[[1]][2]
   x2F <- strsplit(xF, "_")[[1]][1]
   x2F <- as.numeric(x2F)
   return(x2F)
}
mod.results <- read.csv( sptl.files[1] ) 
mod.results$DThresh <- get.DThresh( sptl.files[1] )
for(ii in 2:length(sptl.files)){ 
   results.ii <- read.csv( sptl.files[ii]) 
   results.ii$DThresh <- get.DThresh( sptl.files[ii] )
   mod.results <- rbind( mod.results, results.ii)
}
mod.results$DThresh <- -1 *  mod.results$DThresh #To get the negative back
mod.results <- mod.results[ order(mod.results$DThresh, decreasing=T) ,]

mod.results <- mod.results[ which(mod.results$yvar=="Resistance_z") , ]
mod.strDiv <- mod.results[ which(mod.results$rhs==str.metric) , ]
mod.spRich <- mod.results[ which(mod.results$rhs==sp.metric) , ]
mod.strEven <- mod.results[ which(mod.results$rhs==str.even.metric) , ]
mod.spEven <- mod.results[ which(mod.results$rhs==sp.even.metric) , ]

#Structural Richness 
str.main <- mod.strDiv[ , "Estimate"]
str.up95 <- str.main + 1.96 * (mod.strDiv[ , "Std.err"] )  #upper 95% confidience limits
str.lo95 <- str.main - 1.96 * (mod.strDiv[ , "Std.err"] )
y.max.str <- ceiling( max(str.up95) *20)/20
y.min.str <- floor( min(str.lo95) *20)/20
#Species Richness
spRich.main <- mod.spRich[ , "Estimate"]  
spRich.up95 <- spRich.main + 1.96 * (mod.spRich[ , "Std.err"] )
spRich.lo95 <- spRich.main - 1.96 * (mod.spRich[ , "Std.err"] )
y.max.sp <- ceiling( max(spRich.up95) *20)/20
y.min.sp <- floor( min(spRich.lo95) *20)/20
#Species Evenness
spEven.main <- mod.spEven[ , "Estimate"]  
spEven.up95 <- spEven.main + 1.96 * (mod.spEven[ , "Std.err"] )
spEven.lo95 <- spEven.main - 1.96 * (mod.spEven[ , "Std.err"] )
y.max.speven <- ceiling( max(spEven.up95) *20)/20
y.min.speven <- floor( min(spEven.lo95) *20)/20

threshs <- mod.strDiv$DThresh 

#Create Plot
setwd(base.folder.figs)
png( "Fig4_StressGradient.png" , width=8, height=3.5 , units="in" , res=900)

   par(mfrow=c(1,3))
   par( oma = c(0, 4, 0, 0))
   par( mar=c(5,2,3,2) )

   #Str Plot
   plot( x=threshs, y=str.main, pch=15, cex=2 , cex.lab=axis.size , 
         xlim=c(-1+side.spacing, -2-side.spacing) , ylim=c(y.min.str ,  y.max.str), axes=FALSE, xlab="SPEI Drought Threshold" , ylab="Model Coefficient")
   mtext("Path Coefficient", side = 2, line = 3 , cex = mtext.size)
   axis(1, at=threshs, labels=mod.strDiv$DThresh , cex.axis = label.size )
   ylab.seq <- seq(y.min.str, y.max.str, 0.05)
   axis(2, at=ylab.seq, labels=as.character(ylab.seq) , cex.axis = label.size  )
   title( "(A) Structural Richness" , line=0.5, adj=0.02 , cex.main=title.size)
   box()
   abline(h=0, lty=5)
   arrows(x0=threshs, y0=unlist(str.lo95) , y1=unlist(str.up95),
       length=0, lwd=3)

   #Sp Plot
   plot( x=threshs, y=spRich.main, pch=15, cex=2 , cex.lab=axis.size , 
         xlim=c(-1+side.spacing, -2-side.spacing), ylim=c( y.min.sp ,  y.max.sp ), axes=FALSE, xlab="SPEI Drought Threshold" , ylab="") 
   axis(1, at=threshs, labels=mod.strDiv$DThresh, cex.axis = label.size )
   ylab.seq <- seq(y.min.sp, y.max.sp, 0.05)
   axis(2, at=ylab.seq, labels=as.character(ylab.seq) , cex.axis = label.size  )
   title( "(B) Species Richness" , line=0.5, adj=0.02, cex.main=title.size)
   box()
   abline(h=0, lty=5)
   arrows(x0=threshs, y0=unlist(spRich.lo95) , y1=unlist(spRich.up95),
          length=0, lwd=3)

   #Sp Evenness
   plot( x=threshs, y=spEven.main, pch=15,  cex=2 , cex.lab=axis.size , 
         xlim=c(-1+side.spacing, -2-side.spacing), ylim=c(y.min.speven,  y.max.speven), axes=FALSE, xlab="SPEI Drought Threshold" , ylab="") 
   axis(1, at=threshs, labels=mod.strDiv$DThresh , cex.axis = label.size)
   ylab.seq <- seq(y.min.speven, y.max.speven, 0.05)
   axis(2, at=ylab.seq, labels=as.character(ylab.seq) , cex.axis = label.size  )
   title( "(C) Species Evenness" , line=0.5, adj=0.02, cex.main=title.size)
   box()
   abline(h=0, lty=5)
   arrows(x0=threshs, y0=unlist(spEven.lo95) , y1=unlist(spEven.up95),
          length=0, lwd=3)
dev.off()



##(3.5) Supplemental Figures ---------------------------------------------------

## Fig S1 - Map created in ArcGIS using the EPA Level Ecoregions file (dissolved to level 1)
#  and the plot coordinates

## Fig S2 - In the pre-processing code, the thresholds of 500MgC/ha/yr and 100MgC/ha/yr were not applied.
#  Otherwise this figure was created using the same process described for Fig 1 above.

## Fig S3 - See pre-processing code, end of section #8 for code to create this figure.

## Fig S4 - Map created in ArcGIS using the EPA Level Ecoregions file (dissolved to level 1)
#  and the plot coordinates for plots that experienced any of the 5 drought thresholds.
#  (i.e. from columns LATx and LONx in the ForestDroughtData.csv file created in the processing steps)

## Fig S5 - Created with the same process as Fig 1, 
#  but uses model m11 and files within the m11.folder

## Fig S6 - Created with the same process as Fig 1, 
#  but adds a line in the model 11 formula for Biomass_log_z as a predictor variable

## Fig S7 - Created with the same process as Fig 1, 
#  but uses model m15 and files within the m15.folder

## Fig S8 - Created with the same process as Fig 1, 
#  but uses model m17 and files within the m17.folder

## Fig S9 - Created with the same process as Fig 1, 
#  but uses model m19 and files within the m19.folder

## Fig S10 - Multiple Drought Categories - created with the same process as Fig 4,
#  but in the pre-processing steps, end of section 9, follow steps to extract plots
#  that have experienced multiple drought events and then repeate the analyses above

## Fig S11 - Multiple Regression 
setwd(m12.multi.folder)
my.results <- read.csv("MultiReg_coeffs_BWsel.csv")
my.results$X <- c("Intercept", "Structural Richness", "Species Richness", 
                  "Structural Evenness","Species Evenness",
                  "Stand Age", "Avg 20yr Temp", "Bulk Density")
colnames(my.results)[1] <- "CNames"
my.results <- my.results[-1, ] #Remove intercept
my.results$Ord <- 1:nrow(my.results)
my.results <- my.results[ order(my.results$Ord, decreasing = T), ]
my.results$lower <- my.results$Estimate - (1.96*my.results$Std..Error)
my.results$upper <- my.results$Estimate + (1.96*my.results$Std..Error)
my.results$CNames <- factor(my.results$CNames, levels=my.results$CNames )
#Create Plot
setwd(base.folder.figs)
png( "Fig_MultiReg.png" , width=6.5, height=6, res=900 , units="in")
   par( mar=c(5,9,3,2) )
   plot( x=my.results$Estimate , 1:nrow(my.results) ,
      xlim=c(-0.5 , 0.25) , pch=15 , cex=1.5 ,
      xlab="Standardized Slope Estimate" , yaxt="n", ylab="" )
   axis(2, at=1:nrow(my.results) , labels=my.results$CNames , las=1  )
   arrows(x0=my.results$lower , y0=1:nrow(my.results) , 
          x1=my.results$upper , y1=1:nrow(my.results) , length=0, lwd=3)
   abline(v=0, lty=5)
dev.off()


## Fig S12  - this follows the same approach as m11.formula, but replaces StrRich_z with 
#  other structural diversity variables

str.vars <- c("Biomass_log_z", "Height_mean_z",  "BasalArea_log_z", 
              "SDI_log_z", "nStems_log_z")
str.names <- c("Biomass", "Avg Tree Height", "Basal Area",
               "Stand Density Index", "Number of Stems")
#Create New Folders
dd.fol <- c("T-2.00")
mm.fol <- c("m41_Biomass_log_z", "m42_Height_mean_z",  "m43_BasalArea_log_z", 
            "m44_SDI_log_z", "m45_nStems_log_z")
ss.fol <- c("SEM_reg", "SEM_SptlCoeffs", "SEM_EPA_reg", "SEM_EPA_SptlCoeffs",
            "SEM_ForestType_reg","SEM_ForestType_SptlCoeffs",
            "SEM_MultiThresh", "SEM_MultiThresh_SptlCoeffs")
for(dd in 1:length(dd.fol) ){
   dir.x1 <- paste0( base.folder.stats, "/", dd.fol[dd] )
   if( !dir.exists(dir.x1) ){ dir.create(dir.x1) }
   for(mm in 1:length(mm.fol) ){
      dir.x2 <- paste0( dir.x1, "/", mm.fol[mm] )
      if( !dir.exists(dir.x2) ){ dir.create(dir.x2) }
      for(ss in 1:length(ss.fol) ){
         dir.x3 <- paste0( dir.x2, "/", ss.fol[ss] )
         if( !dir.exists(dir.x3) ){ dir.create(dir.x3) }
      }#ss
   }#mm
}#dd
getAIC.value <- function( dirF, rowF=ss , aic.dfF=aic.df ){
   setwd(dirF)
   ss.filesF <- list.files()
   ss.filesF <- ss.filesF[grep("csv", ss.filesF)]
   str.csvF <- read.csv( ss.filesF[ length(ss.filesF) ] )
   aic.dfF[rowF,"AIC"] <- str.csvF[ which(str.csvF$lhs=="AIC"), 4]  
   return(aic.dfF)
}
aic.df <- data.frame( StrVar=str.vars, AIC=NA, StrName=str.names)

for(ss in 1:length(str.vars) ){
   m.ss.folder <- paste0( base.folder.stats,"/", full.drought.name,"/m4",ss,"_",str.vars[ss])
   run.SEM( model.nameF="m11.basic" , model.formulaF=m11.formula , 
         str.varF = str.vars[ss]  , sp.varF = "SpRich_log_z"  ,
         base.folderF=m.ss.folder , groupF=NULL  )
   aic.df <- getAIC.value( paste0(m.ss.folder, "/SEM_reg"), rowF=ss )
}
#Add Structural Richness and Structural Shannon
aic.df[6:7,1] <- c("StrRich_z","StrShannon_z")
aic.df[6:7,3] <- c("Str Diversity: Richness","Str Diversity: Shannon")
aic.df <- getAIC.value( paste0(m11.folder, "/SEM_reg"), rowF=6)
aic.df <- getAIC.value( paste0(m13.folder, "/SEM_reg"), rowF=7)
#Create Barplot
aic.df <- aic.df[ order(aic.df$AIC), ]
aic.df$StrName <- factor( aic.df$StrName, levels=aic.df$StrName)
setwd( "F:/DataFiles" )
png( "FigS12_AICBarplot2.png", width=14, height=12, units="cm", res=300 )
   par( mar=c(10,4,1,1) )
   xx <- barplot( aic.df$AIC ~ aic.df$StrName , 
               xlab="", ylab="AIC Value", las=2,
               ylim=c(0, max(aic.df$AIC)+1500) , col="steelblue3")  
   vals.round <- round(aic.df$AIC) #to whole number
   text(x = xx, y=vals.round-200, label = vals.round, pos=3, cex=0.7, col="grey60")
   box()
dev.off()   


## Fig S13 - Created with the same process as Fig 1, 
#  but uses model m13 and files within the m13.folder


##Fig S14 - Species Richness and Juniperus Abundance
#Note: user must run the preprocessing code to obtain FIA data & metrics and generate a new ForestDroughtData.csv
setwd( "F:/DataFiles" )
data.df.full <- read.csv( "ForestDroughtData.csv")
data.df <- data.df.full[ which(data.df.full$SPEIcat == -2), ]
states <- c("MO","KS","NE","SD","ND","NM","AZ","NV",
            "CO","UT","WY","MT","ID","WA","OR","CA",
            "TN","TX","KY","VA","OK",
            "WV","MD","DE","NJ","PA","NY",
            "CT","RI","MA","VT","NH","ME",
            "OH","MI","IN","IL","WI","MN","IA",
            "SC","GA","FL","AL","MS","AR","LA","NC") 
##Create Template File
setwd( paste0( "F:/FIA_Public_Data/", states[1] ) )
tree.df <- read.csv( paste0( states[1], "_TREE.csv") )
tree.df <- tree.df[ which(tree.df$PLT_CN %in% data.df$Before_PLT_CN), ]
##Run Other States
for( ii in 2:length( states ) ){
   setwd( paste0( "F:/FIA_Public_Data/", states[ii] ) )
   fia.tree.ii <- read.csv( paste0( states[ii], "_TREE.csv") )
   fia.tree.ii <- fia.tree.ii[ which(fia.tree.ii$PLT_CN %in% data.df$Before_PLT_CN), ]
   tree.df <- rbind( tree.df, fia.tree.ii )
   print(ii)
}
colnames(tree.df)[2] <- "Before_PLT_CN"
colnames(tree.df)[16] <- "SPECIESCode"
#Calc abundance of each type of tree
tree.df$Abundance <- 1
comm.mat <- reshape::cast(tree.df[ ,c("Before_PLT_CN","SPECIESCode","Abundance")], 
                             Before_PLT_CN ~ SPECIESCode, sum)  #PLT_CN is unique plot id
rownames(comm.mat) <- comm.mat$Before_PLT_CN 	    #Make Plot ID number the rownames
comm.mat$Before_PLT_CN <- NULL			             #Remove this to have a clean community matrix
comm.mat <- as.data.frame(comm.mat)	       
#Reorder
comm.mat <- comm.mat[ as.character( data.df$Before_PLT_CN ), ]
all.equal( rownames(comm.mat), as.character( data.df$Before_PLT_CN ) ) 
#Plot
setwd(base.folder.figs)
png( "FigS14_Juniperus.png" , width=6.5, height=6, res=900 , units="in")
   plot( comm.mat$`61` , data.df$q0,  pch=16 ,
      xlab= expression( paste("Number of ", italic("Juniperus ashei"), " trees" ) ) , 
      ylab="Species Richness")
dev.off()


##Fig S15 - Map: SPEI-NPP
#Figure created in ArcGIS using EPA Ecoregions Level 1 as the polygon boundaries
#See end of Pre-processing section 8 for code to extract median R2 values
#from Productivity-SPEI relationships.


##End##

