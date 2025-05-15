# This R code provides functions to run the core analyses for the paper: 

# Crockett ETH, Qingfeng G, Atkins JE, Sun G, Potter KM, Coztanza J,
# Ollinger S, Woodall C, McNulty  S, Trettin C, Holgerson, J, and Xiao J.
# Influences of structural and species diversity on forest resistance to drought.

# (c) Erin Crockett, 2025
# erin.crockett@unbc.ca


## Arguments for functions below:
#model.nameF: what to call when save to file
#model.formulaF: the formula to read into lavaan path model
#y.varF: the y variable (resistance)
#str.varF: the structural variable
#sp.varF: the species varaible, eg sp richness
#full.summaryF: df, sheet from the all.mod.summary file
#base.folderF: where to read in/save files -- including initial mod formula
#restrictPathsF: used to indicate shared covariance between richness and evenness
#dataF: df, which data to use to create the models
#groupF: chr, if want to do a multigruop model
#n.grpF: num, min number of points within each group
#add.sptlCorF: if should modify the coeffs for spatial correatlion
#sptl.ptsF: sptl df: with coordinates of locations


##[1] Run path analysis with lavaan with str div and sp div OR str even and sp even
run.SEM <- function( model.nameF , model.formulaF, 
                      y.varF="Resistance" , str.varF, sp.varF ,
                      full.summaryF = all.mod.summary ,  
                      base.folderF=base.save.folder , 
                      dataF=data.df , groupF=NULL , n.grpF=20 , 
                      add.sptlCorF=TRUE , sptl.ptsF=mod.prj ){
   
   run.mod.againF <- TRUE
   mod.numberF <- 101      
   if( sp.varF %in% c("Sp_Even_z", "Str_Even_z") ){   
      model.formulaF[1] <- paste0( y.varF, "~", str.varF, "+", sp.varF, "+", model.formulaF[1] )
      model.formulaF[2] <- paste0( str.varF, "~", model.formulaF[2] )
      model.formulaF[3] <- paste0( sp.varF, "~", model.formulaF[3] )
   }else{
      model.formulaF[1] <- paste0( y.varF, "~", str.varF, "+", sp.varF, "+", model.formulaF[1] )
      model.formulaF[2] <- paste0( str.varF, "~", sp.varF, "+", model.formulaF[2] )
      model.formulaF[3] <- paste0( sp.varF, "~", model.formulaF[3] )
   }
   model.formulaF <- gsub( " ", "", model.formulaF)
   all.models.listF <- list()
   all.models.listF[[mod.numberF-100]] <- model.formulaF 
   all.namesF <- unlist( strsplit( model.formulaF, split="\\+") )
   all.namesF <- all.namesF[ -grep( "~" , all.namesF ) ]
   all.namesF <- unique(all.namesF)
   all.namesF <- c(y.varF, str.varF, all.namesF)
   cnames.sumryF <- c("ModName","ModNum","AIC","ModPval","y.R2","str.R2","sp.R2")
   cnames.sumryF <- c(cnames.sumryF,
                      paste0(all.namesF,"_P"   ) ,
                      paste0(all.namesF,"_Est" ) , 
                      paste0(all.namesF,"_SE"  )  )
   sumry.regF <- as.data.frame( matrix( NA, nrow=0, ncol=length(cnames.sumryF)))
   colnames(sumry.regF) <- cnames.sumryF
   sumry.SptlCorF <- sumry.regF 
   sumry.psemF <- sumry.regF 
   
   sem.modF <- lavaan::sem( model.formulaF, data=dataF  )
   name.of.modF <- paste(model.nameF, mod.numberF , sep=".")
   setwd( paste0( base.folderF, "/SEM_reg") )
   writeLines( model.formulaF, con= paste0( name.of.modF, ".txt"))
   full.summaryF <- EC_addInfo(  modF=sem.modF, name.of.modF, str.varF, sp.varF, 
                                 y.varF, mod.sumryF=full.summaryF )
   write.csv(full.summaryF, file= paste0( "_Summary_", model.nameF, ".csv") , row.names=FALSE )
   save.SEMinfo( sem.modF , name.of.modF , dirF=paste0( base.folderF, "/SEM_reg") )
   print( paste( "Mod Pvalue ", round( unlist(summary(sem.modF)$test[[1]][6] ), digits=4 ) ) )
   if( add.sptlCorF ){
      sem.mod.sptlCorF <- lavSpatialCorrect( sem.modF, sptl.ptsF@coords[ ,1], sptl.ptsF@coords[ ,2])
      sptl.coeffs <- all.SptlCoeffs( sem.mod.sptlCorF )
      setwd( paste0( base.folderF, "/SEM_SptlCoeffs") )
      write.csv(sptl.coeffs, file= paste0( name.of.modF, "SptlCoeffs.csv") , row.names=FALSE)
   }
   if( !is.null(groupF) ){
      group.nF <- table( dataF[ ,groupF] )
      group.dfF <- as.data.frame(group.nF)  #use later to save the table
      group.nF <- group.nF[ group.nF >= n.grpF ]  #min number of pts per group (eg 20)
      data.groupF <- dataF[ which(dataF[ ,groupF] %in% names(group.nF)), ]
      sem.mod.groupF <- lavaan::sem( model.formulaF, data=data.groupF , group=groupF )
      if( groupF == "EPA_L1CODE"){ This.FolderF <- paste0( base.folderF, "/SEM_EPA_reg") }
      if( groupF == "ForestType"){ This.FolderF <- paste0( base.folderF, "/SEM_ForestType_reg") }
      save.SEMinfo( sem.mod.groupF , name.of.modF, groupF=TRUE,  dirF=This.FolderF)
      for(ff in 1:length(group.nF) ){
         group.idF <- names(group.nF)[ff]
         data.ffF <- data.groupF[ which(data.groupF[ ,groupF] == group.idF ), ]
         sptl.pts.ffF <- sptl.ptsF[ which(sptl.ptsF@data[ ,groupF] == group.idF ) , ]
         sem.mod.group.ffF <- lavaan::sem( model.formulaF, data=data.ffF  )
         sem.mod.sptlCor.ffF <- lavSpatialCorrect( sem.mod.group.ffF, sptl.pts.ffF@coords[ ,1], sptl.pts.ffF@coords[ ,2])
         sptl.coeffs.ff <- all.SptlCoeffs( sem.mod.sptlCor.ffF )
         if( groupF == "EPA_L1CODE"){
            setwd( paste0( base.folderF, "/SEM_EPA_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
         if( groupF == "ForestType"){
            setwd( paste0( base.folderF, "/SEM_ForestType_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
         if( groupF == "SPEI_factor"){
            setwd( paste0( base.folderF, "/SEM_MultiThresh_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_Thresh",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
      }
      colnames(group.dfF)[1] <- groupF
      write.csv( group.dfF, file="_nPts_perGroup.csv"  , row.names=FALSE )
   }#close if add Groups (eg EPA)
   
   while( run.mod.againF ){
      mod.numberF <- mod.numberF + 1
      print( paste( "----- Starting " , mod.numberF, "-----"))
      ##Check if Should dd a new variable to the model
      x <- subset(modindices(sem.modF), mi > 2 )
      if( any(which(x$rhs %in% c(y.varF , str.varF ) ) ) ){
         x <- x[ - which(x$rhs %in% c(y.varF , str.varF ) ) , ]
      }
      #If still have vars to add
      if( nrow(x) > 0 ){
         x2 <- x[order(x$mi, decreasing=TRUE), ]
         #Add the top one (but not directions that dont make sense)
         for(ff in 1:nrow(x2)){
            if( x2[ff,"lhs"] %in% c(y.varF, str.varF, sp.varF)  ){
               if( x2[ff,"lhs"] == y.varF   ){ model.formulaF[1] <- paste( model.formulaF[1] , x2[ff,"rhs"] , sep="+" )  } 
               if( x2[ff,"lhs"] == str.varF ){ model.formulaF[2] <- paste( model.formulaF[2] , x2[ff,"rhs"] , sep="+" )  } 
               if( x2[ff,"lhs"] == sp.varF  ){ model.formulaF[3] <- paste( model.formulaF[3] , x2[ff,"rhs"] , sep="+" )  } 
               break() #ff loop, so only add one at a time
            }
         }#close ff loop
      }else{
         coeff.tabF <- summary(sem.modF)$pe
         coeff.tabF <- coeff.tabF[ which(coeff.tabF$op == "~"), ]
         coeff.tabF <- coeff.tabF[ order(coeff.tabF$pvalue, decreasing=TRUE), ]
         if( coeff.tabF[1,"pvalue"] > 0.05 ){
            lhs.remF <- coeff.tabF[1,"lhs"]
            rhs.remF <- coeff.tabF[1,"rhs"]
            if( lhs.remF == y.varF ){   rem.row <- 1 }
            if( lhs.remF == str.varF ){ rem.row <- 2 }
            if( lhs.remF == sp.varF ){  rem.row <- 3 }
            new.fileF <- strsplit( model.formulaF[rem.row] , split="\\+" )[[1]]  #Split the list
            new.file.part1F <- strsplit( new.fileF, split="~")[[1]]   #Split to get the y var
            new.fileF[1] <- new.file.part1F[2]                        #Re-Add the first variable to the explanatory variable string
            new.fileF <- new.fileF[ -which(new.fileF == rhs.remF)]    #Remove the least important variable
            new.formF <- paste(new.fileF , collapse="+")              #Paste together to create RHS of formula
            new.formF <- paste(new.file.part1F[1], new.formF, sep="~") #Create full New formula
            model.formulaF[rem.row] <- new.formF           #Replace Line in existing formula
         }else{
            run.mod.againF <- FALSE
         }
      }#close else
      
      model.current.listF <- list()
      model.current.listF[[1]] <- model.formulaF
      if( model.current.listF %in% all.models.listF ){
         run.mod.againF <- FALSE
      }else{
         all.models.listF[[mod.numberF-100]] <- model.formulaF
      }
      if( run.mod.againF ){
         sem.modF <- lavaan::sem( model.formulaF, data=dataF  )
         name.of.modF <- paste(model.nameF, mod.numberF , sep=".")
         setwd( paste0( base.folderF, "/SEM_reg") )
         writeLines( model.formulaF, con= paste0( name.of.modF, ".txt"))
         full.summaryF <- EC_addInfo(  modF=sem.modF, name.of.modF, str.varF, sp.varF, 
                                       y.varF, mod.sumryF=full.summaryF )
         setwd( paste0( base.folderF, "/SEM_reg") )
         write.csv(full.summaryF, file= paste0( "_Summary_", model.nameF, ".csv") , row.names=FALSE )
         save.SEMinfo( sem.modF , name.of.modF , dirF=paste0( base.folderF, "/SEM_reg") )
         print( paste( "Mod Pvalue ", round( unlist(summary(sem.modF)$test[[1]][6] ), digits=4 ) ) )
         if( add.sptlCorF ){
            sem.mod.sptlCorF <- lavSpatialCorrect( sem.modF, sptl.ptsF@coords[ ,1], sptl.ptsF@coords[ ,2])
            sptl.coeffs <- all.SptlCoeffs( sem.mod.sptlCorF )
            setwd( paste0( base.folderF, "/SEM_SptlCoeffs") )
            write.csv(sptl.coeffs, file= paste0( name.of.modF, "SptlCoeffs.csv") , row.names=FALSE)
         }
         #SEM by GROUP
         if( !is.null(groupF) ){
            group.nF <- table( dataF[ ,groupF] )
            group.nF <- group.nF[ group.nF >= n.grpF ]  #min number of pts per group (eg 20)
            data.groupF <- dataF[ which(dataF[ ,groupF] %in% names(group.nF)), ]
            sem.mod.groupF <- lavaan::sem( model.formulaF, data=data.groupF , group=groupF )
            if( groupF == "EPA_L1CODE"){ This.FolderF <- paste0( base.folderF, "/SEM_EPA_reg") }
            if( groupF == "ForestType"){ This.FolderF <- paste0( base.folderF, "/SEM_ForestType_reg") }
            save.SEMinfo( sem.mod.groupF , name.of.modF, groupF=TRUE,  dirF=This.FolderF )
            for(ff in 1:length(group.nF) ){
               group.idF <- names(group.nF)[ff]
               data.ffF <- data.groupF[ which(data.groupF[ ,groupF] == group.idF ), ]
               sptl.pts.ffF <- sptl.ptsF[ which(sptl.ptsF@data[ ,groupF] == group.idF ) , ]
               sem.mod.group.ffF <- lavaan::sem( model.formulaF, data=data.ffF  )
               sem.mod.sptlCor.ffF <- lavSpatialCorrect( sem.mod.group.ffF, sptl.pts.ffF@coords[ ,1], sptl.pts.ffF@coords[ ,2])
               sptl.coeffs.ff <- all.SptlCoeffs( sem.mod.sptlCor.ffF )
               if( groupF == "EPA_L1CODE"){
                  setwd( paste0( base.folderF, "/SEM_EPA_SptlCoeffs") )
                  write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
               }
               if( groupF == "ForestType"){
                  setwd( paste0( base.folderF, "/SEM_ForestType_SptlCoeffs") )
                  write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
               }
            }
         }#close if add Groups (eg EPA)
      }#Close if run model again
      if( mod.numberF >= 150 ){
         run.mod.againF <- FALSE
      }
   }#Close while loop
} #Close function


##[2] Run path analysis with Multiple SPEI Thresholds
SEM.MultiThresh <- function( model.nameF , model.formulaF , 
                             y.varF="Resistance" , str.varF, sp.varF ,
                             full.summaryF=all.mod.summary  ,  
                             base.folderF=base.save.folder , dataF=data.df , 
                             groupF=NULL , n.grpF=20 , 
                             add.sptlCorF=TRUE , sptl.ptsF=mod.prj ){
   
   name.of.modF <- model.nameF 
   if( !is.null(groupF) ){
      group.nF <- table( dataF[ ,groupF] )
      group.dfF <- as.data.frame(group.nF)  #use later to save the table
      group.nF <- group.nF[ group.nF >= n.grpF ]  #min number of pts per group (eg 20)
      data.groupF <- dataF[ which(dataF[ ,groupF] %in% names(group.nF)), ]
      sem.mod.groupF <- lavaan::sem( model.formulaF, data=data.groupF , group=groupF )
      if( groupF == "EPA_L1CODE"){ save.folderF <- paste0( base.folderF, "/SEM_EPA_reg") }
      if( groupF == "SPEI_factor"){ save.folderF <- paste0( base.folderF, "/SEM_MultiThresh") }      
      save.SEMinfo( sem.mod.groupF , name.of.modF, groupF=TRUE,  dirF=save.folderF )
      for(ff in 1:length(group.nF) ){
         group.idF <- names(group.nF)[ff]
         data.ffF <- data.groupF[ which(data.groupF[ ,groupF] == group.idF ), ]
         sptl.pts.ffF <- sptl.ptsF[ which(sptl.ptsF@data[ ,groupF] == group.idF ) , ]
         sem.mod.group.ffF <- lavaan::sem( model.formulaF, data=data.ffF  )
         sem.mod.sptlCor.ffF <- lavSpatialCorrect( sem.mod.group.ffF, sptl.pts.ffF@coords[ ,1], sptl.pts.ffF@coords[ ,2])
         sptl.coeffs.ff <- all.SptlCoeffs( sem.mod.sptlCor.ffF )
         if( groupF == "EPA_L1CODE" ){
            setwd( paste0( base.folderF, "/SEM_EPA_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
         if( groupF == "SPEI_factor" ){
            setwd( paste0( base.folderF, "/SEM_MultiThresh_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_Thresh",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
      }
      colnames(group.dfF)[1] <- groupF
      write.csv( group.dfF, file="_nPts_perGroup.csv"  , row.names=FALSE )
   }
}


##[3] Run Path Analysis with Div and Even
run.SEM.wEven <- function( model.nameF , model.formulaF, 
                           y.varF="Resistance", str.varF, sp.varF ,
                           str.evenF = NULL , sp.evenF = NULL ,
                           restrictPathsF = "ForceCovary",
                           full.summaryF = all.mod.summary , base.folderF=base.save.folder , 
                           dataF=data.df , groupF=NULL , n.grpF=20 , 
                           add.sptlCorF=TRUE , sptl.ptsF=mod.prj ){

   run.mod.againF <- TRUE
   mod.numberF <- 101      
   if( is.null(str.evenF) ){
      if( sp.varF %in% c("Sp_Even_z", "Str_Even_z") ){   #then dont have sp-->str
         model.formulaF[1] <- paste0( y.varF, "~", str.varF, "+", sp.varF, "+", model.formulaF[1] )
         model.formulaF[2] <- paste0( str.varF, "~", model.formulaF[2] )
         model.formulaF[3] <- paste0( sp.varF, "~", model.formulaF[3] )
      }else{
         model.formulaF[1] <- paste0( y.varF, "~", str.varF, "+", sp.varF, "+", model.formulaF[1] )
         model.formulaF[2] <- paste0( str.varF, "~", sp.varF, "+", model.formulaF[2] )
         model.formulaF[3] <- paste0( sp.varF, "~", model.formulaF[3] )
      }
   #With all four: str+sp rich+even  
   }else{
      model.formulaF[1] <- paste0( y.varF,    "~", str.varF, "+", sp.varF, "+", 
                                   str.evenF, "+", sp.evenF, "+", model.formulaF[1] )
      model.formulaF[2] <- paste0( str.varF,  "~", sp.varF, "+", model.formulaF[2] )
      model.formulaF[3] <- paste0( sp.varF,   "~", model.formulaF[3] )
      model.formulaF[4] <- paste0( str.evenF, "~", model.formulaF[4] )
      model.formulaF[5] <- paste0( sp.evenF,  "~", model.formulaF[5] )
   }
   if( restrictPathsF == "ForceCovary" ){
      model.formulaF[6] <- paste0( str.varF, "~~" ,str.evenF )
      model.formulaF[7] <- paste0( sp.varF, "~~" ,sp.evenF )
   }
   
   model.formulaF <- gsub( " ", "", model.formulaF)
   all.models.listF <- list()
   all.models.listF[[mod.numberF-100]] <- model.formulaF 
   all.namesF <- unlist( strsplit( model.formulaF, split="\\+") )
   all.namesF <- all.namesF[ -grep( "~" , all.namesF ) ]
   all.namesF <- unique(all.namesF)
   all.namesF <- c(y.varF, str.varF, all.namesF)
   cnames.sumryF <- c("ModName","ModNum","AIC","ModPval","y.R2","str.R2","sp.R2")
   cnames.sumryF <- c(cnames.sumryF,
                      paste0(all.namesF,"_P"   ) ,
                      paste0(all.namesF,"_Est" ) , 
                      paste0(all.namesF,"_SE"  )  )
   sumry.regF <- as.data.frame( matrix( NA, nrow=0, ncol=length(cnames.sumryF)))
   colnames(sumry.regF) <- cnames.sumryF
   sumry.SptlCorF <- sumry.regF 
   sumry.psemF <- sumry.regF 

   sem.modF <- lavaan::sem( model.formulaF, data=dataF  )
   name.of.modF <- paste(model.nameF, mod.numberF , sep=".")
   setwd( paste0( base.folderF, "/SEM_reg") )
   writeLines( model.formulaF, con= paste0( name.of.modF, ".txt"))
   full.summaryF <- EC_addInfo(  modF=sem.modF, name.of.modF, str.varF, sp.varF, 
                                 y.varF, mod.sumryF=full.summaryF )
   write.csv(full.summaryF, file= paste0( "_Summary_", model.nameF, ".csv") , row.names=FALSE )
   save.SEMinfo( sem.modF , name.of.modF , dirF=paste0( base.folderF, "/SEM_reg") )
   print( paste( "Mod Pvalue ", round( unlist(summary(sem.modF)$test[[1]][6] ), digits=4 ) ) )
   if( add.sptlCorF ){
      sem.mod.sptlCorF <- lavSpatialCorrect( sem.modF, sptl.ptsF@coords[ ,1], sptl.ptsF@coords[ ,2])
      sptl.coeffs <- all.SptlCoeffs( sem.mod.sptlCorF )
      setwd( paste0( base.folderF, "/SEM_SptlCoeffs") )
      write.csv(sptl.coeffs, file= paste0( name.of.modF, "SptlCoeffs.csv") , row.names=FALSE)
   }
   #SEM by GROUP
   if( !is.null(groupF) ){
      group.nF <- table( dataF[ ,groupF] )
      group.dfF <- as.data.frame(group.nF)  #use later to save the table
      group.nF <- group.nF[ group.nF >= n.grpF ]  #min number of pts per group (eg 20)
      data.groupF <- dataF[ which(dataF[ ,groupF] %in% names(group.nF)), ]
      sem.mod.groupF <- lavaan::sem( model.formulaF, data=data.groupF , group=groupF )
      if( groupF == "EPA_L1CODE"){ This.FolderF <- paste0( base.folderF, "/SEM_EPA_reg") }
      if( groupF == "ForestType"){ This.FolderF <- paste0( base.folderF, "/SEM_ForestType_reg") }
      save.SEMinfo( sem.mod.groupF , name.of.modF, groupF=TRUE,  dirF=This.FolderF)
      for(ff in 1:length(group.nF) ){
         group.idF <- names(group.nF)[ff]
         data.ffF <- data.groupF[ which(data.groupF[ ,groupF] == group.idF ), ]
         sptl.pts.ffF <- sptl.ptsF[ which(sptl.ptsF@data[ ,groupF] == group.idF ) , ]
         sem.mod.group.ffF <- lavaan::sem( model.formulaF, data=data.ffF  )
         sem.mod.sptlCor.ffF <- lavSpatialCorrect( sem.mod.group.ffF, sptl.pts.ffF@coords[ ,1], sptl.pts.ffF@coords[ ,2])
         sptl.coeffs.ff <- all.SptlCoeffs( sem.mod.sptlCor.ffF )
         if( groupF == "EPA_L1CODE"){
            setwd( paste0( base.folderF, "/SEM_EPA_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
         if( groupF == "ForestType"){
            setwd( paste0( base.folderF, "/SEM_ForestType_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
         if( groupF == "SPEI_factor"){
            setwd( paste0( base.folderF, "/SEM_MultiThresh_SptlCoeffs") )
            write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_Thresh",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
         }
      }
      colnames(group.dfF)[1] <- groupF
      write.csv( group.dfF, file="_nPts_perGroup.csv"  , row.names=FALSE )
   }#close if add Groups (eg EPA)
   
   while( run.mod.againF ){
      mod.numberF <- mod.numberF + 1
      print( paste( "----- Starting " , mod.numberF, "-----"))
      #Check if should add a new variable to the model
      x <- subset(modindices(sem.modF), mi > 2 )
      if( any(which(x$rhs %in% c(y.varF , str.varF ) ) ) ){
         x <- x[ - which(x$rhs %in% c(y.varF , str.varF  ) ) , ]
      }
      my.explan.vars <- c("FLDAGE_log_z", "Tavg_z", "PrcpAvg_z" , 
                          "bdod_0_200_avg_z", "nitrogen_0_200_avg_log_z")
      if( any( which(x$lhs %in% my.explan.vars) ) ){
         x <- x[ - which(x$lhs %in% my.explan.vars) , ]
      }
      if( restrictPathsF == "ForceCovary" ){
         #Rem StrEven --> StrRich
         x.formulas <- paste0( x$lhs,x$op,  x$rhs)
         if( any(x.formulas == paste0(str.varF, "~", str.evenF) ) ){
            x <- x[ - which(x.formulas==paste0(str.varF, "~", str.evenF) ) , ]
         }
         #Rem SpEven --> SpRich
         x.formulas <- paste0( x$lhs,x$op,  x$rhs)
         if( any(x.formulas == paste0(sp.varF, "~", sp.evenF) ) ){
            x <- x[ - which(x.formulas==paste0(sp.varF, "~", sp.evenF) ) , ]
         }
         #Rem SpRich --> SpEven
         x.formulas <- paste0( x$lhs,x$op,  x$rhs)
         if( any(x.formulas == paste0(sp.evenF, "~", sp.varF) ) ){
            x <- x[ - which(x.formulas==paste0(sp.evenF, "~", sp.varF) ) , ]
         }
      }#close if      
      #If still have vars to add
      if( nrow(x) > 0 ){
         x2 <- x[order(x$mi, decreasing=TRUE), ]
         for(ff in 1:nrow(x2)){
            if( x2[ff,"lhs"] %in% c(y.varF, str.varF, sp.varF, str.evenF, sp.evenF)  ){
               if( x2[ff,"lhs"] == y.varF   ){ model.formulaF[1] <- paste( model.formulaF[1] , x2[ff,"rhs"] , sep="+" )  } 
               if( x2[ff,"lhs"] == str.varF ){ model.formulaF[2] <- paste( model.formulaF[2] , x2[ff,"rhs"] , sep="+" )  } 
               if( x2[ff,"lhs"] == sp.varF  ){ model.formulaF[3] <- paste( model.formulaF[3] , x2[ff,"rhs"] , sep="+" )  } 
               if( x2[ff,"lhs"] == str.evenF ){ model.formulaF[4] <- paste( model.formulaF[4] , x2[ff,"rhs"] , sep="+" )  }
               if( x2[ff,"lhs"] == sp.evenF  ){ model.formulaF[5] <- paste( model.formulaF[5] , x2[ff,"rhs"] , sep="+" )  }
               break() 
            }
         }#close ff loop
      }else{
         coeff.tabF <- summary(sem.modF)$pe
         coeff.tabF <- coeff.tabF[ which(coeff.tabF$op == "~"), ]
         coeff.tabF <- coeff.tabF[ order(coeff.tabF$pvalue, decreasing=TRUE), ]
         if( coeff.tabF[1,"pvalue"] > 0.05 ){
            lhs.remF <- coeff.tabF[1,"lhs"]
            rhs.remF <- coeff.tabF[1,"rhs"]
            print( paste( "Remove:", lhs.remF , "~", rhs.remF))
            if( lhs.remF == y.varF ){    rem.row <- 1 }
            if( lhs.remF == str.varF ){  rem.row <- 2 }
            if( lhs.remF == sp.varF ){   rem.row <- 3 }
            if( lhs.remF == str.evenF ){ rem.row <- 4 }
            if( lhs.remF == sp.evenF ){  rem.row <- 5 }
            new.fileF <- strsplit( model.formulaF[rem.row] , split="\\+" )[[1]]  #Split the list
            new.file.part1F <- strsplit( new.fileF, split="~")[[1]]   #Split to get the y var
            new.fileF[1] <- new.file.part1F[2]                        #Re-Add the first variable to the explanatory variable string
            new.fileF <- new.fileF[ -which(new.fileF == rhs.remF)]    #Remove the least significnat variable
            new.formF <- paste(new.fileF , collapse="+")              #Paste together to create RHS of formula
            new.formF <- paste(new.file.part1F[1], new.formF, sep="~") #Create full New formula
            model.formulaF[rem.row] <- new.formF           #Replace Line in existing formula
         }else{
            run.mod.againF <- FALSE
         }
      }

      model.current.listF <- list()
      model.current.listF[[1]] <- model.formulaF
      if( model.current.listF %in% all.models.listF ){
         run.mod.againF <- FALSE
      }else{
         all.models.listF[[mod.numberF-100]] <- model.formulaF
      }
      if( run.mod.againF ){
         sem.modF <- lavaan::sem( model.formulaF, data=dataF  )
         name.of.modF <- paste(model.nameF, mod.numberF , sep=".")
         setwd( paste0( base.folderF, "/SEM_reg") )
         writeLines( model.formulaF, con= paste0( name.of.modF, ".txt"))
         full.summaryF <- EC_addInfo(  modF=sem.modF, name.of.modF, str.varF, sp.varF, 
                                       y.varF, mod.sumryF=full.summaryF )
         setwd( paste0( base.folderF, "/SEM_reg") )
         write.csv(full.summaryF, file= paste0( "_Summary_", model.nameF, ".csv") , row.names=FALSE )
         save.SEMinfo( sem.modF , name.of.modF , dirF=paste0( base.folderF, "/SEM_reg") )
         print( paste( "Mod Pvalue ", round( unlist(summary(sem.modF)$test[[1]][6] ), digits=4 ) ) )
         if( add.sptlCorF ){
            sem.mod.sptlCorF <- lavSpatialCorrect( sem.modF, sptl.ptsF@coords[ ,1], sptl.ptsF@coords[ ,2])
            sptl.coeffs <- all.SptlCoeffs( sem.mod.sptlCorF )
            setwd( paste0( base.folderF, "/SEM_SptlCoeffs") )
            write.csv(sptl.coeffs, file= paste0( name.of.modF, "SptlCoeffs.csv") , row.names=FALSE)
         }
         #SEM by GROUP
         if( !is.null(groupF) ){
            group.nF <- table( dataF[ ,groupF] )
            group.nF <- group.nF[ group.nF >= n.grpF ]  #min number of pts per group (eg 20)
            data.groupF <- dataF[ which(dataF[ ,groupF] %in% names(group.nF)), ]
            sem.mod.groupF <- lavaan::sem( model.formulaF, data=data.groupF , group=groupF )
            if( groupF == "EPA_L1CODE"){ This.FolderF <- paste0( base.folderF, "/SEM_EPA_reg") }
            if( groupF == "ForestType"){ This.FolderF <- paste0( base.folderF, "/SEM_ForestType_reg") }
            save.SEMinfo( sem.mod.groupF , name.of.modF, groupF=TRUE,  dirF=This.FolderF)
            for(ff in 1:length(group.nF) ){
               group.idF <- names(group.nF)[ff]
               data.ffF <- data.groupF[ which(data.groupF[ ,groupF] == group.idF ), ]
               sptl.pts.ffF <- sptl.ptsF[ which(sptl.ptsF@data[ ,groupF] == group.idF ) , ]
               sem.mod.group.ffF <- lavaan::sem( model.formulaF, data=data.ffF  )
               sem.mod.sptlCor.ffF <- lavSpatialCorrect( sem.mod.group.ffF, sptl.pts.ffF@coords[ ,1], sptl.pts.ffF@coords[ ,2])
               sptl.coeffs.ff <- all.SptlCoeffs( sem.mod.sptlCor.ffF )
               if( groupF == "EPA_L1CODE"){
                  setwd( paste0( base.folderF, "/SEM_EPA_SptlCoeffs") )
                  write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
               }
               if( groupF == "ForestType"){
                  setwd( paste0( base.folderF, "/SEM_ForestType_SptlCoeffs") )
                  write.csv(sptl.coeffs.ff, file= paste0( name.of.modF,"_EPA",group.idF,"_SptlCoeffs.csv") , row.names=FALSE)
               }
            }
         }#close if add Groups (eg EPA)
       }
      if( mod.numberF >= 150 ){
         run.mod.againF <- FALSE
      }
    }#Close while loop
} #Close function


##[4] Remove a variable from the model formula (text file) ---------------------
remove.var.from.file <- function( fileF, rowF, rhs.remF ){
   new.fileF <-    strsplit( fileF[rowF] , split="\\+" )[[1]]
   new.fileF <- new.fileF[ -which(new.fileF == rhs.remF)]
   new.formF <- paste(new.fileF , collapse="+")
   fileF[rowF] <- new.formF
   return(fileF)
}   


##[5] Create Summary File - all mods together rather than individual summaries
mnames <- c("ModName","AIC","Mod_Pval", "StrVar","StrCoeff","StrSE","StrPval",  "SpVar" ,"SpCoeff" , "SpSE","SpPval" )
all.mod.summary <- as.data.frame( matrix(NA, nrow=0, ncol=length(mnames)))
colnames(all.mod.summary) <- mnames

EC_addInfo <- function(  modF , mod.nameF , 
                         str.varF , sp.varF , y.varF , 
                         mod.sumryF=all.mod.summary ){
   
   coeff.tabF <- summary(modF)$pe 
   coeff.tabF <- coeff.tabF[ which(coeff.tabF$lhs == y.varF), ]
   str.coeff <- coeff.tabF[ which(coeff.tabF$rhs == str.varF),  ]
   sp.coeff <- coeff.tabF[ which(coeff.tabF$rhs == sp.varF),  ]
   add.infoF <- data.frame( ModName = mod.nameF,
                            AIC = AIC(modF) ,
                            Mod_Pval = unlist( summary(modF)$test[[1]][6] ) ,
                             StrVar = str.varF , 
                            StrCoeff = str.coeff[1,"est"] ,
                            StrSE = str.coeff[1,"se"] ,
                            StrPval = str.coeff[1,"pvalue"] ,
                            SpVar = sp.varF  ,
                            SpCoeff = sp.coeff[1,"est"] ,   
                            SpSE = sp.coeff[1,"se"] ,
                            SpPval = sp.coeff[1,"pvalue"] 
   )
   mod.sumryF <- rbind(mod.sumryF, add.infoF)
   return(mod.sumryF)
}


##[6] Save Main SEM Coefficients to File (from reg model, without sptl correction)
save.SEMinfo <- function( sem.modelF , file.nameF, groupF=FALSE,  dirF=sem.save.folder){
   
   summry.modF <- summary(sem.modelF)
   my.model.infoF <- summry.modF$pe
   n.row <- nrow(my.model.infoF)
   my.model.infoF[ n.row+2, "lhs" ] <- "Pvalue"
   my.model.infoF[ n.row+2, "exo" ] <- summry.modF$test$standard$pvalue
   my.model.infoF[ n.row+3, "lhs" ] <- "n.df"
   my.model.infoF[ n.row+3, "exo" ] <- summry.modF$test$standard$df
   my.model.infoF[ n.row+4, "lhs" ] <- "n.stat"
   my.model.infoF[ n.row+4, "exo" ] <- summry.modF$test$standard$stat
   my.model.infoF[ n.row+5, "lhs" ] <- "AIC"
   my.model.infoF[ n.row+5, "exo" ] <- AIC(sem.modelF)

   if( groupF ){
      group.infoF <-  attributes(summry.modF$test) 
      grp.namesF <- unlist( group.infoF[[2]][2] )
      my.model.infoF[ n.row+6, 1] <- "GroupPosition"
      my.model.infoF[ n.row+6, 2:(1+length(grp.namesF)) ] <- 1:length(grp.namesF)
      my.model.infoF[ n.row+7, 1] <- "GroupName"
      my.model.infoF[ n.row+7, 2:(1+length(grp.namesF)) ] <- grp.namesF
      my.model.infoF[ n.row+8, 1] <- "nData"
      my.model.infoF[ n.row+8, 2:(1+length(grp.namesF)) ] <- summry.modF$data[[2]]
      my.model.infoF[ n.row+9, 1] <- "TestStatistic"
      my.model.infoF[ n.row+9, 2:(1+length(grp.namesF)) ] <- summry.modF$test[[1]][[3]]
   }
   setwd(dirF)
   write.csv( my.model.infoF, file= paste0(file.nameF, ".csv") , row.names=FALSE)
}


##[7] Function to Write final/Best Model to the Main Folder Menu
write.mod.formula <- function( folderF=paste0(m11.folder, "/SEM_reg") ,
                               newNameF="m11.MultiT", 
                               newSaveFolF = base.folder.stats ){
   setwd(folderF) 
   files.xF <- list.files()
   files.xF <- files.xF[ -grep("Summary", files.xF) ]
   files.xF <- files.xF[ -grep("csv", files.xF) ]
   if( any( grep("Output", files.xF)) ){
      files.xF <- files.xF[ -grep("Output", files.xF) ]
   }
   last.xF<- tail(files.xF, 1)
   x.fileF <- readLines(last.xF)
   setwd( newSaveFolF )
   writeLines( x.fileF, con= paste0( newNameF, ".txt") )
}

##[8] Get Sptl Adjusted Coefficients 
each.SptlCoeff <- function(modF, param.numF){
   my.paramsF <- modF$parameters[[ param.numF ]]
   my.params.listF <- strsplit(my.paramsF$Parameter, "~")
   for(ii in 1:length(my.params.listF)){
      rhs <- my.params.listF[[ii]][2]
      my.paramsF[ii,"rhs"] <- rhs
      lhs <- my.params.listF[[ii]][1]
      my.paramsF[ii,"yvar"] <- lhs
   }
   my.paramsF$Parameter <- NULL
   rownames(my.paramsF) <- 1:nrow(my.paramsF)
   return(my.paramsF)
}


##[9] Extract Spatial Coefficients to a df 
all.SptlCoeffs <- function( modF ){
   sptl.coeffs.yvarF <- each.SptlCoeff(modF,  1) #yvar
   sptl.coeffs.strdivF <- each.SptlCoeff(modF,  2) #strDiv
   sptl.coeffs.sprichF <- each.SptlCoeff(modF,  3) #spRich
   sptl.coeffsF <- rbind(sptl.coeffs.yvarF, sptl.coeffs.strdivF , sptl.coeffs.sprichF)
   return(sptl.coeffsF)
}


##[10] Extract Morans I values 
get.morans <- function(modF, param.numF ){
   moran.dfF <- modF$Morans_I[[ param.numF]]
   moran.dfF$yvar <- names( modF$Morans_I )[ param.numF]
   return(moran.dfF)
}


##[11] Correct for Spatial Correlation
#Functions Written by Jarret Byrnes 
#https://github.com/jebyrnes/spatial_correction_lavaan

lavResidualsY <- function(object,
                          ynames = lavNames(object, "ov.nox"),
                          xnames = lavNames(object, "ov.x")){
   
   pred <- lavPredictY(object,
                       ynames = ynames,
                       xnames = xnames) |> 
      as.data.frame()
   d <- inspect(object, "data") |> 
      as.data.frame()
   r <- lapply(names(pred), function(x){
      pred[[x]] - d[[x]]
   })
   res <- do.call(cbind, r) |> 
      as.data.frame()
   names(res) <- names(pred)
   res
}
lavSpatialCorrect <- function(obj, xvar, yvar, alpha=0.05){
   require(lavaan)
   require(ape)
   
   # first, get the residuals from the model
   resids <- lavResidualsY(obj)
   # get only endogenous variables
   # no longer need
   #resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
   # make a distance matrix
   distMat <- as.matrix(dist(cbind(xvar, yvar)))
   # invert this matrix for weights
   distsInv <- 1/distMat
   diag(distsInv) <- 0
   morans_i <- lapply(resids,  function(x){
      mi <- Moran.I(c(x), distsInv)
      if(mi$p.value>alpha){
         mi$n.eff <- nrow(resids) #don't correct sample size
      }else{
         #large sample size approximation
         mi$n.eff <- nrow(resids)*(1-mi$observed)/(1+mi$observed)
      }
      #return the list
      mi
   })

   # get the vcov matrix
   v <- diag(vcov(obj))
   n <- nrow(resids)
   # using new sample sizes, for each variable, calculate new Z-scores
   params <- lapply(names(morans_i), function(acol){
      idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
      idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
      v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
      ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                         n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
      ret[["Z-value"]] <- ret$Estimate/ret$Std.err
      ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
      ret
   })
   names(params) <- names(morans_i)
   mi <- lapply(morans_i, function(m) {
      data.frame(observed=m$observed, expected=m$expected, 
                 sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
   })
   return(list(Morans_I = mi, parameters=params)) 
}


