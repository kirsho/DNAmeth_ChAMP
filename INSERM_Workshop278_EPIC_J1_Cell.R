##############################################################################
# INSERM278 Workshop DNA methylation analysis workflow with ChAMP package
##############################################################################
# 2024 May 14th-16th
# Teachers :
# Guillaume Velasco : guillaume.velasco@u-paris.fr
# Olivier Kirsh : olivier.kirsh@u-paris.fr

# This is a R script to perform Illumina Infinium DNA methylationEPICv1 analysis
# with R-4.1.x and ChAMP 1.29.1 & ChAMPdata 1.31.1 packages
# Add your comment and notes after a # mark
# Unmarked lines are R commands 
# Can be executed with Run button or with Ctrl + ENTER 
# Always check that the previous command is correctly executed and the prompt ( >) is ready

# A reminder to R 
# https://intro2r.com
#

##############################################################################
## IFB Rstudio session
##############################################################################

## Open Ondemand
## https://ondemand.cluster.france-bioinformatique.fr/pun/sys/dashboard

## Start an interactive Rsudio session
### R version: 4.1.1
### Account: tp_tps_epic_dname_165656
### Number of CPUs: 6
### Amount of memory: 36G
### Number of GPUs: 0
### Number of hours: 5

## Connect to Rstudio session

## Open an new Rscript

##############################################################################
# Set Path2Files & Directories
##############################################################################

## Workshop shared folder should be
## /shared/projects/tp_tps_epic_dname_165656/LearnerName/

## Example on IFB rstudio
### Create a learner folder in tp_tps_epic_dname_165656 project
Path2TP <- "/shared/projects/tp_tps_epic_dname_165656"
learner.folder <- "Mag_beta" # put your own name! 

Path2dir <- file.path(Path2TP, learner.folder)

print(Path2dir)

dir.create(Path2dir)

### Copy idat files to learner folder
Path2iDat_tocopy <- file.path(Path2TP,"all_idat_atelier")
print(Path2iDat_tocopy)

### create idat folder
Path2iDat <- file.path(Path2dir, "idat_files")
dir.create(Path2iDat)
print(Path2iDat)

### list files to copy
idat_files <- list.files(Path2iDat_tocopy)
idat_files

file.copy(from = paste0(Path2iDat_tocopy,"/", idat_files),   # Copy files
          to = paste0(Path2iDat,"/", idat_files))

list.files(Path2iDat)


### Path to metadata file
PD <- file.path(Path2TP,"Pheno_WS.txt")



## Exemple on IFB Cluster with a private project folder already existing
### Path to directories
Path2dir <- "/shared/projects/remeth/AtelierEPIC"
Path2iDat <- "/shared/projects/remeth/AtelierEPIC/all_idat_atelier/"

### Path to metadata file
PD <- "/shared/projects/remeth/AtelierEPIC/Pheno_WS.txt"


#############################################################################
# Set working directory
############################################################################# 

# Check actual WD
getwd()

# Set WD if needed
setwd(Path2dir)
getwd()

# Set Explorer panel to WD
rstudioapi::filesPaneNavigate(Path2dir)
################################################


##############################################################################
## Check ChAMP & ChAMPdata packages installation !!! MUST READ !!!
##############################################################################

## Check that # ChAMP 2.29.1 & ChAMPdata 2.31.1 are installed
### Explore Packages Panel (lower rightpanel)
### If not Follow instructions on ChAMP_Install.R script


## Rename `probe.features.epic.rda` to `probe.features.epicv1.rda`
### This can be done "by hand"

### go in "/shared/home/<projectname>/R/x86_64-conda-linux-gnu-library/4.1/ChAMPdata/data"
### rename probe.features.epic.rda to probe.features.epicv1.rda

### This can be done  with CLI with R

rstudioapi::filesPaneNavigate("~/R/x86_64-conda-linux-gnu-library/4.1/ChAMPdata/data") # Go to the folder

file.rename("~/R/x86_64-conda-linux-gnu-library/4.1/ChAMPdata/data/probe.features.epic.rda",
            "~/R/x86_64-conda-linux-gnu-library/4.1/ChAMPdata/data/probe.features.epicv1.rda")

rstudioapi::filesPaneNavigate(Path2dir) # return to default folder



##############################################################################
## ChAMP & ChAMPdata packages Documentation
##############################################################################

# https://bioconductor.org/packages/release/bioc/html/ChAMP.html
# https://github.com/YuanTian1991/ChAMP


## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()
#champ.norm()
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()
#champ.DMP()
#champ.DMR()


##############################################################################
# R session set UP
##############################################################################

## Check you library paths (You should get 2)
.libPaths()

## Load required libraries
library(dplyr) # Load dplyr package (object handling & mods)

library(ChAMP ,                                                # Load ChAMP
        lib.loc = "~/R/x86_64-conda-linux-gnu-library/4.1/")   # Specify the path

library(ChAMPdata,
        lib.loc = "~/R/x86_64-conda-linux-gnu-library/4.1/") 



## Libraries can be loaded manually on Packages panel





##############################################################################
## Prepare & reformat metadata file
##############################################################################

###################
### Read "Pheno_WS.txt" file
pD_EPIC_WS <- read.table(PD,  
                         header=T,
                         sep="\t",
                         stringsAsFactors=F) 

#### Check Size
dim(pD_EPIC_WS)

#### Check import (first 6 rows)
head(pD_EPIC_WS)

#### Check entire object
pD_EPIC_WS

#### Sample Names
pD_EPIC_WS$Sample_Name

#### Colnames
colnames(pD_EPIC_WS)


###################
### Format pD_files EPIC for ChAMP (change col names)

#### Extract GSM (Sentrix )ID
Sentrix_ID <- unlist(lapply(basename(as.character(pD_EPIC_WS$Basename)),
                         function(x) strsplit(gsub("_2",":2",x),":")[[1]][[2]]))

Sentrix_ID <- unlist(lapply(basename(as.character(Sentrix_ID)),
              function(x) strsplit(gsub("_R",":R",x),":")[[1]][[1]]))

#### Extract (Sentrix) Position
Sentrix_Position = unlist(lapply(basename(as.character(pD_EPIC_WS$Basename)),
                               function(x) strsplit(gsub("_R",":R",x),":")[[1]][[2]]))



#### Reformat pD_WS, Change Colnames and some infos
ChAMP_pD_WS <- data.frame(Sample_Name = pD_EPIC_WS$Sample_Name, 
                          Sample_Plate = pD_EPIC_WS$Chips, 
                          Sample_Group = pD_EPIC_WS$mutation,
                          Sample_Group_2 = pD_EPIC_WS$Sex,
                          Sample_GSE_ID = "",
                          Sample_Status = pD_EPIC_WS$Sample_status,
                          Mutation_Status = pD_EPIC_WS$Mutation_status,
                          Pool_ID = "",
                          Project = "",
                          Sample_Well = "",
                          Sentrix_ID = Sentrix_ID,  
                          Sentrix_Position = Sentrix_Position,
                          Basename = pD_EPIC_WS$Basename,
                          Tissue = pD_EPIC_WS$tissue,
                          stringsAsFactors = F)

#### Check entire updated object
ChAMP_pD_WS



#############################################################################
# champ.load()
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load()  # >>>>>>>>>>>>>
#champ.QC()
#champ.norm()
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()
#champ.DMP()
#champ.DMR()


## DNA methylation analysis J1  : cell line -----------------------------------

### Filters ChAMP_CSV_WS
### keeping all Tissue == Cell

ChAMP_f_pD_WS <- subset(ChAMP_pD_WS, 
                       Tissue =="Cell")

dim(ChAMP_f_pD_WS)
ChAMP_f_pD_WS



### Create csv file of pD define in Samples Selections, in directory of .idat files

write.table(ChAMP_f_pD_WS,
            file = paste0(Path2iDat,"/ChAMP_f_pD_WS.csv"),
            row.names = F,
            quote = F, 
            sep=",")

### Check "idat" directory
dir(Path2iDat)


### Load idat files
myLoad_Cell <- champ.load(Path2iDat,
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,         
                     force=FALSE,
                     arraytype="EPICv1")    



# Explore myLoad_Cell
myLoad_Cell
is(myLoad_Cell)
is(myLoad_Cell$beta)

# View myLoad_Cell data
myLoad_Cell$beta
myLoad_Cell$pd$Sample_Status

myLoad_Cell$intensity
myLoad_Cell$pd



### Remove .
file.remove(paste0(Path2iDat, "/ChAMP_f_pD_WS.csv"))


#############################################################################
# champ.QC()
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   # >>>>>>>>>>>>>
#champ.norm()
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()
#champ.DMP()
#champ.DMR()

### Quality controls
champ.QC(beta = myLoad_Cell$beta,
         pheno = myLoad_Cell$pd$Sample_Status,
         mdsPlot = TRUE,
         densityPlot = TRUE,
         dendrogram = TRUE,
         PDFplot = TRUE,
         Rplot = TRUE,
         Feature.sel = "None",
         resultsDir = "./CHAMP_QCimages_Cell/")


### Explore graphs in Plots panel &/or in CHAMP_QCimages/ folder


### reset Plots panel
dev.off()


#############################################################################
# champ.norm()
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   
#champ.norm() # >>>>>>>>>>>>>
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()
#champ.DMP()
#champ.DMR()


## Normalization
myNorm_Cell <-champ.norm(beta = myLoad_Cell$beta,
                         rgSet = myLoad$rgSet,
                         mset = myLoad$mset,
                         resultsDir = "./CHAMP_Normalization_Cell/",
                         method = "BMIQ",  # ("PBC","BMIQ","SWAN")                  
                         plotBMIQ = TRUE,
                         arraytype = "EPICv1",
                         cores = 5)   # use n-1 core (with n = reserved core)

# It takes some time (5 to 10')... annotate your code.

# Check myNorm_Cell object
is(myNorm_Cell)
head(myNorm_Cell)
dim(myNorm_Cell)



## You can save it
# write.table(myNorm_Cell,file="beta_EPIC_Cell.tsv", row.names=T, quote=F, sep="\t")



## Methylation Data Visualization 

### With Rbase graphs
#### Methylation Data extraction as dataframe
methdata = data.frame(myNorm_Cell)

head(methdata)
summary(methdata)

#### Methylation Data boxplot
boxplot(methdata,cex.axis = 0.8)

#### Reset Plots panel
dev.off()



### With GGPLOT2 graph
library(ggplot2)
library(reshape2)

#### Methylation Data extraction as vector
methdata <- setNames(reshape2::melt(methdata), c("Group","Methylation"))
head(methdata)

p <- ggplot(methdata,   
            aes(x=Group, 
                y=Methylation,
                fill=Group, 
                color=Group)) + 
    geom_violin(trim=FALSE) 

p 


# Add layers and group colors
p +
  stat_summary(fun=mean, geom="crossbar", width=0.2, color="black") +
  scale_fill_manual(values=c("#00FF00","#00FF00", "#00FF00","#FF0000", "#FF0000","#FF0000")) +
  scale_color_manual(values=c("#00FF00","#00FF00", "#00FF00","#FF0000", "#FF0000","#FF0000")) +
  theme_classic() + 
  theme(axis.text=element_text(size=6))

rm(methdata) 
dev.off()


#############################################################################
# champ.SVD() / COMBAT  / REFBASE
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   
#champ.norm() 
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect # >>>>>>>>>>>>>
#champ.QC()
#champ.DMP()
#champ.DMR()

## Check metadata and $slide column
myLoad_Cell$pd

## SVD
### rename myNorm_Cell
myNorm <- myNorm_Cell

champ.SVD(beta = myNorm,
          pd = myLoad_Cell$pd,
          resultsDir = "./CHAMP_SVDimages_Cell/")

### Explore graph

### reset Plots panel
rm(myNorm)
dev.off()



## ChAMP  COMBAT

### If necessary, Combat correction for batch effect (450K vs EPIC Samples, samples prep, date etc...).
### Not always applicable
## Only one slide, not applicable here

# myNorm_combat <- champ.runCombat(beta = myNorm_Cell,
#                    pd=myLoad_Cell$pd,
#                    variablename = "Sample_Group",
#                    batchname = c("Array"),
#                    logitTrans = TRUE)





#############################################################################
# champ.QC() 
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   
#champ.norm() 
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()   # >>>>>>>>>>>>>
#champ.DMP()
#champ.DMR()

champ.QC(beta = myNorm_Cell,
         pheno = myLoad_Cell$pd$Sample_Status,
         mdsPlot = TRUE,
         densityPlot = TRUE,
         dendrogram = TRUE,
         PDFplot = TRUE,
         Rplot = TRUE,
         Feature.sel = "None",
         resultsDir = "./CHAMP_QCimages_Norm_Cell/")


### Explore graphs in Plots panel &/or in CHAMP_QCimages_Norm folder

### reset Plots panel
dev.off()




#############################################################################
# champ.DMP() 
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   
#champ.norm() 
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()   
#champ.DMP()  # >>>>>>>>>>>>>
#champ.DMR()

## Check Status
unique(myLoad_Cell$pd$Sample_Status)

## Set as factor, not mandatory if you specify compare.group
### ORDER MATTERS

# myLoad_Cell$pd$Sample_Status <- factor(myLoad_Cell$pd$Sample_Status,
#                                       levels = c("CCL_WT",  "CCL_DKO"))

## Call DMPs
myDMP_CCL <- champ.DMP(beta = myNorm_Cell,
                       pheno = myLoad_Cell$pd$Sample_Status,
                       compare.group = c("CCL_WT","CCL_DKO"), 
                       adjPVal = 0.05,
                       adjust.method = "BH",
                       arraytype = "EPICv1")



## Explore outputs
head(myDMP_CCL$CCL_WT_to_CCL_DKO)
dim(myDMP_CCL$CCL_WT_to_CCL_DKO)


## write.table(myDMP_CCL$CCL_WT_to_CCL_DKO,"myDMP_CCL.tsv", row.names=T, quote=F, sep="\t")

## Searching specific gene in DMP list or probes

myDMP_CCL$CCL_WT_to_CCL_DKO[myDMP_CCL$CCL_WT_to_CCL_DKO$gene=="MAEL",]
myDMP_CCL$CCL_WT_to_CCL_DKO[myDMP_CCL$CCL_WT_to_CCL_DKO$gene=="SYCE1",]
myDMP_CCL$CCL_WT_to_CCL_DKO["cg16208053",]


## Filter DMPs with thresholds
head(myDMP_CCL$CCL_WT_to_CCL_DKO)
dim(myDMP_CCL$CCL_WT_to_CCL_DKO)

### All_DM Delta 20% & FDR < 1%                                
FDR <- 0.01
DELTAB <- 0.2

myDMP20_CCL <- subset(myDMP_CCL$CCL_WT_to_CCL_DKO,
                    abs(myDMP_CCL$CCL_WT_to_CCL_DKO$deltaBeta)> DELTAB & myDMP_CCL$CCL_WT_to_CCL_DKO$adj.P.Val < FDR)

### All_DM Delta 60% & FDR < 1%
FDR <- 0.01
DELTAB <- 0.6

myDMP60_CCL <- subset(myDMP_CCL$CCL_WT_to_CCL_DKO,
                    abs(myDMP_CCL$CCL_WT_to_CCL_DKO$deltaBeta)> DELTAB & myDMP_CCL$CCL_WT_to_CCL_DKO$adj.P.Val < FDR)


### Comapre table sizes
dim(myDMP_CCL$CCL_WT_to_CCL_DKO)
dim(myDMP20_CCL)
dim(myDMP60_CCL)

# Exercice : Count Hyper and hypo methylated probes



## Volcano plot
library(EnhancedVolcano)

png("volcano_DMP_Cell.png",width = 250, height = 200, units='mm',res = 600)
p_volc <- EnhancedVolcano(myDMP_CCL$CCL_WT_to_CCL_DKO,
                          lab = "",
                          x = "deltaBeta",
                          y = "adj.P.Val",
                          pCutoff = 0.01,
                          FCcutoff = 0.20,
                          xlim = c(-1,1), 
                          xlab = bquote(~Delta~"Beta"),
                          title = "DMP_CCL",
                          legendLabels = c("NS","dBETA","p.adj < 0.05","p.adj < 0.01 & abs(dBeta) > 0.2"))
p_volc
dev.off()


## Look at the graph, comments.


#############################################################################
# champ.DMR()
############################################################################# 

## Champ Pipeline in a nutshell
#champ.load() 
#champ.QC()   
#champ.norm() 
#champ.SVD() champ.runCombat() champ.refbase() # In case of batch effect
#champ.QC()   
#champ.DMP()  
#champ.DMR() # >>>>>>>>>>>>>

## Does it make sens?

## Call DMR
myDMR_CCL <- champ.DMR(beta = myNorm_Cell,
                         pheno = myLoad_Cell$pd$Sample_Status,
                         compare.group = c("CCL_WT","CCL_DKO"),
                         arraytype = "EPICv1",
                         method = "Bumphunter",
                         minProbes=7,
                         adjPvalDmr=0.05,
                         cores = 5)   # use n-1 core (with n = reserved core)


myDMR_CCL %>% head()
myDMR_CCL$BumphunterDMR %>% dim()

## Exercice : Explore output as in the DMP section








