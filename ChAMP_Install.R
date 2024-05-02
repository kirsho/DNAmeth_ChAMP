##############################################################################
# INSERM278 Workshop DNA methylation analysis workflow with ChAMP package
##############################################################################
# 2024 May 14th-16th
# Teachers :
# Guillaume Velasco : guillaume.velasco@u-paris.fr
# Olivier Kirsh : olivier.kirsh@u-paris.fr

##############################################################################
# R package Documentation
##############################################################################

## Intro 2 R
## https://intro2r.com

## R
## https://cran.r-project.org/

## Bioconductor
## https://bioconductor.org/packages/release/bioc/html/ChAMP.html

## GitHub
## https://github.com/YuanTian1991/ChAMP



##############################################################################
## Software requirement
##############################################################################

# R version 4.1. to 4.1.3

# ChAMP 1.29.1
# ChAMPdata 1.31.1


##############################################################################
## IFB Rstudio session
##############################################################################

## Open Ondemand
## https://ondemand.cluster.france-bioinformatique.fr/pun/sys/dashboard

## Start an interactive Rsudio session
### R version: 4.1.1
### Number of CPUs: 6
### Amount of memory: 36G
### Number of GPUs: 0

## Connect to Rstudio session

## Install Last ChAMP package version (Do not update anything)
### Install "ChAMP" # version 2.29.1
devtools::install_github("YuanTian1991/ChAMP") # add , force = T if it's not working
# Click in the console & select nothing, press ENTER (no updates)

# run if required & rerun previous line
install.packages("mime")

### Install "ChAMPdata" # version 2.31.1
devtools::install_github("YuanTian1991/ChAMPdata")
# Click in the console & select nothing, press ENTER (no updates)




##############################################################################
# from scratch on your computer
##############################################################################
## Install R4.1.x
## Install Rtools
###https://cran.r-project.org/bin/windows/Rtools/rtools40.html



### Install "Devtools" package 
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools", dependencies = T)


### Install this packages packages if devtools does not install properly (non zero status exit message)
if (!requireNamespace("usethis", quietly=TRUE))
  install.packages("usethis", dependencies = T)

if (!requireNamespace("shiny", quietly=TRUE))
  install.packages("shiny", dependencies = T)

if (!requireNamespace("htmlwidgets", quietly=TRUE))
  install.packages("htmlwidgets", dependencies = T)


### Install Last ChAMP package version (Do not update anything)
### Install "ChAMP" # version 2.29.1
devtools::install_github("YuanTian1991/ChAMP")
# select nothing, press ENTER (no updates)


### Install "ChAMPdata" # version 2.31.1
devtools::install_github("YuanTian1991/ChAMPdata")
# select nothing, press ENTER (no updates)


