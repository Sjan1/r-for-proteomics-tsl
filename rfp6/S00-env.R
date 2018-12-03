#source("http://www.bioconductor.org/biocLite.R")
#library("BiocInstaller")
#biocLite("MSnbase")
#biocLite("limma")
#install.packages("devtools")

# update all packages
#BiocManager::install()

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("MSnbase")

library("MSnbase")
library("ggplot2")
library("dplyr")
library("magrittr")
library("readr")
library("MSnID")
