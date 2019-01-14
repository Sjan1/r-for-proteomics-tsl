#source("http://www.bioconductor.org/biocLite.R")
#library("BiocInstaller")
#biocLite("MSnbase")
#biocLite("limma")
#install.packages("devtools")

# update all packages
#BiocManager::install()

#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("MSnbase")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")

library("MSnbase")
library("ggplot2")
#library("dplyr")
library("tidyverse")
library("magrittr")
library("readr")
library("MSnID")