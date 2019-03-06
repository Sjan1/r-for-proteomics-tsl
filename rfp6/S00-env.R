#source("http://www.bioconductor.org/biocLite.R")
#library("BiocInstaller")
#biocLite("MSnbase")
#biocLite("limma")
#install.packages("devtools")

# update all packages
#BiocManager::install()

## to install packages
#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("MSnbase")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biobroom", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MSnID", version = "3.8")

#install.packages("tidyverse")

library("MSnbase")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("magrittr")
library("readr")
library("MSnID")
library("broom")
library("biobroom")
