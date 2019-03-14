##  install package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msdata", version = "3.8")

library("MSnbase")
#library("msdata")
library("magrittr")

#fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)[2]
#basename(fl)
#data_prof <- readMSData(fl, mode = "onDisk", centroided = FALSE)

## for the data on the disc
mzml <- "data/mzml"
filenames <- c("js190313_P6-pep12_hcd-sid",
               "js190313_pep10_hcd-sid",
               "js190313_P6-pep12_cid",
               "js190313_pep11_cid",
               "js190313_pep10_cid",
               "js190313_pep10_hcd",
               "js190313_P6-pep12_hcd")
afile <- NULL
mzml_files <- NULL
for (i in filenames){
    afile <- file.path(mzml, paste(i, "mzML", sep = "."))
    mzml_files <- c(mzml_files,afile)
}
mzml_files
fl <- mzml_files[1]
fl
data_prof <- readMSData(fl, mode = "onDisk", msLevel=1, centroided = FALSE)


## We next extract the profile MS data for the [M+H]+ adduct of serine
## with the expected m/z of 106.049871. 
## We thus filter the data_prof object using an m/z range
## containing the signal for the metabolite and a retention time window
## from 175 to 187 seconds corresponding to the time
## when the analyte elutes from the LC.


## Define the mz and retention time ranges
serine_mz <- 928.9807 #968.9591
mzr <- c(serine_mz - 0.02, serine_mz + 0.02)
rtr <- c(5000, 8000)

## Filtering the object
serine <- data_prof %>%
  filterRt(rtr) %>%
  filterMz(mzr)
## We can now plot the profile MS data for serine.

plot(serine, type = "XIC")
abline(h = serine_mz, col = "red", lty = 2)
