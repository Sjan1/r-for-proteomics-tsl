rm(list = ls());
source("S00-env.R")
library("rtslprot")
library("dplyr")
library("msmsTests")

## READ THE EXPERIMENT TABLE that describes:
## (i) samples
## (ii) measurements
## (iii) searches
## (iv) experimental factors

mzid <- "mascot_fdr1pc"
mzid <- "mascot"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)
## APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()

## DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples) that will be combined
## all this have to be described in the experiment table
fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)
etab

## READ INTO MSnSets
## This chunk reads each entry of the experimental design etab
## and creates the corresponding MSnSet from the mzid files.
## The MSnSets are returned as a list of MSnSets.
msnl <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(biorep == .etab[["biorep"]],
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
           select(name)
  mzid_files <- file.path(mzid, paste(filenames[[1]], "mzid",
                                      sep = "."))
  ## make MSnSet (choose peptide or PSMs - two functions in rtslprot)
    e <- rtslprot:::make_pep_MSnSet(mzid_files,
                                  fdr = 0.05,
                                  level = "peptide",
                                  fcol = "pepSeq")
  e@processingData@files <- mzid_files
  sampleNames(e) <- paste(.etab, collapse = "_")
  e <- updateFvarLabels(e, sampleNames(e))
  if (validObject(e))
    return(e)
})

## COMBINE MSnSets IN ONE
## Convert the list of MSnSets to a single MSnSet where
## each element of the list becomes a columns (sample)
## and update the experimental design
e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl)) {
  e <- MSnbase::combine(e, msnl[[i]])
}
rownames(etab) <- sampleNames(e)
pData(e) <- etab


## DEAL WITH NAs ON PEPIDE LEVEL
e <- impute(e, method = "zero")


## keep the original
saveRDS(e,"e.rds")
## open when needed
e <- readRDS("e.Rds")
## save as
saveRDS(e,"e_mascot_fdr1pc.rds")




## 20-03-2019: COMBINE PEPTIDES INTO PROTEINS
## Make vector with all accessions matching every peptide
## Append the vector to the feature data
## Step through all the features (rows and columns)
## and copy accessions when found.

e <- readRDS("e_mascot_fdr1pc.rds")

## concatenate all  accessions - NEW
i <- grep("accession\\.", fvarLabels(e))    # e is a peptide-level MSnSet
k <- apply(fData(e)[, i], 1,
           function(x) unique(na.omit(as.character(x))))
fData(e)$nprots <- lengths(k)
#fData(e)$accession <- sapply(k, paste, collapse = ";") # Laurent's suggestion
fData(e)$accession <- sapply(k, paste) # But when not collapsed, we then easily make a list
l <- as.list(fData(e)$accession)        # the list is nedded for combineFeatures

#save modified MSnSet
saveRDS(e,"e.rds")

eprot_m <- combineFeatures(e, groupBy = l,
                           fun = "sum",redundancy.handler = "multiple")

eprot_u <- combineFeatures(e, groupBy = fData(e)$accession,
                           fun = "sum", redundancy.handler = "unique")



## save results as RDS
saveRDS(eprot,"eprot.rds")
saveRDS(eprot_m,"eprot.rds")
saveRDS(eprot_u,"eprot.rds")
## read saved results
eprot <- readRDS("eprot.rds")


stop("never mind the error, execution stops here")
##################################################



## A NEW FUNCTION TEST
## WHY: examples of two different results from mzid vignette

## read msnid as data.table
td <- as(msnid,"data.table")
colnames(td)
td[,c(24,26)]
td[td$accession=="04266",26]

## read msnid as MSnSet
ts <- as(msnid,"MSnSet")
class(ts)
colnames(ts)
fvarLabels(ts)
head(fData(ts),20)
fData(ts)$accession


## THE SOLUTION
## 12th Apr 2019 from Lauent in email
## Test data, to be run in rfp6
msnid <- MSnID()
msnid <- read_mzIDs(msnid,"mascot_fdr1pc/PDLP1_C1_2_1_1_120416.mzid")
msnid

## Better to test samples (merged fractions)
## We can use more files from experiment design table - mzid_files
## filtered and grouped mzid_files from the experimental design table
## using only the firsst part of the S01_code.R, without building MSnSet

mzid <- "mascot_fdr1pc"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)
## APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()

## DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples) that will be combined
## all this have to be described in the experiment table
fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)
etab


mzid_files <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(biorep == .etab[["biorep"]],
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
    select(name)
  mzid_files <- file.path(mzid,paste(filenames[[1]], "mzid",
                                     sep = "."))
  ## make a list of files
  return(mzid_files)
})
## add names
names(mzid_files) <- apply(etab, 1,
                           function(.etab) paste(.etab, collapse = "_"))

## then from 'etab' and 'exp' we get
mzid_files

##  Now we can compare the proteins we had problem earlier -
##  combined from several mzids
##  samples P_I_3 (list #23) and P_I_1 (list #32)
msnid <- MSnID()
#msnid <- read_mzIDs(msnid,mzid_files[["P_I_3"]])
msnid <- read_mzIDs(msnid,mzid_files[["P_I_1"]])
msnid

## The new function
#as_MSnSet <- function(x, fcol = NULL) {
#  td <- as(x, "data.table")
#  td$e <- 1
#  ## Create a PSM-level MSnSet
#  x <- readMSnSet2(td, ecol = which(colnames(td) == "e"))
#  if (!is.null(fcol)) {
#    ## If there's an fcol, combine at that level by summing PSM counts
#    stopifnot(fcol %in% fvarLabels(x))
#    x <- combineFeatures(x, fcol = fcol, fun = sum)
#  }
#  return(x)
#}

## Counts at the PSM level
psm <- rtslprot::as_MSnSet(msnid)

## Counts at the peptide level - not here that peptide is K.ASVGFK.A while
## pepSeq is ASVGFK in the data table
pep2 <- rtslprot::as_MSnSet(msnid, fcol = "peptide")
pep1 <- rtslprot::as_MSnSet(msnid, fcol = "pepSeq")

## Counts at the protein (accession) level
prot <- rtslprot::as_MSnSet(msnid, fcol = "accession")

dim(fData(psm))
dim(fData(pep1))
dim(fData(pep2))
dim(fData(prot))

length(unique(fData(pep1)$pepSeq))
length(unique(fData(pep2)$pepSeq))
length(unique(fData(pep2)$peptide))

dim(fData(pep2)[,c("pepSeq","peptide")])



## check the accesions == "04266"

table(fData(pep1)$accession=="04266")
rownames(fData(psm)[fData(psm)$accession=="04266",])
unique(fData(psm)[fData(psm)$accession=="04266",c(23,25)])

## WHY THE TWO FOLLOWING LINES PRODUCE DIFFERENT RESULTS?
unique(fData(psm)[fData(psm)$accession=="04266",25])
rownames(fData(pep1)[fData(pep1)$accession=="04266",])
rownames(fData(pep2)[fData(pep2)$accession=="04266",])

## check the accesions == "20144"