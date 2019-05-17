rm(list = ls());
source("S00-env.R")
library("rtslprot")
library("dplyr")
#library("msmsTests")

## 1 - read THE EXPERIMENT TABLE describing the experimental design:
## (i) samples
## (ii) measurements
## (iii) searches
## (iv) experimental factors

mzid <- "mascot_fdr1pc"
#mzid <- "mascot"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)

## 2 - APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()


## 3 - DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples) that will be combined
## all this have to be described in the experiment table
fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)
etab


## 4.1 - READ PATHS TO SEARCH DATA
## This chunk reads each row of the experimental design
## and creates a list of paths to mzid files.
## useful for re-process only subset of samples
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
## thus from 'exp' and 'etab' we get
mzid_files


## 4.2 - READ PATHS TO SEARCH DATA
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

# ## 5 - MAKE MSnSets
# ## make MSnSet (choose peptide or PSMs - two functions in rtslprot)
#   e <- rtslprot:::make_pep_MSnSet(mzid_files,
#                                   fdr = 0.05,
#                                   level = "peptide",
#                                   fcol = "pepSeq")

## 6 - USING as_MSnSet (choose protein, peptide or PSM)
## READ SEARCH DATA
  msnid <- MSnID()
  msnid <- read_mzIDs(msnid, mzid_files)
  
  e <- rtslprot::as_MSnSet(msnid)
  #e <- rtslprot::as_MSnSet(msnid, fcol = "pepSeq")
  #e <- rtslprot::as_MSnSet(msnid, fcol = "accession")
  
  e@processingData@files <- mzid_files
  sampleNames(e) <- paste(.etab, collapse = "_")
  e <- updateFvarLabels(e, sampleNames(e))
  if (validObject(e))
    return(e)
})


## 7 - COMBINE LIST OF MSnSets INTO ONE
## Convert the list of MSnSets to a single MSnSet where
## each element of the list becomes a columns (sample)
## and update the experimental design
e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl)) {
  e <- MSnbase::combine(e, msnl[[i]])
}
rownames(etab) <- sampleNames(e)
pData(e) <- etab


## 8 - DEAL WITH NAs ON PEPIDE LEVEL
e <- impute(e, method = "zero")

## 9 - SAVE IT
## keep the original
saveRDS(e,"e.rds")
## open when needed
e <- readRDS("e.Rds")
## save as
saveRDS(e,"e_mascot_fdr1pc.rds")
e <- readRDS("e_mascot_fdr1pc.rds")

# ## 20-03-2019: COMBINE PEPTIDES INTO PROTEINS
# ## Make vector with all accessions matching every peptide
# ## Append the vector to the feature data
# ## Step through all the features (rows and columns)
# ## and copy accessions when found.
#
# e <- readRDS("e_mascot_fdr1pc.rds")
#
# ## concatenate all  accessions - NEW
# i <- grep("accession\\.", fvarLabels(e))    # e is a peptide-level MSnSet
# j <- grep("pepSeq\\.", fvarLabels(e))    # e is a peptide-level MSnSet
# k <- apply(fData(e)[, i], 1,
#            function(x) unique(na.omit(as.character(x))))
# fData(e)$nprots <- lengths(k)
# #fData(e)$accession <- sapply(k, paste, collapse = ";") # Laurent's suggestion
# fData(e)$accession <- sapply(k, paste) # But when not collapsed, we then easily make a list
# l <- as.list(fData(e)$accession)        # the list is nedded for combineFeatures
# l 
# #save modified MSnSet
# saveRDS(e,"e.rds")
# 
# eprot_m <- combineFeatures(e, groupBy = l,
#                            fun = "sum",redundancy.handler = "multiple")
# 
# eprot_u <- combineFeatures(e, groupBy = fData(e)$accession,
#                            fun = "sum", redundancy.handler = "unique")

# ## save results as RDS
# saveRDS(eprot,"eprot.rds")
# saveRDS(eprot_m,"eprot.rds")
# saveRDS(eprot_u,"eprot.rds")
# ## read saved results
# eprot <- readRDS("eprot.rds")


stop("Never mind the error, execution stops here!")
##################################################




## THE SOLUTION TO THE PROBLEM of different outputs of as_MSnSet
## 12th Apr 2019 from Lauent in email
## Test data, one LC-MS/MS run
# msnid <- MSnID()
# msnid <- read_mzIDs(msnid,"mascot_fdr1pc/PDLP1_C1_2_1_1_120416.mzid")
# msnid

## Better to test samples that are to be merged (fractions)
## The corresponding file names (mzid_files) are recorded
## in the experiment (design) table

## Here we filter and group mzid_files from the experimental design table
## using only the first part of the S01_code.R, without building MSnSet
mzid <- "mascot_fdr1pc"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)
## APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()

## DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples)
## that will be combined
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

## thus from 'exp' and 'etab' we get
mzid_files
## that is a useful list
## to reprocess only subset of samples

##  Now we can compare the proteins we had problem earlier -
##  combined from several mzids
##  samples P_I_3 (list #23) and P_I_1 (list #32)
msnid <- MSnID()
#msnid <- read_mzIDs(msnid,mzid_files[["P_I_3"]])
msnid <- read_mzIDs(msnid,mzid_files[["P_I_1"]])
msnid


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

length(unique(fData(psm)$pepSeq))
length(unique(fData(pep1)$pepSeq))
length(unique(fData(pep2)$pepSeq))
length(unique(fData(pep2)$peptide))
dim(fData(pep2)[,c("pepSeq","peptide")])

## check the accesions == "20144"
## check the accesions == "04266"

table(fData(pep1)$accession=="04266")
rownames(fData(psm)[fData(psm)$accession=="04266",])
unique(fData(psm)[fData(psm)$accession=="04266",c(23,25)])
## WHY THE TWO FOLLOWING LINES PRODUCE DIFFERENT RESULTS?
unique(fData(psm)[fData(psm)$accession=="04266",25])
rownames(fData(pep1)[fData(pep1)$accession=="04266",])
rownames(fData(pep2)[fData(pep2)$accession=="04266",])
## both peptides should be in present in "04266" 
fData(pep1)[fData(pep1)$pepSeq=="FICTTGK",][c(25,23)]
fData(pep1)[fData(pep1)$pepSeq=="YPDHMK",][c(25,23)]

## protein MSnSet shoudl be correct in terms of unique accessions
unique(fData(psm)$accession)
length(unique(fData(psm)$accession))
## which is the same as prtoein MSnset
length(unique(fData(prot)$accession))



## HERE IS AN IDEA... what to do with missing information
## peptides and accessions in peptide and protein MSnSets
df <- fData(psm)
#saveRDS(df,"df.rds")
df <- fData(psm)[,c(25,23)]
dim(df)
df

## to get total spectral count (SPC) for peptides
spc <- df%>%group_by(pepSeq) %>% summarise(n=n())
dim(spc)
View(spc)
## total SPC per protein
pc <- df%>%group_by(accession) %>% summarise(n=n())
dim(pc)
View(pc)

## remove clear duplicates - the same accession and sequence
## count here is still not total SPC
df <- df%>%group_by(pepSeq,accession) %>% summarise(n=n())
dim(df)
View(df)

## Looking good, let's start again to get UPC
df <- fData(psm)
dim(df)
upc <- df%>%group_by(pepSeq, accession) %>% distinct() %>% summarise(n=n())
dim(upc)
View(head(upc))

## aggregate protein accessions
#df <- aggregate(pepSeq~accession, df, paste, collapse=",")
df <- fData(psm)
dim(df)
pe <- aggregate(accession~pepSeq, df, paste, collapse=",")
dim(pe)
head(pe)
View(pe)

## here we make list of peptide sequences,
## where each can have multiple accesions
pel <- as.list(pe$accession)
pel
names(pel) <- pe$pepSeq
pel
length(pel)
View(pel)

## aggregate peptide sequences 
df <- fData(psm)
dim(df)
pr <- aggregate(pepSeq~accession, df, paste, collapse=",")
dim(pr)
head(pr)
View(head(pr))

## here we make list of accession,
## where each can have multiple peptide sequences
prl <- as.list(pr$pepSeq)
prl
names(prl) <- pr$accession
prl
length(prl)
View(prl)

## protein list with unique peptides
uprl <- lapply(prl,function(x) strsplit(x,",")) 
uprl
uprl <- lapply(uprl,function(x) unique(unlist(x)))
uprl
uprlcnt <- lapply(uprl,function(x) length(x))
uprlcnt
## append the lists to peptide and protein MSnSets the list

## check it
tail(exprs(prot))
tail(prl)
tail(uprlcnt)

fvarLabels(prot)
fData(prot)$upepSeq <- prl
rownames(fData(prot))

## reordering uprl list in the same order as prot MSnSet
uprl[rownames(fData(prot))]
fData(prot)$upepSeq <- uprl[rownames(fData(prot))]
fData(prot)$upepSeq


## Now we could use the peptide-accession list for
## the redundancy handler in combineFearutes
## to obtain protein MSnSet
prottest <- combineFeatures(pep1,
                        groupBy = pel,
                        fun="sum",
                        redundancy.handler = "multiple")




