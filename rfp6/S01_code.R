source("S00-env.R")
library("rtslprot")
library("dplyr")
library("msmsTests")

## READ THE EXPERIMENT TABLE that describes:
## (i) samples
## (ii) measurements
## (iii) searches
## (iv) experimental factors
mzid <- "msgf"
mzid <- "mascot"
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
e0 <- e
saveRDS(e0,"e0.rds")
## save others
saveRDS(e,"e_msgf.rds")
saveRDS(e,"e_mascot.rds")
saveRDS(e,"e_mascot_fdr1pc.rds")

## open if needed
e <- readRDS("e0.Rds")


## 20-03-2019: COMBINE PEPTIDES INTO PROTEINS
## Make vector with all accessions matching every peptide
## Append the vector to the feature data
## Step through all the features (rows and columns)
## and copy accessions when found.

load("e.rda") #MSnSet with peptides from the script code.R
e <- readRDS("e.rds")
e <- readRDS("e_mascot.rds")
e <- readRDS("e_msgf.rds")
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
saveRDS(eprot0,"eprot0.rds")
saveRDS(eprot,"eprot_msgf.rds")
saveRDS(eprot,"eprot_mascot.rds")
saveRDS(eprot,"eprot_mascot_fdr1pc.rds")
saveRDS(eprot,"eprot_m.rds")
saveRDS(eprot0,"eprot_u.rds")
## read the original
eprot <- readRDS("eprot.rds")
eprot <- readRDS("eprot0.rds")
eprot <- readRDS("eprot_msgf.rds")
eprot <- readRDS("eprot_mascot.rds")
eprot <- readRDS("eprot_mascot_fdr1pc.rds")
eprot <- readRDS("eprot_m")
eprot <- readRDS("eprot_u")






## statistical tests

## null_TR__alt_PH+TR 
null.f <- "y~treatment"
alt.f <- "y~phenotype+treatment"

## Can variance be explained by phenotype? 
####null_1__alt_PH
null.f <- "y~1"
alt.f <- "y~phenotype"

## is treatment doing anything?   
## null_PH+TR__alt_PH*TR
null.f <- "y~phenotype+treatment"
alt.f <- "y~phenotype*treatment"


#e.notTreatedSamples <- e[e$treatment=="H",]
e <- rtslprot:::msms_edgeR_test(e, 
                                null.f = null.f, 
                                alt.f = alt.f, 
                                fnm = "phenotype",
                                test_name = "null_TR__alt_PH+TR")
e <- rtslprot:::msms_edgeR_test(e, 
                                null.f = null.f, 
                                alt.f = alt.f, 
                                fnm = "phenotype",
                                test_name = "null_1__alt_PH")

e <- rtslprot:::msms_edgeR_test(e, 
                           null.f = null.f, 
                           alt.f = alt.f, 
                           fnm = "phenotype",
                           test_name = "null_PH+TR__alt_PH*TR")



## null_TR__alt_PH+TR
plot(fData(e)$`LogFC_null_TR__alt_PH+TR`,-log10(fData(e)$`p.value_null_TR__alt_PH+TR`))
hist(fData(e)$`p.value_null_TR__alt_PH+TR`)
hist(fData(e)$`adj.p.values_null_TR__alt_PH+TR`)
## null_1__alt_PH
plot(fData(e)$`LogFC_null_1__alt_PH`,-log10(fData(e)$`p.value_null_1__alt_PH`))
hist(fData(e)$`p.value_null_1__alt_PH`)
hist(fData(e)$`adj.p.values_null_1__alt_PH`)
## null_PH+TR__alt_PH*TR
plot(fData(e)$`LogFC_null_PH+TR__alt_PH*TR`,-log10(fData(e)$`p.value_null_PH+TR__alt_PH*TR`))
hist(fData(e)$`p.value_null_PH+TR__alt_PH*TR`)
hist(fData(e)$`adj.p.values_null_PH+TR__alt_PH*TR`)

length(unique(fData(e)$`p.value_null_TR__alt_PH+TR`))
length(unique(fData(e)$`p.value_null_1__alt_PH`))
length(unique(fData(e)$`p.value_null_PH+TR__alt_PH*TR`))


saveRDS(e,"e.rds")
##
head(exprs(e))
head(fData(e))[,1:3]
pData(e)
