source("S00-env.R")
library(rtslprot)
library(dplyr)

mzid <- "msgf"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)

exp <- exp %>% 
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()


fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)

msnl1 <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(fcols[[]] == .etab[[fcols[3]]], 
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
    select(name)
  mzid_files <- file.path(mzid, paste(filenames[[1]], "mzid", 
                                      sep = "."))
  #e <- D_make_pep_MSnSet(mzid_files, fcol = "pepSeq")
  e <- D_make_psm_MSnSet(mzid_files)
    e@processingData@files <- mzid_files
  sampleNames(e) <- paste(.etab, collapse = "_")
  e <- updateFvarLabels(e, sampleNames(e))
  if (validObject(e))
    return(e)
})

e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl))
  e <- MSnbase::combine(e, msnl[[i]])
rownames(etab) <- sampleNames(e)
pData(e) <- etab

save(e, file = "e.rda")

e0 <- e
e <- impute(e, method = "zero")

library("msmsTests")

null.f <- "y~treatment"
alt.f <- "y~phenotype+treatment"

e <- rtslprot:::msms_edgeR_test(e, 
                           null.f = null.f, 
                           alt.f = alt.f, 
                           fnm = "phenotype",
                           test_name = "phenotype")





D_make_psm_MSnSet <- function(mzid_files) {
  msnid <- MSnID()
  msnid <- read_mzIDs(msnid, mzid_files)
  ## TODO - filter MSnID and remove non-tryptic peptides
  td <- as(msnid, "data.table")
  td$count <- 1
  output <- readMSnSet2(td, ecol = ncol(td))
  return(output)
  }

D_make_pep_MSnSet <- function(mzid_files, fcol = "pepSeq") {
  e <- D_make_psm_MSnSet(mzid_files)
  output2 <- combineFeatures(e, groupBy = fData(msnset)$accession, fun = sum)
  return(output2)
}
