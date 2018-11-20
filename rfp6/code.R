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

etab <- experimentHierarchy(exp, fcol)

msnl <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(biorep == .etab[["biorep"]], 
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
    select(name)
  mzid_files <- file.path(mzid, paste(filenames[[1]], "mzid", sep = "."))
  e <- rtslprot:::make_pep_MSnSet(mzid_files, fcol = "pepSeq")
  e <- combineFeatures(e, fcol = "pepSeq", fun = sum)
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


epep <- combineFeatures(e, 
                        fcol = "pepSeq", 
                        fun = sum)
