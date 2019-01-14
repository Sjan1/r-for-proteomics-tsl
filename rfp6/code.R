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

msnl <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(biorep == .etab[["biorep"]], 
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
    select(name)
  mzid_files <- file.path(mzid, paste(filenames[[1]], "mzid", 
                                      sep = "."))
## make MSnSet (choose peptide or PSMs - two functions in rtslprot)
    e <- rtslprot:::make_pep_MSnSet(mzid_files, fcol = "pepSeq")
    e@processingData@files <- mzid_files
  sampleNames(e) <- paste(.etab, collapse = "_")
  e <- updateFvarLabels(e, sampleNames(e))
  if (validObject(e))
    return(e)
})

e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl)){
  e <- MSnbase::combine(e, msnl[[i]])}
rownames(etab) <- sampleNames(e)
pData(e) <- etab

save(e, file = "e.rda")
saveRDS(e,"e.rds")
##keep the original
e0 <- e
e <- e0
## open if needed
e <- readRDS("e.Rds")

## deal with NAs
e <- impute(e, method = "zero")

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

plot(fData(e)$LogFC_phenotype,-log10(fData(e)$adj.p.values_phenotype))

## null_TR__alt_PH+TR
hist(fData(e)$`p.value_null_TR__alt_PH+TR`)
hist(fData(e)$`adj.p.values_null_TR__alt_PH+TR`)
## null_1__alt_PH
hist(fData(e)$`p.value_null_1__alt_PH`)
hist(fData(e)$`adj.p.values_null_1__alt_PH`)
## null_PH+TR__alt_PH*TR
hist(fData(e)$`p.value_null_PH+TR__alt_PH*TR`)
hist(fData(e)$`adj.p.values_null_PH+TR__alt_PH*TR`)

length(unique(fData(e)$`p.value_null_TR__alt_PH+TR`))
length(unique(fData(e)$`p.value_null_1__alt_PH`))
length(unique(fData(e)$`p.value_null_PH+TR__alt_PH*TR`))


save(e, file = "e.rda")
##
head(exprs(e))
head(fData(e))[,1:3]
pData(e)
