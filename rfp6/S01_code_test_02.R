source("S00-env.R")
library("rtslprot")
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
       e <- make_pep_MSnSet(mzid_files, fcol = "pepSeq")
  #    e <- make_psm_MSnSet(mzid_files)
  e@processingData@files <- mzid_files
  sampleNames(e) <- paste(.etab, collapse = "_")
  e <- updateFvarLabels(e, sampleNames(e))
  if (validObject(e))
    return(e)
})
#save msnl based on PSMs
#msnl_psms <- msnl

e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl)){
  e <- MSnbase::combine(e, msnl[[i]])}
rownames(etab) <- sampleNames(e)
pData(e) <- etab

## TROUBLESHOOTING START ######################################################################
## we had a suspicion exprs data exists as multiples of basic integers...

#MSnSets, msnl, based on PSMs was saved here: 
msnl_psms <- msnl

## get peptide seq from the 1st MSnSet of PSMs
pep_msnl_psm1 <- fData(msnl_psms[[1]])[26]
## coutn them
PepFromPsm <- pep_msnl_psm1 %>% group_by(pepSeq.C_H_2) %>% summarise(Unq_pep = n())
PepFromPsm
## comparison wih peptides from make_pep_MSnSet
## should give the equal numbers, and it gives:
## SPC sum of the whole sample
sum(exprs(msnl[[1]]))     # peptides combined by readMSnSet2
sum(PepFromPsm$Unq_pep)   # peptides combined here with dplyr
## The problem must be in MSnbase::combine 

## 2nd sample
pep_msnl_psm2 <- fData(msnl_psms[[2]])[26]
PepFromPsm2 <- pep_msnl_psm2 %>% group_by(pepSeq.C_H_3) %>% summarise(Unq_pep2 = n())
## SPC sum of the whole sample
sum(exprs(msnl[[2]]))  
sum(PepFromPsm2$Unq_pep2)

## 3nd sample
pep_msnl_psm3 <- fData(msnl_psms[[3]])[26]
PepFromPsm3 <- pep_msnl_psm3 %>% group_by(pepSeq.C_H_4) %>% summarise(Unq_pep3 = n())
## SPC sum of the whole sample
sum(exprs(msnl[[3]]))  
sum(PepFromPsm3$Unq_pep3)

## to visualize a single peptide in msnl[[1]]
## PSMs just need to be just summarized
pep_msnl_psm1 %>% filter(pepSeq.C_H_2 =="AAALNIVPTSTGAAK") %>% 
  group_by(pepSeq.C_H_2) %>% 
  summarise(Unq_pep = n())

## peptides in the resulting MSnSet "e" (all "msnl" combined)
d <- as.data.frame(exprs(e))
d <- rownames_to_column(d, var = "rowname") #rownames go to the first column
rownames(d) <- d$rowname
d %>% filter(rowname =="AAALNIVPTSTGAAK") %>% 
  group_by(rowname)
summarise(Unq_pep = n())

## the calculation seems to be correct!!

## view only prtoeins without NAs (complete cases)
enona <- omit.na(exprs(e))
View(enona)
## also this look good to me

## TROUBLESHOOTING END ###################################################################

## save the result
save(e, file = "e.rda")
saveRDS(e,"e.rds")
## open it
e <- readRDS("e.rds")
##keep the original
e0 <- e
e <- e0

## deal with NAs
e <- impute(e, method = "zero")

## STATISTICAL TESTS

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
