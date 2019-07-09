rm(list = ls());
source("S00-env.R")
library("rtslprot")
library("dplyr")
library("MSnID")
#library("msmsTests")

####################################################################
## 1 - read THE EXPERIMENT TABLE describing the experimental design:
## (i) samples
## (ii) measurements
## (iii) searches
## (iv) experimental factors

## mzid <- "mascot_fdr1pc"
mzid <- "msgf"
## mzid <- "msgf2"
exp <- readSampleExperimentTable("SampleExperimentTable.csv",
                                 mzid = mzid)
################################
## 2 - APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()

############################
## 3 - DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples) that will be combined
## all this have to be described in the experiment table
fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)
etab

#################################
## 4 - READ PATHS TO SEARCH DATA
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
mzid_files
## add names
names(mzid_files) <- apply(etab, 1,
                           function(.etab) paste(.etab, collapse = "_"))
#rename also etab
rownames(etab) <- names(mzid_files)
etab

## let us consider only subset od data
## thus from 'exp' and 'etab' we get
mzid_files
#x <- mzid_files["P_I_1"]
#x
y <- mzid_files[c("P_I_1","P_I_4","P_I_5" )]
y
mzid_files <- y
## and update mzid_files with it
mzid_files

######################################################
## 5 - continuing from the list of paths (mzid_files)
## advantage is that if we want, we can subset the list
index <- names(mzid_files)
names(index) <- names(mzid_files) 
index

msnl <- lapply(index, function(.x){
    msnid <- MSnID()
    msnid <- read_mzIDs(msnid, mzid_files[[.x]])
    
    #################################
    ## HERE ADJUST FDRs
    ## call another script 
    #source("S02-2_FDR.R")
    #################################
    
    psm <- rtslprot::as_MSnSet(msnid)
    pep <- rtslprot::as_MSnSet(msnid, "pepSeq")
    prot <- rtslprot::as_MSnSet(msnid, "accession")

## TO ADD MISSING INFO WE HAVE TO GO BACK TO PSM MSnSet
## Here, at the level of 'list of MSnSets',
## we modify all the items (MSnSets of individual samples)

    ## aggregate protein accessions
    #df <- aggregate(pepSeq~accession, df, paste, collapse=",")
    df <- fData(psm)
    #dim(df)
    pe <- aggregate(accession~pepSeq, df, paste, collapse=",")
    ## here we make list of peptide sequences,
    ## where each can have multiple accesions
    pel <- as.list(pe$accession)
    names(pel) <- pe$pepSeq
    
    ## lets do the same for protein
    ## aggregate peptide sequences 
    df <- fData(psm)
    pr <- aggregate(pepSeq~accession, df, paste, collapse=",")
    ## here we make list of accession,
    ## where each can have multiple peptide sequences
    prl <- as.list(pr$pepSeq)
    names(prl) <- pr$accession
    
    ## lists need to have unique peptides/accessions
    ## protein list with unique peptides
    uprl <- lapply(prl,function(x) strsplit(x,",")) 
    uprl <- lapply(uprl,function(x) unique(unlist(x)))
    uprlcnt <- lapply(uprl,function(x) length(x))
    #uprlcnt
    
    ## protein list with unique peptides
    upel <- lapply(pel,function(x) strsplit(x,",")) 
    upel <- lapply(upel,function(x) unique(unlist(x)))
    upelcnt <- lapply(upel,function(x) length(x))
    #head(upelcnt)
    
    
    ## the list to append to protein MSnSets
    fvarLabels(prot)
    fData(prot)$pepSeq_smpl <- uprl
    #rownames(fData(prot))
    
    ## the list to append to peptide MSnSets
    fvarLabels(pep)
    fData(pep)$accession_smpl <- upel
    #rownames(fData(prot))
    ## reordering uprl list in the same order as prot MSnSet
    uprl[rownames(fData(prot))]
    fData(prot)$pepSeq_smpl <- uprl[rownames(fData(prot))]
    #fData(prot)$upepSeq
    
    ## Now we could use the peptide-accession list for
    ## the redundancy handler in combineFearutes
    ## to obtain (hopefuly better)protein MSnSet
    ## Aternativelly, could we just merge protein MSnSets?
    
    e <- combineFeatures(pep,
                           groupBy = upel,
                           fun="sum",
                           redundancy.handler = "multiple")
    ## This new proteins need to complete peptide sequence info,
    ## we have them in uprl.
    fData(e)$pepSeq_smpl <- uprl
    
    ## update samples names
    sampleNames(e) <- .x
    e <- updateFvarLabels(e, sampleNames(e))
    if (validObject(e))
      return(e)
})
## list of MSnSets 'msnl' from all samples completed

## check the result
msnl
names(msnl)
sampleNames(msnl[[1]])
head(exprs(msnl[[1]]))
head(fData(msnl[[1]]))
for(x in msnl){print(length(fvarLabels(x)))}
length(fvarLabels(msnl[[1]]))
length(fvarLabels(msnl[[2]]))
## why some feature data are missing one column?
## which one is missing?
fvarLabels(msnl[[1]])
fvarLabels(msnl[[2]])
## to get rid of .SampleName
one <- fvarLabels(msnl[[1]])
two <- fvarLabels(msnl[[2]])
one <- gsub("\\..+","", one)
two <- gsub("\\..+","", two)
one
two
## any difference?
setdiff(one,two)
setdiff(two,one)
## WHY THE NUMBERS OF FEATURE-DATA COLUMNS DIFFER? 
#View(fData(msnl[[1]])[,c(26,24,34,25)])

#######################################
## 7 - COMBINE LIST OF MSnSets INTO ONE
## Convert the list of MSnSets to a single MSnSet where
## each element of the list becomes a columns (sample)
## and update the experimental design
e <- MSnbase::combine(msnl[[1]], msnl[[2]])
for (i in 3:length(msnl)) {
  e <- MSnbase::combine(e, msnl[[i]])
}
## 'msnl' list of MSnSets from all samples completed


## add sample meta data
pData(e) <- etab[index,]
pData(e)

######################################
## 8 - DEAL WITH NAs ON PEPIDE LEVEL
e <- impute(e, method = "zero")

######################################
## 9 - SAVE IT
## keep the original
saveRDS(e,"e.rds")
## open when needed
e <- readRDS("e.Rds")
## save as
saveRDS(e,"e_mascot_fdr1pc.rds")
e <- readRDS("e_mascot_fdr1pc.rds")


## check it
i <- grep(".?accession.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))
i <- grep(".?pepSeq.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))
#fData(e)[,i][1,4]

#######################################
## 10 - Global unique accession/peptide lists for combined MSnSet
## PROBLEM - the lists do not work well
## c(\"acc\",\"acc\")

## concatenate all  accessions
i <- grep("accession.?", fvarLabels(e))    # e is a peptide-level MSnSet
i
k <- apply(fData(e)[, i], 1,
            function(x) unique(na.omit(as.character(x))))
fData(e)$accession_all <- lengths(k)
fData(e)$accession <- sapply(k, paste, collapse = ";") # Laurent's suggestion
fData(e)$accession_all <- sapply(k, paste, collapse=NULL) # But when not collapsed, we then easily make a list
l <- as.list(fData(e)$accession_all)        # the list is nedded for combineFeatures
 
fData(e)$accession_all <- l

## concatenate all unique peptides
#i <- grep("AllAccession\\.", fvarLabels(e))    # e is a peptide-level MSnSet
j <- grep("pepSeq_smpl\\.?", fvarLabels(e))    # e is a peptide-level MSnSet
j
k <- apply(fData(e)[, j], 1,
           function(x) unique(na.omit(as.character(x))))
fData(e)$NpepSeq_all <- lengths(k)
#fData(e)$accession <- sapply(k, paste, collapse = ";") # Laurent's suggestion
fData(e)$pepSeq_all <- sapply(k, paste0, collapse=NULL) # But when not collapsed, we then easily make a list
l <- as.list(fData(e)$pepSeq_all)        # the list is nedded for combineFeatures
l 
fData(e)$pepSeq_all <- l

## check it
i <- grep(".?accession.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))
i <- grep(".?pepSeq.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))

## PROBLEMS
## It should not be needed to combine features again.
## 1) 'e' should have unique proteins
## 2) it should have peptides added to every sample
## 3) both peptide sequence and accessions were concatenated in *_all

## but how keep the best peptide feature data?
## by combineFeatures further up  
#save modified MSnSet
saveRDS(e,"e.rds")

######################################################
######################################################
stop("Never mind the error, execution stops here!")###
######################################################
######################################################







## tests, previous versions
 
eprot_m <- combineFeatures(e, groupBy = lpr,
                            fun = "sum",redundancy.handler = "multiple")
 
eprot_u <- combineFeatures(e, groupBy = fData(e)$accession,
                            fun = "sum", redundancy.handler = "unique")

## save results as RDS
 saveRDS(eprot,"eprot.rds")
 saveRDS(eprot_m,"eprot.rds")
 saveRDS(eprot_u,"eprot.rds")
# ## read saved results
# eprot <- readRDS("eprot.rds")


 ## check it
 i <- grep("accession.?", fvarLabels(eprot_m))
 i
 colnames(fData(eprot_m)[,i])
 View(head(fData(eprot_m)[,i],100))
 i
 
 fData(eprot_m)[,487]
 
 
 dim(fData(eprot_m))

## check the accesions == "20144"
## check the accesions == "04266"
table(fData(eprot_m)$accession=="04266")
rownames(fData(eprot_m)[fData(eprot_m)$accession=="04266",])
unique(fData(eprot_m)[fData(eprot_m)$accession=="04266",c(23,25)])
## WHY THE TWO FOLLOWING LINES PRODUCE DIFFERENT RESULTS?
unique(fData(eprot_m)[fData(eprot_m)$accession=="04266",25])
rownames(fData(eprot_m)[fData(eprot_m)$accession=="04266",])
rownames(fData(eprot_m)[fData(eprot_m)$accession=="04266",])
## both peptides should be in present in "04266" 
fData(eprot_m)[fData(eprot_m)$pepSeq=="FICTTGK",][c(25,23)]
fData(eprot_m)[fData(eprot_m)$pepSeq=="YPDHMK",][c(25,23)]

## protein MSnSet shoudl be correct in terms of unique accessions
unique(fData(eprot_m)$accession)
length(unique(fData(eprot_m)$accession))
## which is the same as prtoein MSnset
length(unique(fData(eprot_m)$accession))



##################################################
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
head(fData(pep2)[,c("pepSeq","peptide")])
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
head(df)

## to get total spectral count (SPC) for peptides
spc <- df%>%group_by(pepSeq) %>% summarise(n=n())
dim(spc)
head(spc)

## total SPC per protein
pc <- df%>%group_by(accession) %>% summarise(n=n())
dim(pc)
head(pc)

## remove clear duplicates - the same accession and sequence
## count here is still not total SPC
df <- df%>%group_by(accession,pepSeq) %>% summarise(n=n())
dim(df)
head(df)

## Looking good, let's start again to get UPC
df <- fData(psm)
dim(df)
upc <- df%>%group_by(pepSeq,accession) %>% summarise(n=n())
dim(upc)
head(upc)

## aggregate protein accessions
#df <- aggregate(pepSeq~accession, df, paste, collapse=",")
df <- fData(psm)
dim(df)
pe <- aggregate(accession~pepSeq, df, paste, collapse=",")
dim(pe)
head(pe)

## here we make list of peptide sequences,
## where each can have multiple accesions
pel <- as.list(pe$accession)
head(pel)
names(pel) <- pe$pepSeq
head(pel)
length(pel)
head(pel)

## protein list with unique peptides
head(pel)
upel <- lapply(pel,function(x) strsplit(x,",")) 
head(upel)
upel <- lapply(upel,function(x) unique(unlist(x)))
head(upel)
upelcnt <- lapply(upel,function(x) length(x))
head(upelcnt)
## append the lists to peptide and protein MSnSets the list


## aggregate peptide sequences 
df <- fData(psm)
dim(df)
pr <- aggregate(pepSeq~accession, df, paste, collapse=",")
dim(pr)
head(pr)
head(pr)

## here we make list of accession,
## where each can have multiple peptide sequences
prl <- as.list(pr$pepSeq)
head(prl)
names(prl) <- pr$accession
head(prl)
length(prl)
head(prl)

## protein list with unique peptides
head(prl)
uprl <- lapply(prl,function(x) strsplit(x,",")) 
head(uprl)
uprl <- lapply(uprl,function(x) unique(unlist(x)))
head(uprl)
uprlcnt <- lapply(uprl,function(x) length(x))
head(uprlcnt)
## append the lists to peptide and protein MSnSets the list

## check it
tail(exprs(prot))
tail(prl)
tail(uprl)
tail(uprlcnt)

fvarLabels(prot)
#fData(prot)$upepSeq <- prl
rownames(fData(prot))
## reordering uprl list in the same order as prot MSnSet
uprl[rownames(fData(prot))]
fData(prot)$upepSeq <- uprl[rownames(fData(prot))]
fData(prot)$upepSeq
#View(head(fData(prot)[,c(32,34)]))

## Now we could use the peptide-accession list for
## the redundancy handler in combineFearutes
## to obtain protein MSnSet
prottest <- combineFeatures(pep1,
                        groupBy = upel,
                        fun="sum",
                        redundancy.handler = "multiple")



## check it 
i <- grep(".?ccession.?",fvarLabels(prottest))
i
View(head(fData(prottest)[,i],100))
rownames(fData(prottest))
