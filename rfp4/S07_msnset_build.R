## Date: 12th-November-2018
## This script needs rdf from S06_table.R 
## S07_msnset_build.R creates MSnSets from the files in from DOE table.
## 
##    
##rm(list=ls())
source("S00-env.R")

head(rdf)

# Here we neeed to bring in the data (like in the scripts S02, S03)

fn <- dir("MSGF", full.names = TRUE)
## fn <- dir("Mascot", full.names = TRUE)

msnid <- MSnID()
msnid <- read_mzIDs(msnid, fn)

labs <- rdf$labs
names(labs) <- rdf$file
labs["PDLP1_wt-1t_1_1_1_111206"]

labs_sub <- labs[sub("\\.mzML$","",msnid$spectrumFile)]
# new column called 'label' to put our sample-experiment description in
msnid$label <- labs[sub("\\.mzML$","",msnid$spectrumFile)]

## different fractions, same lane (all same sample)
td <- as(msnid, "data.table")

# check for missing labels 
# due to samples not being measured or inconsistency in the Smpl-Exp table
x <- unique(td$label)
length(x)
length(rdf)
rdf
# idealy comparison of rdf and unique(td$label)
#sel <- is.na(td$label)
#sum(sel)
#table(sel)
#View(td[sel,])

## 180716 
## all one 
## here we count spectra, different spectral counts in the new column
sel <- !duplicated(td$pepSeq)
count <- table(td$pepSeq)
x <- td[sel, ]
x$count <- as.vector(count[x$pepSeq])
#View(x)

## read all (samples, fractions) into one msnset object
i <- which(names(x) == "count")
e <- readMSnSet2(x, i)
featureNames(e) <- fData(e)$pepSeq


## keeping samples separately
list_msnsets <- list()
reps <- unique(sub("\\d+$","",rdf$labs,perl = TRUE))

## reps contains unique biosamples (fractions removed)
## the counter is needed for the loop
counter <- 1
for (biorep in reps) {
  #grepl(reps[1],td$label)
  tds <- td[grepl(biorep,td$label),]

#merging fraction, making msnset 
  sel <- !duplicated(tds$pepSeq)
  count <- table(tds$pepSeq)
  x <- tds[sel, ]
  x$count <- as.vector(count[x$pepSeq])
  #View(x)

## read into msnset object
  i <- which(names(x) == "count")
  e <- readMSnSet2(x, i)
  featureNames(e) <- fData(e)$pepSeq

  ## update sample names and Feature Labels
  sampleNames(e) <- biorep
  #e <- updateSampleNames(e,biorep)
  e <- updateFvarLabels(e,biorep)
  
  ## First MSnSet
  if(counter==1){
    combined_e <- e
  }
  ## Second MSnSet and up -> here we combine generated MSnSets 
  if(counter>1){
    combined_e <- BiocGenerics::combine(combined_e,e)
  }
  
  #initial idea that did not work  
  #list_msnsets[[biorep]] <- e
  
    
  ##to stop it
  if(counter>999){break}
  
  ## visualize the progress
    print(paste(counter, biorep, sep=" - "))
    Sys.sleep(0.01)
    flush.console()
    
    counter <- counter+1  
}
saveRDS(combined_e,"combined_e.rds")
combined_e <- readRDS("combined_e.rds")

##saveRDS(combined_e,"combined_e_Mascot.rds")
##combined_e <- readRDS("combined_e_Mascot.rds")

## sums of peptides in individual samples
sample_peptide_sums <- colSums(exprs(combined_e),na.rm = TRUE)
temp <- names(sample_peptide_sums)
otemp <- order(temp)
## plot total SPCs for samples
barplot(sample_peptide_sums[otemp])

## add some more QC here...



## 2-10-2018: ASSIGN PEPTIDES TO PROTEINS ###################################
## Let us create a vector with all accessions matching every peptide matches.
## Mind the possibility one peptide can match more accessions
## This vector we then append to feature data
## In a loop we will step through features and remember accession when found
accvec <- c()
acc <- grep("accession",fvarLabels(combined_e),value = TRUE)
df_temp <- (fData(combined_e)[,acc])
for (rn in rownames(df_temp)) {
  onerow <- df_temp[rn,]
  for (cn in 1:length(onerow)) {
    oneacc <- onerow[1,cn]
    if (is.na(oneacc)){}else{accvec <- c(accvec,oneacc);
      break}    
  }
  
}

## lenght of accvec should be the same as the length of featurenames
length(accvec)
length(featureNames(combined_e))
## The lenght of accvec should be the number of non redundant protein
## after features are combined this will be number of combined features
## length(featureNames(comb1))
length(unique(accvec))
## Are there more accesions inthe same row?
grep(",",accvec)
grep(",",fData(combined_e)$acc)
## append combined accessions to feature data
fData(combined_e)$acc <- accvec

## remove NAs before combining features
head(exprs(combined_e))
combined_e <- impute(combined_e, method = "zero")
head(exprs(combined_e))

## PROTEOTYPIC PEPTIDES ########################################
## combining proteotypic peptides to the corresponding proteins
comb1 <- combineFeatures(combined_e, groupBy = fData(combined_e)$acc, 
                                               fun = "sum")

## sanity checks
## number of combined features
length(featureNames(comb1))
##number of features in the original combined_e
length(fData(combined_e)$acc)

## Let's check the abundant proteins present in more than one sample
## were combined correctly
test <- exprs(comb1)
#vector to order test
ordervec <- base::rowSums(test,na.rm=TRUE)
#order rows of test using ordervec
head(test[order(ordervec,decreasing=TRUE),])
## Let's pick one of the accessions
## and find which peptides it containts 
one <- fData(combined_e)
pep <- rownames(one[one$acc=="AT1G23410.1",])
pep
## This is how the peptides were distribute among smaples
exprs(combined_e)[pep,]
## This is after features were combineed
ex1 <- exprs(comb1)
s <- rownames(ex1)=="AT1G23410.1"
ex1[s,,drop=FALSE]


## save results as RDS
saveRDS(comb1,"comb1.rds")
comb1 <- readRDS("comb1.rds")


## ALL PEPTIDES ###########################################
## combining ALL peptides to  to the corresponding proteins
comb2 <- combineFeatures(combined_e, groupBy = accvec, 
                      redundancy.handler = "multiple",
                      fun = "sum")

## new sanity checks
## number of combined features
length(featureNames(comb2))
##number of features in the original combined_e
length(fData(combined_e)$acc)

## Let's check the abundant proteins present in more than one sample
## were combined correctly
test <- exprs(comb2)
#vector to order test
ordervec <- base::rowSums(test,na.rm=TRUE)
#order rows of test using ordervec
head(test[order(ordervec,decreasing=TRUE),])
## Let's pick one of the accessions
## and find which peptides it containts 
one <- fData(combined_e)
pep <- rownames(one[one$acc=="AT1G23410.1",])
pep
## This is how the peptides were distribute among smaples
exprs(combined_e)[pep,]
## This is after features were combineed
ex2 <- exprs(comb2)
s <- rownames(ex1)=="AT1G23410.1"
ex2[s,,drop=FALSE]

## save results as RDS
saveRDS(comb2,"comb2.rds")
comb2 <- readRDS("comb2.rds")


## END END END #########################################

## 6-10-2018
## check all peptides are unique
pepvec <- c()
pep <- grep("peptide",fvarLabels(combined_e),value = TRUE)
df_temp <- (fData(combined_e)[,pep])
for (rn in rownames(df_temp)) {
  onerow <- df_temp[rn,]
  for (cn in 1:length(onerow)) {
    onepep <- onerow[1,cn]
    if (is.na(onepep)){}else{pepvec <- c(pepvec,onepep);
    break}    
  }
}
length(pepvec)
length(unique(pepvec))























########################################################
## some older testing
## combine all into one msnset
msnset <-  do.call(BiocGenerics::combine,list_msnsets)

msnset <-  BiocGenerics::combine(list_msnsets[[1]],list_msnsets[[2]],list_msnsets[[3]])

msnset <-  BiocGenerics::combine(list_msnsets[[1:3]])


#testing on two msnsets
e1 <- list_msnsets[[1]]
e2 <- list_msnsets[[2]]
sampleNames(e1) <- c("sample1")
sampleNames(e2) <- c("sample2")
e1 <- updateSampleNames(e1)
e2 <- updateSampleNames(e2)
e1 <- updateFvarLabels(e1)
e2 <- updateFvarLabels(e2)
c <- BiocGenerics::combine(e1,e2)

#e1_1 <- list_msnsets[[1]][1:10]
#e1_2 <- list_msnsets[[1]][11:21]
#c <- do.call(combine,list(e1_1,e1_2))

e11 <- e[1:10,]
e22 <- e[11:21,]
sampleNames(e11) <- c("sample11")
sampleNames(e22) <- c("sample22")

e11 <- updateSampleNames(e11)
e22 <- updateSampleNames(e22)

#e11 <- updateFeatureNames(e11)
#e22 <- updateFeatureNames(e22)

e11 <- updateFvarLabels(e11)
e22 <- updateFvarLabels(e22)

c <- BiocGenerics::combine(e11,e22)


sampleNames(comb2)
tab[1:3, ]

pd <- data.frame(tab[1:3, ])
sampleNames(msnset)[1:3]
rownames(pd) <- sampleNames(msnset)[1:3]
pData(msnset) <- pd

msnset$fraction == 1
