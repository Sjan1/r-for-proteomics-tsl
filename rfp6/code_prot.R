## 25-11-2018: ASSIGN PEPTIDES TO PROTEINS ##################

## Make vector with all accessions matching every peptide
## Append the vector to the feature data
## Step through all the features (rows and columns)
## and copy accessions when found.
## Problem:
## For one peptide can match several protein accessions
## and in some mzid files (perhaps exported from Scaffold?),
## the accessions were on the same one row, separaded by comma.
## If this is the case, We read them all in a vectake the first accession.

load("e.rda") #MSnSet with peptides from the script code.R

accvec <- c() # vector with the acctual Accession numbers 
acc <- grep("accession",fvarLabels(e),value = TRUE) # headers (featureNames) with protein accessions
df_temp <- (fData(e)[,acc])
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
length(featureNames(e))
## length of MSnSet
length(exprs(e))
length(fData(e))
fData(e)
dim(fData(e))
## The lenght of accvec should be the number of non redundant protein
## after features are combined this will be number of combined features
## length(featureNames(comb1))
length(unique(accvec))
## Are there comma separated Accession numbers on the same row?
grep(",",accvec)
## Are there comma separated Accession numbers in the MSnSet?
grep(",",fData(e)$acc)
## append combined accessions to feature data
fData(e)$acc <- accvec

## remove NAs before combining features
head(exprs(e))
e <- impute(e, method = "zero")
head(exprs(e))

## how data looks like?
table(exprs(e))
## plot historam
hist(exprs(e))

## PROTEOTYPIC PEPTIDES ########################################
## combining proteotypic peptides to the corresponding proteins
eprot <- combineFeatures(e, groupBy = fData(e)$acc, 
                         fun = "sum")

## sanity checks
## number of combined features
length(featureNames(eprot))
##number of features in the original combined_e
length(fData(e)$acc)

## Let's check the abundant proteins present in more than one sample
## were combined correctly
test <- exprs(eprot)
#vector to order test
ordervec <- base::rowSums(test,na.rm=TRUE)
#order rows of test using ordervec
head(test[order(ordervec,decreasing=TRUE),])
## Let's pick one of the accessions
## and find which peptides it containts 
one <- fData(e)
pep <- rownames(one[one$acc=="AT1G23410.1",])
pep
## This is how the peptides were distribute among smaples
exprs(e)[pep,]
## This is after features were combineed
ex1 <- exprs(eprot)
s <- rownames(ex1)=="AT1G23410.1"
ex1[s,,drop=FALSE]
## how data looks like?
table(exprs(eprot))
## plot historam
hist(exprs(eprot))



## save results as RDS
saveRDS(eprot,"eprot.rds")
eprot <- readRDS("eprot.rds")

#hist(exprs(eprot))
#hist(log(exprs(eprot),2)[exprs(eprot)>1])
#hist(log(exprs(eprot),2))

eprot0 <- eprot
eprot <- impute(eprot, method = "zero")


## phenotype
null.f <- "y~1"
alt.f <- "y~phenotype"

## treatment
null.f <- "y~1"
alt.f <- "y~treatment"

#phenotype ph0-treatment
null.f <- "y~phenotype"
alt.f <- "y~phenotype+treatment"

#block treat0_phenotype
null.f <- "y~treatment"
alt.f <- "y~phenotype+treatment"



eprot <- rtslprot:::msms_edgeR_test(eprot, 
                                null.f = null.f, 
                                alt.f = alt.f, 
                                fnm = "phenotype",
                                test_name = "treat0_phenotype")


## volcano plot
plot(fData(eprot)$LogFC_treat0_phenotype,-log10(fData(eprot)$adj.p.values_treat0_phenotype))

plot(fData(eprot)$LogFC_treatment,-log10(fData(eprot)$adj.p.values_treatment))
plot(fData(eprot)$LogFC_block_treat,-log10(fData(eprot)$adj.p.values_block_treat))
