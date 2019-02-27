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
e <- readRDS("e_mascot.Rds")
e <- readRDS("e_msgf.Rds")


uniqacc <- function(e) {
accvec <- c() # vector with the acctual Accession numbers 
acc <- grep("^[Aa]ccession",fvarLabels(e),value = TRUE) # headers (featureNames) with protein accessions
df_temp <- (fData(e)[,acc])
for (rn in rownames(df_temp)) {
  rowvec <- c();
  onerow <- df_temp[rn,]
  for (cn in 1:length(onerow)) {
    oneacc <- onerow[1,cn]
    if (is.na(oneacc)){}else{rowvec <- c(rowvec,oneacc);
    break}    
  }
  acc.str <- paste(unique(rowvec), collapse = ",")
  accvec <- c(accvec, acc.str)
    if(length(unique(rowvec)) != 1) {
  message(paste("There are several accessions for a peptide on row No.",rn, 
                       "of MSnSet!"))
}}
return(accvec);
}

rm(accvec)
accvec <- uniqacc(e);
length(accvec)

fData(e)$Prot_Acc <- accvec
fvarLabels(e)
saveRDS(e,"e.rds")



length(accvec)
## lenght of accvec should be the same as the length of featurenames
length(accvec)
length(featureNames(e))
## length of MSnSet
length(exprs(e))
length(fData(e))
head(fData(e))
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

saveRDS(e,"e.rds")

## remove NAs before combining features
head(exprs(e))
e <- impute(e, method = "QRILC")
e <- impute(e, method = "zero")
head(exprs(e))

saveRDS(e,"e.rds")

## how data looks like?
table(exprs(e))
## plot historam
hist(exprs(e))


head(exprs(eprot))
## PROTEOTYPIC PEPTIDES ########################################
## combining proteotypic peptides to the corresponding proteins
eprot <- combineFeatures(e, groupBy = fData(e)$Prot_Acc, 
                         fun = "sum")

## save results as RDS
saveRDS(eprot,"eprot.rds")
eprot <- readRDS("eprot.rds")
## save as the original
eprot0 <- eprot
saveRDS(eprot0,"eprot0.rds")
saveRDS(eprot,"eprot_msgf.rds")
saveRDS(eprot,"eprot_mascot.rds")
## read the original
eprot <- readRDS("eprot0.rds")
eprot <- readRDS("eprot_msgf.rds")
eprot <- readRDS("eprot_mascot.rds")
eprot.msgf <- eprot



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




## statistical tests

## 01
## Can variance be explained by phenotype and treatment? 
## null_TR__alt_PH+TR 
null.f <- "y~treatment"
alt.f <- "y~phenotype+treatment"

eprot <- rtslprot:::msms_edgeR_test(eprot, 
                                null.f = null.f, 
                                alt.f = alt.f, 
                                fnm = "phenotype",
                                test_name = "null_TR__alt_PH+TR")

## 02
## Can variance be explained by phenotype? 
## null_1__alt_PH
null.f <- "y~1"
alt.f <- "y~phenotype"

eprot <- rtslprot:::msms_edgeR_test(eprot, 
                                    null.f = null.f, 
                                    alt.f = alt.f, 
                                    fnm = "phenotype",
                                    test_name = "null_1__alt_PH")

## 03
## Is ohenotype an treatment influence each other?   
## null_PH+TR__alt_PH*TR
null.f <- "y~phenotype+treatment"
alt.f <- "y~phenotype*treatment"

eprot <- rtslprot:::msms_edgeR_test(eprot, 
                                    null.f = null.f, 
                                    alt.f = alt.f, 
                                    fnm = "phenotype",
                                    test_name = "null_PH+TR__alt_PH*TR")


## histogram
    ## not adjisted p-val
    ## adjusted p-val 
## null_TR__alt_PH+TR
hist(fData(eprot)$`p.value_null_TR__alt_PH+TR`)
hist(fData(eprot)$`adj.p.values_null_TR__alt_PH+TR`)
## null_1__alt_PH
hist(fData(eprot)$`p.value_null_1__alt_PH`)
hist(fData(eprot)$`adj.p.values_null_1__alt_PH`)
## null_PH+TR__alt_PH*TR
hist(fData(eprot)$`p.value_null_PH+TR__alt_PH*TR`)
hist(fData(eprot)$`adj.p.values_null_PH+TR__alt_PH*TR`)

## volcano plot
# non adjusted p-values
# null_TR__alt_PH+TR
plot(fData(eprot)$`LogFC_null_TR__alt_PH+TR`,
     -log10(fData(eprot)$`p.value_null_TR__alt_PH+TR`))
# null_1__alt_PH
plot(fData(eprot)$`LogFC_null_1__alt_PH`,
     -log10(fData(eprot)$`p.value_null_1__alt_PH`))
# null_PH+TR__alt_PH*TR
plot(fData(eprot)$`LogFC_null_PH+TR__alt_PH*TR`,
     -log10(fData(eprot)$`p.value_null_PH+TR__alt_PH*TR`), ylim=c(0,10))

## p-value distribution interpretation
## http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

## what next -- distribusion of p-values seems to be influenced by imputation of the missing values.
## what imputiaton to use?
## at least remove the hits with SPC=<1, scoring in only one replicate of the treatment or phenotype
