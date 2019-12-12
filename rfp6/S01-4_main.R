rm(list = ls());
source("../rfp4/S00-env.R")
source("run_msnid.R")
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
y <- mzid_files[c("P_H_2","P_H_3","P_H_4","P_H_5","C_H_2","C_H_3","C_H_4")]
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
    .files <- mzid_files[[.x]]
    msnid <- run_msnid(.files)
    psm <- rtslprot::as_MSnSet(msnid)
    df <- fData(psm)
    
    ## TO ADD MISSING INFO WE HAVE TO GO BACK TO PSM MSnSet
    ## Here, at the level of 'list of MSnSets',
    ## we modify all the items (MSnSets of individual samples)
    
    ## aggregate protein accessions
    pe <- aggregate(accession~pepSeq, df, c)
    ## here we make list of peptide sequences,
    ## where each can have multiple accesions
    pel <- as.list(pe$accession)
    names(pel) <- pe$pepSeq
    pel <- lapply(pel, unique)
    gb <- pel[df$pepSeq]
    names(gb) <- featureNames(psm)
    
    pe2 <- aggregate(pepSeq~accession, df, c)
    pe2l <- as.list(pe2$pepSeq)
    names(pe2l) <- pe2$accession
    pep_count <- sapply(pe2l, function(x) length(unique(x)))
   
    ## new for unique peptides
    pe2lu <- sapply(pe2l, function(x) unique(x))
    pe2ludf <- pe2lu[df$accession]
    ## a check - should give the same result from pe2lu list and df
    # unique(df[,c(24,32)][which(df$accession=="04266"),])
    # pe2lu[which(names(pe2lu)=="04266")]
    # pe2ludf[which(names(pe2ludf)=="04266")]
    
    ## re-ordered protein group
    ## protein with the highest spectral count => master protein
    gb2 <- gb_count <- gb
    for (i in seq_along(gb)) {
        gb2[[i]] <- gb[[i]][order(pep_count[gb[[i]]], decreasing = TRUE)]
        gb_count[[i]] <- paste(pep_count[gb2[[i]]], collapse = ";")
    }
    ## put the groups into feature data of PSM MSnSet
    fData(psm)$protein_groups <- sapply(gb2, paste, collapse = ";")
    fData(psm)$group_size <- lengths(gb2)
    fData(psm)$pep_counts <- unlist(gb_count)
    fData(psm)$master_prot <- sapply(gb2, "[[", 1)
    #unique peptides list
    fData(psm)$pepSeqUnq <- sapply(pe2ludf, paste, collapse = ";")
    
    for (k in unique(fData(psm)$master_prot)) {
        i <- which(fData(psm)$master_prot == k)
        j <- which.max(fData(psm)[i, "group_size"])
        ii <- i[j]
        fData(psm)[i, c("protein_groups",
                        "group_size",
                        "pep_counts")] <-
            fData(psm)[ii, c("protein_groups",
                             "group_size",
                             "pep_counts")]
    }
    

    ## comine features based on master proteins    
        e <- combineFeatures(psm, fcol = "master_prot",
                         fun = "sum", cv = FALSE)

    ## update samples names
    sampleNames(e) <- .x
    e <- updateFvarLabels(e, sampleNames(e))
    if (validObject(e))
      return(e)
})
## list of MSnSets 'msnl' from all samples completed


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
## 8 - DEAL WITH heNAs ON PEPIDE LEVEL
e <- impute(e, method = "zero")



#######################################
## 9 - Global unique accession/peptide lists for combined MSnSet
## PROBLEM - the lists do not work well
## c(\"acc\",\"acc\")

# ## concatenate all  accessions --- THIS is useless, accessions, PepSeq. and peptide. 
## are both "first from a lot". 
## To get all accessions we have to use protein groups,
## for peptides we need pepSeqUnq.


## concatenate all unique peptides in master proteins in all samples
j <- grep("pepSeqUnq\\.?", fvarLabels(e))    # e is a peptide-level MSnSet
j
k <- apply(fData(e)[, j], 1,
           function(x) unique(na.omit(as.character(x))))

k <- sapply(k, paste0,collapse=";")
k <- sapply(k, function(x) strsplit(x,";"))
k <- sapply(k, function(x) unique(x))

fData(e)$pepSeqUnq_counts <- lengths(k)
fData(e)$pepSeqUnq_all <- sapply(k, paste, collapse = ";")

## if this is neede in a form of a list for combineFeatures...
#fData(e)$accession <- sapply(k, paste, collapse = ";") # Laurent's suggestion
l <- as.list(sapply(k, paste0, collapse=NULL)) # But when not collapsed, we then easily make a list
# fData(e)$pepSeqUnq_all_list <- l        # the list is nedded for combineFeatures
head(l)                       



## global protein groups
## concatenate all  accessions in all protein groups in all samples
i <- grep("protein_groups\\.?", fvarLabels(e))    # e is a peptide-level MSnSet
i
k <- apply(fData(e)[, i], 1,
           function(x) unique(na.omit(as.character(x))))
k
k <- sapply(k, paste0,collapse=";")
k <- sapply(k, function(x) strsplit(x,";"))
k <- sapply(k, function(x) unique(x))
k
fData(e)$group_size_global <- lengths(k)
fData(e)$protein_group_global <- sapply(k, paste, collapse = ";")

## and the list if needed
pggl <- as.list(sapply(k, paste, collapse=NULL)) # But when not collapsed, we then easily make a list
#fData(e)$protein_group_global_list <- pggl      # the list is nedded for combineFeatures







########################################
## 10 - USEFUL
a <- grep(".?accession.?", fvarLabels(e))
a
#View(head(fData(e)[,i],100))
m <- grep(".?master.?", fvarLabels(e))
m
#View(head(fData(e)[,y],100))
pg <- grep(".?protein_groups.?", fvarLabels(e))
pg
#View(head(fData(e)[,pg],100))
s <- grep(".?pepSeq.?", fvarLabels(e))
s
su <- grep(".?pepSeqUnq.?", fvarLabels(e))
su
gs <- grep("group_size\\.?",fvarLabels(e))
gs
pc <- grep("pep_counts\\.?",fvarLabels(e))
pc
pgg <- grep("protein_group_global.?",fvarLabels(e))
pgg
#View(head(fData(e)[,pgg],100))
gsg <- grep("group_size_global\\.?",fvarLabels(e))
gsg

fData(e)[,c(m,a,pg,s,su,gs,pc,pgg,gsg)]
View(head(fData(e)[,c(m,pg,su,gs,pc,pgg,gsg)],100))


######################################
## 11 - SAVE IT
## keep the original
saveRDS(e,"e.rds")
## open when needed
e <- readRDS("e.Rds")
## save as
saveRDS(e,"e_mascot_fdr1pc.rds")
e <- readRDS("e_mascot_fdr1pc.rds")




## CHECKS OF RESULTS INTEGRITY
## Check I ############################
## featureNames should be master proteins.
## Master proteins can be either the same or NAs
## for all the samples. Are they?
View(fData(e)[,m])
c1 <- featureNames(e) == fData(e)[,m[1]]
c2 <- featureNames(e) == fData(e)[,m[2]]
c3 <- featureNames(e) == fData(e)[,m[3]]
#c
#c[is.na(c)] <- TRUE
#c
#any(isFALSE(c))
#any(isTRUE(c))
#c <- c(c,FALSE)
#table(c)
#table(!is.na(c))

## using table
table(c1, useNA="always")
table(c2, useNA="always")
table(c3, useNA="always")
## using dplyr
group_by(as.data.frame(c1),c1)%>%summarise(count=n())
group_by(as.data.frame(c2),c2)%>%summarise(count=n())
group_by(as.data.frame(c3),c3)%>%summarise(count=n())
## tally can be used instead of summarise(count=n())

## Check II ############################



## to view a protein group - select a master protein and show all
featureNames(e)[featureNames(e) == "AT2G39730.1"]
featureNames(e)[featureNames(e) == "35a12"]
featureNames(e)[featureNames(e) == "AT5G64570.1"]

## find an accession in rownames (should be master protein) and filter feature data
##DF.new -> DF %>% filter(row.names(DF) %in% c("12a","13a"))
View(fData(e)[,c(m,a,pg,s,su,gs,pc,pgg,gsg)] %>%
       filter(row.names(fData(e)) %in% featureNames(e)[featureNames(e) == "35a12"]))

## seaching for an accession in global group in one item  
pgg[1]
fData(e)[1:2,pgg[1]]
l <- fData(e)[1:3,pgg[1]]
class(l)
l
l <- unlist(strsplit(l,";"))
grep("77",l)
l[2]


k <- sapply(fData(e)[, pgg[1]], paste, collapse=";")
k[2]
class(k[1])
ka <- t(apply(k, 1, function(x) unique(na.omit(as.character(x)))))
ka
apply(fData(e)[1, pgg[1]], 1,
           function(x) unique(unlist(strsplit(na.omit(as.character(x)),";"))))

## find all GFP proteins
find_string <- "35a12"
test <- fData(e)[,c(m,su,gs,pc,pg,pgg[1:2])][grep(find_string,fData(e)[,pgg[1]]),]
#test <- fData(e)[,c(m,su,gs,pc,pg,pgg[1:2])][grep(find_string,fData(e)[,a]),]
dim(test)
test
View(test)

## test of GFP proteins
l3 <- lapply(as.list(test[3,grep("pepSeqUnq",colnames(test),value = FALSE)[1:3]]),function(x) strsplit(unique(x),";"))
l2 <- lapply(as.list(test[2,grep("pepSeqUnq",colnames(test),value = FALSE)[1:3]]),function(x) strsplit(unique(x),";"))
l1 <- lapply(as.list(test[1,grep("pepSeqUnq",colnames(test),value = FALSE)[1:3]]),function(x) strsplit(unique(x),";"))

l3
l2
l1

l11 <- unlist(l1[[1]])
l21 <- unlist(l2[[1]])
l23 <- unlist(l2[[3]])
l32 <- unlist(l3[[2]])

l11
l32
l213 <- unique(c(l21,l23))

#unique(c(l11,l21,l23,l32))
l11
l32
setdiff(l11,l32)
setdiff(l32,l11)    # nothing unique in 32... it is subset of 11

setdiff(l11,l213)
setdiff(l213,l11)

setdiff(l11,c(l32,l213))


setdiff(l213,l32)
setdiff(l32,l213)  # nothing unique in 32... it is subset of 213

## sample peptide, FEGDTLVNR that should be unique, it is?

x <- grep(".?pepSeq.?", colnames(df))
df[1:10,26]
df$pepSeq[1:10]
df$pepSeq=="FEGDTLVNR"
table(df$pepSeq=="FEGDTLVNR")
df$accession[df$pepSeq=="FEGDTLVNR"]          


length(unlist(pel[names(pel)=="FEGDTLVNR"]))





## Check III ############################ 


View(fData(e)[1:500,c(24,66,108,42,84,126,39,81,123,26,68,110)])
test <- fData(e)[,c(24,66,108,42,84,126,39,81,123,26,68,110)]
test <- fData(e)[,c(24,66,108,150,42,84,126,168,39,81,123,165,26,68,110,152)]
View(head(test))
## accesion and master should not be the same
identical(fData(e)[,a[1]],fData(e)[,m[1]])
ni <- which(fData(e)[,a[1]] != fData(e)[,m[1]])
## rownames and masters should be the same
ac <- fData(e)[,c(a,m,g)]
head(ac)
table(rownames(ac)==ac[,4])
table(rownames(ac)==ac[,5])
table(rownames(ac)==ac[,6])
table(rownames(ac)==ac[,7])
head(ac[,7])
## rownames and first accession in protein groups shoudl be the same 
l <- strsplit(ac[,7],";")
sl <- sapply(l, "[[", 1)
table(rownames(ac)==sl)

l <- strsplit(ac[,8],";")
sl <- sapply(l, "[[", 1)
table(rownames(ac)==sl)

l <- strsplit(ac[,9],";")
sl <- sapply(l, "[[", 1)
table(rownames(ac)==sl)

## where is accession "04266"?
unique(df[,c(24,32)][which(df$accession=="04266"),])






## Check IV ############################
## just one replicate of a sample
fData(e)[1,grep("pep_counts.C_H_4",fvarLabels(e))]
strsplit(fData(e)[1,grep("pep_counts.C_H_4",fvarLabels(e))],";")
x <- as.data.frame(strsplit(fData(e)[1,grep("pep_counts.C_H_4",fvarLabels(e))],";"))
xn <- as.numeric(x[,1])
xn
sum(xn)
## how expression values look?
exprs(e)[1,"C_H_4"]
## too bad!!!


## Check V ############################
## protein groups
length(unlist(strsplit(fData(e)[1,"protein_group_global"],";")))
## OK

## Check VI ###########################
de <- grep("isDecoy.\\.?",fvarLabels(e))
de
table(fData(e)[,de])




## how keep the best peptide feature data?
## by combineFeatures further up  

#save modified MSnSet
saveRDS(e,"e.rds")

######################################################
######################################################
######################################################
stop("Never mind the error, execution stops here!")###
######################################################
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




## 2 Dec 2019
## PROBLEM: sum of pep_counts.* do not match to exprs(e) values for proteins
## The numbers in exprs are higher than in Scaffold
## e.g
fData(e)["35a12","pep_counts.C_H_3"]
exprs(e)["35a12",]


