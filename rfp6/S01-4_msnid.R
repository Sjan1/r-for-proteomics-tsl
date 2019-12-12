rm(list = ls());
source("S00-env.R")
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

psml <- lapply(index, function(.x){
    .files <- mzid_files[[.x]]
    msnid <- MSnID()
    msnid <- read_mzIDs(msnid, .files)
    #msnid <- run_msnid(.files)
    
    psm.dt <- as(msnid, "data.table")
    
    # psm <- rtslprot::as_MSnSet(msnid)
    # df <- fData(psm)
    # 
    # ## TO ADD MISSING INFO WE HAVE TO GO BACK TO PSM MSnSet
    # ## Here, at the level of 'list of MSnSets',
    # ## we modify all the items (MSnSets of individual samples)
    # 
    # ## aggregate protein accessions
    # pe <- aggregate(accession~pepSeq, df, c)
    # ## here we make list of peptide sequences,
    # ## where each can have multiple accesions
    # pel <- as.list(pe$accession)
    # names(pel) <- pe$pepSeq
    # pel <- lapply(pel, unique)
    # gb <- pel[df$pepSeq]
    # names(gb) <- featureNames(psm)
    # 
    # pe2 <- aggregate(pepSeq~accession, df, c)
    # pe2l <- as.list(pe2$pepSeq)
    # names(pe2l) <- pe2$accession
    # pep_count <- sapply(pe2l, function(x) length(unique(x)))
    # 
    # ## new re-ordered protein group
    # gb2 <- gb_count <- gb
    # for (i in seq_along(gb)) {
    #     gb2[[i]] <- gb[[i]][order(pep_count[gb[[i]]], decreasing = TRUE)]
    #     gb_count[[i]] <- paste(pep_count[gb2[[i]]], collapse = ";")
    # }
    # 
    # fData(psm)$protein_groups <- sapply(gb2, paste, collapse = ";")
    # fData(psm)$group_size <- lengths(gb2)
    # fData(psm)$pep_counts <- unlist(gb_count)
    # fData(psm)$master_prot <- sapply(gb2, "[[", 1)
    # 
    # for (k in unique(fData(psm)$master_prot)) {
    #     i <- which(fData(psm)$master_prot == k)
    #     j <- which.max(fData(psm)[i, "group_size"])
    #     ii <- i[j]
    #     fData(psm)[i, c("protein_groups",
    #                     "group_size",
    #                     "pep_counts")] <-
    #         fData(psm)[ii, c("protein_groups",
    #                          "group_size",
    #                          "pep_counts")]
    # }
    # 
    # # e <- combineFeatures(psm, fcol = "master_prot",
    # #                      method = "sum", cv = FALSE)
    # e <- psm
    
    
    ## update samples names
    psm.dt$sample <- .x
    #e <- updateFvarLabels(e, sampleNames(e))
    #if (validObject(e))
      return(psm.dt)
    
})

## mzid tables are now in msnl list

####################
## END ##
####################


## list of MSnSets 'msnl' from all samples completed

#######################################
## 7 - COMBINE all PSM from all samples INTO ONE table
x <- do.call(rbind, psml)
dim(x)
unique(x$sample)
colnames(x)

## create category from sample name
table(gsub("_","",strsplit(x$sample,".\\d"), ))
x$categ <- gsub("_","",strsplit(x$sample,".\\d"), )
colnames(x)
head(x$categ)

saveRDS(x,"PDLP1_PSM_NoFDRfilter.rds")
## open when needed
x <- readRDS("PDLP1_PSM_NoFDRfilter.Rds")
#x <- readRDS("PDLP1_PSM_FDR1.rds")
## create category

## decoys
table(x$isDecoy)
table(x$sample)
table(x[x$isDecoy==FALSE,]$sample)
table(x[x$isDecoy==TRUE,]$sample)
#plot(table(x[x$isDecoy==FALSE,]$sample))

## remove decoys
x <- x[x$isDecoy==FALSE,]

## count samples within categories
df <- x %>% group_by(sample,categ) %>% 
  summarise(totalSPC=n())

## better plot
ggplot(data=df) +aes(x=sample) + 
  geom_bar(aes(weight=totalSPC, fill=categ)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
colnames(x)
y <- as.data.frame(x[,c(14,39,40)])
class(y[,1])
## score distributions
ggplot(y) + 
  aes(y[,1]) + 
  geom_histogram(aes(fill=categ)) +
  facet_wrap(~categ, ncol = 2)


## if we want them separatelly...  
a1 <- psml[[1]][,c("accession","pepSeq","sample")]
a2 <- psml[[2]][,c("accession","pepSeq","sample")]
a3 <- psml[[3]][,c("accession","pepSeq","sample")]
dim(a1)
dim(a2)
dim(a3)
#number of rows
dim(a1)[1] + dim(a2)[1] + dim(a3)[1]

aps <- rbind(a1,a2,a3)
dim(aps)
class(aps)
colnames(aps)
head(aps)

a_p <- aggregate(pepSeq~accession,aps,paste,
                 collapse = ",")

p_a <- aggregate(accession~pepSeq,aps,paste,
                 collapse = ",")
#View(head(a_p))
dim(a_p)
colnames(a_p)
length(a_p[,1])
a_p[[1]][500]
a_p[[2]][500]

dim(p_a)
colnames(p_a)
length(p_a[,1])
p_a[[1]][500]

a_pl <- as.list(a_p$pepSeq)
names(a_pl) <- a_p$accession

p_al <- as.list(p_a$accession)
names(p_al) <- p_a$pepSeq


p_al[500:502]
a_pl[500:502]
length(p_al)
length(a_pl)

p_al <- sapply(p_al, function(x) strsplit(x,","))
p_al[502:504]
p_ual <- sapply(p_al, function(x) unique(x))
p_ual[502:504]

p_al_count <- sapply(p_al, function(x) length(x))
p_al_count[502:504]
p_ual_count <- sapply(p_ual, function(x) length(unique(x)))
p_ual_count[502:504]



a_pl <- sapply(a_pl, function(x) strsplit(x,","))
a_pl[500:502]
a_upl <- sapply(a_pl, function(x) unique(x))
length(a_upl)
# now count them
a_pl_count <- sapply(a_pl, function(x) length(unique(x)))
a_pl_count[500:502]
a_pl[500:502]

a_upl_count <- sapply(a_upl, function(x) length(unique(x)))
a_upl_count[500:502]
a_upl[500:502]

## BOLD RED <=> Exclusively unique peptides
# show accession with count = 1
# working with unique list, repetitive strings must be removed
length(p_ual)
table(p_ual_count==1)
length(p_ual[p_ual_count == 1])
up <- (p_ual[p_ual_count == 1])
up <- unlist(up)
up <- as.data.frame(sapply(up,"[[",1))
class(up)
names(up)
updf <- data.frame("pepSeq" = names(up),"Accession" = up)
updf[1:3,]
dim(updf)
length(unique(updf$pepSeq))
length(unique(updf$Accession))

## GROUPS
a_p <- aggregate(pepSeq~accession,aps,paste,collapse=",")
head(a_p)
# accession with multiple peptides
a_p[grepl(",",a_p$pepSeq),]$pepSeq
# accession with one peptide only
a_p[!grepl(",",a_p$pepSeq),]$pepSeq


## BOLD RED <=> EXCLUSIVE UNIQUE PEPTIDES
p_a <- aggregate(accession~pepSeq,up,paste,collapse=",")
head(p_a)
# accession with multiple peptides
p_a[grepl(",",p_a$pepSeq),]$accession
# accession with one peptide only
View(p_a[!grepl(",",p_a$pepSeq),])
test <- p_a[!grepl(",",p_a$pepSeq),]
dim(test)
length(unique(test$accession))


a_upl[["35a12"]]
p_a[grep("35a12",p_a$accession),]
#  35a12
#  AT5G64570.1
## maybe we should rename variables to be clear how they are aggregated
## all variant of  many Peptide ~ unique Accession lists 
AccAll_PepAll    # all data in to count
AccAll_PepUnq    # to see coverage
AccAll_PepUnqEx    # to see what is specific
AccEx_PepUnq    #  the best proteins
AccEx_PepUnqEx    # the best protein with the best peptides 

# Acc1  Pep1, Pep2
# Acc2  Pep3, Pep2
# Acc3  Pep1
# Acc4  Pep4, Pep1
#
# Pep1  Acc1, Acc3, Acc4
# Pep2  Acc1, Acc3
# Pep3  Acc2
# Pep4  Acc4
#
## Problems
## 1- label and count unique exclusive (bold red) peptides to select the best accessions.
## 2- use (count) low scoring hits to help decide protein inference problem.
## 3- we select master protein based on highest number of unique peptides, but
## maybe we should look at number of bold reds as well and perhaps even at number of low scoring hits

## list opf groups for each accession
# 1) generate manyAcc ~ unqPep aggregated list aggregate(acc~pep,...)
# 2) remove duplicity in Accessio1n  
# 3) replace peptides in manyPep ~ unqAcc list aggregate(pep~acc,...)
#    with accessions from 2)
# 4) remove peptide duplicity from 3)
# 5) the result is the list of protein_groups for each unique Accession

## how to make BOLD REDs
# 1) count number of accessions in manyAcc~unqPep aggregated list
# 2) where count = 1, the peptide hits only one accession <=> BOLD RED

## 10;9;9;9;8;8;8;2;1
## acc = Bold red count


######################################
## 9 - SAVE IT
## keep the original
saveRDS(e,"e.rds")
## open when needed
e <- readRDS("e.Rds")
## save as
saveRDS(e,"psm.rds")
psm <- readRDS("psm.rds")


## check it
i <- grep(".?accession.?", fvarLabels(e))
i
t <- fData(e)[,i]
class(t)
head(t)
i <- grep(".?master.?", fvarLabels(e))
i
u <- fData(e)[,i]

#how to join u and t?

View(head(t))

i <- grep(".?protein_groups.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))
i <- grep(".?pepSeq.?", fvarLabels(e))
i
View(head(fData(e)[,i],100))
#fData(e)[,i][1,4]
## check all manually...
View(fData(e)[1:500,c(24,66,108,42,84,126,39,81,123,26,68,110)])


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
j <- grep("pepSeq\\.?", fvarLabels(e))    # e is a peptide-level MSnSet
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
