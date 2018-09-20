## Date: 16th-July-2018
## Parsing SmplExp Table
## We are creating a unique stting that link raw file names with experiments, 
## then we shall cut this string to find which groups of files form samples, 
## replicates and how they should be combined and how to assemle the proteins.   

source("S00-env.R")

tab <- readr::read_csv("data/SampleExperimentTable_fixed.csv")
names(tab)[1:2] <- c("file", "raw")
tab
rownames(tab) <- tab$file

fns <- c()
labs <- c()
line_cnt <- 0
for(file in rownames(tab)){
    line_cnt <- line_cnt +1  
    row <- tab[file,]
    new_lable <- paste0(row$sample,row$treatment,row$biorep,row$cleavage,row$fraction)
    cat(line_cnt, file, new_lable, "\n")
    fns <- c(fns,file)
    labs <- c(labs,new_lable)
}
ndf <- data.frame(file=fns, labs=labs, stringsAsFactors = FALSE)
saveRDS(ndf, "data/tab_unq")
write.table(ndf,"data/tab_unq.txt", sep="\t", row.names = TRUE, col.names = TRUE,
            quote=FALSE)
rm(ndf)

rdf <- readRDS("data/tab_unq")
head(rdf)


# Here we neeed to bring in the data (like in the scripts S02, S03)

fn <- dir("MSGF", full.names = TRUE)

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

  list_msnsets[[biorep]] <- e
}

## combine all into one msnset
msnset = BiocGenerics::do.call(combine,list_msnsets)

#testing
e1 <- list_msnsets[[1]]
e2 <- list_msnsets[[2]]
c <- do.call(combine,list(e1,e2))

e1_1 <- list_msnsets[[1]][1:10]
e1_2 <- list_msnsets[[1]][11:21]
c <- do.call(combine,list(e1_1,e1_2))

e11 <- e[1:10,]
e22 <- e[11:21,]
sampleNames(e11) <- c("sample11")
sampleNames(e22) <- c("sample22")

e11 <- updateSampleNames(e11)
e22 <- updateSampleNames(e22)

e11 <- updateFeatureNames(e11)
e22 <- updateFeatureNames(e22)

e11 <- updateFvarLabels(e11)
e22 <- updateFvarLabels(e22)

c <- do.call(combine,list(e11,e22))

