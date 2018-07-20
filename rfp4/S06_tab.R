## Date: 16th-July-2018
## Parsing SmplExp Table
## We are creating a unique stting that link raw file names with experiments, 
## then we shall cut this string to find which groups of files form samples, 
## replicates and how they should be combined and how to assemle the proteins.   

source("S00-env.R")

tab <- readr::read_csv("data/SampleExperimentTable.csv")
names(tab)[1:2] <- c("file", "raw")
tab
rownames(tab) <- tab$filea

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
rm(ndf)

rdf <- readRDS("data/tab_unq")
head(rdf)


# Here we connec to scripts S02, S03
msnid <- MSnID()
msnid <- read_mzIDs(msnid, fn)

labs <- rdf$labs
names(labs) <- rdf$file
labs["PDLP1_wt-1t_1_1_1_111206"]
labs_sub <- labs[sub("\\.mzML$","",msnid$spectrumFile)]

msnid$spectrumFile <- labs[sub("\\.mzML$","",msnid$spectrumFile)]
# we should add a new column with our special string that describes smpl-exp
