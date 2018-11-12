## Date: 16th-July-2018
## Parsing SmplExp Table
## We are creating a unique stting that link raw file names with experiments, 
## then we shall cut this string to find which groups of files form samples, 
## replicates and how they should be combined and how to assemle the proteins.   
rm(list=ls())
source("S00-env.R")
## this creates tibble
##tab <- readr::read_csv("data/SampleExperimentTable_fixed_Mascot.csv", comment = "#",
##                       col_names = TRUE)
##this creates dataframe
tab <- read.table("data/SampleExperimentTable_fixed_Mascot.csv", comment = "#",
                       header = TRUE, sep=",")
tab <- tab[tab[[1]]!= "not_used",]
head(tab)

if(colnames(tab)[1]=="not_used"){
names(tab)[2:3] <- c("file", "raw")
} else{
names(tab)[1:2] <- c("file", "raw")}

temp <- colnames(tab)
if(any(grepl("uniq",temp))){
  print("Uniq already exists!")
} else{
  
  head(tab)
  rownames(tab) <- tab$file
  ## here we combine columns to create unique sample labels
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
  ## to attach unique string to the original table
  tab$uniq <- labs
}
write.table(tab,"data/SamplesExperimentsTable_edited.csv",col.names = FALSE)



ndf <- data.frame(file=fns, labs=labs, stringsAsFactors = FALSE)
saveRDS(ndf, "data/tab_unq")
write.table(ndf,"data/tab_unq.txt", sep="\t", row.names = TRUE, col.names = TRUE,
            quote=FALSE)
rm(ndf)

rdf <- readRDS("data/tab_unq")
head(rdf)

