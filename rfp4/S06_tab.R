## Date: 16th-July-2018
## Parsing SmplExp Table
## We are creating a unique stting that link raw file names with experiments, 
## then we shall cut this string to find which groups of files form samples, 
## replicates and how they should be combined and how to assemle the proteins.   
rm(list=ls())
source("S00-env.R")

## this creates tibble -  we has a problem here: 
## empty spaces in "not_used" were not read/saved correcly.
#tab <- readr::read_csv("data/SampleExperimentTable_fixed_Mascot.csv", comment = "#",
#                       col_names = TRUE)
## this creates dataframe
tab <- read.table("data/SampleExperimentTable_fixed.csv",
                  header = TRUE, sep=",")
## only the records without "not_used" are allowed  
tab <- tab[tab[[1]]!= "not_used",]
head(tab)

## first three columns must be "not_used", "file", "raw"
## I think we should check the other header exist in the form we expect them
## In the current form this script is tailored to PDLP1
## We should aim for variable DOE into account...
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
## Save the tab that will be used.
#View(tab)
write.table(tab,"data/SamplesExperimentsTable_used.csv", col.names = TRUE, row.names=FALSE)


ndf <- data.frame(file=fns, labs=labs, stringsAsFactors = FALSE)
saveRDS(ndf, "data/tab_unq")
write.table(ndf,"data/tab_unq.txt", sep="\t", row.names = TRUE, col.names = TRUE,
            quote=FALSE)
rm(ndf)

rdf <- readRDS("data/tab_unq")
head(rdf)

## ...conntinue with script S07_msnset_build.R

