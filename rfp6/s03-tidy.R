source("https://bioconductor.org/biocLite.R")
biocLite("biobroom")
library(biobroom)
library("ggplot2")

if (require("MSnbase")) {
  library(MSnbase)
  # import MSnSet object

  #already imported...
  
  # Use tidy to extract genes, sample ids and measured value
  tidy(eprot)
  # add phenoType data
  et <- tidy(eprot, addPheno=TRUE)
}
saveRDS(et,"et.rds")
et <- readRDS("et.rds")

et_c <- et %>% filter(phenotype == "C")
et_p <- et %>% filter(phenotype == "P")
et_h <- et %>% filter(treatment == "H")
et_i <- et %>% filter(treatment == "I")

#ggplot(et_p) + aes(x=sample) + geom_bar()

ggplot(et_c) + aes(x=sample) + geom_bar(aes(weight=value))
ggplot(et_p) + aes(x=sample) + geom_bar(aes(weight=value))
ggplot(et_h) + aes(x=sample) + geom_bar(aes(weight=value))
ggplot(et_i) + aes(x=sample) + geom_bar(aes(weight=value))