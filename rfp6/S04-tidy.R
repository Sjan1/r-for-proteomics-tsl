source("https://bioconductor.org/biocLite.R")
#biocLite("biobroom")
BiocManager::install("biobroom")
library(biobroom)
library("ggplot2")
library("tidyverse")
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

## split data (for stats and plots)
etg <- et %>% group_by(phenotype)
head(etg)
dim(etg)
dim(et)
etg <- et %>% group_by(phenotype) %>% summarise(value = mean(value, na.rm = TRUE))
etg <- et %>% group_by(phenotype) %>% summarise(value = n())
etg <- et %>% group_by(treatment) %>% summarise(value = n())
etg <- et %>% group_by(phenotype,treatment) %>% summarise(value = n())
# average spectral count
etg <- et %>% group_by(protein,phenotype,treatment) %>% summarise(value = mean(value, na.rm = TRUE))
etg <- et %>% group_by(protein,phenotype,treatment) %>% summarise(value = mean(value, na.rm = TRUE))


# wider format for scatter plot
ec  =  mutate(et,categ = paste(phenotype, treatment, sep = '_'))
ecg <- ec %>% group_by(protein,categ) %>% summarise(value = mean(value, na.rm = TRUE))
ecgw <- ecg %>% spread(key = "categ", value = "value", fill = NA)
ecgw <- ecgw %>% replace(is.na(.),0)
ggplot(data = ecgw, mapping = aes(x=C_H,y=P_H)) + geom_point()
head(ecgw)
ecgw <- as.data.frame(ecgw)
# remove protein column, keeping protein info as rownames
rownames(ecgw) <- ecgw$protein
ecgwm <- ecgw[,2:ncol(ecgw)]

# scatter plot
head(ecgw[,2:ncol(ecgw)])
plot(ecgw[,2:ncol(ecgw)])

# correaltion
cr <- cor(ecgwm)
cr

library(corrplot)
corrplot(cr, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

library("PerformanceAnalytics")
# my_data <- mtcars[, c(1,3,4,5,6,7)]
chart.Correlation(cr, histogram=TRUE, pch=19)

# Get some colors
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cr, col = col, symm = TRUE)

# heatmap...

