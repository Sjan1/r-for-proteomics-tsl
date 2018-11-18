## what is in the data
library("tidyverse")

## transform expression values into long format
expval <- as.data.frame(exprs(combined_e))
dim(expval)
head(rownames(expval))
head(expval[1])


expval <- data.frame(seq = row.names(expval), expval)
dim(expval)
head(rownames(expval))
head(expval[1])
#expval <- as.data.frame(c(n,expval))

expval.long <- gather(expval, key = "samples", value = "SPC", -seq)
head(expval.long)
#exprs.long <- select(exprs.long, -samples) 

## adding factors from Exp-Sample table

f <- data.frame(tab$uniq,tab$sample,tab$treatment)
names(f)<-gsub("tab.","",names(f))
f$uniq <- gsub("\\d+$","",f$uniq,perl=TRUE)
f <- unique(f)
expval.long.f <- left_join(expval.long, f,by=c("samples"="uniq"))
dim(expval.long.f)

## bar plot1
ggplot(expval.long.f) + aes(x=samples, y=SPC, fill=sample) + 
  geom_bar(stat = "identity") + 
  labs(title = "Spectral counts per sample") + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5))

## bar plot2
ggplot(expval.long.f) + aes(x=samples, y=SPC, fill=treatment) + 
  geom_bar(stat = "identity") + 
  labs(title = "Spectral counts per sample") + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5))
