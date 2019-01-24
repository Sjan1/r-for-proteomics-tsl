##
## Fold change analysis ##
## ==================== ##
##
rm(list = ls());
## protein data
eprot <- readRDS("eprot.rds")
df <- exprs(eprot)
head(df)

## column names - sample repicates
allcol <- colnames(df);
## for the results
sellist <- list();

## which pair of sample types to compare  
interested <- c("C_I", "C_H");
cols.done <- c();

for(sample in interested) {
sre <- paste("^", sample, sep = "");
  
col.lgl <- grepl(sre, allcol); 

just.this <- as.matrix(df[, col.lgl])
cols.done <- c(cols.done, colnames(just.this))

## filter - criterium of reproducibility; 
## protein must be identified in more than 50% of samples within replicates 
## e.g. at least 2 out of 3 replicated must hit a protein
to.qualify <- ceiling(ncol(just.this) / 2);
to.qualify

# head(just.this)



tempv <- c();
for(rn in seq(1, nrow(just.this))) {
  ro <- just.this[rn,]
  lgl <- ro > 1;
  mto <- table(lgl)["TRUE"]
  if(is.na(mto)) { mto <- 0; }
  if(mto >= to.qualify) {
    tempv <- c(tempv, rn)
  }
  # cat(mto)
#  if(rn > 30) {
#    break()
#  }
}
sellist[[sample]] <- tempv;
}
## list of row numbers list = all from [1] or [2]
sellist
## unnique
all <- unique(c(sellist[[1]], sellist[[2]]));

## selected proteins hits with at least two peptides in one group 
seldf <- df[all, cols.done]
seldf

## calculate means of selected groups
## what
cic <- grepl("^C_I", colnames(seldf))
cic <- grepl(paste0("^",interested[1]), colnames(seldf))  # beter?
ci <- seldf[, cic]

chc <- grepl("^C_H", colnames(seldf))
chc <- grepl(paste0("^",interested[2]), colnames(seldf))  # beter?
ch <- seldf[, chc]

rownames(ci) <- rownames(seldf)
rownames(ch) <- rownames(seldf)
cim <- rowMeans(ci)
chm <- rowMeans(ch)
## means of both samples we are comparing
mdf <- data.frame(ci = cim, ch = chm)
rownames(mdf) <- rownames(seldf)
head(mdf)


