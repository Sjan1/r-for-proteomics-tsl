##
## Fold change analysis ##
## ==================== ##
## 1. set a filter to select for quality data
## 2. calculte mean value and fold change
## 3. set threshold and color S-shaped bar plot
## 4. provide list of proteins
## 5. link to relevant peptides

rm(list = ls());
## protein data
eprot <- readRDS("eprot_mascot.rds")

df <- exprs(eprot)
head(df)

## column names - sample repicates
allcol <- colnames(df);
## for the results
sellist <- list();

## which pair of sample types to compare  
interested <- c("C_H", "P_H");
names(interested) <- c("C_H", "P_H")
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
all
## selected proteins hits with at least two peptides in one group 
seldf <- df[all, cols.done]
seldf

## calculate means of selected groups
## what
phc <- grepl("^P_H", colnames(seldf))
phc <- grepl(paste0("^",interested["P_H"]), colnames(seldf))  # beter?
phc
ph <- seldf[, phc]

chc <- grepl("^C_H", colnames(seldf))
chc <- grepl(paste0("^",interested["C_H"]), colnames(seldf))  # beter?
ch <- seldf[, chc]

rownames(ph) <- rownames(seldf)
rownames(ch) <- rownames(seldf)
phm <- rowMeans(ph)
chm <- rowMeans(ch)
## means of both samples we are comparing
mdf <- data.frame(ph = phm, ch = chm)
rownames(mdf) <- rownames(seldf)
head(mdf)

x <- mdf;
head(mdf)
nrow(mdf)
#nrow(y)

x[x$ph == 0, "ph"] <- 0.001;
x[x$ch == 0, "ch"] <- 0.001;

x$ch <- log2(x$ch)
x$ph <- log2(x$ph)
head(x)
x$lfc <- x$ph - x$ch;

ordvec <- order(x$lfc,decreasing = TRUE)
y <- x[ordvec,]
head(y)
# y is in the right order (reordered x) for plotting... lfc small to high
#windows()
barplot(y$lfc)

#adding lfc to mdf
mdf$lfc <- y[rownames(mdf), "lfc"]
head(mdf)
## test of an accession number
y["AT1G02500.1",]

## filter to visualize fold changes
nrow(mdf[abs(mdf$lfc) >= 2, ])
mdf[abs(mdf$lfc) >= 2,]

## color the fold changes
cols <- c("blue", "red")[(abs(y$lfc) >= 2) + 1]
barplot(y$lfc, col=cols)

## more exxamples on how to produce colored bar plot here:
## https://stackoverflow.com/questions/13112974/change-colours-of-particular-bars-in-a-bar-chart

## check for decoys
temp <- grepl("^XXX_", rownames(mdf))
table(temp)

## go back to mgf, subscript to get differentialy abundant proteins mdf[mdf$lfc>1...], 
## lits the proteins for user to evaluate,
## link it with raw spectral counts
## link with petides

## get the list of ordered LFC with the original SPC data
m <- merge(mdf,seldf,by="row.names")
rownames(m) <- m$Row.names;
m <- m[-c(1)]

ordvec <- order(m$lfc, decreasing = TRUE)
m <- m[ordvec,]
head(m)
#swap some columns
m <- m[,c(2,1,3,4:ncol(m))]
head(m)


## Calculate unique peptide couints - UPC
## Read MSnset with peptides in it
## It needs to have Prot_Acc column (calculated in S03_code_prot.R)
e <- readRDS("e.Rds")
f <- fData(e)
f[,grepl("pepSeq",colnames(f))]

colnames(mdf)

m2 <- m;
for(acc in rownames(m)) {
for(cn in colnames(m[4:ncol(m)])) {
  re <- paste0("pepseq\\.", cn);
  chcol <- colnames(f)[grepl(re, colnames(f), ignore.case = T)]
  ch_pep <- f[f$Prot_Acc == acc, c(chcol)]
  ch_pep2 <- !is.na(ch_pep)
  upc <- sum(ch_pep2);
  m2[acc, cn] <- upc;
#break;
  }
#  break;
}

## filter accessinos with UPC < 2 - we want at least 2 unique peptides per protein"
# function
upcfun <- function(x) {
  lgl <- x > 1;
  rowlgl <- any(lgl)
  return(rowlgl);
}
# apply function
m2$accept <- apply(as.matrix(m2[, 4:ncol(m2)]), 1, upcfun)
head(m2)
m2[m2$accept == FALSE,]


#the same plot after filtering
#windows()
barplot(m2$lfc)
cols <- c("blue", "red")[(abs(m2$lfc) >= 2) + 1]
barplot(m2$lfc, col=cols)


## last step
## plot selected accessions

res <- which(rownames(m2) == "AT1G02500.1")
abline(v=res)
text(x=res, y=11, labels="AT1G02500.1")


