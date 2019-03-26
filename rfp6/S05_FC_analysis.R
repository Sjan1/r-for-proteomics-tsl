##
## Fold change analysis ##
## ==================== ##
## 1. set a filter to select for quality data
## 2. calculte mean value and fold change
## 3. set threshold and color S-shaped bar plot
## 4. provide list of proteins
## 5. link to relevant peptides
rm(list = ls());
if(!exists("foo", mode="function")) source("footnote.R")
## load MSnSets
## protein data
#eprot <- readRDS("eprot_mascot.rds")
#eprot <- readRDS("eprot_mascot_fdr1pc.rds")
eprot <- readRDS("eprot_m.rds")
## peptide data
## It needs to have Prot_Acc column (calculated in S03_code_prot.R)
#e <- readRDS("e_mascot.rds")
e <- readRDS("e_mascot_fdr1pc.rds")


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

## FOR FUTURE IMPROVEMENTS
## create labels for the data compared and link them to the pattern used above 
treat <- gsub("[_]","",names(interested[2]))
ctrl <- gsub("[_]","",names(interested[1]))
treat
ctrl
  
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
## selected protein hits with at least two peptides in one group 
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

## replacing zeros
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
# swap some columns
m <- m[,c(2,1,3,4:ncol(m))]
head(m)


## Calculate unique peptide couints - UPC
## Read MSnset with peptides in it
## It needs to have Prot_Acc column (calculated in S03_code_prot.R)
#e <- readRDS("e_mascot.rds")
#e <- readRDS("e_mascot_fdr1pc.rds")

f <- fData(e)
f[,grepl("pepSeq",colnames(f))]
f
colnames(mdf)

m2 <- m;
for(acc in rownames(m)) {
for(cn in colnames(m[4:ncol(m)])) {
  re <- paste0("pepseq\\.", cn);
  chcol <- colnames(f)[grepl(re, colnames(f), ignore.case = T)]
  #ch_pep <- f[f$Prot_Acc == acc, c(chcol)]
  ch_pep <- f[f$accession == acc, c(chcol)]
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
# apply the function
m2$accept <- apply(as.matrix(m2[, 4:ncol(m2)]), 1, upcfun)
head(m2)
m2[m2$accept == FALSE,]
m2[m2$accept == TRUE,]
## filtered upc - unique peptide count 
upc <- m2[m2$accept == TRUE,]
# adjust spc in m in the same way - filter m to have the same accession as upc

if (identical(rownames(m), rownames(m2))){
m$accept <- m2$accept  
} else{
stop("m and m2 are not identical")  
}
## filtered spc - total spectral count
spc <- m[m$accept == TRUE,]
spc[,"accept"] <- NULL ## now remove accept column when no longer needed 
head(spc)
## lists of all enriched (up/down-regulated) proteins
res1 <- spc[spc$lfc >= 0,]
res2 <- spc[spc$lfc < 0,]
# re-order res2 to have the lowes values on top
ordvec <- order(res2$lfc, decreasing = FALSE)
res2 <- res2[ordvec,]
nrow(res1)
nrow(res2)
head(res1)
head(res2)
#format(round(res1, 2), nsmall = 2)
#format(round(res2, 2), nsmall = 2)
res1 <- format(round(res1, 2))
res2 <- format(round(res2, 2))

## resulting list -  IT WOULD BE NICE TO ADD DESCRIPTION TO THE 'res' TABLES 
head(res1,20)
head(res2,20)
write.csv(res1, file = paste0("res1_mascot",treat,".csv", collapse = NULL))
write.csv(res2, file = paste0("res2_mascot",ctrl,".csv", collapse = NULL))

## remains to do - add description to the result table
head(fData(eprot)[,c("description.C_H_2", "description.C_H_3")],10)

## venn diagram - table
uch <- length(cols[spc$ph == 0])
uph <- length(cols[spc$ch == 0])
common <- length(spc$lfc)-uph-uch
uph
common
uch
## IT WOULD BE NICE TO PLOT VENN DIAGRAM HERE

## plot
#windows()
#png(filename='Figs/eprot_mascot_S-curve.png')
barplot(spc$lfc)
#cols <- c("blue", "red")[(abs(spc$lfc) >= 2) + 1]
cols <- rep("red", nrow(spc))
cols[abs(spc$lfc) >= 2] <- "yellow"
cols[spc$ch == 0] <- "green"
cols[spc$ph == 0] <- "green"

barplot(spc$lfc,
        col=cols,
        space = c(0,0),
        ylim = c(min(spc$lfc),
        max(spc$lfc)),
        main=paste("Spectral count fold change (",treat,"/",ctrl,")"),
        xlab="protein hits",
        ylab=paste(c("log2 (",treat,"/",ctrl,")"),collapse = "")
        )

## plot labels
for(rn in c(1, nrow(spc))) {
temp <- as.matrix(spc[rn,1:3])[1,]
cn1 <- names(temp)[1];
cn2 <- names(temp)[2];
md <- ceiling(nrow(spc[spc$lfc > 0,]) + abs(nrow(spc[spc$lfc > 0,])-nrow(spc[spc$lfc < 0,]))/2)

if(temp[cn1] > temp[cn2]) {
    labpos = cn1;
} else {
    labpos = cn2; 
}
if (temp["lfc"] > 0) {
  text(md, 5, label = labpos, pos=4, offset = -5);
  } else {
  text(md, -5, label = labpos, pos = 4, offset = 5);
}
}
## unique hit numbers (like in Venn diagram)
text(md, -10, label = paste0(treat,"-unique =", uph,
                             ",  common = ",common,",  ",
                             ctrl,"-unique =", uch), pos = 1, offset = 0);
makeFootnote(footnote)
makeFootnote(footnote)

#dev.copy(png,'LFC_mascot_fdr1pc.png')
#dev.copy(png,'LFC_mascot.png')
#dev.off()


## Visualize input from a collaborator (cherry-picked proteins)

## PDLP1
res <- which(rownames(m2) == "AT5G43980.1")
#res <- c(150,200,220)  #test
#abline(v=res)
segments(x0=res,y0=-4,x1=res,y1=0,lwd=1,lty="dotted", col="black")
#arrows(x0=res,y0=-10,x1=res,y1=0,lwd=5,lty="dotted", col="grey")
text(x=res, y=-5, labels="AT5G43980")

## RUBISCO
res <- which(rownames(m2) == "ATCG00490.1")
segments(x0=res,y0=-6,x1=res,y1=0,lwd=1,lty="dotted", col="black")
text(x=res, y=-7, labels="ATCG00490")

## VAMP721
res <- which(rownames(m2) == "AT1G04750.1")
segments(x0=res,y0=-6,x1=res,y1=0,lwd=1,lty="dotted", col="black")
text(x=res, y=-7, labels="AT1G04750")

## print to file
dev.copy(png,filename = "LFC_mascot.png",
         width = 680,
         height = 480,
         units = "px",
         bg = "white")
dev.off()

stop("never mind the error, execution stops here")



## 2-plot several accesion stored in 'cherry' 
cherry <- c("AT1G02500.1","AT5G43980.1","AT1G23410.1")
res <- NULL
for(i in 1:length(cherry)) {
   res <- which(rownames(m2) == cherry[i])
   abline(v=res)
   text(x=res, y=11, labels=cherry[i])
}



match(0,spc$lfc)  # a test



## FURTHER IMPROVEMENTS FOR CONSIDERATION
## ======================================
## - DISTRIBUTE PROTEIN LABELS WITHOUT OVERLAPPING
## - USE ctrl, treat VARIABLES CONSISTENTLY THROUGHOUT THE SCRIPT 
## - MAKE IT UNIVERSAL FOR ANY PAIRWISE COMPARISON
## - WOULD SCATTEP PLOT HAVE SOME ADVANTAGE?
## - FILL THE BARS PROPORTINALLY TO A NUMBER OF REPLICATES
## - add protein description to 'cherry_desc' and swap it with cherry (accessions)


## MA-plot
## better to 

ph <- grep("P_H", colnames(df))
ch <- grep("C_H", colnames(df))


rsch <- rowSums(df[, ch])
rsph <- rowSums(df[, ph])

cls <- rep("black", nrow(df))
cls[rsch == 0] <- "red"
cls[rsph == 0] <- "blue"

## na values were removed by addition of 1 to enable fold change calculation
plot((rsph + rsch + 2)/2, log2((rsph + 1)/(rsch + 1)), 
     col = cls,
     log = "x")
grid()



