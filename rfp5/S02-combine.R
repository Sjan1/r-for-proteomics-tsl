library("MSnbase")
library("BiocInstaller")
biocLite("rlang")

x1 = msnset[1:10, ]
x2 = msnset[11:30, ]
x3 = msnset[31:55, ]

l = list(x1, x2, x3)

res = do.call(combine, l)

## combine(combine(combine(x1), x2), x3)

## EXAMPLE OF WEIGTED SPECTRAL COUNT
## Jan's weighing function - sums*counts, more is better + a big reward for reproducibility
## Laurent suggested  better to weight peptides than proteins

## let's generate some data
## dataset a
set.seed(11)
m <- matrix(sample(0:5, 12, replace = TRUE), nrow = 4)
m[4, 1] <- 2
m
colnames(m) <- LETTERS[1:3]
rownames(m) <- letters[1:4]
m
a <- readMSnSet2(as.data.frame(m), eco = 1:3)
exprs(a)
fData(a)$notZero <- rowSums(exprs(a) != 0)
fData(a)$rws <- rowSums(exprs(a))
fData(a)$fac <- fData(a)$rws * fData(a)$notZero
fData(a)
## dataset b
set.seed(1)
m <- matrix(sample(0:5, 12, replace = TRUE), nrow = 4)
colnames(m) <- LETTERS[10:12]
rownames(m) <- letters[2:5]
b <- readMSnSet2(as.data.frame(m), eco = 1:3)
fData(b)$notZero <- rowSums(exprs(b) != 0)
fData(b)$rws <- rowSums(exprs(b))
fData(b)$fac <- fData(b)$rws * fData(b)$notZero

## The following changes expression values.
## Do we need to mutiply peptide spectral counts and then use a mean value to compare sample groups?
## It is certainly a possibility.
## Jan  have been using calculated "fac" as it is to compare sample groups.
exprs(a) <- exprs(a) * fData(a)$fac
exprs(b) <- exprs(b) * fData(b)$fac

## this resolves combining the feature variables that should not merge
## e.g. score.A, score.B
a <- updateFvarLabels(a, "A")
b <- updateFvarLabels(b, "B")
fData(a)
## This combines peptides into unique list based on the sequence.
## ???Coud we combine based on sequence and PTMs? Or even Seq-PTM-charge???
x <- combine(a, b)
exprs(x)
fData(x)
pData(x)
## ading some proteins - this will be with the real data by the search engine
## 1-1:  for every PROTEOTYPIC peptide matching only one protein 
## This works for single protein match per peptide
fData(x)$acc <- c("P1", "P1", "P2", "P2", "P3")
fData(x)
## 1-2:  for combining peptides to multiple proteins  <=> works for multipe proteins
## This works for multiple  proteins matches per peptide
## ???How to make following list from search, comma separated, output???
gbl <- list(c("P1"), c("P1", "P4"), c("P2"),
            c("P2"), c("P3", "P4"))
## to add accession to the msnset 
fData(x)$acc <- list(c("P1"), c("P1", "P4"), c("P2"),
            c("P2"), c("P3", "P4"))
fData(x)
## This replaces NA with zeros
## ???How to replace with e.g. 0.001???
exprs(x)
x <- impute(x, method = "zero")
exprs(x)
## combineFeatures(x, fcol = "acc", method = "sum")
## 1-1: combining PROTEOTYPIC peptides to the corresponding proteins
y1 <- combineFeatures(x, groupBy = fData(x)$acc, 
                     fun = "sum")
exprs(y1)
## 1-2: combining peptides that can match multiple proteins
y2 <- combineFeatures(x, groupBy = gbl, 
                     redundancy.handler = "multiple",
                     fun = "sum")

exprs(y2)

exprs(combineFeatures(x, fcol = "acc", fun = "median"))
exprs(combineFeatures(x, fcol = "acc", fun = "mean"))

## !!!A problem - combineFeatures needs one accession per peptide sequence!!!
## A peptide with several protein accessions on one row, as a vector,
## this metod fails
exprs(combineFeatures(x, groupBy = fData(x)$acc, fun = "median"))


?iPQF
