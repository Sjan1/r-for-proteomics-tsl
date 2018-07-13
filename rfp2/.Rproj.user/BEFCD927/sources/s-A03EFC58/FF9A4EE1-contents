## vector

x <- 1:10
letters
y <- c(TRUE, FALSE, TRUE)

x[3]
letters[1:4]

## data.frame
df <- 
  read.csv("./data/mulvey2015.csv",
           header = TRUE,
           stringsAsFactors = FALSE)

class(df)
dim(df)
ncol(df)
nrow(df)

dfdim <- dim(df)
dfdim[2]

dim(df)[2]


dim(df[1:4, ])
dim(df[, 1:4])
dim(df[1:10, 1:4])

colnames(df)
df[, 1]
df[, "FeatureNames"]
df$FeatureNames

## installing packages
# install.packages("pkg")

source("http://www.bioconductor.org/biocLite.R")

library("BiocInstaller")
biocLite("limma")

## MSnSet
library("MSnbase")
?readMSnSet2

f <- "./data/mulvey2015.csv"
e <- 2:19

x <- readMSnSet2(f, e)
class(x)
dim(x)

exprs(x)[1:4, 1:3]
fData(x)[1:4, 1:2]
pData(x)

sampleNames(x)
fvarLabels(x)
featureNames(x)

x <- readMSnSet2(f, e, 
                 fnames = 1)
featureNames(x)[1:5]

x <- readMSnSet2(f, e)
featureNames(x) <- 
  fData(x)$FeatureNames

## saving and loading
save(x, 
     file = "./data/msnsets.rda")

ls()
rm(x)
ls()

load("./data/msnsets.rda")
ls()

## plotting

set.seed(123)
x <- 1:10
y <- rnorm(10)

plot(x, y)

## produce a scatter plot
## of replicates 1 and 2 
## of timepoint 72

load("./data/msnsets.rda")
ms <- x

sampleNames(ms)

i <- 5
j <- 11
plot(exprs(ms)[, i],
     exprs(ms)[, j])


plot(exprs(ms)[, i],
     exprs(ms)[, j], 
     main = "72 hr",
     xlab = "rep1",
     ylab = "rep2")

log10(1:10)

plot(log10(exprs(ms)[, i]),
     log10(exprs(ms)[, j]))

i <- exprs(ms)[, 5]
j <- exprs(ms)[, 11]

plot(i, j, 
     col = c("red", "blue"))

View(fData(ms))

m <- fData(ms)$markers
table(m)

cls <- rep("black", 
           nrow(ms))

plot(i, j, col = cls)
cls <- rep("black", 
           nrow(ms))

cls[m == "pluri"] <- "red"
cls[m == "diff"] <- "blue"
table(cls)

cls <- as.character(fData(ms)$markers)
cls[cls == "pluri"] <- "red"
cls[cls == "diff"] <- "blue"
cls[cls == "unknown"] <- "black"

plot(i, j, col = cls)

hist(exprs(ms)[, 5])

boxplot(exprs(ms)[, 5])

boxplot(exprs(ms))

