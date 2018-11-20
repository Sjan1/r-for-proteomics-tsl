library("MSnbase")

source("https://bioconductor.org/biocLite.R")
biocLite("MSnbase")

f <- "./data/mulvey2015.csv"
e <- 2:19

m <- readMSnSet2("./data/mulvey2015.csv", 
                 ecol = 2:19)
m

m2 <- m[1:10, 1:5]

exprs(m2)
fData(m2)

sampleNames(m)
featureNames(m) <- fData(m)$FeatureNames
exprs(m)[1:5, 1:3]

m <- readMSnSet2("./data/mulvey2015.csv", 
                 ecol = 2:19, 
                 fnames = 1)
exprs(m)[1:5, 1:3]


