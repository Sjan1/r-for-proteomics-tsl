library("pRolocdata")
f0 <- dir(system.file("extdata", 
                      package = "pRolocdata"),
          full.names = TRUE,
          pattern = "Dunkley2006")

## f0 is a files exported by ProteomeDiscoverer

getEcols(f0, split = ",")
i <- grepEcols(f0, "M", split = ",")
## i <- 5:20

x <- readMSnSet2(f0, i)

## From raw
?quantify

## raw data example, can be read in with readMSData
data("itraqdata")
itraqdata ## MSnExp (raw data)
e <- quantify(itraqdata, 
              reporters = iTRAQ4, 
              method = "max")
## e is an MSnSet (quantitative data)

## combineFeatures(e, fcol = "accession")

fData(msexp)$DatabaseAccess

e2 <- quantify(msexp, method = "count")
fData(e2)$DatabaseAccess[2] <- fData(e2)$DatabaseAccess[1] 

e3 <- combineFeatures(e2, fcol = "DatabaseAccess", fun = "sum")
MSnbase::exprs(e3)
