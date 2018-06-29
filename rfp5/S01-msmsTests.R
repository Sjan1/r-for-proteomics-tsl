library("msmsTests")

vignette(package = "msmsTests")

vignette("msmsTests-Vignette", 
         package = "msmsTests")

## to select the vignette of interest
openVignette() 

data(msms.dataset) 

dim(msms.dataset)
pData(msms.dataset)

table(msms.dataset$treat, msms.dataset$batch)

library(msmsEDA)
counts.pca(msms.dataset)

null.f <- "y ~ batch"
## alt.f <- "y ~ treat"
alt.f <- "y ~ treat + batch"
### Normalizing condition
div <- apply(exprs(msms.dataset),2,sum)

## run the negative binomial test for 
## count data
res <- msms.edgeR(
  msms.dataset,
  alt.f,
  null.f,
  div = div,
  fnm = "treat")

## adjust the p-values using the 
## Benjamini Hochberg method
res$padj <- p.adjust(res$p.value, method = "BH")

## a data.frame
head(res)

## check that the order has been maintained
identical(featureNames(msms.dataset), 
          rownames(res))

## add the stats results to the MSnSet's 
## feature data
fData(msms.dataset) <- 
  cbind(fData(msms.dataset), res)
