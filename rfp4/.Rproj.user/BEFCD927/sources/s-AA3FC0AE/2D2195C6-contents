fn <- tab %>% 
  filter(biorep == 1 & techrep == 1) %>% 
  filter(category == "CI" & cleavage == "tryp") %>% 
  select(file) 

## one bio-sample/lane
fn <- dir("MSGF2", pattern = "lp", full.names = TRUE)
## these files will be processed by MSnID, will 
## be exported as a single table with as(., "data.table").
## That table will be processed (duplicates, counts, ...)
## to generate a single MSnSet for that sample. Let call it 
## sample 1, and we will want to rename it accordingly with 
## sampleNames(.) <- "sample1"


## another bio-sample/lane
fn2 <- dir("MSGF2", pattern = "wt", full.names = TRUE)
## these files will be processed by MSnID, will 
## be exported as a single table with as(., "data.table").
## That table will be processed (duplicates, counts, ...)
## to generate a single MSnSet for that sample. Let call it 
## sample 2, and we will want to rename it accordingly with 
## sampleNames(.) <- "sample2"

msnid <- MSnID()
msnid <- read_mzIDs(msnid, fn)

msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))

msnid$numCys <- sapply(lapply(strsplit(msnid$peptide,''),'==','C'), sum)
msnid$PepLength <- nchar(msnid$peptide) - 4
pepLen <- unique(psms(msnid)[,c("PepLength", "isDecoy", "peptide")])
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
ppm <- mass_measurement_error(msnid)

## create filter
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
filtObj$msmsScore <- list(comparison=">", threshold=10.0)
show(filtObj)


## optimise filter
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)
show(filtObj.grid)

msnid <- apply_filter(msnid, filtObj.grid)

## different fractions, same lane (all same sample)
td <- as(msnid, "data.table")
## View(td)

sel <- !duplicated(td$pepSeq)
count <- table(td$pepSeq)
x <- td[sel, ]
x$count <- as.vector(count[x$pepSeq])
## View(x)

i <- which(names(x) == "count")
e <- readMSnSet2(x, i)
featureNames(e) <- fData(e)$pepSeq

fvarLabels(e)
eprots <- combineFeatures(e, 
                          fcol = "accession", 
                          fun = "sum")
eprots

sampleNames(eprots) <- "sample1"

## eprots1, named sample1
## eprots2, named sample2
## ...

## efinal <- combine(eprots1, eprots2)

## use tab to populate pData(efinal)
## for example

sampleNames(msnset)
tab[1:3, 1]

pd <- data.frame(tab[1:3, ])
rownames(pd) <- sampleNames(msnset)
pData(msnset) <- pd

msnset$fraction == 1
