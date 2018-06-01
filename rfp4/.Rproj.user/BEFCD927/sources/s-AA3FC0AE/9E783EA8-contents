library("MSnID")

## To read the vignette:openVignette("MSnID")

# id1 <- readMzIdData("MSGF/PDLP1_lp-1t_1_1_1_111206.mzid")
# View(id1)


## Create the commands to run MSGF+
for (f in raw_files) {
  cmd <- paste("java -jar ~/bin/MSGFPlus.jar -s", f ,
               "-d TAIR_contaminants.fasta -t 10ppm -tda 1 -inst 0 -thread 48")
  cat(cmd, "\n", file = "msgf.sh", append = TRUE)
}

## New files with decoy data
msgf2_files <- dir("MSGF2", full.names = TRUE)


## Running MSnID

msnid <- MSnID()
msnid <- read_mzIDs(msnid, msgf2_files)
msnid

msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) +
  geom_bar(stat='identity', position='dodge') +
  ggtitle("Number of Missed Cleavages")

msnid$numCys <- sapply(lapply(strsplit(msnid$peptide,''),'==','C'), sum)
msnid$PepLength <- nchar(msnid$peptide) - 4
pepLen <- unique(psms(msnid)[,c("PepLength", "isDecoy", "peptide")])
ggplot(pepLen, aes(x=PepLength, fill=isDecoy)) +
  geom_histogram(position='dodge', binwidth=3) +
  ggtitle("Distribution on of Peptide Lengths")

msnid
msnid <- apply_filter(msnid, "numIrregCleavages == 0")
msnid

ppm <- mass_measurement_error(msnid)
ggplot(as.data.frame(ppm), aes(x=ppm)) +
  geom_histogram(binwidth=100)

msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)

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

msnid0 <- msnid
msnid <- apply_filter(msnid, filtObj.grid)

## make spectral counts table

## 1 id file per sample
msnset <- as(msnid, "MSnSet")

prots <- combineFeatures(msnset, 
                         fcol = "accession",
                         fun = "sum")

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


