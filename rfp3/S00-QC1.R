getwd()

library(readr)
sequence_export <- read_csv("data/sequence_export.csv")
View(sequence_export)

i <- c(1, 47, 190)
fn <- sequence_export[i, 1]
fn


library("MSnbase")
id1 <- readMzIdData("./data/HeLa_180116_m36_r1.mzid")
id2 <- readMzIdData("./data/HeLa_180123_m43_r2.mzid")
id3 <- readMzIdData("./data/HeLa_180124_m43_r1.mzid")


id1 <- filterIdentificationDataFrame(id1)
id2 <- filterIdentificationDataFrame(id2)
id3 <- filterIdentificationDataFrame(id3)


sum(id1$MS.GF.QValue < 0.05)
sum(id2$MS.GF.QValue < 0.05)
sum(id3$MS.GF.QValue < 0.05)

summary(id1$experimentalMassToCharge - id1$calculatedMassToCharge)
summary(id2$experimentalMassToCharge - id2$calculatedMassToCharge)
summary(id3$experimentalMassToCharge - id3$calculatedMassToCharge)



rw <- readMSData("data/HeLa_180116_m36_r1.mzML", mode = "onDisk")
rw <- addIdentificationData(rw, "./data/HeLa_180116_m36_r1.mzid")

table(fData(rw)$isDecoy)
fData(rw)$MS.GF.QValue <- as.numeric(fData(rw)$MS.GF.QValue)

which(fData(rw)$MS.GF.QValue < 0.00001)

sp4027 <- rw[[4027]]
plot(sp4027, centroided = TRUE, full = TRUE)

s <- fData(rw)[4027, "sequence"]

centroided(rw, msLevel = 1) <- FALSE
centroided(rw, msLevel = 2) <- TRUE

plot(rw[[4027]], s, centroided = TRUE)

calculateFragments(s, z = 2)
?calculateFragments

k <- 11075
spk <- rw[[k]]
plot(spk, full = TRUE)

s <- fData(rw)[k, "sequence"]

plot(rw[[k]], s, main = s)

calculateFragments(s, z = 2)


fls <- dir("data", pattern = "mzid$", full.names = TRUE)
fls


rw <- readMSData("data/HeLa_180116_m36_r1.mzML", 
                 mode = "onDisk", msLeve = 2)
rw <- addIdentificationData(rw, "./data/HeLa_180116_m36_r1.mzid")


## unique peptide count
## total spectral count
## protein count

