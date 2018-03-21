library("msdata")
f <- proteomics(full.names = TRUE)
f2 <- f[2]

rw <- readMSData(f2, mode = "onDisk")

fvarLabels(rw)
length(rw)

msLevel(rw)

rw[[1]]
rw[[2]]
rw[[6]]

centroided(rw)

cnt <- isCentroided(rw)
isCentroidedFromFile(rw)

table(cnt, msLevel(rw))

centroided(rw) <- cnt

table(centroided(rw), 
      msLevel(rw))


rw <- readMSData(f2, 
                 mode = "onDisk",
                 centroided = c(FALSE, TRUE, FALSE))

?addIdentificationData

sp1 <- rw[[1]]
plot(sp1)

plot(rw[[2]], full = TRUE)
plot(rw[[2]], full = TRUE)

plot(rw[[6]], 
     reporters = TMT10,
     full = TRUE)

fileNames(rw)
