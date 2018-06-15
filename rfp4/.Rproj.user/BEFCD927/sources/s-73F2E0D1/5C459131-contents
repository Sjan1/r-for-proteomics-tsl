td ## from MSnID - see S02-id.R or S03-id.R
View(td)

td %>% 
  filter(accession == "AT3G44310.1")

ti <- td[1, ]
ti

f <- paste0("raw/", ti$spectrumFile)
i <- 3378

rw <- readMSData(f, mode = "onDisk")
sp <- rw[[i]]
plot(sp, full = TRUE)


tj <- td %>% 
  filter(accession == "AT3G44310.1") %>% 
  as_tibble

tj$spectrumFile
tj$`scan number(s)`

for (k in seq_len(nrow(tj))) {
  f <- paste0("raw/", tj$spectrumFile[k])
  i <- tj$`scan number(s)`[k]
  rw <- readMSData(f, mode = "onDisk")
  sp <- rw[[i]]
  pdfn <- paste0("./figs/spectrum_", k, ".pdf")
  pdf(pdfn)
  plot(sp, full = TRUE)
  dev.off()
}


pl <- vector("list", 13)

for (k in seq_len(nrow(tj))) {
  f <- paste0("raw/", tj$spectrumFile[k])
  i <- tj$`scan number(s)`[k]
  rw <- readMSData(f, mode = "onDisk")
  sp <- rw[[i]]
  pl[[k]] <- plot(sp, full = TRUE)
}
