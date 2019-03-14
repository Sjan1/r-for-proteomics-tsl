i <- grep("accession\\.", fvarLabels(e))
#View(fData(e)[,i])
k <- apply(fData(e)[, i], 1, 
      function(x) unique(na.omit(as.character(x))))
#View(k)
fData(e)$nprots <- lengths(k)
fData(e)$accession <- sapply(k, paste, collapse = ";")


