##-----------------------------------------------------
## from Laurent to replace accvec 
## get accession from different samples in feature data
i <- grep("accession\\.", fvarLabels(e))
#View(fData(e)[,i])
k <- apply(fData(e)[, i], 1, 
      function(x) unique(na.omit(as.character(x))))
#View(k)
fData(e)$nprots <- lengths(k)
fData(e)$accession <- sapply(k, paste, collapse = ";")

##------------------------------------------------------
## modified for peptides
## combine peptide sequences in one column
i <- grep("pepSeq\\.", fvarLabels(e))
k <- apply(fData(e)[, i], 1, 
           function(x) unique(na.omit(as.character(x))))
fData(e)$npepseq <- lengths(k)
fData(e)$peptide_sequence <- sapply(k, paste, collapse = ";")


