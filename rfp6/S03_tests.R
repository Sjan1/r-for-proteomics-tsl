## some tests of and eprot (Jan)
# all feature data columns
cn <- colnames(fData(eprot))
cn
# select just peptide sequences columns
cn_pep <- cn[grepl("^peptide.",cn)]
cn_pep
head(fData(e)[cn_pep])

# select one protein row - logical
pr <- rownames(fData(eprot))=="35a12"
head(pr)
# features of this protein

# only sequences
d <- fData(eprot)[pr,cn_pep]
View(d)

# it does not make much sense with proteins do it all again with peptides

cn <- colnames(fData(e))
cn
# select just peptide sequences columns
cn_pep <- cn[grepl("^peptide.",cn)]
cn_pep
head(fData(e)[cn_pep])


# select all accesion numbers - logical
acc <- grepl("^accession.",colnames(fData(e)))
head(acc)
acc <- fData(e)[,acc]
head(acc)

# select only one protein from accessions
pr <- grepl("35a12",acc[,2])
pr
table(pr)
apr <- fData(eprot)[,pr]
apr <- acc[,2]
apr[apr==pr]
table(apr)

apr[apr=="35a12"]

head(apr)
View(apr)
# features of this protein

# only sequences
d <- fData(eprot)[pr,cn_pep]
View(d)     
