## some tests of and eprot (Jan)
eprot <- readRDS("eprot_msgf.Rds")
eprot <- readRDS("eprot_mascot.Rds")
e <- readRDS("e_msgf.Rds")
e <- readRDS("e_mascot.Rds")
# all feature data columns
cn <- colnames(fData(eprot))
cn
# select just peptide sequences columns
cn_pep <- cn[grepl("^pepSeq.",cn)]
cn_pep
head(fData(eprot)[cn_pep])

# select one protein row - logical
pr <- rownames(fData(eprot))=="ART4-eGFP"
head(pr)
table(pr)
# features of this protein

# only the sequences
d <- fData(eprot)[pr,cn_pep]
View(d)

# But this does not make much senseto do with the proteins
# do it all again with the peptides
# input above
cn <- colnames(fData(e))
cn
# select just peptide sequences columns
cn_pep <- cn[grepl("^pepSeq.",cn)]
cn_pep
head(fData(e)[cn_pep])


# select all accesion numbers - logical
ac <- grepl("^accession.",colnames(fData(e)))
head(ac)
table(ac)
acc <- fData(e)[,ac]
head(acc)
dim(acc)
# select only the samples that have a protein from accessions
pr <- grepl("35a12",acc)
pr
table(pr)
#apr <- fData(eprot)[acc==pr]
apr <- acc[,pr]
head(apr)
dim(apr)
apr[apr==pr]
table(apr)
#
which()


# features of this protein

# only sequences
d <- fData(eprot)[pr,cn_pep]
View(d)     
