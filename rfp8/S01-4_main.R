get_proteins_in_group <- function(x, group, i) {
    x2 <- x[group == i, , drop = FALSE]
    x2 <- x2[, colSums(x2, na.rm = TRUE) > 0, drop = FALSE]
    n_zeros <- apply(x2, 1, function(xx) sum(xx == 0, na.rm = TRUE))
    m_zeros <- apply(x2, 2, function(xx) sum(xx == 0, na.rm = TRUE))
    x2[order(n_zeros), order(m_zeros), drop = FALSE]
}

plot_proteins_in_group <- function(x, group, i,
                                   yticks = 10,
                                   x.cex.axis = 0.75,
                                   y.cex.axis = 0.75,
                                   ...) {
    x <- get_proteins_in_group(x, group, i)
    x[x == 0] <- NA
    nc <- ncol(x)
    nr <- nrow(x)
    par(mar = c(4, 7, 1, 1))
    graphics::image(t(x),
                    xaxt = "n",
                    yaxt = "n", ...)
    xticks <- seq(0, 1, length.out = ncol(x))
    axis(1, xticks,
         labels = colnames(x),
         cex.axis = x.cex.axis)
    yticks <- seq(0, 1, length.out = nrow(x))
    axis(2, yticks,
         labels = rownames(x),
         cex.axis = y.cex.axis,
         las = 2)
    invisible(NULL)
}


library("MSnbase")
library("ggplot2")
library("dplyr")
library("magrittr")
library("readr")
library("MSnID")
## BiocManager::install("lgatto/rtslprot")
library("rtslprot")


source("../rfp6/run_msnid.R")


####################################################################
## 1 - read THE EXPERIMENT TABLE describing the experimental design:
## (i) samples
## (ii) measurements
## (iii) searches
## (iv) experimental factors

## mzid <- "mascot_fdr1pc"
## needs simlink to ../rfp6/msgf
mzid <- "msgf"
## mzid <- "msgf2"
exp <- readSampleExperimentTable("../rfp6/SampleExperimentTable.csv",
                                 mzid = mzid)
################################
## 2 - APPLY FILTER if necessary
exp <- exp %>%
  filter(cleavage == "tryp") %>%
  select(-category) %>%
  dropConstantVariables()

############################
## 3 - DEFINE UNIQUE SAMPLES
## lower in hierarchy are only fractions (of samples) that will be combined
## all this have to be described in the experiment table
fcols <- c("phenotype", "treatment", "biorep")
etab <- experimentHierarchy(exp, fcols)
etab

#################################
## 4 - READ PATHS TO SEARCH DATA
## This chunk reads each row of the experimental design
## and creates a list of paths to mzid files.
## useful for re-process only subset of samples
mzid_files <- apply(etab, 1, function(.etab) {
  filenames <- exp %>%
    filter(biorep == .etab[["biorep"]],
           phenotype == .etab[["phenotype"]],
           treatment == .etab[["treatment"]]) %>%
    select(name)
  mzid_files <- file.path(mzid,paste(filenames[[1]], "mzid",
                                     sep = "."))
  
## make a list of files
  return(mzid_files)
})
mzid_files
## add names
names(mzid_files) <- apply(etab, 1,
                           function(.etab) paste(.etab, collapse = "_"))
#rename also etab
rownames(etab) <- names(mzid_files)
etab

## let us consider only subset od data
## thus from 'exp' and 'etab' we get
mzid_files
(mzid_files <- mzid_files[c("P_H_2","P_H_3","P_H_4","P_H_5","C_H_2","C_H_3","C_H_4")])


######################################################
## 5 - continuing from the list of paths (mzid_files)
## advantage is that if we want, we can subset the list
index <- names(mzid_files)
names(index) <- names(mzid_files)
index

msnl <- lapply(index, function(.x){
    .files <- mzid_files[[.x]]
    msnid <- run_msnid(.files)
    psm <- rtslprot::as_MSnSet(msnid)
    ## create protein by peptide matrix
    prot_by_pep <- matrix(NA,
           nrow = length(unique(fData(psm)$accession)),
           ncol = length(unique(fData(psm)$pepSeq)))
    rownames(prot_by_pep) <- sort(unique(fData(psm)$accession))
    colnames(prot_by_pep) <- unique(fData(psm)$pepSeq)
    for (.prot in rownames(prot_by_pep)) {
        for (.pep in colnames(prot_by_pep)) {
            prot_by_pep[.prot, .pep] <-
                sum(fData(psm)$accession == .prot & fData(psm)$pepSeq == .pep)
        }
    }
    psm@experimentData@other$prot_by_pep <- prot_by_pep
    ## calculate protein groups based on Jaccard distance
    d <- proxy::dist(prot_by_pep, method = "Jaccard")
    cl <- hclust(d)
    ## group as soon as one shared peptide
    pgroup <- cutree(cl, h = 0.25)
    names(pgroup) <- names(d)
    ## create protein group identifers
    fData(psm)$pgroup_id <- pgroup[fData(psm)$accession]
    fData(psm)$pgroup <- NA
    for (i in pgroup) {
        pg_i <- unique(fData(psm)[fData(psm)$pgroup_id == i, "accession"])
        pg_i <- paste(pg_i, collapse = ";")
        fData(psm)[fData(psm)$pgroup_id == i, "pgroup"] <- pg_i
    }
    psm
})


## Define protein groups accross all samples
all_prots <- unique(unlist(lapply(msnl, function(x) fData(x)$accession)))
all_peps <- unique(unlist(lapply(msnl, function(x) fData(x)$pepSeq)))

## Create global proteins by peptides matrix
prots_by_peps <- matrix(0,
                        nrow = length(all_prots),
                        ncol = length(all_peps))
rownames(prots_by_peps) <- all_prots
colnames(prots_by_peps) <- all_peps

pb <- progress::progress_bar$new(total = 7)
for (x in msnl) {
    pb$tick()
    df <- expand.grid(unique(fData(x)$accession),
                      unique(fData(x)$pepSeq))
    for (i in 1:nrow(df)) {
        .prot <- as.character(df[i, 1])
        .pep <- as.character(df[i, 2])
            prots_by_peps[.prot, .pep] <-
            sum(fData(x)$accession == .prot & fData(x)$pepSeq == .pep)/length(.prot)
    }
}

## clustering and cutting at 0.25
d <- proxy::dist(prots_by_peps, method = "Jaccard")
cl <- hclust(d)
plot(cl)
## group as soon as one shared peptide
pgroup <- cutree(cl, h = 0.5)
abline(h = 0.5, col = "red")
names(pgroup) <- names(d)

## define protein groups
pgroups <- tapply(names(pgroup), pgroup, c)

res <- matrix(0, ncol = length(msnl),
              nrow = length(pgroups))
colnames(res) <- names(msnl)
rownames(res) <- names(pgroups)

for (i in 1:length(pgroups)) {
    for (j in 1:length(msnl)) {
        .pg <- pgroups[[i]]
        .mset <- msnl[[j]]
        res[i, j] <- sum(exprs(.mset)[fData(.mset)$accession %in% .pg ,])/length(.pg)
    }
}
fd <- data.frame("pgroups"=sapply(pgroups, paste, collapse = ";"))
res <- MSnSet(exprs = res, fData = fd, pData = etab[index, ])
#fData(res)$pgroups <- sapply(pgroups, paste, collapse = ";")
#pData(res) <- etab[index, ]

for (i in 1:max(pgroup)) {
    plot_proteins_in_group(prots_by_peps, pgroup, i)
    print(get_proteins_in_group(prots_by_peps, pgroup, i))
    scan(n = 1)
}

## looking for one protein in particular
k <- which(sapply(pgroups, function(x) "AT5G43980.1" %in% x))
exprs(res[k, ])
get_proteins_in_group(prots_by_peps, pgroup, k)
plot_proteins_in_group(prots_by_peps, pgroup, k)

fData(res)$lfc <- log2(rowMeans(exprs(res)[, 1:4]+1) / rowMeans(exprs(res)[, 5:7]+1))
fData(res)$meanInt <- rowMeans(exprs(res))
with(fData(res), plot(meanInt, lfc))
abline(h = 0)

ms2df(res) %>%
    arrange(desc(lfc)) %>%
    DT::datatable()




#############
## The End ##
#############

## SAVE IT
saveRDS(res,"res.rds")
saveRDS(msnl,"msnl.rds")
saveRDS(pgroups,"pgroups.rds")
## open when needed
res <- readRDS("res.Rds")
