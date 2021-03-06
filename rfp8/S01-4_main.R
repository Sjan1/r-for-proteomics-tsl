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
    par(mar = c(10, 7, 1, 1))
    graphics::image(t(x),
                    xaxt = "n",
                    yaxt = "n",...)
    xticks <- seq(0, 1, length.out = ncol(x))
    axis(1, xticks,
         labels = colnames(x),
         cex.axis = x.cex.axis,
         las=2)
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
pgroup <- cutree(cl, h = 0.99)
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
  writeLines(c(".............",paste("|||||",i,"|||||")))
      plot_proteins_in_group(prots_by_peps, pgroup, i)
    print(get_proteins_in_group(prots_by_peps, pgroup, i))
    #scan(n = 1)
    title(main = paste("pgroup: ", i),
                    cex.main = 1,
                    font.main= 2,
                    col.main= "navyblue")
scan(n=1)
}


## looking for one protein in particular
k <- which(sapply(pgroups, function(x) "AT3G01500.1" %in% x))
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




#######################
## The End  Dec 2019 ##
######################

## SAVE WHAT WE MIGHT NEED
saveRDS(res,"res.rds")
saveRDS(msnl,"msnl.rds")
saveRDS(pgroup,"pgroup.rds")
saveRDS(pgroups,"pgroups.rds")
saveRDS(prots_by_peps,"prots_by_peps.rds")
## open when needed
res <- readRDS("res.Rds")
msnl <- readRDS("msnl.rds")
pgroup <- readRDS("pgroup.Rds")
pgroups <- readRDS("pgroups.Rds")
prots_by_peps <- readRDS("prots_by_peps.Rds")


## PEPTIDE SCATTER PLOTS WITHIN PGROUPS
tmp <- list()
GroupSamplePep <- list()
for (g in 1:length(pgroups)){
  #for (p in 1:length(all_peps)){
  for (m in 1:length(msnl)){
    .pg <- pgroups[[g]]
    #.peps <- all_peps[p]
    .mset <- msnl[[m]]
  
    #print(c(g,"-",m))
    tmp[[m]] <- fData(.mset)$pepSeq[fData(.mset)$accession %in% .pg]
  }
  names(tmp) <- names(index)  
  GroupSamplePep[[g]] <- tmp
}  
names(GroupSamplePep) <- names(pgroups)

## Plot scatter plots for each pgroup
for (g in 1:length(pgroups)){
  GroupSamplePep[[g]][1:4]
  x <- unlist(GroupSamplePep[[g]][1:4])
  y <- unlist(GroupSamplePep[[g]][5:7])
  x <- as.data.frame(table(x))
  y <- as.data.frame(table(y))
  if (isEmpty(x[1])) {x <- data.frame(PepSeq=factor(),Freq.PH=integer())}
  if (isEmpty(y[1])) {y <- data.frame(PepSeq=factor(),Freq.CH=integer())}
  colnames(x) <- c("pepSeq","Freq.PH")
  colnames(y) <- c("pepSeq","Freq.CH")
  df.pept <- merge(x,y, by="pepSeq",all=TRUE)
  df.pept[is.na(df.pept)] <- 0
 
  plot(x = df.pept$Freq.PH, y = df.pept$Freq.CH, pch = 16, cex = 1, 
       main = "Peptide scatterplot for selected samples", 
       sub = paste("pgroup= ",g),
       col.main="black", col.lab="blue", col.sub="red",
       col = ifelse(df.pept$Freq.PH == 0 | df.pept$Freq.CH == 0, "red", "blue"))
      abline(a=0, b=1, lty = 5)
  text(x = df.pept$Freq.PH,
       y = df.pept$Freq.CH,
       labels=df.pept$pepSeq, cex=0.65, pos=3, col="brown")
  
  writeLines(c(".............",paste("|||||",g,"|||||")))
  scan(n=1)

  #require(ggplot2
  #ggplot(df.pept, aes(x,y)) + geom_point() + geom_text(aes(label=pepSeq))
  #scan()  
  }


        