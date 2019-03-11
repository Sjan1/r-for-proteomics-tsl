## rtslprot functions
## meant for local use and fast testing without the need of changing the package 
make_psm_MSnSet <- function(mzid_files) {
  msnid <- MSnID()
  msnid <- read_mzIDs(msnid, mzid_files)
  ## TODO - filter MSnID and remove non-tryptic peptides
  td <- as(msnid, "data.table")
  td$count <- 1
  readMSnSet2(td, ecol = ncol(td))
}

make_pep_MSnSet <- function(mzid_files, fcol = "pepSeq") {
  e <- make_psm_MSnSet(mzid_files)
  combineFeatures(e, fcol = fcol, fun = sum)
}
