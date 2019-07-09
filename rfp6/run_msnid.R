run_msnid <- function(f, fdr = 0.01) {
    if (file.exists(".Rcache"))
        unlink(".Rcache", recursive = TRUE)
    suppressMessages(msnid <- MSnID("."))
    on.exit(unlink(".Rcache", recursive = TRUE))
    suppressMessages(msnid <- read_mzIDs(msnid, f))

    ## add columns
    msnid <- assess_termini(msnid, validCleavagePattern = "[KR]\\.[^P]") ## adds numIrregCleavages
    msnid <- assess_missed_cleavages(msnid, missedCleavagePattern = "[KR](?=[^P$])") ## adds numMissCleavages
    msnid$numCys <- sapply(lapply(strsplit(msnid$peptide,''),'==','C'), sum)
    msnid$PepLength <- nchar(msnid$peptide) - 4
    msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
    msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

    ## Hard filters
    msnid <- apply_filter(msnid, "numIrregCleavages == 0")
    msnid <- apply_filter(msnid, "numMissCleavages <= 2")

    ## Recalibrate
    msnid.chopped <- apply_filter(msnid, "abs(mass_measurement_error(msnid)) < 20")
    msnid <- recalibrate(msnid.chopped)

    ## set initial filter
    filtObj <- MSnIDFilter(msnid)
    filtObj$absParentMassErrorPPM <- list(comparison = "<", threshold = 10.0)
    filtObj$msmsScore <- list(comparison = ">", threshold = 10.0)


    ## optimise filter
    filtObj <- optimize_filter(filtObj, msnid,
                               fdr.max = fdr, level = "accession",
                               method = "Nelder-Mead", n.iter = 500)
    print(evaluate_filter(msnid, filtObj))

    ## apply filter
    msnid <- apply_filter(msnid, filtObj)
    ## remove contaminants
    msnid <- apply_filter(msnid, "!grepl('cont', accession)")
    ## remove decoys
    apply_filter(msnid, "!grepl('XXX', accession)")
}
