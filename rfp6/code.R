#' Read the sample experiment table
#' 
#' This function read the sample experiment table 
#' and checks that all the mzML and mzid files exist.
#'
#' @param f `character(1)` defining the name of the 
#' sample experiment table in comma separated values 
#' (csv) format. Default is `"SampleExperimentTable.csv`.
#' @param mzml `character(1)` defining the directory name
#' where the mzML files can be found. Default is 
#' `"mzml"`.
#' @param mzid `character(1)` defining the directory name
#' where the mzid files can be found. Default is 
#' `"msgf"`.
#'
#' @return `data.frame` with the sample experiment table.
#' @export 
#' @md
readSampleExperimentTable <- function(f = "SampleExperimentTable.csv",
                                      mzml = "mzml", 
                                      mzid = "msgf") {
    exp <- read.csv(f, stringsAsFactors = FALSE)
    ## check that mzml files exist
    mzml_files <- file.path(mzml, paste(exp$name, "mzML", sep = "."))
    mzml_exist <- file.exists(mzml_files)
    if (!all(mzml_exist))
      stop(sum(!mzml_exist), "mzML files don't exist.")
    ## check that mzid files exist
    mzid_files <- file.path(mzid, paste(exp$name, "mzid", sep = "."))
    mzid_exist <- file.exists(mzid_files)
    if (!all(mzid_exist))
      stop(sum(!mzid_exist), "mzid files don't exist.")
    exp
}
