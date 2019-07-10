##  install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msdata", version = "3.8")

library("MSnbase")
library("magrittr")
library("rtslprot")

##############
## INPUT data: 
## -1 directory with the data
## -2 filenames without the extension

# ##PRM - MS2 only data
# mzml <- "data/mzml_long_grad_PRM"
# filenames <- c(##"js190313_P6-pep10-11-12_cid",
#   "js190313_P6-pep10-11-12_hcd",
#   "js190313_P6-pep12_hcd_r2",
#   "js190313_pep11_hcd")
#   #"js190313_P6-pep12_hcd-sid",
#     #"js190313_pep11_hcd-sid",
#   #"js190313_pep10_hcd",
#  #"js190313_pep10_hcd-sid")

## MS1-MS2 data
mzml <- "data/190430_LongGrad_MS1MS2"
filenames <- c("js190501_P6-pep10-11-12chs",
               "js190501_P6-pep10-11-12ch")

############################
## INPUT RT: retention times in seconds
rts <- 100*60
rte <- 117*60
####################
## INPUT m/z: masses of which to extract chromatograms 
mz <- c(928.9807, 968.9591, 776.9290)
################
## INPUT: colors
col <-  c("red","blue","black")
#################
## INPUT; margins
marg <- 0.01
mzs <- mz - marg
mze <- mz + marg

## XIC input table
xic_tab <- data.frame(mzs, mze, rts, rte, col)
xic_tab

## to restart RStudio plot area
graphics.off()
par("mar")
#par(mar=c(1,1,1,1))
par(mar=c(2,2,2,2))

length(filenames)
par(mfcol=c(length(filenames),1))
#par(mfcol=c(1,1))

#######
## main 
for (i in filenames){
  afile <- file.path(mzml, paste(i, "mzML", sep = "."))
  mzml_files <- c(mzml_files,afile)
  mzml_files
  x <- readMSData(afile, mode = "onDisk", msLevel = 2, centroided = FALSE)

  afile <- NULL
  mzml_files <- NULL
  res <- NULL
  resl <- NULL
    ## calculate XICs
    for (r in rownames(xic_tab)){
    mzr <- c(xic_tab[r,]$mzs, xic_tab[r,]$mze)
    rtr <- c(xic_tab[r,]$rts, xic_tab[r,]$rte)
    
    ## Filtering the object with raw data
    xic <- x %>% filterRt(rtr) %>% filterMz(mzr)
    res <- rtslprot:::plotxic(xic, rtr, mzr, plot=FALSE) ## produces the plot
    
    ## add a name to the result (source filename) and recalculate TIC to %
    res$pc <- (res$tic-min(res$tic))/(max(res$tic)-min(res$tic))
    res$nm <- i
    
    ## save all XIC results to a list
    res <<- list(res)       # results for one plot
    resl <<- c(resl,res)    # results for all plots
    }

  ###############################  
  ## plot it - scaled to 100% TIC
  # test with colors and legend
  plot(resl[[1]]$rt, resl[[1]]$pc, type="l", 
       ylim=c(range(resl[[1]]$pc, resl[[2]]$pc, resl[[3]]$pc)),
       xlim=c(xic_tab$rts[1],xic_tab$rte[1]),
       col=col[1], main=i, pch = 16, cex = .7)
  for (j in 2:length(mz)){
    myplot <-  lines(resl[[j]]$rt, resl[[j]]$pc, type="l",
                     col=col[j], main=i, pch = 16, cex = .7)}
  ## legend
  legend("topleft", 
         legend=c(paste0(mz,"+/-",marg)), col=col, lty=1, cex=1)
    
}

#############
## THE END ##
#############




