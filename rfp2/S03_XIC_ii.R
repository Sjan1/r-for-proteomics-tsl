##  install package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msdata", version = "3.8")

library("MSnbase")
#library("msdata")
library("magrittr")
library("rtslprot")

## to restart RStudio plot area
graphics.off()
par("mar")
#par(mar=c(1,1,1,1))
par(mar=c(2,2,2,2))

## for the data on the disc
## directory with the data
mzml <- "data/mzml"

## filenames without the extension
filenames <- c("js190313_P6-pep10-11-12_cid",
"js190313_P6-pep10-11-12_hcd",
"js190313_P6-pep12_hcd_r2",
"js190313_P6-pep12_hcd-sid",
"js190313_P6-pep12_cid",
"js190313_pep11_hcd",
"js190313_pep11_hcd-sid",
"js190313_pep11_cid",
"js190313_pep10_cid",
"js190313_pep10_hcd",
"js190313_pep10_hcd-sid")

#"js190313_P4_hcd" # acetylated V-R peptide
#"js190313_P5_hcd"  # acetylated V-R peptide
#"js190313_P6-pep12_hcd"  # bad run

afile <- NULL
mzml_files <- NULL

par(mfcol=c(6,2))

for (i in filenames){
    afile <- file.path(mzml, paste(i, "mzML", sep = "."))
    
    mzml_files <- c(mzml_files,afile)
    mzml_files
    #fl <- mzml_files[i]
    x <- readMSData(afile, mode = "onDisk", msLevel = 2, centroided = FALSE)

## Define the mz and retention time ranges
#mz <- 968.9591    #phospho
#mz <- 928.9807    #non phospho
    rtr <- c(6100, 7000)
    mz1 <- 928.9807
    mz2 <- 968.9591
    #mz3 <- 636.869
    mz4 <- 776.9298
    
    mzr1 <- c(mz1 - 0.01, mz1 + 0.01)
    mzr2 <- c(mz2 - 0.01, mz2 + 0.01)
    #mzr3 <- c(mz3 - 0.01, mz3 + 0.01)
    mzr4 <- c(mz4 - 0.001, mz4 + 0.001)
## Filtering the object
    xic1 <- x %>% filterRt(rtr) %>% filterMz(mzr1)
    xic2 <- x %>% filterRt(rtr) %>% filterMz(mzr2)
    #xic3 <- x %>% filterRt(rtr) %>% filterMz(mzr3)
    xic4 <- x %>% filterRt(rtr) %>% filterMz(mzr4)
  
    
    res1 <- rtslprot:::plotxic(xic1, rtr, mzr1) ## produces the plot
    res2 <- rtslprot:::plotxic(xic2, rtr, mzr2) ## produces the plot
    #res3 <- rtslprot:::plotxic(xic3, rtr, mzr3) ## produces the plot
    res4 <- rtslprot:::plotxic(xic4, rtr, mzr4) ## produces the plot
    
        
        
        ## plot MS2 XIC
        #par(mfrow=c(2,2))
        #dev.off()
        #window()
        ## superimposed plots
plot(res1$rt,res1$tic,type="o", ylim=c(range(res1$tic)), col="red", main=i,pch = 16, cex = .7)
par(new=TRUE)
plot(res2$rt,res2$tic,type="o", ylim=c(range(res2$tic)), col="blue",axes=FALSE, pch = 16, cex = .7)
par(new=TRUE)
#plot(res3$rt,res3$tic,type="l", ylim=c(range(res3$tic)), col="green",axes=FALSE, pch = 16, cex = .7)
#par(new=TRUE)
plot(res4$rt,res4$tic,type="o", ylim=c(range(res4$tic)), col="black",axes=FALSE, pch = 16, cex = .7)
axis(4)

}

#########################################################
stop("NEVER MIND THE ERROR, THIS EXECUTION STOPS HERE")
#########################################################

## print to file
dev.copy(png,filename = "XIC_01.png",
         width = 1440,
         height = 2160,
         units = "px",
         bg = "white")
dev.off()

plot(mpg ~ hp, data = mtcars, pch = 16, cex = .9)

mzr <- c(928.9807-0.01, 928.9807+0.01) ## 968.9591 +/- 0.01
rtr <- c(6300, 7000)
#window()
res_plot2 <- rtslprot:::plotxic(x, rtr, mzr) ## produces the plot
res ## contains the data used to plot
dev.off()

## superimposed plots
plot(res_plot1$rt,res_plot1$tic,type="l", ylim=c(range(res_plot1$tic)), col="red")
par(new=TRUE)
#lines(res_plot2$rt,res_plot2$tic,col="blue", axes=F)#ylim = c(range(res_plot2$tic))) 
plot(res_plot2$rt,res_plot2$tic,type="l", ylim=c(range(res_plot2$tic)), col="blue",axes=FALSE)
axis(4)

### END ###
##example

x <- 1:20
y1 <- x
y2 <- (x)*x

par(mfrow=c(2,1))
for (i in x){
  
}

plot(x,y1,ylim=c(0,25), col="blue")
par(new=TRUE)
plot(x,y2,ylim=c(0,100), col="red",axes=FALSE)
axis(4)




