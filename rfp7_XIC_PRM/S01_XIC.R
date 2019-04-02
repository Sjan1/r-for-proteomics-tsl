##  install package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("msdata", version = "3.8")

library("MSnbase")
library("magrittr")
library("rtslprot")

## to restart RStudio plot area
graphics.off()
par("mar")
#par(mar=c(1,1,1,1))
par(mar=c(2,2,2,2))

## for the data on the disc
## directory with the data
mzml <- "data/mzml_long_grad_PRM"

## INPUT: filenames without the extension
filenames <- c(#"js190313_P6-pep10-11-12_cid",
               "js190313_P6-pep10-11-12_hcd",
"js190313_P6-pep12_hcd_r2",
#"js190313_P6-pep12_hcd-sid",
#"js190313_P6-pep12_cid",
"js190313_pep11_hcd",
#"js190313_pep11_hcd-sid",
#"js190313_pep11_cid",
#"js190313_pep10_cid",
"js190313_pep10_hcd")#,
#"js190313_pep10_hcd-sid")

#"js190313_P4_hcd" # acetylated V-R peptide
#"js190313_P5_hcd"  # acetylated V-R peptide
#"js190313_P6-pep12_hcd"  # bad run

afile <- NULL
mzml_files <- NULL
res <- NULL
resl <- NULL
length(filenames)
par(mfcol=c(length(filenames),1))
for (i in filenames){
  
  afile <- file.path(mzml, paste(i, "mzML", sep = "."))
  mzml_files <- c(mzml_files,afile)
  mzml_files
  #fl <- mzml_files[i]
  x <- readMSData(afile, mode = "onDisk", msLevel = 2, centroided = FALSE)
  
  ## Define the mz and retention time ranges
  #mz <- 968.9591    #phospho #mz <- 928.9807    #non phospho
  
  ## INPUT: retention times in seconds
  rts <- 6000
  rte <- 7000
  
  ## INPUT: m/z
  mz1 <- 928.9807
  mz2 <- 968.9591
  mz3 <- 636.869
  mz4 <- 776.929
  
  ## INPUT: select m/z to extract chromatogrmas from 
  mz <- c(mz1, mz2, mz4)
  
  ## INPUT; margins
  marg <- 0.01
  mzs <- mz - marg
  mze <- mz + marg
  
  ## XIC input table
  xic_tab <- data.frame(mzs, mze, rts, rte)
  
  resl <- NULL
  for (r in rownames(xic_tab)){
    mzr <- c(xic_tab[r,]$mzs, xic_tab[r,]$mze)
    rtr <- c(xic_tab[r,]$rts, xic_tab[r,]$rte)
    
    ## Filtering the object with raw data
    xic <- x %>% filterRt(rtr) %>% filterMz(mzr)
    res <- rtslprot:::plotxic(xic, rtr, mzr) ## produces the plot
    
    ## INPUT: colors
    col <-  c("red","blue","black")  
    
    ## add a name to the result (source filename) and recalculate TIC to %
    res$pc <- (res$tic-min(res$tic))/(max(res$tic)-min(res$tic))
    res$nm <- i
  
    ## save all XIC results to a list
    res <- list(res)       # results for one plot
    resl <- c(resl,res)    # results for all plots
    
  }
  ## plot it -> not scaled
  #              plot(resl[[1]]$rt,resl[[1]]$tic,type="o", ylim=c(range(resl[[1]]$tic,resl[[2]]$tic,resl[[3]]$tic)), 
  #                col=col[1], main=i, pch = 16, cex = .7)
  #            for (j in 2:length(mz)){
  #              lines(resl[[j]]$rt, resl[[j]]$tic, type="o",col=col[j], main=i, pch = 16, cex = .7)
  #            }
  
  ## plot it - scaled to 100% TIC
  plot(resl[[1]]$rt,resl[[1]]$pc,type="o", ylim=c(range(resl[[1]]$pc,resl[[2]]$pc,resl[[3]]$pc)), 
       col=col[1], main=i, pch = 16, cex = .7)
  for (j in 2:length(mz)){
    lines(resl[[j]]$rt, resl[[j]]$pc, type="o",
          col=col[j], main=i, pch = 16, cex = .7)
  }
  ## legend
  legend("topleft", legend=c(paste0(mz[1],"+/-",marg), paste0(mz[2],"+/-",marg), paste0(mz[3],"+/-",marg)),
         col=c(col[1], col[2], col[3]), lty=1, cex=1)
    
}




#########################################################
stop("NEVER MIND THE ERROR, THE EXECUTION STOPS HERE")
#########################################################


## EXAMPLES

##---------------------------------------------
## print to file
dev.copy(png,filename = "XIC_01.png",
         width = 1440,
         height = 2160,
         units = "px",
         bg = "white")
dev.off()
graphics.off()
##---------------------------------------------
## two plots on top of each other
x <- 1:20
y1 <- x
y2 <- (x)*x

par(mfrow=c(2,1))
for (i in x){
  
}

plot(x,y1,ylim=c(0,25), col="blue")
par(new=TRUE)
plot(x,y2,ylim=c(0,100), col="red",axes=FALSE, xaxt = "n", yaxt = "n",
     ylab = "", xlab = "")
axis(4)


##---------------------------------------------
##multiple series in one plot 
y_vals<-rnorm(100) 
x_vals<-1:100
ii=1
plot(y_vals, col=(257+ii*10)) 
for(ii in 1:3) 
{ 
  y_vals<-rnorm(100) 
  points(x_vals, y_vals, col=(200+ii*10)) 
} 


##---------------------------------------------
## adding lines to a plot
#windows()
par(mfrow=c(2,1))

x<-matrix(rnorm(20000,5,3), nrow=200, ncol=100)
y<-matrix(0, nrow=200, ncol=100)

for (i in 1:200) {
  for (j in 1:100) {
    y[i,j] <- mean(x[i,1:j])
  } 
}

#png(filename="./a1.png")

#here is the ugly bit
#plot(1:100,y[1,1:100],type="l", ylim=range(c(10,0)))
#par(new = TRUE)
#for (j in 2:200) {
#  plot(1:100,y[j,1:100],type="l", ylim=range(c(10,0)), xaxt='n', yaxt='n', ann=FALSE)
#  par(new = TRUE)
#}
#
#graphics.off()

#better
for (p in 1:2) {
    plot(1:100, y[1,1:100], type="l", ylim=range(c(10,0)))
        for (j in 2:200) {
          lines(1:100, y[j,1:100], col=(200+j*10))
        }
}


#graphics.off()


##---------------------------------------------
## from pirate's guide to R
library("yarrr")
par(mfrow = c(2, 2))  # Set up a 2 x 2 plotting space

# Create the loop.vector (all the columns)
loop.vector <- 1:4

for (i in loop.vector) { # Loop over loop.vector
  
  # store data in column.i as x
  x <- examscores[,i]
  
  # Plot histogram of x
  hist(x,
       main = paste("Question", i),
       xlab = "Scores",
       xlim = c(0, 100))
}
##--------------------------------------------- 
##lists

l1 <- NULL
l2 <- NULL
l <- NULL
a <- c("a","b","c")
b <- c(1,2,3)
d <- c("d","e","f")
e <- c(4,5,6)
l1 <- data.frame(a,b)
l2 <- data.frame(d,e)
l1
l2
l <- c(l1,l2)
l <- list(l1,l2)
l


##---------------------------------------------
## list maximum
L1 <- list(a = c(3.4, 5.6, -2.1, -7.8), b = c(2.1, 6.7), c = c(-6.7, 0.001, 8.9))
aa <-  c(4.4, 56.6, -21.1, -7.9)
bb <-  c(1,2,3,4)
aabb <- data.frame(aa,bb)
aabb
L2 <- c(L1,aabb)
L2
m <- lapply(L1, function(x) x[which.max(abs(x))])
class(m)
range(m)

##---------------------------------------------
## legend and main title
##You can use the oma parameter to increase the outer margins, then add the main title with mtext, and try to position the legend by hand.

op <- par(
  oma=c(5,1,5,1),# Room for the title and legend
  mfrow=c(2,2)
)
for(i in 1:4) {
  plot( cumsum(rnorm(100)), type="l", lwd=3,
        col=c("navy","orange")[ 1+i%%2 ], 
        las=1, ylab="Value",
        main=paste("Random data", i))
  
  }
par(op) # Leave the last plot
mtext("Main title", line=2, font=2, cex=1.5)
op <- par(usr=c(0,1,0,0), # Reset the coordinates
          xpd=NA)         # Allow plotting outside the plot region
legend(5,3, # Find suitable coordinates by trial and error
       c("one", "two"), lty=1, lwd=3, col=c("navy", "orange"), box.col=NA)


legend(35,3, legend = c("0", "1"), 
       lty = c(1,1), lwd = c(2,2), col=c("navy", "orange"),
       title = "Subgroup", horiz=TRUE, xpd=NA)

##---------------------------------------------
## multiple scales
samp <- rnorm(50)
d1 <- density(samp)
d2 <- density(samp, adj = 0.5)
d3 <- density(samp, adj = 2.0)
plot( c(d1$x, d2$x, d3$x), c(d1$y, d2$y, d3$y), type = "n")
lines( d1$x, d1$y, lty = 2 )  # maybe lines(d1, lty = 2) would work
lines( d2$x, d2$y, lty = 3 )
lines( d3$x, d3$y, lty = 4 )

##---------------------------------------------
## legend 
# http://www.sthda.com/english/wiki/add-legends-to-plots-in-r-software-the-easiest-way



