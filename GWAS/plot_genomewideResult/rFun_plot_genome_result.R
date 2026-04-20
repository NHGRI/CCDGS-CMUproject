###############
# Functions defined under the current directory
# Last updated:: Time-stamp: <2022-08-29 11:33:16 (qli)>
## REQUIRE!!!###!!!
# Function to plot the result for whole genome, including 22 chromosomes
# plotGenomeResult(result.fileId, outDir, fid, plottype="pdf", my.ylim=NULL, my.hline=c(5), my.col=c("blue", "darkgreen"))
###############


## the file of the result file will look like following
## col1: chromosome (integer); col2: position; col3: results
##$$ 2   1251111        1.76032   
##$$ 2   1252790        1.76056   
##$$ 2   1253590        1.77219   
##$$ 2   1254487        1.77309  

##################################
##$$parameters for the function
##################################
##$$ result.fileId: uniq identifier of the genomewide result files
##$$ 
##$$ outDir:   output directory for plot
##$$ plottype: pdf, ps, tiff
##$$ fid:      uniq identifier of the output plot file
##$$ 
##$$ my.ylim:  range of y axis
##$$ my.hline: horizontal line for threshold
##$$ my.col:   an array of color choices to be used repeatedly to plot each chromosome

plotGenomeResult = function(result.fileId, withHeader=T, outDir, fid, plottype="pdf", my.title="", my.ylab="log10(p-val)", my.ylim=NULL, my.hline=c(5), my.col=c("blue", "darkgreen")){
    
## fid ="plotGWAS.genome"
## result.fileId = ""
## plottype="pdf"
## 
## outfile=""
## my.ylim = c(0, 5)
## my.hline = c(5,6)
## my.col = c("blue", "darkgreen")

## plot different types
if (plottype=="pdf"){
  outfile = file.path(outDir, paste(fid, ".pdf", sep=""))
  pdf(outfile)
} else if (plottype=="ps") {
  outfile = file.path(outDir, paste(fid, ".ps", sep=""))
  ps(outfile)
} else if (plottype=="tiff") {
  outfile = file.path(outDir, paste(fid, ".tiff", sep=""))
  tiff(filename = outfile, width = 10, height = 8, units = "in", pointsize = 16,
       compression = "lzw", bg = "white",res=500)
}


par(las=2, cex.axis=.8)

x.interGap = 10
x.gap = NULL
## if the files are in 22 files, identified by the same file id
for (ch.num in 1:22){
    fn = paste(result.fileId, ch.num, ".txt", sep="")
     
    ch.dat = read.csv(fn, sep="", header=withHeader, as.is=T)
    print(str(ch.dat))

    x.gap = c(x.gap, ch.dat[ nrow(ch.dat), 2 ]+x.interGap)

    my.ylim.t = min(ch.dat[,3][is.finite(ch.dat[,3])])
    my.ylim.t2 = max(ch.dat[,3][is.finite(ch.dat[,3])])
    my.ylim = c( min(my.ylim.t, my.ylim[1]), max(my.ylim.t2, my.ylim[2]))
}

x.gap.end = cumsum(x.gap)
x.gap.start = x.gap.end
x.gap.start = c( 0, x.gap.start[-22] )
#print(cbind(x.gap.start, x.gap.end))


leg.col = rep(my.col, times=ceiling(22/length(my.col)))
## plot the range of the graph
plot(0, 0, lty="dashed", xlim=c(0, x.gap.end[22]*1.02), ylim=my.ylim, main=my.title, ylab=my.ylab, xlab="", type="n", axes=F)
abline(h=my.hline, lty="dashed", col="grey")
abline(v=x.gap.start, lty="solid", col="grey")  ## plot the divider bw chromosomes

## plot actual scatterplot
for (ch.num in 1:22){
    fn = paste(result.fileId, ch.num, ".txt", sep="")
    ch.dat = read.csv(fn, sep="", header=withHeader, as.is=T)
    points(x.gap.start[ch.num]+ch.dat[,2], ch.dat[,3], pch=16, col=leg.col[ch.num])
}

axis(2)
axis(1, at=(x.gap.start+x.gap.end)/2, labels=paste(1:22, "", sep=""))

box()
dev.off()

}


plotGenomeResult_oneInput= function(result.fileName, withHeader=T, outDir, fid, plottype="pdf", my.title="", my.ylab="log10(p-val)", my.ylim=NULL, my.hline=c(5), my.col=c("blue", "darkgreen")){
    
## plot different types
if (plottype=="pdf"){
  outfile = file.path(outDir, paste(fid, ".pdf", sep=""))
  pdf(outfile)
} else if (plottype=="ps") {
  outfile = file.path(outDir, paste(fid, ".ps", sep=""))
  ps(outfile)
} else if (plottype=="tiff") {
  outfile = file.path(outDir, paste(fid, ".tiff", sep=""))
  tiff(filename = outfile, width = 10, height = 8, units = "in", pointsize = 16,
       compression = "lzw", bg = "white",res=500)
}


par(las=2, cex.axis=.8)

x.interGap = 10
x.gap = NULL
## if the files are in 22 files, identified by the same file id
     
dd = read.csv(result.fileName, sep="", header=withHeader, as.is=T)
print(str(dd))
    
for (ch.num in 1:22){
    ch.dat = dd[ dd[,1]==ch.num, ]
    #print(str(ch.dat))

    x.gap = c(x.gap, ch.dat[ nrow(ch.dat), 2 ]+x.interGap)

    my.ylim.t = min(ch.dat[,3][is.finite(ch.dat[,3])])
    my.ylim.t2 = max(ch.dat[,3][is.finite(ch.dat[,3])])
    my.ylim = c( min(my.ylim.t, my.ylim[1]), max(my.ylim.t2, my.ylim[2]))
}

x.gap.end = cumsum(x.gap)
x.gap.start = x.gap.end
x.gap.start = c( 0, x.gap.start[-22] )
#print(cbind(x.gap.start, x.gap.end))


leg.col = rep(my.col, times=ceiling(22/length(my.col)))
## plot the range of the graph
plot(0, 0, lty="dashed", xlim=c(0, x.gap.end[22]*1.02), ylim=my.ylim, main=my.title, ylab=my.ylab, xlab="", type="n", axes=F)
abline(h=my.hline, lty="dashed", col="grey")
abline(v=x.gap.start, lty="solid", col="grey")  ## plot the divider bw chromosomes

## plot actual scatterplot
for (ch.num in 1:22){
    ch.dat = dd[ dd[,1]==ch.num, ]
    points(x.gap.start[ch.num]+ch.dat[,2], ch.dat[,3], pch=16, col=leg.col[ch.num])
}

axis(2)
axis(1, at=(x.gap.start+x.gap.end)/2, labels=paste(1:22, "", sep=""))

box()
dev.off()

}
## EOF
