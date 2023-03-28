###############
# Functions defined under the current directory
# Last updated:: Time-stamp: <2022-08-29 11:17:06 (qli)>
## REQUIRE!!!###!!!
## procedure to prepare the dataset to plot manhatton plots for genomewide result
###############


source("rFun_plot_genome_result.R")
## plot genomewide result, based on 22 files, each from one chromosomes
plotGenomeResult(result.fileId="qingExample.plot_genome_result.chr", outDir="testDir", fid="testGenomeResultPlot1", my.title="Example of Plot Genomewide Result", my.ylab="log10(p-val)", my.ylim=c(0,7), my.hline=c(5,6), my.col=c("blue", "darkgreen", "black"))


## plot genomewide result, based on one file
plotGenomeResult_oneInput(result.fileName="qingExample.plot_genome_result.all.txt", outDir="testDir", fid="testGenomeResultPlot2", my.title="Example of Plot Genomewide Result", my.ylab="log10(p-val)", my.ylim=c(0,7), my.hline=c(5,6), my.col=c("blue", "darkgreen", "black"))

## plots are saved under the sub-directory testDir/testGenomeResultPlot1.pdf
## testDir/testGenomeResultPlot2.pdf
