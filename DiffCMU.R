####  Set Directory ######
# to path/to/CMU 

#### Set Library and Packages ######
library(psych)   # stats  test
library(reticulate)  # to use python
library(karyoploteR)  # for plots
library(matrixStats)
library(GenomicRanges)
library(ggplot2)
library(gg.gap)
library(ggsignif)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(annotate)
library(gg.gap)

#setwd(paste0('~/Downloads//CodeforCU/'))
setwd("~/Library/CloudStorage/Box-Box/Methylcorrelation_MS/CodePreperationForReleasing/CodeForCMUFreezeVer01/")

source_python('PythonIndexing/ReadCpGfromBigFileVer01.py')
source('Residuals.R')
source('pvalForDiffCMU.R')
source('getClusDetails.R')
source('grangeForCMU.R')
source('manhattanCaseControl.R')


qqnormx <- function(x) {qqnorm(x,plot.it = FALSE)$x}
corForCluster <- function(XP)
{
  XP[abs(XP)<cutoffCor] <- 0
  XP
}


###### Parameters  #####
pvalCthresh <- .05
pvalWthresh <-  .05
NCpGThresh <- 9 # > ie. set one less 


####  Get P values for comparison  ####
###### Cont and Non Cont seperately #####
# set the following variables 
# data_file methylation tsv file prepared  in Step 1 
# index_file corresponding index file for  methylation data prepared in Step 2
# idx_ctrl_X      index for controls/normal 
# idx_case_X      index for case/Disease subjects  
# covX covariates to be used, in same order as it is in data_file. idx_ctrl_X and idx_case_X would be used to subset 
#                this matrix for diseaes and controls. 
# CMU_file_for_ctrl path_to_CMU_ctrl_Output/ContM0.6.txt
# CMU_file_for_case path_to_CMU_case_Output/ContM0.6.txt
# pvalFilenameForNormalCMU Provide filename with path and extension Rdata
# pvalFilenameForDiseaseCMU Provide filename with path and extebsuib Rdata

        pvalForDiffCor(filename1 = data_file , IndexFile1 = index_file,
                       indexControl = idx_ctrl_X,
                       indexCase = idx_case_X,
                       X = covX,
                       ResultFileControl = CMU_file_for_ctrl,
                       ResultFileCase = CMU_file_for_case,
                       pvalFileName = pvalFilenameForNormalCMU)
        pvalForDiffCor(filename1 = data_file , IndexFile1 = index_file,
               indexControl = idx_case_X,
               indexCase = idx_ctrl_X,
               X = covX,
               ResultFileControl = CMU_file_for_case,
               ResultFileCase = CMU_file_for_ctrl,
               pvalFileName = pvalFilenameForDiseaseCMU)

      
        
        
        ds1 <- grangeForCMU(pvalFilenameForNormalCMU,
                            CMU_file_for_ctrl,
                            pvalCthresh,pvalWthresh,NCpGThresh)
        ds1$pval <- ds1$pvalWNorm
        ds1$grtVal <- ds1$grtValNorm
        pCut <- rbind(pCut,c(paste0(itype,iloc,icor,"N"),max(ds1$pval[ds1$grtVal])))
        
        #                  #
        # Check Get P values section to generate p-values 
        #                  #
        
        ds2 <- grangeForCMU(pvalFilenameForDiseaseCMU,
                            CMU_file_for_case,
                            pvalCthresh,pvalWthresh,NCpGThresh)
        ds2$pval <- ds2$pvalWNorm
        ds2$grtVal <- ds2$grtValNorm
        if(icor == "0.6") pCut <- rbind(pCut,c(paste0(itype,iloc,icor,"E"),max(ds1$pval[ds1$grtVal])))
        
        figFileName = 'pathToManhattan.pdf'        
        png(file=figFileName,width = 3000,height = 1500,pointsize = 30)  
        manhattanCaseControl(ds1,ds2,.05)
        dev.off()  
        figFileName = 'pathToVolcano.pdf'        
        ds2$effectSizeNorm <- -ds2$effectSizeNorm
        volData <- rbind(cbind(ds1$effectSizeNorm,-log10(ds1$pvalWNorm),ds1$grtValNorm),cbind(ds2$effectSizeNorm,-log10(ds2$pvalWNorm),ds2$grtValNorm))
        volData <- data.frame(volData)
        png(file=figFileName,width = 1000,height = 1000,pointsize = 30)  
        print(ggplot(volData,aes(x=X1,y=X2,col=X3))+geom_point(aes(size=10))+theme_classic())+ylim(c(0,14))
        dev.off()  
        