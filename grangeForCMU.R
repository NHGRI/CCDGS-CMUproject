grangeForCMU <- function(pvalFileName = paste0('~/Work/DNAMethylationPaper/CodeforCU/SCM/Results/',itype,iloc,icor,"ESCMControl.Rdata"),
                             ResultFileControl = paste0('~/Work/DNAMethylationPaper/CodeforCU/SCM/Results/',itype,iloc,"ESCM/ContM",icor,".txt"),
                         pvalCthresh=.01,pvalWthresh=.01,Ncpg=3)
{
  clusDetailsControl <- read.table(ResultFileControl,stringsAsFactors = FALSE)
  load(pvalFileName)
  clusDetailsControl[,3] <- getStart(clusDetailsControl)
  # keeping start and end same to be able to plot manhattan plot if range given it plots bar instead of dots
  CMUIndex <- getNCpG(clusDetailsControl)>Ncpg
  clusDetailsControl <- clusDetailsControl[CMUIndex,]
  pvalCcortest <- pvalCcortest[CMUIndex]
  pvalCwilcox <- pvalCwilcox[CMUIndex]
  pvalCcortestNorm <- pvalCcortestNorm[CMUIndex]
  pvalCwilcoxNorm <- pvalCwilcoxNorm[CMUIndex]
  effectSizeNorm <- effectSizeNorm[CMUIndex]
  ds <- cbind(clusDetailsControl[,2:3],clusDetailsControl[,3])
  colnames(ds) <- c("chr","start","end")
  ds$chr <- paste0("chr",ds$chr)
  ds <- makeGRangesFromDataFrame(ds)
  pvalCcortest[pvalCcortest==0] <-  min(pvalCcortest[pvalCcortest>0])
  grtVal <- p.adjust(pvalCwilcox,"bonferroni")<pvalWthresh&p.adjust(pvalCcortest,"bonferroni")<pvalCthresh
  pvalCcortestNorm[pvalCcortestNorm==0] <-  min(pvalCcortestNorm[pvalCcortestNorm>0])
  grtValNorm <- p.adjust(pvalCwilcoxNorm,"bonferroni")<pvalWthresh&p.adjust(pvalCcortestNorm,"bonferroni")<pvalCthresh
  values(ds) <- data.frame(effectSizeNorm=effectSizeNorm,pval=pvalCcortest,pvalW=pvalCwilcox,grtVal=grtVal,pvalNorm=pvalCcortestNorm,pvalWNorm=pvalCwilcoxNorm,grtValNorm=grtValNorm,End=getEnd(clusDetailsControl),cpgNames=clusDetailsControl$V4,cpgLocation=clusDetailsControl$V5)
  ds
}
