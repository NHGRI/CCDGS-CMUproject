pvalForDiffCor <- function(filename1,IndexFile1,indexControl,indexCase,X=NULL,ResultFileControl,ResultFileCase,pvalFileName)
{
  library(psych)
  source_python('PythonIndexing/ReadCpGfromBigFileVer01.py')
  library(matrixStats)
  qqnormx <- function(x) {qqnorm(x,plot.it = FALSE)$x}
  
  clusDetailsControl <- read.table(ResultFileControl,stringsAsFactors = FALSE)
  
  NClusControl <- dim(clusDetailsControl)[1]
  
  offsetIndex <- 0 
  
  pvalCcortest <- matrix(0,NClusControl) 
  pvalCwilcox <- matrix(0,NClusControl) 
  pvalCcortestNorm <- matrix(0,NClusControl) 
  pvalCwilcoxNorm <- matrix(0,NClusControl) 
  effectSizeNorm <- matrix(0,NClusControl) 
  psudoEffectSizeNorm <- matrix(0,NClusControl) 
  
  presentWindow <- clusDetailsControl[1,3]
  
  for (i in 1:NClusControl)
  {
    if (presentWindow!=clusDetailsControl[i,3])
    {presentWindow=clusDetailsControl[i,3]
    offsetIndex = offsetIndexNew
    }
    
    cpgclusnames <- unlist(strsplit(clusDetailsControl[i,4],split = ";"))  
    write.table(cpgclusnames,file = 'trash.txt',col.names = FALSE,row.names = FALSE,quote = FALSE)
    # input 1 can be changed to 0 to speed up the process but all the cpg must be in linear order 
    offsetIndexNew = ReadCpGfromBigFile(filename1,IndexFile1,'trash.txt','MethTrash.txt',1,offsetIndex)
    # set new offset after window changed, within window things may not be linear 
    
    M=read.table('MethTrash.txt',stringsAsFactors = FALSE,header = FALSE,row.names = 1)
    if (!is.null(X))
    {MControl <- Residuals(X[indexControl,],M[,indexControl])
    MCase <- Residuals(X[indexCase,],M[,indexCase])} else
    {
      MControl <-M[,indexControl]
      MCase <- M[,indexCase]
    }
    conCor <- cor(t(MControl))
    caseCor <- cor(t(MCase))
    pvalCcortest[i] <- cortest(conCor,caseCor,dim(MControl)[2],dim(MCase)[2],fisher=TRUE)$prob
    pvalCwilcox[i] <- wilcox.test( conCor[lower.tri(conCor)], caseCor[lower.tri(caseCor)], paired=FALSE,alternative = "greater",mu=0.1)[3]$p.value
    
    M[,indexControl] <- t(apply(M[,indexControl],1,qqnormx))
    M[,indexCase] <- t(apply(M[,indexCase],1,qqnormx))
    
    if(!is.null(X))
    {MControl <- Residuals(X[indexControl,],M[,indexControl])
    MCase <- Residuals(X[indexCase,],M[,indexCase])}    
    else
    {
      MControl <-M[,indexControl]
      MCase <- M[,indexCase]
    }
    conCor <- cor(t(MControl))
    caseCor <- cor(t(MCase))
    pvalCcortestNorm[i] <- cortest(conCor,caseCor,dim(MControl)[2],dim(MCase)[2],fisher=TRUE)$prob
    wilcoxResult <- wilcox.test( conCor[lower.tri(conCor)], caseCor[lower.tri(caseCor)], paired=FALSE,alternative = "greater",mu=0.1)
    pvalCwilcoxNorm[i] <- wilcoxResult[3]$p.value
    #  WilcoxDiff <- conCor[lower.tri(conCor)] - caseCor[lower.tri(caseCor)]
    #WilcoxDiff[abs(WilcoxDiff)<.1] <- 0
    effectSizeNorm[i] <- wilcoxResult[1]$statistic/(sum(lower.tri(conCor)))^2
    psudoEffectSizeNorm[i] <- median(conCor[lower.tri(conCor)]) - median(caseCor[lower.tri(caseCor)])
  }
  save(file=pvalFileName,list=c("pvalCcortest","pvalCwilcox","pvalCcortestNorm","pvalCwilcoxNorm","effectSizeNorm","psudoEffectSizeNorm"))
}

