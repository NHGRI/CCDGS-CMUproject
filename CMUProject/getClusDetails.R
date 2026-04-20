getEnd <- function(clusDetailForPosition)
{
  unlist(lapply(clusDetailForPosition$V5,a <- function(x) {stenCpG  <- strsplit(x,split=";"); as.numeric(stenCpG[[1]][length(stenCpG[[1]])])}))
}

getStart <- function(clusDetailForPosition)
{
  unlist(lapply(clusDetailForPosition$V5,a <- function(x) {as.numeric(strsplit(x,split=";")[[1]][1])}))
}


getNCpG <- function(clusDetail)
{
  unlist(lapply(strsplit(clusDetail$V5,split=";"), function(x) {length(x)}))
}



createContPlusNonCont <- function(ResultFileCont,ResultFileControl,pvalFileNameCont,
                                  pvalFilenameControl,ResultFileContPlusNonCont,pvalFileNameContPlusNonCont)
{ library(data.table)
  # get clus for continuous 
  clusDetailsCont <- read.table(ResultFileCont,stringsAsFactors = FALSE)          
  # get non cont clus 
  clusDetailsControl <-  read.table(ResultFileControl,stringsAsFactors = FALSE)          
  # Find cont clus not in non cont 
  keep <- apply(clusDetailsCont,1, function(x) sum((clusDetailsControl$V2==as.numeric(x[2]) &
              clusDetailsControl$V3==as.numeric(x[3]))[clusDetailsControl$V4[(clusDetailsControl$V2==as.numeric(x[2]) &
              clusDetailsControl$V3==as.numeric(x[3]))] %like% x[4]])==0)
 # sub-setting everything    clusDetailsCont <- clusDetailsCont[keep,]
  load(pvalFileNameCont)
  pvalCcortesttemp <- pvalCcortest[keep]
  pvalCcortestNormtemp <- pvalCcortestNorm[keep]
  pvalCwilcoxNormtemp <- pvalCwilcoxNorm[keep]
  pvalCwilcoxtemp <- pvalCwilcox[keep]
  effectSizeNormtemp <- effectSizeNorm[keep]
  psudoEffectSizeNormtemp <- psudoEffectSizeNorm[keep]
  
  # merge with non cont 
  load(pvalFilenameControl)
  clusDetailContNonCont <- rbind(clusDetailsControl,clusDetailsCont[keep,])
  orderChange <- order(as.numeric(clusDetailContNonCont[,2]),as.numeric(clusDetailContNonCont[,3]))
  clusDetailContNonCont <- clusDetailContNonCont[orderChange,]
  write.table(clusDetailContNonCont,file=ResultFileContPlusNonCont,row.names = FALSE,col.names = FALSE)
  pvalCcortest <- c(pvalCcortest,pvalCcortesttemp)[orderChange] 
  pvalCcortestNorm <- c(pvalCcortestNorm,pvalCcortestNormtemp)[orderChange] 
  pvalCwilcoxNorm<- c(pvalCwilcoxNorm,pvalCwilcoxNormtemp)[orderChange] 
  pvalCwilcox<- c(pvalCwilcox,pvalCwilcoxtemp)[orderChange] 
  effectSizeNorm<- c(effectSizeNorm,effectSizeNormtemp)[orderChange] 
  psudoEffectSizeNorm<- c(psudoEffectSizeNorm,psudoEffectSizeNormtemp)[orderChange] 
  save(file=pvalFileNameContPlusNonCont,list=c("pvalCcortest","pvalCwilcox","pvalCcortestNorm","pvalCwilcoxNorm","effectSizeNorm","psudoEffectSizeNorm"))
  
}


RatioContThatMakesNonCont <- function(ResultFileCont,ResultFileControl)
{ library(data.table)
  # get clus for continuous 
  clusDetailsCont <- read.table(ResultFileCont,stringsAsFactors = FALSE)          
  # get non cont clus 
  clusDetailsControl <-  read.table(ResultFileControl,stringsAsFactors = FALSE)          
  # Find cont clus not in non cont 
 # keep <- apply(clusDetailsCont,1, function(x) sum((clusDetailsControl$V2==as.numeric(x[2]) &
 #        clusDetailsControl$V3==as.numeric(x[3]))[clusDetailsControl$V4[(clusDetailsControl$V2==as.numeric(x[2]) &
 #        clusDetailsControl$V3==as.numeric(x[3]))] %like% x[4]])==0)
keep <- apply(clusDetailsCont,1, function(x) length(grep( x[4],clusDetailsControl$V4[(clusDetailsControl$V2==as.numeric(x[2]) &
                                                                                        clusDetailsControl$V3==as.numeric(x[3]))],perl = TRUE))==0)
  
  table(keep)
}

NumberOfContInNonCont <- function(ResultFileCont,ResultFileControl)
{ library(data.table)
  # get clus for continuous 
  clusDetailsCont <- read.table(ResultFileCont,stringsAsFactors = FALSE)          
  # get non cont clus 
  clusDetailsControl <-  read.table(ResultFileControl,stringsAsFactors = FALSE)          
  # Find cont clusin non cont 
  
  keep <- apply(clusDetailsControl,1, function(x) sum(unlist(sapply(clusDetailsCont$V4[(clusDetailsCont$V2==as.numeric(x[2]) &
            clusDetailsCont$V3==as.numeric(x[3]))],function(y) grep(y,x[4],perl = TRUE)))))
  print(sum(keep==1))
  keep <- keep[keep!=1]
  data.frame(table(keep))
}

CMUCpGByChrWin <- function(ResultFileCont)
{
  clusDetailsCont <- read.table(ResultFileCont,stringsAsFactors = FALSE)
  clusDetailsCont$chrWin <- paste(clusDetailsCont$V2,clusDetailsCont$V3,sep="_")
  temp <- aggregate(clusDetailsCont[,c("V4")],by=list(clusDetailsCont$chrWin),paste,collapse="")
  temp1 <- aggregate(clusDetailsCont[,c("V3")],by=list(clusDetailsCont$chrWin),length)
  colnames(temp1) <- c("X1","X3")
  full_join(data.frame(cbind(temp$Group.1,unlist(lapply(strsplit(temp$x,split=";"), function(x) {length(unique(x))})))),
            data.frame(temp1))
  
}