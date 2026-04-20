# INPUT 
#    PathToFig    folder to store data, use \ or / accordingly based on the OS 
#    data   in BigMemory desc format, rows CpG sites, col subjects 
#   probenameFile  desc format not able to keep rowname (some error), order of CpG site 
#                  names must be same as it is for data file 
CUsAcrossGenome <- function(annotationFile=paste0(mountPath,"/DNAMethylationPaper/CodeforCU/Illumina_probe_annotation.txt"),PathToFig = 'Results/',
                                                 data="/Volumes/Toshiba/DNAMethylationPaper/TCGAMissingData/GBM/AllInOneGene202108.txt",
                            probenameFile ="/Volumes/Toshiba//DNAMethylationPaper/TCGAMissingData/GBM/probeTCGA-GBM.txt" ,
                            indexFile = paste0(mountPath,'DNAMethylationPaper/TCGAMissingData/GBM/DataPrimary/Index.txt'),
                            GroupIndex,confounderX)
{
library(matrixStats)
library(dplyr)
source('GetProbeNamesForARange.R')
source_python('PythonIndexing//ReadCpGfromBigFileVer01.py')
#source('RBigMemoryLoadLibraries.R')
# Clustering Method 
source('ImgeEdgeDetection_10.R')
source('Residuals.R')
#####  Illumina 450k Data ##################
IllProbeAnnotation <- read.table(annotationFile,fill = TRUE,colClasses = 'character',header = TRUE)
IllProbeAnnotation <- IllProbeAnnotation[c("TargetID","CHR","MAPINFO")]
IllProbeAnnotation$TargetID <- as.character(IllProbeAnnotation$TargetID)
IllProbeAnnotation$CHR <- as.character(IllProbeAnnotation$CHR)
IllProbeAnnotation$MAPINFO <- as.numeric(IllProbeAnnotation$MAPINFO)


#####  Data ##################
# getting probe names 
probeNames  <- unlist(read.table(file=probenameFile,stringsAsFactors = FALSE))

#-----------------------------------------------------------------------------------------------#
windowSize = 250000; minProbeinWindow <- 5; ProbeNearestNbrDist <- 20000; NotLeaveProbebpAway <- 1000; cutoffCor <- .6
qqnormx <- function(x) {qqnorm(x,plot.it = FALSE)$x}


#####################################
########### Functions #####################################
#####################################
allCor <- function(X) 
{        
  aa <- t(apply(X,1,qqnormx))
  #aa <- X
  # missing values deleted case by case 
  XP <- cor(t(aa),method="pearson",use="complete.obs")
  XP
}

corForCluster <- function(XP)
{
  XP[abs(XP)<cutoffCor] <- 0
  XP
}


clusSaving <- function(FileName,aaaB,corControl)
{
  for (iclusplot in 1:length(clus))
  {
    localclusnames <- rownames(aaaB)[clus[[iclusplot]]]
    if (length(localclusnames)<1) { next}
     count <- count + 1    
     textFileDetails <- filter(IllProbeAnnotation,TargetID %in% localclusnames) %>% arrange(MAPINFO) %>% mutate(Distance = MAPINFO-MAPINFO[1])
     textFileDetails <- cbind(FileName,ichr,ProbeStartMAPINFO,paste0((textFileDetails[,1]),sep=";",collapse ="" ),paste0((textFileDetails[,3]),sep=";",collapse ="" ))
     write.table(textFileDetails,FileName,append = TRUE,row.names = FALSE,col.names = FALSE)
      }
}

#------------------------------#####################################
#-----------------------------########### Functions End ###########
#------------------------------#####################################

####################################################
########### Scanning Across Genome #####################################
####################################################
count <- 0 
for (ichr in 1:22)
{ 
  print(paste('Scanning Chr',ichr))
  localProbeChr <- IllProbeAnnotation[IllProbeAnnotation$CHR==ichr,]
  rownames(localProbeChr) <- localProbeChr$TargetID
  
  # Find probes that are available in data 
  localProbeChr <- localProbeChr[intersect(localProbeChr$TargetID,probeNames),]
  localProbeChr$MAPINFO<- as.numeric(localProbeChr$MAPINFO)
  localProbeChr <- localProbeChr[order(localProbeChr$MAPINFO),]
  # clean memory 
  gc()
  nlocalProbe  <- dim(localProbeChr)[1]
  localProbebpdist <- localProbeChr$MAPINFO[2]-localProbeChr$MAPINFO[1]
  localProbebpdist[2:(nlocalProbe -1)] <- pmin( -localProbeChr$MAPINFO[2:(nlocalProbe -1)]+localProbeChr$MAPINFO[3:nlocalProbe ],-localProbeChr$MAPINFO[1:(nlocalProbe -2)]+localProbeChr$MAPINFO[2:(nlocalProbe -1)])
  localProbebpdist[nlocalProbe ] <- localProbeChr$MAPINFO[nlocalProbe ]-localProbeChr$MAPINFO[nlocalProbe -1]
  ProbeStartMAPINFO <- localProbeChr$MAPINFO[1]
  ProbeEndMAPINFO <- localProbeChr$MAPINFO[nlocalProbe ]
  # probes whos nearest neighbour is more than 20k away are removed 
  localProbeChr <- localProbeChr[localProbebpdist<ProbeNearestNbrDist,]
  nlocalProbe  <- dim(localProbeChr)[1]
  localProbebpdist <- localProbeChr$MAPINFO[2]-localProbeChr$MAPINFO[1]
  localProbebpdist[2:(nlocalProbe -1)] <- pmin( -localProbeChr$MAPINFO[2:(nlocalProbe -1)]+localProbeChr$MAPINFO[3:nlocalProbe ],-localProbeChr$MAPINFO[1:(nlocalProbe -2)]+localProbeChr$MAPINFO[2:(nlocalProbe -1)])
  localProbebpdist[nlocalProbe ] <- localProbeChr$MAPINFO[nlocalProbe ]-localProbeChr$MAPINFO[nlocalProbe -1]
  ProbeStartMAPINFO <- localProbeChr$MAPINFO[1]
  ProbeEndMAPINFO <- localProbeChr$MAPINFO[nlocalProbe ]
  
  Nwin = ceiling((ProbeEndMAPINFO - ProbeStartMAPINFO)/windowSize)
  # scaning within a window 
  Scanned <- FALSE
  iwin <- 0
  print("Starting")
  while(ProbeStartMAPINFO<ProbeEndMAPINFO)
  { CloseNbr <- 0
  print(ProbeStartMAPINFO)
    #if (!iwin %%  500) { print(paste('Win',iwin))}
    MInd <- localProbeChr[localProbeChr$MAPINFO >= ProbeStartMAPINFO & localProbeChr$MAPINFO <= ProbeStartMAPINFO + windowSize,"TargetID"]
    localLastProbe <- which(localProbeChr$TargetID %in% MInd[length(MInd)])
    # Adding probes outside windowsize but close to the last probe within the window 
    trash <- diff(localProbeChr$MAPINFO[localLastProbe:nlocalProbe])
    CloseNbr <- Position(function(trash) trash > NotLeaveProbebpAway,trash,nomatch = 0)
    winmapinfo <- ProbeStartMAPINFO
    if  (CloseNbr>1)
    { 
        MInd <- c(MInd,localProbeChr$TargetID[seq((localLastProbe+1),(localLastProbe+CloseNbr-1),1)])    
    } 
    ProbeStartMAPINFO <- localProbeChr$MAPINFO[which(localProbeChr$TargetID%in%MInd[length(MInd)])+CloseNbr]  
   
    if (length(MInd)>minProbeinWindow)
    { 
    cpgclusnames <- MInd[ MInd %in% probeNames]
    write.table(cpgclusnames,file = paste0(PathToFig,'trash.txt'),col.names = FALSE,row.names = FALSE,quote = FALSE)
    offsetIndex <- 0
    offsetIndexNew = ReadCpGfromBigFile(data,indexFile,paste0(PathToFig,'trash.txt'),paste0(PathToFig,'MethTrash.txt'),0,offsetIndex)
    # might now work SCM data 
    M=read.table(paste0(PathToFig,'MethTrash.txt'),stringsAsFactors = FALSE,header = FALSE,row.names = 1)    
    rownames(M) <- cpgclusnames
    # atleast 30 subjects with non-zero meth value  
    M <- M[rowSums(!is.na(M))>10,]
    if (nrow(M)<minProbeinWindow) next 
    chrBVal <- M 
    if (hasArg(GroupIndex))
      {chrBVal <-  M[,GroupIndex]}
    if (hasArg(confounderX))
    { if (hasArg(GroupIndex)) 
    {chrBVal <- Residuals(confounderX[GroupIndex,],chrBVal)}
      else 
      {chrBVal <- Residuals(confounderX,chrBVal)} 
    }
     M <- chrBVal
     
      cutoffCor <- .4 
      corM <- allCor(M)
      ################## 
      aaaB <- corForCluster(corM)
      if (length(aaaB)>=minProbeinWindow)
      {
        
      Bothclus <- ImageEdgeDetection(aaaB)
      if (length(Bothclus)!=0)
      {  
      clus <- Bothclus[[1]]
      
      if (length(clus)!=0) {
      FileName <- paste(PathToFig,'NonContM',cutoffCor,'.txt',sep = "")
      clusSaving(FileName,aaaB,corM)
      }
      
      clus <- Bothclus[[2]]
      
      if (length(clus)!=0) {
        FileName <- paste(PathToFig,'ContM',cutoffCor,'.txt',sep = "")
        clusSaving(FileName,aaaB,corM)
      }
      
      } # if BothClus lenght bigger than 1
      }  #  if (length(aaaB)>=16)
      
     
      cutoffCor <- .6
      corM <- allCor(M)
      ################## 
      aaaB <- corForCluster(corM)
      if (length(aaaB)>=minProbeinWindow)
      {
        
        Bothclus <- ImageEdgeDetection(aaaB)
        if (length(Bothclus)!=0)
        {  
          clus <- Bothclus[[1]]
          
          if (length(clus)!=0) {
            FileName <- paste(PathToFig,'NonContM',cutoffCor,'.txt',sep = "")
            clusSaving(FileName,aaaB,corM)
          }
          
          clus <- Bothclus[[2]]
          
          if (length(clus)!=0) {
            FileName <- paste(PathToFig,'ContM',cutoffCor,'.txt',sep = "")
            clusSaving(FileName,aaaB,corM)
          }
          
        } # if BothClus lenght bigger than 1
      }  #  if (length(aaaB)>=16)
      
     
      cutoffCor <- .8 
      corM <- allCor(M)
      ################## 
      aaaB <- corForCluster(corM)
      if (length(aaaB)>=minProbeinWindow)
      {
        Bothclus <- ImageEdgeDetection(aaaB)
        if (length(Bothclus)!=0)
        {  
          clus <- Bothclus[[1]]
          
          if (length(clus)!=0) {
            FileName <- paste(PathToFig,'NonContM',cutoffCor,'.txt',sep = "")
            clusSaving(FileName,aaaB,corM)
          }
          
          clus <- Bothclus[[2]]
          
          if (length(clus)!=0) {
            FileName <- paste(PathToFig,'ContM',cutoffCor,'.txt',sep = "")
            clusSaving(FileName,aaaB,corM)
          }
          
        } # if BothClus lenght bigger than 1
      }  #  if (length(aaaB)>=16)
      
      
      
    }
  } # while scanning a chromosome 
} # while scannning across chromosome 
} # end of function 
