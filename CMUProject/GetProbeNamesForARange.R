GetProbeNamesForARange <- function(annotationFile="~/Work/DNAMethylationPaper/CodeforCU/Illumina_probe_annotation.txt",chr=1,startRange=1000,endRange=250000)
  
{
  IllProbeAnnotation <- read.table(annotationFile,fill = TRUE,colClasses = 'character',header = TRUE)
  IllProbeAnnotation <- IllProbeAnnotation[c("TargetID","CHR","MAPINFO","UCSC_REFGENE_NAME")]
  IllProbeAnnotation$TargetID <- as.character(IllProbeAnnotation$TargetID)
  IllProbeAnnotation$CHR <- as.character(IllProbeAnnotation$CHR)
  IllProbeAnnotation$MAPINFO <- as.numeric(IllProbeAnnotation$MAPINFO)
  
  # Get probes of chr 
  chrIll <- IllProbeAnnotation[IllProbeAnnotation$CHR==chr,]
  probeIll <- chrIll[chrIll$MAPINFO< endRange & chrIll$MAPINFO> startRange,]
  probeIll <- probeIll[order(probeIll$MAPINFO),]
  probeIll$TargetID <- noquote(probeIll$TargetID)
  probeIll
}