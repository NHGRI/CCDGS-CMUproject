# same as Ver01 removed bigmemory related creations 
GenerateDataFormat <- function(filename = "Mimic_EPIC_Exact.txt",
                               PathToResult = "./CMU_EPIC",
                               SortedProbeFilename = "./PythonIndexing/SortedCpGNames450k.txt")
{
  
  # create folder for PathToResult
  system(paste('mkdir',PathToResult))
  # create Results inside  
  ResultPath <- paste(PathToResult,'/Results/',sep="")
  system(paste('mkdir',ResultPath))
  # create Data to PathToResults     
  DataPath <- paste(PathToResult,'/Data/',sep="")
  system(paste('mkdir',DataPath))
  
  # package to use python function 
  library(reticulate)
  # create directory Data under PathToResults 
  # give the path for python that you want to use, needs pandas and numpy installed
  # you can check which python on terminal if you do not know the version 
  # Also, you can install packages using reticulate if you do not have things already installed 
  # e.g. use_python('/home/ajajoo/anaconda3/bin/python3.7m') 
  # py_install("pandas")
  source_python('./PythonIndexing/IndexingFiles.py')
  # Index the data so a given line can be quickly accessed without going over all end characteres before it 
  IndexingFiles(filename,paste0(DataPath,'Index.txt'),SortedProbeFilename)
}

