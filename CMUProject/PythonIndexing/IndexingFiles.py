def IndexingFiles(filename,IndexFilename,Sorted450K):
  # delimeter is changed \t 
  import pandas as pd 
  import numpy as np  
  fo = open(filename)
  Index = open(IndexFilename,"w+")
  # get the offset for first line 
  lineOffset = fo.tell()
  line = fo.readline()
  while line:
  # SPLIT BY TAB explicitly
    parts = line.split('\t')
    if len(parts) > 0:
        cgId = parts[0]
        cgId = cgId.replace('\"', '')
        if ("cg" in cgId ) or ("ch" in cgId):
            Index.write("{} \t {} \n".format(cgId,lineOffset))
    lineOffset = fo.tell()      
    line = fo.readline()
  fo.close()
  Index.close()
  CpGf = open(Sorted450K)
  Dict = {}
  count=0
  for line in CpGf:
     CpGId = line.rstrip().split()[0]
     Dict[CpGId] = count
     count=count+1   
  Index= pd.read_csv(IndexFilename,sep='\t',header=None)
  Index[0] = Index[0].map(str.strip)
  a = Index[0].apply(lambda x: Dict[x])
  b=a.sort_values()
  Index = Index.reindex(b.index)
  np.savetxt(IndexFilename,Index.values,fmt=('%s','%d'),delimiter='\t')
