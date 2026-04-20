def ReadCpGfromBigFile(MethFile,IndexFile,CpGWindowFile,SubsetFile,randomflag,offsetForIndex):
   ''' Ver01  Does not add Result prefix for SubsetFile
   Extract 450K data for specific CpG sites 
   INPUT
   MethFile  450K data with First column CpG IDs as used in 450KIllumina Chip 
   Indexfile Indexed filing pointing to offset for each CpG site, use Indexing450KFiles.py to generate the file
   CpGWindowFile list of CpG IDs as used in 450Illumina, one ID each line
   SubsetFile  filename where extracted data will be stored
   randomflag  if 0, CpG would be assumed to be in genomic order and search would be linear ''' 
   fo = open(MethFile,"r")
   Index = open(IndexFile,"r")
   CpGWindow = open(CpGWindowFile,"r")
   Subsetfo = open(SubsetFile,"w+")
   if not offsetForIndex:
      offsetForIndex=0  
   Index.seek(offsetForIndex) 
   for line in CpGWindow:
      if randomflag:
            Index.seek(offsetForIndex)
      notfound = 1
      while notfound:
        IndexLine = Index.readline()
        if IndexLine=='':
            notfound=0
        if line.split()[0] in IndexLine:
            fo.seek(int(IndexLine.split()[1]))
            Subsetfo.write("{}".format(fo.readline()))
            notfound = 0
   offsetForIndex = Index.tell()
   fo.close()
   Index.close()
   CpGWindow.close()
   Subsetfo.close()
   return offsetForIndex

