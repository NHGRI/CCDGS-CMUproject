#!/bin/bash

##This is a template script to extract the viral contigs detected from all samples after running VirusFinder2

WORKDIR=/path/to/working/directory
## head -1 ${WORKDIR}/C6Y23ACXX-2-ID02/virus-list.txt >> ${WORKDIR}/CombinedVirusList.txt
## edit the above line to extract header info from virus-list.txt
for dir in ${WORKDIR}/*
do
   cd ${dir}
   OUTPUTFILENAME=${dir}
   awk "{if(NR!=1) print \$0\"\t${OUTPUTFILENAME}\"}" virus-list.txt >>  ${WORKDIR}/CombinedVirusList.txt
   cd ../
done

