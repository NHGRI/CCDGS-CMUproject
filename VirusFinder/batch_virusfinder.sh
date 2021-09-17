#!/bin/bash
## This script will run all the fastq files in the path BAMDIR


WORKDIR=/mnt/hanchardlab/savannah/Software/VirusFinder2.0
BAMDIR=/mnt/hanchardlab/savannah/CAFGEN_BAM

cd ${WORKDIR}

for file in ${BAMDIR}/*_R1_001.fastq.gz
do
   OUTPUTFILENAME=`basename ${file} _R1_001.fastq.gz`
   FILEPATH1=/mnt/hanchardlab/savannah/CAFGEN_BAM/${OUTPUTFILENAME}_R1_001.fastq.gz
   FILEPATH2=/mnt/hanchardlab/savannah/CAFGEN_BAM/${OUTPUTFILENAME}_R2_001.fastq.gz
   echo "Starting ${file}"
   echo "Basename of ${file} is ${OUTPUTFILENAME}"
   mkdir ./${OUTPUTFILENAME}
   cd ${OUTPUTFILENAME}
   cp /path/to/config.txt ./config2.txt
   echo "${FILEPATH1} ${FILEPATH2}"
   sed -e 's/\/scratch\/kingw\/virusFinder\/simulation\/seq_1.fastq.gz/\/path\/to\/$OUTPUTFILENAME\_R1_001.fastq.gz/g' -e 's/\scratch\/kingw\/virusFinder\/simulation\/seq_2.fastq.gz/\/path\/to\/$OUTPUTFILENAME\_R2_001.fastq.gz/g' ./config2.txt>config.txt
   rm config2.txt
  perl /path/to/VirusFinder.pl -c ./config.txt -o ${WORKDIR}/${OUTPUTFILENAME}
   echo "Finished ${file}. Created file ${OUTPUTFILENAME}"
   cd ../
done

