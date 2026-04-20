#!/bin/bash

### Script extracts sample ID and renames the file with the sample ID


picard=/home/bin/picard/picard.jar

INPUTDIR=/home/lesediw/gVCFs/WES_august_mapped
OUTPUTDIR=/home/lesediw/gVCFs/SampleRenamedGvcfs/August

cd $INPUTDIR

#set up naming loop

for infile in *.snps.indels.g.vcf.gz
do
    #obtain sample id prefix
    sname=$(basename $infile | cut -d\. -f1)
    echo "Started renaming process for $sname at $(date)"
    echo " "
    java -Xmx150G -jar  $picard RenameSampleInVcf I=$infile O=$OUTPUTDIR/$sname.raw.snps.indels.g.vcf.gz NEW_SAMPLE_NAME=$sname
    echo "Run for $sname comlpleted at: $(date)"
    echo " "
done