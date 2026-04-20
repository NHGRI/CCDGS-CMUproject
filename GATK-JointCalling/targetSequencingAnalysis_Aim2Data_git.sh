#! /bin/bash

# directories and paths
WORKDIR=/data/Hanserv/NMAMP
cd $WORKDIR

#### mode and modules
sinteractive --mem=100g --cpus-per-task=8 --time 24:00:00

#### 1) alignment
mkdir bwa_output_aim2
OUTDIR=$WORKDIR/bwa_output_aim2

# read ID format that meet the GATK4 requirements
ID:NB551314_126_H5H3JAFX5_1 SM:NB551314_126 PU:H5H3JAFX5.1.CTAGCGCT LB:NB551314_126.lib1 PL:ILLUMINA
ID:NB551314_126_H5H3JAFX5_2 SM:NB551314_126 PU:H5H3JAFX5.2.CTAGCGCT LB:NB551314_126.lib1 PL:ILLUMINA
ID:NB551314_126_H5H3JAFX5_3 SM:NB551314_126 PU:H5H3JAFX5.3.CTAGCGCT LB:NB551314_126.lib1 PL:ILLUMINA
ID:NB551314_126_H5H3JAFX5_4 SM:NB551314_126 PU:H5H3JAFX5.4.CTAGCGCT LB:NB551314_126.lib1 PL:ILLUMINA

BWAINDEX=/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa
FASTQDIR=$WORKDIR/NMAMP2_Aim2AmpliseqData/NMAMP2-385496369/FASTQ_Generation_2023-04-04_20_35_17Z-663216733

for i in $(ls $FASTQDIR); do
FILE1=$(ls $FASTQDIR/$i | head -n 1)
FILE2=$(ls $FASTQDIR/$i | tail -n 1)
NAME=$(echo $FILE1 | cut -d _ -f 1-3)
SAMPLEID=$(echo $FILE1 | cut -d _ -f 1)
header=$(zcat $FASTQDIR/$i/$FILE1 | head -n 1)
id=$(echo $header | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
#sm=$(echo $header | cut -f 1-2 -d":" | sed 's/@//' | sed 's/:/_/g')
pu=$(echo $header | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
idx=$(echo $header | grep -Eo "[ATGCN]+$")
#echo \"@RG\\tID:$id\\tSM:$sm\\tPU:$pu"."$idx\\tLB:$sm".lib1"\\tPL:ILLUMINA\"
echo "bwa mem -M -O 30 -E 4 -T 20 -v 3 -t 36 -R $(echo \"@RG\\tID:$id\\tSM:$SAMPLEID\\tPU:$pu"."$idx\\tLB:$sm".lib1"\\tPL:ILLUMINA\") \
$BWAINDEX $FASTQDIR/$i/$FILE1 $FASTQDIR/$i/$FILE2 \
| samblaster -M | samtools fixmate - - | samtools sort -t 36 -O bam -o $OUTDIR/${NAME}_sorted.bam"
done >> ./scripts/01_bwa_aim2.swarm

swarm --file ./scripts/01_bwa_aim2.swarm --module bwa,samblaster,samtools --gb-per-process 32 -t 8 --verbose 1 --time 8:00:00 

#### 2) merge the 4 replicates bam files by samples
mkdir bam_aim2
BAMFILES=$WORKDIR/bam_aim2

# 2.1 get the unique sample list
for i in $(ls $OUTDIR); do 
SAMPLEID=$(echo $i | cut -d _ -f 1)
echo $SAMPLEID
done >> smapleList_aim2.txt

cat smapleList_aim2.txt | uniq > sample_list_aim2.txt
rm smapleList_aim2.txt

# 2.2 get the merge swarm file
for s in $(cat sample_list_aim2.txt); do 
echo "samtools merge $BAMFILES/${s}_merged.bam `find $OUTDIR -name ${s}_*_sorted.bam | xargs`"
done >> ./scripts/02_merge_aim2.swarm

swarm --file ./scripts/02_merge_aim2.swarm --module samtools --gb-per-process 32 -t 8 --verbose 1 --time 2:00:00 

#### 3) sort and index the merged bam files
for b in $(ls $BAMFILES); do
SAMPLEID=$(echo $b | cut -d . -f 1)
samtools sort $BAMFILES/$b -t 36 -O bam -o $BAMFILES/${SAMPLEID}_sorted.bam
samtools index $BAMFILES/${SAMPLEID}_sorted.bam
done

#### 4) mark PCR duplicates with GATK4 MarkDuplicatesSpark, Ref: https://hpc.nih.gov/training/gatk_tutorial/markdup.html
for m in $(ls $BAMFILES/*merged_sorted.bam); do
SAMPLEID=$(basename $m | cut -d . -f 1)
echo "gatk --java-options \"-Djava.io.tmpdir=$WORKDIR/tmp -Xms60G -Xmx60G\" MarkDuplicatesSpark \
-I $m \
-O $BAMFILES/${SAMPLEID}_markdup.bam \
--spark-runner LOCAL"
done >> ./scripts/04_markduplicates.swarm 

swarm --file ./scripts/04_markduplicates.swarm --module GATK/4.2.1.0 -t 20 -g 80 --gres=lscratch:500 --time=12:00:00 

#### 5) BaseRecalibrator (model build) and ApplyBQSR, Ref: https://hpc.nih.gov/training/gatk_tutorial/bqsr.html
# 5.1 BaseRecalibrator
for n in $(ls $BAMFILES/*merged_sorted_markdup.bam); do
SAMPLEID=$(basename $n | cut -d . -f 1)
echo "gatk --java-options \"-Djava.io.tmpdir=$WORKDIR/tmp -Xms4G -Xmx4G -XX:ParallelGCThreads=2\" BaseRecalibrator \
-I $n \
-R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
-O $BAMFILES/${SAMPLEID}_bqsr.report \
--known-sites /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
--known-sites /fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
--known-sites /fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
done >> ./scripts/05_baseRecalibrator.swarm 

swarm --file ./scripts/05_baseRecalibrator.swarm --module GATK/4.2.1.0 -t 2 -g 4 --time=24:00:00 --gres=lscratch:400

# 5.2 ApplyBQSR
for n in $(ls $BAMFILES/*merged_sorted_markdup.bam); do
SAMPLEID=$(basename $n | cut -d . -f 1)
echo "gatk --java-options \"-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2\" ApplyBQSR \
-I $n \
-R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
--bqsr-recal-file $BAMFILES/${SAMPLEID}_bqsr.report \
-O $BAMFILES/${SAMPLEID}_bqsr.bam"
done >> ./scripts/06_applyBQSR.swarm 

swarm --file ./scripts/06_applyBQSR.swarm --module GATK/4.2.1.0 -t 2 -g 2 --gres=lscratch:400 --time=24:00:00

#### 6) Joint Calling
mkdir vcf
VCFFILES=./vcf

# 6.1 HaplotypeCaller
# test for one sample
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2" HaplotypeCaller \
  -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
  -I NA12878_markdup_bqsr.bam \
  -O NA12878.g.vcf.gz \
  -ERC GVCF

for i in $(ls $BAMFILES/*merged_sorted_markdup_bqsr.bam); do
SAMPLEID=$(basename $i | cut -d . -f 1)
echo "gatk --java-options \"-Djava.io.tmpdir=$WORKDIR/tmp -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" HaplotypeCaller \
-R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
-I $i \
-O $VCFFILES/${SAMPLEID}.g.vcf.gz \
-ERC GVCF";
done >> ./scripts/08_haplotypeCaller.swarm

swarm --file ./scripts/08_haplotypeCaller.swarm --module GATK/4.2.1.0 -t 2 -g 20 --time=40:00:00

# re-do index files which might cause the GenomicsDBImport issue of "invalid uncompressedLength"
rm vcf/*vcf.gz.tbi

for i in $(ls $VCFFILES/*.g.vcf.gz); do
SAMPLEID=$(basename $i | cut -d . -f 1)
echo "bcftools index -f -t $i";
done >> ./scripts/08_index_vcf.swarm

swarm --file ./scripts/08_index_vcf.swarm --module bcftools -t 2 -g 20 --time=2:00:00

# 6.2 GenomicsDBImport (replaces CombineGVCFs)
mkdir ./vcf/gdb
# testing with fewer samples
for j in {1..22} X Y M; do
gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenomicsDBImport \
--genomicsdb-workspace-path ./vcf/gdb/chr${j}_gdb \
-R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
-V vcf/PA2098_merged_sorted_markdup_bqsr.g.vcf.gz \
-V vcf/PA2100_merged_sorted_markdup_bqsr.g.vcf.gz \
-V vcf/PA2102_merged_sorted_markdup_bqsr.g.vcf.gz \
--tmp-dir "$WORKDIR/tmp" \
--max-num-intervals-to-import-in-parallel 3 \
--intervals chr${j};
done

# for all of the samples
sbatch --cpus-per-task=2 --mem=2g --gres=lscratch:100 --time=12:00:00 ./scripts/09_GATK_GenomicsDBImport.sh

# 6.3 GenotypeGVCFs
mkdir ./vcf/gvcf

# testing with one chr
gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenotypeGVCFs \
  -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
  -V gendb://vcf/gdb/chr10_gdb -O ./vcf/gvcf/chr10.vcf.gz

# for all of the chr
for j in {1..22} X Y M; do
echo "gatk --java-options \"-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2\" GenotypeGVCFs \
-R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
-V gendb://vcf/gdb/chr${j}_gdb -O ./vcf/gvcf/chr${j}.vcf.gz";
done >> ./scripts/10_GATK_GenotypeGVCFs.swarm 

swarm --file ./scripts/10_GATK_GenotypeGVCFs.swarm  --module GATK/4.2.1.0  -t 2 -g 2 --time=2:00:00

# 6.4 The VCF files from each chromosome should be merged to one compressed file before next step, here we use picard module to do the merge:

sbatch --cpus-per-task=2 --mem=2G --time=2:00:00 ./scripts/11_picard_gatherVCFs.sh

# index the merge vcf
module load bcftools
bcftools index -t vcf/gvcf/merged.vcf.gz 

#### 7) VQSR
# 7.1 VariantRecalibrator
cd vcf/gvcf/; \
gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 \
  -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
  -V merged.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /fdb/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz  \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 /fdb/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 /fdb/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP -O merged_SNP1.recal --tranches-file output_SNP1.tranches \
  --rscript-file output_SNP1.plots.R

cd vcf/gvcf/; \
gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R /fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
  -V merged.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 /fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
  --max-gaussians 4 \
  -mode INDEL -O merged_indel1.recal --tranches-file output_indel1.tranches \
  --rscript-file output_indel1.plots.R

swarm --file scripts/12_GATK_VariantRecalibrator.sh --module GATK/4.2.1.0 -t 2 -g 4 --time=2:00:00

# 7.2 ApplyVQSR
#! /bin/bash
# run ApplyVQSR on SNP then INDEL
module load GATK/4.2.1.0
cd vcf/gvcf/;
gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
  -V merged.vcf.gz \
  --recal-file merged_SNP1.recal \
  -mode SNP \
  --tranches-file output_SNP1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O SNP.recalibrated_99.9.vcf.gz

gatk --java-options "-Djava.io.tmpdir=$WORKDIR/tmp -Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
  -V SNP.recalibrated_99.9.vcf.gz \
  -mode INDEL \
  --recal-file merged_indel1.recal \
  --tranches-file output_indel1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O indel.SNP.recalibrated_99.9.vcf.gz

sbatch --cpus-per-task=2 --mem=2g --time=2:00:00 ./scripts/13_GATK_ApplyVQSR.sh

# annotation by ANNOVAR
module load annovar
cd vcf/gvcf/
# check the databases installed 
ls -l $ANNOVAR_DATA/hg38
# annotate the SNP.vcf
table_annovar.pl \
--vcfinput \
SNP.recalibrated_99.9.vcf.gz \
$ANNOVAR_DATA/hg38 \
--tempdir $WORKDIR/tmp \
--thread 4 \
--buildver hg38 \
--outfile SNP.recalibrated_99.9_annovar \
--remove \
--protocol refGene,knownGene,ensGene,gnomad_genome,ALL.sites.2015_08,avsnp147,clinvar_20210501,cosmic70,dbnsfp30a,dbscsnv11,exac03,intervar_20180118,ALL.sites.2014_10,ljb26_all \
--operation g,g,g,f,f,f,f,f,f,f,f,f,f,f \
--nastring '.' \
--polish

# annotate the INDEL.vcf
table_annovar.pl \
--vcfinput \
indel.SNP.recalibrated_99.9.vcf.gz \
$ANNOVAR_DATA/hg38 \
--tempdir $WORKDIR/tmp \
--thread 4 \
--buildver hg38 \
--outfile indel.SNP.recalibrated_99.9_annovar \
--remove \
--protocol refGene,knownGene,ensGene,gnomad_genome,ALL.sites.2015_08,avsnp147,clinvar_20210501,cosmic70,dbnsfp30a,dbscsnv11,exac03,intervar_20180118,ALL.sites.2014_10,ljb26_all \
--operation g,g,g,f,f,f,f,f,f,f,f,f,f,f \
--nastring '.' \
--polish