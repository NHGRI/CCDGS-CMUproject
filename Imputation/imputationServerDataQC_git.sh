#! /bin/bash

# mode and modules
sinteractive --mem=32g 

# 1. directories and paths
IMPDIR=/data/Hanserv/NMGTH3A/Imputation
cd $IMPDIR

# input data were downloaded from MIS and TOPMed imputation server and saved at
$IMPDIR/output/TOPMed
$IMPDIR/output/MIS

# 4. QC after imputation
# Explore the data and set the cut-off points for R2

# 4.1 stats for all of the chr vcf files
# raw data for imputation
awk -F '\t' '{print $1}' RAW_DATA.bim | sort | uniq -c | sort -nr

# data before imputation
for i in $(seq 1 23); do
bcftools stats $IMPDIR/data/${name}-updated-chr${i}_clean.vcf.gz >> $IMPDIR/output/${name}-updated.stats;
done 
grep -e "RAW_DATA-updated-chr" -e "number of records" RAW_DATA-updated.stats

# data imputed by TOPMed server
for i in $(seq 1 22); do
bcftools stats $IMPDIR/output/TOPMed/chr${i}.dose.vcf.gz >> $IMPDIR/output/${name}_TOPMed.stats;
done 
bcftools stats $IMPDIR/output/TOPMed/chrX.dose.vcf.gz >> $IMPDIR/output/${name}_TOPMed.stats
grep -e "dose.vcf.gz" -e "number of records" RAW_DATA_TOPMed.stats

# data imputed by MIS
for i in $(seq 1 22); do
bcftools stats $IMPDIR/output/MIS/chr${i}.dose.vcf.gz >> $IMPDIR/output/${name}_MIS.stats;
done 
bcftools stats $IMPDIR/output/MIS/chrX.dose.vcf.gz >> $IMPDIR/output/${name}_MIS.stats
grep -e "dose.vcf.gz" -e "number of records" RAW_DATA_MIS.stats

# 4.2 variants count on each chromsome by MAF cutoffs
zcat $IMPDIR/output/MIS/chr*.info.gz > $IMPDIR/output/MIS.info.gz
zcat $IMPDIR/output/TOPMed/chr*.info.gz > $IMPDIR/output/TOPMed.info.gz
## above files are too large for R to process and plot

### MIS
#### <=0.01
for i in $(seq 1 22); do 
echo chr${i} >> MIS.chr.count.01.txt
awk '$5<=0.01' $IMPDIR/output/MIS/chr${i}.info.gz | wc -l >> MIS.chr.count.01.txt
done
echo chrX >> MIS.chr.count.01.txt
awk '$5<=0.01' $IMPDIR/output/MIS/chrX.info.gz | wc -l >> MIS.chr.count.01.txt

#### 0.01 - 0.05
for i in $(seq 1 22); do 
echo chr${i} >> MIS.chr.count.05.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | wc -l >> MIS.chr.count.05.txt
done
echo chrX >> MIS.chr.count.05.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | wc -l >> MIS.chr.count.05.txt

#### >0.05
for i in $(seq 1 22); do 
echo chr${i} >> MIS.chr.count.common.txt
awk '$5>0.05' $IMPDIR/output/MIS/chr${i}.info.gz | wc -l >> MIS.chr.count.common.txt
done
echo chrX >> MIS.chr.count.common.txt
awk '$5>0.05' $IMPDIR/output/MIS/chrX.info.gz | wc -l >> MIS.chr.count.common.txt

### TOPMed
#### <=0.01
for i in $(seq 1 22); do 
echo chr${i} >> TOPMed.chr.count.01.txt
awk '$5<=0.01' $IMPDIR/output/TOPMed/chr${i}.info.gz | wc -l >> TOPMed.chr.count.01.txt
done
echo chrX >> TOPMed.chr.count.01.txt
awk '$5<=0.01' $IMPDIR/output/TOPMed/chrX.info.gz | wc -l >> TOPMed.chr.count.01.txt

#### 0.01 - 0.05
for i in $(seq 1 22); do 
echo chr${i} >> TOPMed.chr.count.05.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | wc -l >> TOPMed.chr.count.05.txt
done
echo chrX >> TOPMed.chr.count.05.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | wc -l >> TOPMed.chr.count.05.txt

#### >0.05
for i in $(seq 1 22); do 
echo chr${i} >> TOPMed.chr.count.common.txt
awk '$5>0.05' $IMPDIR/output/TOPMed/chr${i}.info.gz | wc -l >> TOPMed.chr.count.common.txt
done
echo chrX >> TOPMed.chr.count.common.txt
awk '$5>0.05' $IMPDIR/output/TOPMed/chrX.info.gz | wc -l >> TOPMed.chr.count.common.txt

# 4.3 variants count on each chromsome by Rsq cutoffs
## MIS
for i in $(seq 1 22); do 
echo chr${i} >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0 && $7<=0.1) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> MIS.chr.count.rsq.txt
done 
echo chrX >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0 && $7<=0.1) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> MIS.chr.count.rsq.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> MIS.chr.count.rsq.txt

## TOPMed
for i in $(seq 1 22); do 
echo chr${i} >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0 && $7<=0.1) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
done 
echo chrX >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0 && $7<=0.1) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> TOPMed.chr.count.rsq.txt

# 4.4 check median Rsq by MAF cutoff values on each chromosome
### MIS
#### <=0.01
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/MIS.chr.median.01.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($5<=0.01) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.01.txt
done
echo chrX >> $IMPDIR/output/MIS.chr.median.01.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($5<=0.01) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.01.txt

#### 0.01 - 0.05
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/MIS.chr.median.05.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.05.txt
done
echo chrX >> $IMPDIR/output/MIS.chr.median.05.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.05.txt

#### >0.05
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/MIS.chr.median.common.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($5>0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.common.txt
done
echo chrX >> $IMPDIR/output/MIS.chr.median.common.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($5>0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/MIS.chr.median.common.txt

### TOPMed
#### <=0.01
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/TOPMed.chr.median.01.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($5<=0.01) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.01.txt
done
echo chrX >> $IMPDIR/output/TOPMed.chr.median.01.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($5<=0.01) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.01.txt

#### 0.01 - 0.05
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/TOPMed.chr.median.05.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.05.txt
done
echo chrX >> $IMPDIR/output/TOPMed.chr.median.05.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($5>0.01 && $5<=0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.05.txt

#### >0.05
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/TOPMed.chr.median.common.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($5>0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.common.txt
done
echo chrX >> $IMPDIR/output/TOPMed.chr.median.common.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($5>0.05) print $0}' | sort -k7,7 | awk -f ../scripts/median.awk >> $IMPDIR/output/TOPMed.chr.median.common.txt

# 4.5 count the variants by chromosome at MAF>=0.05 and Rsq>=0.3
### MIS
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/MIS.chr.count.maf05.rsq3.txt
zcat $IMPDIR/output/MIS/chr${i}.info.gz | awk '{if ($5>=0.05 && $7>=0.3) print $0}' | wc -l >> $IMPDIR/output/MIS.chr.count.maf05.rsq3.txt
done
echo chrX >> $IMPDIR/output/MIS.chr.count.maf05.rsq3.txt
zcat $IMPDIR/output/MIS/chrX.info.gz | awk '{if ($5>=0.05 && $7>=0.3) print $0}' | wc -l >> $IMPDIR/output/MIS.chr.count.maf05.rsq3.txt

### TOPMed
for i in $(seq 1 22); do 
echo chr${i} >> $IMPDIR/output/TOPMed.chr.count.maf05.rsq3.txt
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | awk '{if ($5>=0.05 && $7>=0.3) print $0}' | wc -l >> $IMPDIR/output/TOPMed.chr.count.maf05.rsq3.txt
done
echo chrX >> $IMPDIR/output/TOPMed.chr.count.maf05.rsq3.txt
zcat $IMPDIR/output/TOPMed/chrX.info.gz | awk '{if ($5>=0.05 && $7>=0.3) print $0}' | wc -l >> $IMPDIR/output/TOPMed.chr.count.maf05.rsq3.txt

# 5. use target sequencing SNPs as gold standard to benchmark the imputed vcf files
# 5.1 get the target sequencing locus file
## target sequencing gvcf file: 
GVCF=/data/Hanserv/NMAMP/vcf/gvcf/indel.SNP.recalibrated_99.9_annovar.hg38_multianno.vcf
## check the format
bcftools view indel.SNP.recalibrated_99.9_annovar.hg38_multianno.vcf | grep -v "##" | head
bcftools view indel.SNP.recalibrated_99.9_annovar.hg38_multianno.vcf | grep -v "##" | wc -l
# 147753
## get the target sequencing positions from the gvcf file
bcftools view $GVCF | grep -v "#" | awk '{print $1 ":" $2}' | sed 's/chr//g' >> $IMPDIR/output/targetSequencingLoc.txt
bcftools view $GVCF | grep -v "#" | awk '{print $1 ":" $2}' >> $IMPDIR/output/targetSequencingLoc_hg38.txt

targetLoc38=$IMPDIR/output/targetSequencingLoc_hg38.txt

# liftover from hg38 to hg19 at https://genome.ucsc.edu/cgi-bin/hgLiftOver
targetLoc19=$IMPDIR/output/targetSequencingLoc_hg19.txt

# 5.2 format the info.gz file with a chr:pos column to prepare for filtration
## MIS
for i in $(seq 1 22); do 
zcat $IMPDIR/output/MIS/chr${i}.info.gz | tail -n +2 | cut -d " " -f 1 \
| awk -F '[:/]' 'BEGIN{printf "CHROM\tPOS\tREF\tALT\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n"} {printf ("%s\t %s\t %s\t %s\n" ,$1, $2, $3, $4)}' \
| awk '{print $1 ":" $2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' >> $IMPDIR/output/MIS/chr${i}.info.format.txt
done

## TOPMed
for i in $(seq 1 22); do 
zcat $IMPDIR/output/TOPMed/chr${i}.info.gz | tail -n +2 | cut -d " " -f 1 \
| awk -F '[:/]' 'BEGIN{printf "CHROM\tPOS\tREF\tALT\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1\n"} {printf ("%s\t %s\t %s\t %s\n" ,$1, $2, $3, $4)}' \
| awk '{print $1 ":" $2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' >> $IMPDIR/output/TOPMed/chr${i}.info.format.txt
done

# 5.3 select the SNPs that also in the target sequencing locus
## MIS
for i in $(seq 1 22); do 
awk 'FNR==NR {a[$1]; next}; $1 in a' $targetLoc19 $IMPDIR/output/MIS/chr${i}.info.format.txt >> $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt
done

## TOPMed
for i in $(seq 1 22); do 
awk 'FNR==NR {a[$1]; next}; $1 in a' $targetLoc38 $IMPDIR/output/TOPMed/chr${i}.info.format.txt >> $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt
done

# 5.4 variants count on each chromsome by Rsq cutoffs
## MIS
for i in $(seq 1 22); do 
echo chr${i} >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.0 && $7<=0.1) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/MIS/chr${i}.info.targetLoc.txt | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> MIS.chr.count.targetLoc.rsq.txt
done

## TOPMed
for i in $(seq 1 22); do 
echo chr${i} >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.0 && $7<=0.1) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.1 && $7<=0.2) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.2 && $7<=0.3) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.3 && $7<=0.4) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.4 && $7<=0.5) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.5 && $7<=0.6) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.6 && $7<=0.7) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.7 && $7<=0.8) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.8 && $7<=0.9) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
cat $IMPDIR/output/TOPMed/chr${i}.info.targetLoc.txt | awk '{if ($7>0.9 && $7<=1.0) print $0}' | wc -l >> TOPMed.chr.count.targetLoc.rsq.txt
done
