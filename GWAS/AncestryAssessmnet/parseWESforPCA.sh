########## Produce PCA bi-plot for 1000 Genomes Phase III with WES data for COEH ########## 
#### Reference: https://www.biostars.org/p/335605/

# mode and modules
sinteractive --mem=32g --cpus-per-task=8 --time=24:00:00

################## This part 1 files only need to be done once and they were prepared. ##################
################## PLEASE DO NOT run the following script aganin in the same folder.   ##################

#### 1. Reference data
REFDIR=/data/Hanserv/Reference/1000Genomes
cd $REFDIR

# 1.1 Download the files as VCF.gz (and tab-indices)
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr";
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";

# for chr in {1..22}; do
#    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi;
# done

# 1.2 Download 1000 Genomes PED file
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped 

# 1.3 GRCh37 / hg19 reference genome
REFG1K=/data/Hanserv/Reference/1000Genomes/human_g1k_v37.fasta

# 1.4 Convert the 1000 Genomes files to BCF
module load bcftools

# for chr in {1..22}; do
#     bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
#       ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
#       bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
#         bcftools norm -Ob --rm-dup both \
#           > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;

#     bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;
# done

# 1.5 Convert the BCF files to PLINK format
# for chr in {1..22}; do
#     plink --noweb \
#       --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
#       --keep-allele-order \
#       --vcf-idspace-to _ \
#       --const-fid \
#       --allow-extra-chr 0 \
#       --split-x b37 no-fail \
#       --make-bed \
#       --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
# done

# 1.7 Prune variants from each chromosome
# This folder /data/Hanserv/Reference/1000Genomes/Pruned_maf10/ is similar pruning process but with the "--maf 0.10 --indep 50 5 1.1"
# mkdir Pruned ;

# for chr in {1..22}; do
#     plink --noweb \
#       --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
#       --maf 0.20 --indep 50 5 1.1 \
#       --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;

#     plink --noweb \
#       --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
#       --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.prune.in \
#       --make-bed \
#       --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
# done

# 1.8 Get a list of all PLINK files

find . -name "*.bim" | grep -e "Pruned" > ForMerge_maf20.list ;
sed -i 's/.bim//g' ForMerge_maf20.list ;

# 1.9 Merge all projects into a single PLINK file
plink --merge-list ForMerge_maf20.list --out Merged1KG_maf20 ;

################## This part 1 files only need to be done once and they were prepared. ##################
################## PLEASE DO NOT run the ABOVE script aganin in the same folder.   ##################

#### 2. Study data
# Preprocessed from WES BAM to g.vcf data by weaAnalysis.sh

WORKDIR=/data/Hanserv/yxhan/COEH
cd $WORKDIR

REFGEN=/data/Hanserv/Reference/hs37d5.fa

# 2.1 Prepare the PED data 

# Original file: /Users/hany4/Documents/Projects/Hypertension/COEH/COEH_Dictionary.xlsx PED sheet
# /data/Hanserv/yxhan/COEH/supporting_files/COEH_20220620.ped


# 2.2 Convert the vcf files to BCF
module load bcftools

for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f $REFGEN \
      $WORKDIR/vcf/gvcf/chr"${chr}".vcf.gz | \
      bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
          > $WORKDIR/vcf/gvcf/chr"${chr}".bcf ;

    bcftools index $WORKDIR/vcf/gvcf/chr"${chr}".bcf ;
done

# 2.3 Convert the BCF files to PLINK format
module load plink

for chr in {1..22}; do
    plink --noweb \
      --bcf $WORKDIR/vcf/gvcf/chr"${chr}".bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out $WORKDIR/PLINK_output/coeh.chr"${chr}" ;
done

# 2.4 Prune variants from each chromosome
mkdir $WORKDIR/PLINK_output/Pruned ;

for chr in {1..22}; do
    plink --noweb \
      --bfile $WORKDIR/PLINK_output/coeh.chr"${chr}" \
      --maf 0.20 --indep 50 5 1.1 \
      --out $WORKDIR/PLINK_output/Pruned/coeh.chr"${chr}" ;

    plink --noweb \
      --bfile $WORKDIR/PLINK_output/coeh.chr"${chr}" \
      --extract $WORKDIR/PLINK_output/Pruned/coeh.chr"${chr}".prune.in \
      --make-bed \
      --out $WORKDIR/PLINK_output/Pruned/coeh.chr"${chr}" ;
done

# 2.5 Get a list of all PLINK files

cd PLINK_output
find . -name "*.bim" | grep -e "Pruned" > ForMerge_COEH_maf20.list ;
sed -i 's/.bim//g' ForMerge_COEH_maf20.list ;

# 2.6 Merge all projects into a single PLINK file

plink --merge-list ForMerge_COEH_maf20.list --out MergedCOEH_maf20 ;


#### 3. Merge study data with 1000Genome data
# 3.1 find common variants between your dataset and the merged 1000 Genomes dataset (and filter both for these common variants)
REFDIR=/data/Hanserv/Reference/1000Genomes
COEHDIR=/data/Hanserv/yxhan/COEH/PLINK_output

name=MergedCOEH_maf20
refname=Merged1KG_maf20


# 3.1.1 Filter reference and study data for non A-T or G-C SNPs
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' \
$COEHDIR/$name.bim  > $name.ac_gt_snps

awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' \
$REFDIR/$refname.bim  > $refname.ac_gt_snps

# make bed
plink \
--bfile  $REFDIR/$refname \
--exclude $refname.ac_gt_snps \
--allow-extra-chr 0 \
--make-bed \
--threads 8 \
--out $refname.no_ac_gt_snps

plink \
--bfile $COEHDIR/$name \
--exclude $name.ac_gt_snps \
--make-bed \
--out $name.no_ac_gt_snps

# 3.1.2 # Filter reference data for the same SNP set as in study
plink \
--bfile  $REFDIR/$refname \
--extract $name.no_ac_gt_snps.bim \
--allow-extra-chr 0 \
--make-bed \
--out $refname.matched

# 3.2 merge the 1000 Genomes data with your own data

plink \
--bfile $name.no_ac_gt_snps \
--bmerge $refname.matched.bed $refname.matched.bim $refname.matched.fam \
--make-bed \
--out $name.$refname.matched.merged

# 3.3 (1.10) Perform PCA

# Pairwise comparisons of samples, relatedness calculated / Pairwise IBD estimation
plink \
--bfile $name.$refname.matched.merged \
--genome \
--threads 10 \
--out $name.$refname.matched.merged

# PCA
plink \
--bfile $name.$refname.matched.merged \
--read-genome $name.$refname.matched.merged.genome \
--cluster \
--mds-plot 6 \
--out $name.$refname.matched.merged

scp hany4@helix.nih.gov:/data/Hanserv/yxhan/COEH/PLINK_output/*.mds .