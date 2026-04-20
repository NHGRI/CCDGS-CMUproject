Identifying Correlated Methylation Units (CMUs)
================
Aarti Jajoo, Owen Hirchi & Neil Hanchard
2026-04-10

## Overview

Correlated Methylation Units (CMUs) are genomic regions defined by
coordinated DNA methylation patterns across neighboring CpG sites. For a
given 450k or EPIC methylation profile across samples, CMUs are
identified within non-overlapping 250,000 base pair windows across
autosomal chromosomes.

This vignette demonstrates how to:

0.  [Setup R environment, dependencies, and load CMU
    functions](#step-0-setup-environment-dependencies-and-load-cmu-functions)
1.  [Prepare Illumina 450k or EPIC DNA methylation
    data](#step-1-methylation-input-file-formating)
2.  [Generate indexing for efficient CMU
    computation](#step-2-generating-cmu-folder-and-indexing-for-beta-file)
3.  [Identify CMUs](#step-3-running-cmu-detection-on-illumina-450k-data)
4.  [Incorporate sample subsets and
    covariates](#step-4-advanced-options-sample-subsets-and-covariates)
5.  [Perform differential CMUs](#step-5-differential-cmus)
6.  [How to use CMU output files](#how-to-use-cmu-output-files)

The provided toy datasets can help illustrate the correct input format
for the dataset to be used.

------------------------------------------------------------------------

## Step 0: Setup environment, dependencies, and load CMU functions

The CMU pipeline is implemented primarily in R, with Python dependency
used for indexing methylation data.

``` r
library(reticulate)
library(matrixStats)
library(spatstat)
library(imager)
#Python dependency installation:
reticulate::py_install("pandas")
```

> **Note for macOS users:** The `imager` package may require XQuartz. If
> pacakge installation fails, install XQuartz from
> <https://www.xquartz.org/> and restart R.
> <p style="color:red; font-weight:bold;">
>
> 🛑 Do not use ~ (tilde) in the filepath i.e. instead of `~/Data` use
> `/Users/ajajoo/Data`. Otherwise, python reticulate related commands
> give error.
> </p>
>
> ------------------------------------------------------------------------

#### Change working directory & Load CMU functions

Set the working directory to the location where the CMU repository has
been downloaded.

``` r
setwd("path/to/CMU")
```

``` r
source("GetIndexing.R")
source("CUsAcrossGenomeVer05ForPyIndex.R")
```

------------------------------------------------------------------------

## Step 1: Methylation input file formating

CMU analysis requires a **tab-separated text file** containing DNA
methylation beta values with the following strict structure. \####
Example (toy data)

Example of a beta file:

``` text
toyData/EPIC_Mimic_Text.tsv
```

``` text
"Sample_1"    "Sample_2"   
"cg26928153"  0.331238528247923   0.869199158390984
"cg16269199"  0.0595261559366040  0.934328192146495
"cg13869341"  0.354899412253872   0.249734462632203
"cg24669183"  0.912160335108638   0.771520412294194
```

Notes:

- **CpG IDs are quoted**
- **CpG column does not have a column name in header**
- **Sample names are quoted**
- **Values are numeric beta values**

Example of a probe file: File location (toy data):

``` text
toyData/probe_Epic.txt
```

``` text
"cg26928153"
"cg16269199"
"cg13869341"
"cg24669183"
```

Notes:

- **One CpG ID per line**
- **No header**
- **CpG IDs must be quoted**
- **Order must exactly match the methylation file**

------------------------------------------------------------------------

> ⚠️ **Input file formatting**  
> If you encounter errors related to file structure, quoting, or probe
> ordering,  
> see the detailed formatting requirements and ways to convert your data
> in the correct format section  
> [Click here for more details on Input data requirements and data
> formatting](#detailed-input-data-requirements-and-data-formatting).

------------------------------------------------------------------------

## Step 2: Generating CMU Folder and Indexing for beta file

``` r
# For EPIC 
GenerateDataFormat(
filename = "toyData/EPIC_Mimic_Text.tsv",
PathToResult = "CMUResults/EPIC/",
SortedProbeFilename = "./PythonIndexing/SortedCpGNamesEPIC.txt"
)
# For 450K 
# change SortedProbeFilename = "./PythonIndexing/SortedCpGNames450k.txt"
```

This step: - Generates indexing files required for quick reading of beta
file
<p style="color:red; font-weight:bold;">

🛑 DO NOT modify your original beta TSV file after this step.  
If you do, you must rerun Step 4.
</p>

------------------------------------------------------------------------

## Step 3: Running CMU detection on Illumina 450k data

``` r
CUsAcrossGenome(
  annotationFile = "Illumina_probe_annotation_Epic.txt",
  PathToFig = "CMUResults/EPIC/Results/",
  data = "toyData/EPIC_Mimic_Text.tsv",
  probenameFile = "toyData/probe_Epic.txt",
  indexFile = "CMUResults/EPIC/Data/Index.txt"
)
# For 450K 
# change annotationFile = "Illumina_probe_annotation_450k.txt"
```

This function identifies genome-wide CMUs.

------------------------------------------------------------------------

## Step 4: Advanced options: sample subsets and covariates

### Group-specific CMU computation

You can either prepare a new methylation input following Step 3, or CMUs
can also be computed on a subset of samples by providing a logical
vector (`GroupIndex`) matching the order of samples in the methylation
matrix.

``` r
GroupIndex <- SubjectNameInGroup %in% SubjectNameInBetaTable
CUsAcrossGenome(
    annotationFile = "Illumina_probe_annotation_Epic.txt",
    PathToFig = "CMUResults/EPIC/Results/",
    data = "toyData/EPIC_Mimic_Text.tsv",
    probenameFile = "toyData/probe_Epic.txt",
    indexFile = "CMUResults/EPIC/Data/Index.txt",
    GroupIndex = GroupIndex
  )
```

### Covariate adjustment

CMUs can be computed after residualizing a choice of user provided
covariates. The matrix must match the sample order of the input data
`toyData/EPIC_Mimic_Text.tsv` and include an `intercept` column (a
column with all 1s).

``` r
confounderX <- as.matrix(
  bb[, c("intercept", "Gender", "Age")]
)
CUsAcrossGenome(
    annotationFile = "Illumina_probe_annotation_Epic.txt",
    PathToFig = "CMUResults/EPIC/Results/",
    data = "toyData/EPIC_Mimic_Text.tsv",
    probenameFile = "toyData/probe_Epic.txt",
    indexFile = "CMUResults/EPIC/Data/Index.txt",
    confounderX = confounderX
  )
```

## Step 5: Differential CMUs

To identify differential CMUs, modify DiffCMU.R file. You will have to
modify the following variables.

``` r
# set the following variables 
data_file <- # methylation tsv file prepared  in Step 1 
index_file <- # corresponding index file for  methylation data prepared in Step 2
idx_ctrl_X  <- # index for controls/normal 
idx_case_X  <- # index for case/Disease subjects  
covX <- # covariates to be used, in same order as it is in data_file. idx_ctrl_X and idx_case_X would  # be used to subset this matrix for diseaes and controls. 
CMU_file_for_ctrl <- # path_to_CMU_ctrl_Output/ContM0.6.txt
CMU_file_for_case <- # path_to_CMU_case_Output/ContM0.6.txt
pvalFilenameForNormalCMU <- # Provide filename with path and extension Rdata
pvalFilenameForDiseaseCMU <- # Provide filename with path and extebsuib Rdata
```

## How to use CMU output files

CMU calling is performed with three correlation soft thresholds:
`0.4, 0.6 and 0.8`. We recommend `0.6` for general purposes and refer to
our paper for detailed understanding of how lower or higher parameters
can impact the CMU calling. For each parameter there are two files,
e.g. for 0.6 we have `ContM0.6.txt` for contiguous CMUs and
`NonContM0.6.txt` for non-contiguous CMUs.

``` r
# 
source('getClusDetails.R')
CMUs <- read.table('pathToCMU/ContM0.6.txt',stringsAsFactors = FALSE)
head(CMUs)


                 V1                 V2 V3                                                                                                             V4                                                                               V5
1 ../Y123/Amy/Results/ContM0.6.txt  1 1086777 cg24437834;cg15565004;cg19999567;cg19856606;cg12743165;cg04065236;cg04315086;cg09358422;cg24523322;cg12612065; 1022530;1022642;1022900;1022932;1022995;1023150;1023163;1023247;1023281;1023319;
2 ../Y123/Amy/Results/ContM0.6.txt  1 1086777                                                        cg21786289;cg13897241;cg13176867;cg15105996;cg02565062;                                              854966;855046;855425;855445;855544;
3 ../Y123/Amy/Results/ContM0.6.txt  1 1086777                                                        cg07390924;cg10110857;cg10375192;cg02056921;cg24076968;                                              901725;901755;901799;901892;902549;
4 ../Y123/Amy/Results/ContM0.6.txt  1 1086777                                                        cg20062691;cg03811829;cg11211792;cg04788999;cg16526047;                                              949392;949449;949634;949850;949893;
5 ../Y123/Amy/Results/ContM0.6.txt  1 1086777                                                                   cg22226438;cg08029603;cg22699361;cg03890988;                                                     854766;854824;854918;854951;
```

Each row in `CMUs` is a CMU. The second and third columns represent the
chromosome and genomic window, the fourth column is the list of CpGs
that are part of the CMU, and the fifth column is the genomic position
of CpGs based on hg19 coordinates. `getClusDetails.R` has many functions
that can summarize properties of CMUs; for example, `getNCpG` counts the
number of CpGs per CMU.

------------------------------------------------------------------------

------------------------------------------------------------------------

## Detailed Input data requirements and data formatting

#### Required structure (strict)

- **Rows:** CpG probe IDs  
- **Columns:** Samples  
- **Values:** Beta values (0–1)

#### Mandatory formatting rules

1.  **Tab-separated format**
    - The file must be tab-delimited (`.tsv`)
    - Other delimeters are not supported without conversion
2.  **First row has one less column**
    - The first column contains CpG IDs and therefore has **no column
      header**
    - So the first row only has sample names which is one column less
      than the remaining files
3.  **Quoted identifiers are REQUIRED**
    - CpG probe IDs **must be quoted**
    - Sample names **must be quoted**
    - Quoting must be consistent across all input files
4.  **CpG IDs are read as row names**
    - CpG probe IDs appear in the first column

------------------------------------------------------------------------

#### Recommended way to get beta file in correct tsv format via R

Load your beta file in R, rownames must be CpGIDs and column names must
be sample names.

``` r
# load your methylation file 
beta_matrix <- read.table(
  "toyData/EPIC_Mimic_Text.tsv",
  sep = "\t",
)
```

Use below command to generate tsv beta file.

``` r
write.table(
  beta_matrix,
  file = "pathtobeta/beta.tsv",
  sep = "\t"
)
```

#### Example of converting CSV to TSV (terminal)

If your methylation data are stored as a comma-separated file (`.csv`),
convert it to tab-separated format using:

``` bash
sed 's/,/\t/g' pathtobeta/beta.csv > pathtobeta/beta.tsv
```

Always verify formatting:

``` bash
head pathtobeta/beta.tsv
```

------------------------------------------------------------------------

### Probe ID file (`probe_Epic.txt`)

A probe ID file is required and must satisfy the following constraints.

#### Strict requirements

⚠️ **Critical constraint**

The probe ID file must match the **exact probe order** of the
methylation file:

- `probe_Epic.txt` corresponds to the **first column** of the
  methylation file
- The methylation file contains **one additional header row**.
  Therefore, `probe_Epic.txt` must contain **one fewer line** compared
  to beta file.

Example `probe_Epic.txt`:

------------------------------------------------------------------------

\### Generating `probe_Epic.txt` from the methylation file using
beta_matrix (recommended)

``` r
write.table(
rownames(beta_matrix),
file = "toyData/probe_test.txt",
quote = TRUE,
row.names = FALSE,
col.names = FALSE
)
```

This guarantees: - Exact probe ordering - Quoted CpG IDs - One-to-one
correspondence with the methylation matrix

------------------------------------------------------------------------

## Citation

If you use our CMU software, please cite:

> *[A genome-wide map of correlated DNA methylation developed from image
> clustering of multi-tissue methylation
> data](https://www.researchsquare.com/article/rs-2852818/v1)*. (2026).

## References of R packages used

*annotate*: Gentry J (2025). *annotate: Annotation for microarrays*.
<doi:10.18129/B9.bioc.annotate>
<https://doi.org/10.18129/B9.bioc.annotate>, R package version 1.86.1,
<https://bioconductor.org/packages/annotate>.

*GenomicRanges*: Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M,
Gentleman R, Morgan M, Carey V (2013). “Software for Computing and
Annotating Genomic Ranges.” *PLoS Computational Biology*, *9*.
<doi:10.1371/journal.pcbi.1003118>
<https://doi.org/10.1371/journal.pcbi.1003118>,
<http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.

*gg.gap*: Lou J, Zhang J, Lvy Y, Jin Z (2019). *gg.gap: Define Segments
in y-Axis for ‘ggplot2’*. <doi:10.32614/CRAN.package.gg.gap>
<https://doi.org/10.32614/CRAN.package.gg.gap>, R package version 1.3,
<https://CRAN.R-project.org/package=gg.gap>.

*ggplot2*: Wickham H (2016). *ggplot2: Elegant Graphics for Data
Analysis*. Springer-Verlag New York. ISBN 978-3-319-24277-4,
<https://ggplot2.tidyverse.org>.

*ggsignif*: Constantin A, Patil I (2021). “ggsignif: R Package for
Displaying Significance Brackets for ‘ggplot2’.” *PsyArxiv*.
<doi:10.31234/osf.io/7awm6> <https://doi.org/10.31234/osf.io/7awm6>,
<https://psyarxiv.com/7awm6>.

*karyoploteR*: Gel B, Serra E (2017). “karyoploteR : an R / Bioconductor
package to plot customizable genomes displaying arbitrary data.”
*Bioinformatics*, *33*(19), 3088-3090.
<doi:10.1093/bioinformatics/btx346>
<https://doi.org/10.1093/bioinformatics/btx346>.

*matrixStats*: Bengtsson H (2025). *matrixStats: Functions that Apply to
Rows and Columns of Matrices (and to Vectors)*.
<doi:10.32614/CRAN.package.matrixStats>
<https://doi.org/10.32614/CRAN.package.matrixStats>, R package version
1.5.0, <https://CRAN.R-project.org/package=matrixStats>.

*org.Hs.eg.db*: Carlson M (2025). *org.Hs.eg.db: Genome wide annotation
for Human*. R package version 3.21.0.

*psych*: William Revelle (2025). *psych: Procedures for Psychological,
Psychometric, and Personality Research*. Northwestern University,
Evanston, Illinois. R package version 2.5.6,
<https://CRAN.R-project.org/package=psych>.

*reticulate*: Ushey K, Allaire J, Tang Y (2025). *reticulate: Interface
to ‘Python’*. <doi:10.32614/CRAN.package.reticulate>
<https://doi.org/10.32614/CRAN.package.reticulate>, R package version
1.44.1, <https://CRAN.R-project.org/package=reticulate>.

*TxDb.Hsapiens.UCSC.hg19.knownGene*: Carlson M, Maintainer BP (2015).
*TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package for TxDb
object(s)*. R package version 3.2.2.

#### pandas

McKinney, W. (2010). Data Structures for Statistical Computing in
Python.
