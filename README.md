# Oligo analysis

Documents and scripts for generating Figures for Filges *et al.*, "Digital Quantification of Chemical Oligonucleotide Synthesis Errors".

Now published in [Clinical Chemistry](https://doi.org/10.1093/clinchem/hvab136).

Processed data for reproducing the figures is part of this repository. The associated
raw sequencing data (fastq) is available at the NCBI Sequence Read Archive (SRA)
under accession [PRJNA507366](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA727098/).

## Backwards compatibility with UMIErrorCorrect

Code used in the Clinical Chemistry paper (master branch) was adapted to raw sequencing data processed with Debarcer. The code is not compatible with data processed with UMIErrorCorrect. 

An updated version is being developed (devel branch), but currently is not  stable.

# Requirements

## Preprocessing

### Installtions
Raw sequencing data (fastq) must be processed using software that generates
error corrected consensus reads and pileup files, such as UMIErrorCorrect.

In order to process data with UMIErrorCorrect, install the following:

- Burrows-Wheeler-Aligner (BWA), used by the pipeline for genome alignment using any one of the following commands

```
# If using conda:
conda install bwa

# if using apt:
sudo apt-get install bwa

# if using brew:
brew install bwa
```

- UMIErrorCorrect, the main analysis pipeline

```
pip install umierrorcorrect
```

### Reference fasta
If the synthetic oligonucleotide contains mixed genomic or alien sequences a suitable
reference genome is required. If the target sequences match a particular genomic region
exactly, e.g. a human reference genome (such as hg38), that may be used instead. However,
a custom reference will speed up analysis since the only a small region needs to be loaded
into memory. 

Thus, if needed, create a reference genome as a fasta file. Each fasta entry 
contains at least two lines, a header beginning with ">" and a sequence. 
Permissible file endings are '.fa' or '.fasta'.

```
>sequence_1
ACTGA
>sequence_2
ATATATA
```

Note that, if two sequences ("chromosomes") are very similar alignment can occur to either
sequence. Thus, it may be advantageous to align a certain set up samples to fasta file #1
containing only sequence_1 and another set of samples to fasta file #2 containing sequence_2,
to prevent a read that could match either reference to be aligned to the wrong sequence.

In order for BWA to process the reference it needs to be indexed using bwa index.

```
bwa index reference.fa
```

This will generate a number of additional files.

### BED file

Also create a suitable bed file for UMIErrorCorrect containing the locations

```
sequence_1       1       100     Insert_1
sequence_2       20       80     Insert_2
```

Where Insert_1 and Insert_2 will be the names shown in the UMIErrorCorrect output
as the assay names. sequence_1 and sequence_2 must *exactly* match the header
name of the reference fasta, i.e. the text immediately after '>'.

### Wrapper script

In order to process multiple fastq files at once, trim unnecessary parts
from the fastq name as well as generating sequencing quality control reports, 
a wrapper script "umi_pipeline.sh" is provided in this repository in the
"src" folder.

In your computing environment with bwa and umierrorcorrect installed, navigate
to the folder containing the fastq files. The following command assumes the wrapper
script has been placed in the same folder as the fastq files.

Use the following command, specifying the input directory ".", the complete path to the
bed file "-b" and reference fasta "-r". Also set the flags "-f" and "-c". This
will not generate any QC reports.

```
./umipipeline.sh -i . -b ../assay_regions.bed -r ../reference.fa -f -c
```

For details regarding usage see:
```
./umipipeline.sh -h
```

#### QC reports

Optionally install [fastp](https://github.com/OpenGene/fastp) from conda:

```
conda install -c bioconda fastp
```

If not already installed, also download multiqc for collating multiple QC files 
into one

```
pip install multiqc
```

and run the pipeline *without* the "-f" flag: 

```
./umipipeline.sh -i . -b ../assay_regions.bed -r ../reference.fa -c
```


## Required packages 

The following packages are required during the course of the analysis. 

The first code block of the R Markdown script should install these  automatically. If this does not work for one or more package, install them using their normal installation repositories.

## Installation

Clone the repository or go to "Code" -> "Download ZIP". Unpack the archive
and place the folder anywhere on your system.

The main script is "oligo_analysis.Rmd". In the first code blocks of the
main script the user must define the paths to all other relevant directories so
that the script can access auxilliary functions and input data. At this stage
the user also set some global analysis parameters, such as the consensus
threshold.

```{R}
# Should individual files with VAF be printed?
print.individual.files = FALSE

# Select consensus threshold to use for analysis
consensus_cutoff <- 10

# Create table of sample directory paths
# First define the directory containing UMIErrrorCorrect output files
sample_dir = '/Users/stefan/Documents/GitHub/OligoAnalysis/data/umierrorcorrect'

# Load helper functions, change location if needed
source('/Users/stefan/Documents/GitHub/OligoAnalysis/src/helper_functions_oligo.R')
```

Next, the script checks if all required R packages are present and installs them
if needed before loading all packages used in later steps.

```{R}
# CRAN packages to be installed with checkPackages
packages <- c('tidyverse', 'gridExtra', 'RColorBrewer', 'stringr', 'tis')

# Check and if needed install CRAN packages, this function is imported via helper_functions_oligo.R
checkPackages(requiredPackages = packages)

# If needed, install BiocManager for downloading packages from Bioconductor
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# If needed, install Rsamtools from Bioconductor
if (!require("Rsamtools", quietly = TRUE)){
  BiocManager::install("Rsamtools")
}

# Load packages
library(tidyverse)
library(Rsamtools)
library(EnvStats)
library(Hmisc)
library(ggpmisc)
```



