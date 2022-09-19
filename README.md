# Oligo analysis

Documents and scripts for generating Figures for Filges \emph{et al.}, "Digital Quantification of Chemical Oligonucleotide Synthesis Errors".

Now published in [Clinical Chemistry](https://doi.org/10.1093/clinchem/hvab136).

## Backwards compatibility with UMIErrorCorrect

Code used in the Clinical Chemistry paper (master branch) was adapted to raw sequencing data processed with Debarcer. The code is not compatible with data processed with UMIErrorCorrect. 

An updated version is being developed (devel branch), but currently is not  stable.

# Requirements

## Preprocessing

Raw sequencing data (fastq) must be processed using software that generates
error corrected consensus reads and pileup files, such as UMIErrorCorrect.

In order to process data with UMIErrorCorrect, install the following

- Burrows-Wheeler-Aligner (BWA), used by the pipeline for genome alignment
- UMIErrorCorrect, the main analysis pipeline

If the synthetic oligonucleotide contains mixed genomic or alien sequences a suitable
reference genome is required. If the target sequences match a particular genomic region
exactly, e.g. a human reference genome (such as hg38) may be used instead. However,
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


## Required packages 

The following packages are required during the course of the analysis. 

The first code block of the R Markdown script should install these  automatically. If this does not work for one or more package, install them using their normal installation repositories.

