# Trimmomatic on Galaxy

## Introduction

After checking your input sequence reads for quality (e.g. using FastQC) it might be necessary to trim the reads. Here we will use the Trimmomatic tool.

## Learning Objectives

At the end of this tutorial you should be able to:

1. input sequence reads to Trimmomatic
2. trim using appropriate parameters, and
3. examine the output trimmed reads.

## Start

- open your Galaxy instance
- find your quality-checked Illumina sequence reads
- e.g. <fn>mutant_R1.fastq</fn> and <fn>mutant_R2.fastq</fn>
- We want to trim the parts of the reads that are of low quality
- based on the FastQC results, say we want to trim the reads like this:
(FIXME: expand detail)
    - trim Illumina adapters
    - leading and trailing bases - trim if quality is below 15
    - sliding window - trim once average quality is below 20

## Run Trimmomatic

FIXME: change these settings if required (examine FastQC reports)

- Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: QC and manipulation &rarr; Trimmomatic</ss>.
- <ss>Input FASTQ file (R1/first of pair)</ss>: <fn>mutant_R1.fastq</fn>
- <ss>Input FASTQ file (R2/second of pair)</ss>: <fn>mutant_R2.fastq</fn>
- <ss>Perform initial ILLUMINACLIP step </ss>: *Yes*
- <ss>Adapter sequences to use</ss>: FIXME
- <ss>How accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment</ss>: 40
- <ss>How accurate the match between any adapter etc. sequence must be against a read</ss>: 15
- leave the first <ss>Trimmomatic Operation</ss> as is
- click on <ss>+ Insert Trimmomatic Operation</ss>
- <ss>Select Trimmomatic operation to perform</ss>: *Cut bases off the start of a read, if below a threshold quality (LEADING)*
- <ss>Minimum quality required to keep a base</ss>: 15
- click on <ss>+ Insert Trimmomatic Operation</ss>
- <ss>Select Trimmomatic operation to perform</ss>: *Cut bases off the end of a read, if below a threshold quality (TRAILING)*
- <ss>Minimum quality required to keep a base</ss>: 15
- click <ss>Execute</ss>

FIXME: screenshot of these trimmomatic options selected

## Examine output

There are four output files, still in FASTQ format:

- R1 reads that have a pair in the R2 file
- R2 reads that have a pair in the R1 file
- R1 reads with no pair (R2 match was low quality: deleted)
- R2 reads with no pair (R1 match was low quality: deleted)

Examine each file with the eye icon.

- Note that some reads are now shorter.
- Note that the files are different sizes.

## What next?

Next: Assembly

Link to the Trimmomatic [webpage](http://www.usadellab.org/cms/index.php?page=trimmomatic) and the [manual.](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
