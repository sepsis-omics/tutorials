# Trimming reads

## Introduction

After checking your input sequence reads for quality (e.g. using FastQC) it might be necessary to trim the reads. Here we will use the Trimmomatic tool. For more inforamtion, see the Trimmomatic [webpage](http://www.usadellab.org/cms/index.php?page=trimmomatic) and the [manual.](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

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
- based on the FastQC results, we might want to trim the reads like this:
    - trim Illumina adapters
    - leading and trailing bases - trim if quality is below 15
    - sliding window - trim once average quality is below 20

## Trimmomatic functions

- Adapter trimming
    - This function trims adapters, barcodes and other contaminants from the reads.
    - You need to supply a FASTA file of possible adapter sequences, barcodes etc to trim. See Trimmomatic website for detailed instructions.
    - The default quality settings are sensible.
    - This should always be the first trimming step if it is used.

- Sliding window trimming
    - This function uses a sliding window to measure average quality and trims accordingly.
    - The default quality parameters are sensible for this step.

- Trailing bases quality trimming
    - This function trims bases from the end of a read if they drop below a quality threshold. e.g. If base 69 of 75 drops below the threshold, the read is cut to 68 bases.
    - Use FastQC report to decide whether this step is warranted and what quality value to use. A quality threshold value of 10-15 is a good starting point.

- Leading bases quality trimming
    - This function works in a similar fashion to trailing bases trimming except it performs it at the start of the reads.
    - Use FastQC report to determine if this step is warranted. If the quality of bases is poor at the beginning of reads it might be necessary.

- Minimum read length
    - Once all trimming steps are complete, this function makes sure that the reads are still longer than this value. If not, the read is removed from the file and its pair is put into the orphan file.
    - The most appropriate value for this parameter will depend on the FastQC report, specifically the length of the high quality section of the Per Base Sequence Quality graph.

- Each read library should be trimmed separately with parameters dependent on their own FastQC reports.

## Run Trimmomatic

<!---
FIXME: change these settings if required (examine FastQC reports)
--->

- Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: QC and manipulation &rarr; Trimmomatic</ss>.
- <ss>Input FASTQ file (R1/first of pair)</ss>: <fn>mutant_R1.fastq</fn>
- <ss>Input FASTQ file (R2/second of pair)</ss>: <fn>mutant_R2.fastq</fn>
- <ss>Perform initial ILLUMINACLIP step </ss>: *Yes*
<!---
- <ss>Adapter sequences to use</ss>: FIXME
--->
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

<!---
FIXME: screenshot of these trimmomatic options selected
--->

## Examine output

Trimmomatic should produce 2 pairs files (1 left and 1 right hand end) and 1 or 2 single “orphaned reads” files.

The output files are the ones you should use for assembly.

There are four output files, still in FASTQ format:

- R1 reads that have a pair in the R2 file
- R2 reads that have a pair in the R1 file
- R1 reads with no pair (R2 match was low quality: deleted)
- R2 reads with no pair (R1 match was low quality: deleted)

Examine each file with the eye icon. Look for:

- Number of reads orphaned by the trimming / cleanup process.
- Number of pairs lost totally.

## What next?

Next: use the output FASTQ files for Assembly, e.g. with [Spades](../spades/index.md)
