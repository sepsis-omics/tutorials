# Snippy

cmdline

This tutorial will teach you how to find variants in a genome compared to a reference genome.

Variant calling: either map reads directly to a reference genome, or make a de novo assembly first and compare that to the reference genome.

FIXME: snippy does which one?

The former can be better because the depth of coverage at a position can inform the call.

Mapping: phred scores for confidence

Mapped file: A sequence alignment map (SAM) and binary version (BAM)

bayesian methods to call variants cf to ref




[Github link](https://github.com/tseemann/snippy/blob/master/README.md#snippy)

##Pre-requisites

- connect to your GVL (commandline)
- background knowledge: variant calling

##Start

- log in to virtual machine via terminal
- type in "snippy" to check it is there

##Input

- Illumina paired-end sequence reads from a bacteria (fasta or fastq)

- same sp. reference genome (fasta or gbk)

- FIXME: get data

- We will use one sample and one reference, and call the variants.

##How it works:

- uses BWA first to align (= mapping to ref?)

- then finds variants using freebayes

##Run snippy

need to enter in:
- number of cpus [or it uses a default]
- output directory (write a folder name and it will be created)
- ref genome filename
- r1 filename
- r2 filename

## Output

- 17 output files
- look at the VCF file (FIXME: add more here)

## Next

-
