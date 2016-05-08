# Snippy

cmdline

This tutorial will teach you how to find variants in a genome compared to a reference genome.
Input: raw Illumina paired-end sequence reads from a bacteria.

one sample: one ref: call the variants


[Github link](https://github.com/tseemann/snippy/blob/master/README.md#snippy)

**Pre-requisites**
- connect to your GVL (commandline)
- background knowledge: variant calling

**Start**
- log in to virtual machine via terminal
- type in "snippy" to check it is there

**Input**

required data:
sequence reads - eg illumina data R1 and R2 reads (fasta or fastq)
same sp. reference genome (fasta or gbk)

in this tutorial:
we will get the data from [here]

**How it works:**

uses BWA first to align
then finds variants using freebayes

**Run snippy**

need to enter in:
- number of cpus [or it uses a default]
- output directory (write a folder name and it will be created)
- ref genome filename
- r1 filename
- r2 filename

(TBC)
