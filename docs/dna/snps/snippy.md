# Snippy

cmdline

This tutorial will demonstrate how to find variants in a bacterial genome using Snippy. Variants are found by comparing to a reference genome of the same species.

[Github link to Snippy](https://github.com/tseemann/snippy/blob/master/README.md#snippy)

## Pre-requisites
- connect to your GVL - cmdline
- background knowledge: variant calling

## Start
- log in to your virtual machine via terminal
- navigate to the place where you want Snippy to run.
- make a folder called snippy - `mkdir snippy`
- move into that folder - `cd snippy`

## Input

### Raw sequence reads
- Illumina paired-end reads from a bacteria in FASTQ format.
- These reads are from *Pasteurella multocida*, from EMBL-EBI ENA. We will use `wget` to download them via FTP.
- in the snippy folder,
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257473/SRR1257473_1.fastq.gz
```
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/003/SRR1257473/SRR1257473_2.fastq.gz
```
- These files are compressed and so end in .gz. Snippy can use these .gz read files directly without unzipping.

### Reference genome
- Reference genome from the same species, *Pasteurella multocida*, from EMBL-EBI Ensembl genomes, in FASTA format.
- in the snippy folder,
```bash
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-31/fasta/bacteria_104_collection/pasteurella_multocida_subsp_multocida_gca_001027695/dna/Pasteurella_multocida_subsp_multocida_gca_001027695.ASM102769v1.31.dna.genome.fa.gz
```
- This file is also compressed into .gz format. This needs to be unzipped:
```bash
gunzip ftp://ftp.ensemblgenomes.org/pub/bacteria/release-31/fasta/bacteria_104_collection/pasteurella_multocida_subsp_multocida_gca_001027695/dna/Pasteurella_multocida_subsp_multocida_gca_001027695.ASM102769v1.31.dna.genome.fa.gz
```
- The file will now end in .fa (which is fasta format, and Snippy can use).

## How it works
- Reads are mapped to the reference genome using BWA
- this makes a BAM file, sent to Freebayes
- Freebayes finds differences between the reads and the reference, and calls the variants.

## Run Snippy
- cpus: choose number of cpus to use [or it uses a default] - here we will use 16
- outdir: choose a name for the output directory, where results will go - here we will use "mysnps"
- ref: the input reference genome filename
- R1: the input R1 reads filename
- R2: the input R2 reads filename
- to run snippy:
```bash
snippy --cpus 16 --outdir mysnps --ref [filename.fa] --R1 [R1.fastq.gz] --R2 [R2.fastq.gz]
```
## Output
- 17 output files
- list all the output files (that were put into the "mysnps" folder):
```bash
ls mysnps
```
- look at the first 10 lines of the snps.tab file
```bash
head -10 mysnps/snps.tab
```
- look at these columns: chromosome (CHROM), genomic position (POS), variant type (TYPE), nucleotide state in the ref (REF), nucleotide state in the input sample (ALT), and the frequency counts of REF and ALT (EVIDENCE).
- FIXME: screenshot with arrows
- FIXME: filter for quality?
- FIXME: load reference and the tabular vcf file into JBrowse/Artemis/IGV to view the genome and the snps.
- FIXME: is there anything we are looking for in particular? e.g. number of variants, existing known variants, variants in particular genes, AMR variants? 

## Advanced options
- FIXME: add more
- See if there are any unmapped reads
- See what they are (e.g. a plasmid)
```bash
snippy --outdir out --unmapped  etc
```
- then how to view the output fastq.gz

## Next

- ?
