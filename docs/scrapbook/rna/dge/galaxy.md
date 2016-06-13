# DGE using Galaxy

No idea - i think there are some workflows but we dont need Cufflinks etc

***

This tutorial demonstrates how to quantify differences in gene expression between two bacterial samples. We will use RNA-seq data from the samples, and then use several tools to test for differential gene expression. We can then examine the results in the program [Degust](./degust.md)

## Pre-requisites
- your own mGVL instance
- knowledge: Galaxy
- knowledge: Differential gene expression

## Start
- open your mGVL galaxy instance.

## Input data
- RNA-seq reads from the samples
    - Bacterial species E. coli from study SRP027344 from EBI ENA
    - We will use wildtype (called WTA and WTB) vs. condition 23 (called 23A and 23B)
    - [what was this condition? grown in different media?]
    - Illumina, single-end, 51-bp
- a reference genome
    - E coli: NCBI Reference Sequence: NC_000913.3 in fasta format
- FIXME: Get data from GenomeSpace/ or saved galaxy history

## Map transcripts
- We need to map the transcripts to a reference genome.
- Bacteria don't need splice-aware mapping (don't have introns)
- Galaxy: tools: <ss>NGS Analysis: NGS Mapping: Map with BWA-MEM</ss>
- <ss>Will you select a reference genome from your history or use a built-in index?: Use a genome from history and build index</ss>
- <ss>Use the following dataset as the reference sequence:</ss> E_coli_ref_genome
- <ss>Single or Paired-end reads: single </ss>
- <ss>Select FASTQ dataset</ss>:
    - click on the multiple files button (image) in centre
    - make sure all 4 FASTQ files are in there  
    - hold down shift to select them all (they turn blue)
    - this will map each set of transcripts to the ref genome, so there will be 4 output files

## Visualize the mapped reads
- The mapped reads are now as .bam files which can't be viewed by just clicking on them.
- [should we use JBrowse?]
- on local computer, install/open IGV
- Genomes->load genome from file->../igv/genomes/NC_000913.2 (Ecoli ref)
- then in galaxy, go to the mapped transcripts and click on IGV local
- in IGV, zoom in (top right)
- the bottom pane is the ref sequence / zoom in and out to see the reads.

## Count reads
- Generate read counts per gene
- For each transcriptome, count the number of transcripts per gene/feature
E.g. with HTSeq

data:

this is part of the Colombus database of bacterial expression
http://nar.oxfordjournals.org/content/suppl/2013/10/29/gkt1086.DC1/nar-02461-data-e-2013-File001.pdf

wildtype and two mutants
(all grown under minimal conditions or just mutants?)
hopefully some expression changes bn condition 23 and wt so that can be viewed in the Degust plots
