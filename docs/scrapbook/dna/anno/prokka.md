# Prokka on Galaxy

## Background
Genome annotation involves finding and describing particular features, such as genes, tRNAs and rRNAs. This tutorial will demonstrate how to annotate an assembled bacterial genome using the tool Prokka. Link to [Prokka on github](https://github.com/tseemann/prokka#prokka-rapid-prokaryotic-genome-annotation); link to [Prokka citation](http://bioinformatics.oxfordjournals.org/content/30/14/2068.full).

## Learning objectives
At the end of this tutorial you should be able to :

1. input files into Prokka
2. change settings
3. run Prokka, and
4. examine the output: annotated genome.


## Pre-requisites
- a mGVL and galaxy instance

## Start
- open your galaxy instance in your mGVL

## Input data

- assembled contigs, e.g. <fn>SPAdes_contigs.fasta</fn>

## How it works
- Prokka compares the input contigs with various databases to identify coding sequences, rRNA genes, tRNA genes, non-coding RNA, and signal leader peptides.

- These databases are maintained by different organisations, and include information about known genomic features and their locations. Prokka includes a local copy.

## Run Prokka

- In Galaxy, go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: Annotation &rarr; Prokka</ss>  
- Set the following parameters (leave everything else unchanged):
    - <ss>Contigs to annotate</ss>: <fn>SPAdes contigs (fasta)</fn>  
    - <ss>Locus tag prefix (--locustag)</ss>: P
    - FIXME: actually we want to have a different locus tag for each sample for later use in Roary. Is there some problem here with the length of the locus tag that can be used?
    - <ss>Force GenBank/ENA/DDJB compliance (--compliant)</ss>: *Yes*
    - <ss>Sequencing Centre ID (--centre)</ss>: V
    - <ss>Genus Name</ss>: *Staphylococcus*  
    - <ss>Species Name</ss>: *aureus*  
    - <ss>Use genus-specific BLAST database</ss> *No*  
    - Click <ss>Execute</ss>. This may take x minutes.

## Examine the output

Once Prokka has finished, examine each of its output files.

- The gff and gbk files contains all of the information about all of the features annotated (in different formats.)

- `summary.gff`: a list of all the features found, listed in order of their location (starting at the start of contig number 1). Each row is a genomic feature and its location. Column 2 is the source - the database used to find the feature. Column 3 is the feature - e.g. CDS, tRNA.
- `summary.gbk`: the contigs listed in order. For each contig, the features are listed (e.g. CDS name and translation), followed by the sequence of the whole contig.

- The txt file contains a summary of the number of features annotated.
- The faa file contains the protein sequences of the genes annotated.
- The ffn file contains the nucleotide sequences of the genes annotated.

- Download the gff file to your local computer: click on the file name with the .gff extension, and then click on the disk icon ![disk icon](./images/image00.png).

![galaxy file](./images/image01.png)


## Annotated features
Now that we have annotated the draft genome sequence, we would like to view the sequence in the Artemis genome viewer.

- Open Artemis and load the downloaded .gff file.
- The top panel shows an overview - here we can see annotated genes and other features.
- The middle panel shows the DNA sequence and amino acid translations in 6 frames.
- The bottom panel shows a text summary of the features.
- Scroll left and right with the horizontal bars under each panel.
- Zoom with the vertical bars to the right.

![Artemis screenshot](./images/image02.png)

## What Next?
- Determine core and pan genomes using [Roary](../pan/roary.md).
