# Prokka

Cmdline and galaxy

Genome annotation involves finding and describing particular features, such as genes, tRNAs and rRNAs. This tutorial will demonstrate how to annotate an assembled bacterial genome using the tool Prokka.

[Prokka on github](https://github.com/tseemann/prokka#prokka-rapid-prokaryotic-genome-annotation)

[Prokka citation](http://bioinformatics.oxfordjournals.org/content/30/14/2068.full)

## Using Prokka in Galaxy

## Pre-requisites
- a mGVL and galaxy instance

## Start
- open your galaxy instance in your mGVL

## Input data

- one FASTA file containing the scaffolds in nucleotide format, e.g. the fasta output from Spades. e.g. in this tutorial: FIXME: data

- Go to galaxy - start a history - get his file into the history

## How it works

- Prokka compares the input contigs with various databases to identify coding sequences, rRNA genes, tRNA genes, non-coding RNA, and signal leader peptides.

- These databases are maintained by different organisations, and include information about known genomic features and their locations. (FIXME: does Prokka include a local copy or does it access the current dbases?)

## Run Prokka

- In Galaxy, go to **NGS Analysis: NGS Annotation: Prokka**.

- Under Contigs to <ss>annotate</ss> click on the single dataset icon (underneath the word "Contigs") and then select the file in the drop down menu (this will be your input .fna file that is in your history).

- Under **Locus tag prefix (--locustag)** type in "L". FIXME: actually we want to have a different locus tag for each sample for later use in Roary. Is there some problem here with the length of the locus tag that can be used?

- Under **Sequencing centre ID (--centre)** type in "C". FIXME: or something else?

- Leave the other settings as they are.

- Click **Execute**. This may take 10 minutes.

## Output

- There are 11 output files. Look at the output files by clicking on the eye icon under each file. For more information about output files, see [Table 2](http://bioinformatics.oxfordjournals.org/content/30/14/2068.full).

- `summary.gff`: a list of all the features found, listed in order of their location (starting at the start of contig number 1). Each row is a genomic feature and its location. Column 2 is the source - the database used to find the feature. Column 3 is the feature - e.g. CDS, tRNA.

- `summary.gbk`: the contigs listed in order. For each contig, the features are listed (e.g. CDS name and translation), followed by the sequence of the whole contig.

FIXME: best way to view / look at features


## Next

- Determine core and pan genomes. Roary can do this using a set of .gff files from different isolates as input.

# Using Prokka in commandline

- ssh to your mGVL
- wget data - scaffolds of .fna for several species
- (alternatively - have already put data on your mGVL)
- type in:
```
prokka --outdir [name of output folder for 1 sample] --locustag [tag eg sample number] [fna filename]
```

- make a new directory for the gff files e.g. "gff_files"
- mv -v */*.gff gff_files/  (means move any .gff files into that folder)

-FIXME: prokka options see manual
