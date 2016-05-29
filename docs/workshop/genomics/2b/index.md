# Genome annotation using Prokka

## Background

In this section we will use a software tool called Prokka to annotate the draft genome sequence produced in Activity 2a. Prokka is a “wrapper”; it collects together several pieces of software (from various authors) - this avoids “re-inventing the wheel”.

Prokka finds and annotates features (both protein coding regions and RNA genes i.e. tRNA, rRNA) present on on a sequence. Note, Prokka uses a two step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://prodigal.ornl.gov/); second, the function of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

## Learning objectives

At the end of this tutorial you should be able to :

1. input file into Prokka
2. change settings
3. run Prokka, and
4. examine the output: annotated genome.

## Run Prokka

- In Galaxy, load the Prokka interface from the tool panel (**NGS: Annotation -> Prokka**)  
- Set the following parameters:
    - **Contigs to annotate:** <fn>Spades contigs</fn>  
    - **Force GenBank/ENA/DDJB compliance (--compliant): Yes**  
    - **Genus Name:** *Staphylococcus*  
    - **Strain Name:** *aureus*  
    - **Use genus-specific BLAST database: No**  
    - **Locus tag prefix: P**  
    - **Sequencing Centre: V**  
- Click **Execute**  

## Examine the output

Once Prokka has finished, examine each of its output files.

- The gff and gbk files contains all of the information about all of the features annotated (in different formats.)
- The txt file contains a summary of the number of features annotated.
- The faa file contains the protein sequences of the genes annotated.
- The ffn file contains the nucleotide sequences of the genes annotated.

- Download the gff file to your local computer

Now that we have annotated the draft genome sequence, we would like to
view the sequence in the Artemis genome viewer.

- The file is downloaded from the history panel. Identify the file with the .gff extension, expand file header to reveal the expanded file header and download the file using the disk icon ![disk icon](./images/image00.png).

![galaxy file](./images/image01.png)

- Open Artemis and load the downloaded .gff file.

## Discussion - a closer look at the annotated features
FIXME: to add

- Open Artemis and load the downloaded .gff file.
