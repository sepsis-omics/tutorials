# Genome annotation using Prokka on the command line

## Background

In this section we will use a software tool called Prokka to annotate a bacterial genome assembly.

Prokka is a “wrapper”; it collects together several pieces of software (from various authors), and so avoids “re-inventing the wheel”.

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence. Note, Prokka uses a two-step process for the annotation of protein coding regions: first, protein coding regions on the genome are identified using [Prodigal](http://prodigal.ornl.gov/); second, the *function* of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found [here](https://github.com/tseemann/prokka).

## Learning objectives

At the end of this tutorial you should be able to:

1. load a genome assembly into Prokka
2. annotate the assembly using Prokka
3. examine the annotated genome using Artemis

## Input data

Prokka requires assembled contigs.

- Download the assembled contigs from [tba]

## Open the GVL command line

- go to GVL dashboard - go to SSH - find name of virtual machine
- on local machine, open terminal (or equivalent)
- ssh to your GVL

```bash
prokka --outdir mydir --prefix staph contigs.fa
```

output:
- 10 files

prokka will save the output in "mydir"

move there:

```
cd mydir
```

look at the text file:

```
cat file.txt   [check name]
```

look at [screenshot]

look at the annotations in artemis:

```
art mydir/staph.gff
```
look at [screenshot]
