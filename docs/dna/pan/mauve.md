# Mauve

This tutorial demonstrates how to use the Mauve software to align genomes. For example, a de novo assembly can be compared against a reference genome to check the assembly. Or, changes between two or more genomes can be examined.

[Link to Mauve](http://darlinglab.org/mauve/mauve.html)

## Pre-requisites
- mGVL instance if using Lubuntu desktop

## Start
- FIXME: use linux version in Lubuntu desktop, or local version e.g. Mac ?
- or you can use in the cmd line

## Input
- 2+ genomes in FASTA(.fasta or .fna)/gbk format
- FIXME: choose data

## How it works
- Mauve finds structural variants in genomes.
- It aligns genoms and finds homologous regions, either from a common ancestor (orthologs) or lateral transfer (xenologs).

## Run
- <ss>File: Align with progressiveMauve</ss>
- <ss>Add sequences</ss>: select 2+ genome FASTA files
- <ss>Output</ss>: provide a name for the output folder
- <ss>Align</ss>
- the Mauve Console will appear while Mauve is running, showing the job status.

## Output
- FIXME: screenshots
- one genome per line
- genomic regions are coloured blocks: locally co-linear blocks (LCB): a block of genome that is unchanged inside, but as a block might have moved/inverted etc.
- you can export a list of annotated homologous features


## Next
