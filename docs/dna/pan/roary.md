# Roary

cmdline

This tutorial demonstrates how to calculate the pan and core genomes of a set of input bacterial species, using Roary.

[Roary code and manual on github](https://sanger-pathogens.github.io/Roary/)

[Roary paper](http://bioinformatics.oxfordjournals.org/content/31/22/3691)

## Pre-requisites

- background: pan genomes
- a mGVL instance

## Start
- **via local Terminal**: in your terminal, ssh into your mGVL, but make sure you put in -X -Y after ssh so that xquartz can view files later. (FIXME: word better)
- **or via virtual desktop**: Go to your mGVL dashboard. Click on the link to the Lubuntu desktop. A virtual desktop will open in a new browser window. Enter username: ubunutu; and your GVL password. Click on terminal in the top left corner.

## Input data
- Roary takes .gff files produced by Prokka. A gff file has sequences and annotations.
- [file formats](../dna/anno/prokka) FIXME: link to proper page
- Get files into mGVL. (FIXME: explain how, wget etc or from GenomeSpace)
- put all gff files into a folder
- FIXME: choose a good sample set - ideally something that usefully shows how AMR genes can be present/absence in a group? for a draft can use the Listeria tutorial at https://github.com/microgenomics/tutorials/blob/master/pangenome.md

## Run

- navigate into the place where the gff folder is.
```bash
roary -f ./results ./gff_files/*.gff
```
- "-f ./results" puts the output into a directory called results

## How it works
- Based on the input genomes, Roary works out which genes are shared between all (core) and which are not (accessory).
- It uses the protein-coding genes from each of the input genomes.
- converts to protein seqs
- similar protein seqs are clustered progressively.
- each sample: will be labelled with presence/absence of orthologous genes.

## Output

### summary statistics:

```
more summary_statistics.txt
```
- you will see the number of core genes, shell genes, etc.
- `q` to exit viewing

### gene presence/absence graphically:
```
roary2svg.pl gene_presence_absence.csv > pan_genome.svg
```
- (if you have logged in with -X -Y)
```
firefox pan_genome.svg &
```
- then `enter`
- the & makes it run in the background
- a firefox window should open with the svg image
- (later: close the firefox window to stop this job running in the background)

### list of genes that are present/absent:
- view the gene_presence_absence.csv by (FIXME)
- lots of information about this file in the roary website (FIXME summarize?)

### query the pan genome:
- copy the input .gff files into the results folder (FIXME: do this earlier)
- cd into this folder
```
query_pan_genome -a intersection *.gff
```
- this finds the core genes
```
more pan_genome_results
```
- shows the list of genes found in the core genome.
- `q` to exit viewing

## Advanced options
FIXME: update firefox on mGVL so can run phandango

### Run roary and create an alignment of core genes:
```bash
roary -f ./results -e -n -p 8 ./gff_files/*.gff
```
- "-f ./results" puts the output into a directory called results
- "-e -n" creates an alignment of core genes using mafft
-  "-p 8" gives 8 threads - optional, if you know how many you have

### Generate a tree based on the presence/absence of core genes:
- navigate into the results folder that you want to use.
```bash
FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick
```
- (By default, roary will also have created a (very quick) tree from the accessory genes.)
- [FastTree information and options](http://meta.microbesonline.org/fasttree/).

### Use roary_plots.py to generate plots:
- navigate into the results folder that you want to use.
```
python roary_plots.py core_gene_alignment.nwk gene_presence_absence.csv
```
- output: pangenome matrix, frequency plot, pie chart.
- view these by typing `firefox [filename]` and a firefox window will open to show the image. You need to close the window before you open the next image.

## What next?
View using Phandango; [tutorial here.](../../viz/phandango/index.md)

## More information
- [another Roary tutorial](https://github.com/microgenomics/tutorials/blob/master/pangenome.md)
