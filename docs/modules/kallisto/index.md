<br>
#Kallisto


This tutorial is about differential gene expression in bacteria, using tools on the command-line tools (kallisto) and the web (Degust).

## Background
Differential Gene Expression (DGE) is the process of determining whether any genes were expressed at a different level between two conditions. For example, the conditions could be wildtype versus mutant, or two growth conditions. Usually multiple biological replicates are done for each condition - these are needed to separate variation within the condition from that between the conditions.

## Learning Objectives

At the end of this tutorial you should be able to:

1. (Pseudo-)align RNA-Seq data to a reference transcriptome and count: kallisto  
2. Perform statistical analysis to obtain a list of differentially expressed genes: Degust
3. Visualize and interpret the results

## Prepare workspace

e.g. log in to the ABRPI mGVL or use your own workspace with kallisto installed.

(check kallisto is on the abrpi)

Make a folder for these analyses:

```bash
mkdir E_coli_DGE
cd E_coli_DGE
```

## Get data


### RNA-Seq reads

A typical experiment will have 2 conditions each with 3 replicates, for a total of 6 samples.

- Our RNA-seq reads are from 6 samples in <fn>FASTQ</fn> format.
    - We have single-end reads; so one file per sample.
    - Data could also be paired-end reads, and there would be two files per sample.
- These have been reduced to 1% of their original size for this tutorial.
- The experiment used the bacteria *E. coli* grown in two conditions.
    - Files labelled "LB" are the wildtype
    - Files labelled "MG" have been exposed to 0.5% &alpha;MG - alpha methyglucoside (a sugar solution).



### Reference transcriptome

We need to count the number of RNA-seq reads (that exist as fragments) that match different transcripts in the genome, including those for protein-coding sequences (such as genes) and RNA sequences (such as tRNA and mRNA).

Therefore, we need a subset of the whole genome. We will use a genbank file for this sample and extract out all features for which a transcript could map to.

Generate the reference transcriptome:

Here we use a custom python script called "genbank_to_kallisto.py". You also need to enter the name of an <fn>input genbank file</fn>, and the names for two output files that will be generated: a <fn>transcripts.fasta</fn> file and a <fn>table.tsv</fn> file. We type in "python3" at the start to run this script using python.

Genbank file:

(put this script on abrpi? genbank file needs to be in this folder?)


```python
python3 genbank_to_kallisto.py Ecoli.gbk Ecoli_transcripts.fasta Ecoli_table.tsv
```
Alternatively, download the already-generated fasta and table.


### Where to get the files for this tutorial

(where to put - eg swift container? -- how to upload)

- 6 x RNA seq reads in fastq.gz
- 1 x Ref transcriptome in .fasta
- 1 x table of features for degust.


## Run Kallisto



Kallisto will count the reads per transcript.

### Index the transcripts file

``` bash
kallisto index -i transcripts.idx [this is the output file] transcripts.ffn
```

### Run kallisto for every read set

```bash
kallisto quant -i transcripts.idx -o output_LB1 -b 100 --single -l 100 -s 10 Ecoli1.fastq.gz
```

** estimated length (l) and s (SD) . need to check this **


### Combine the counts

We need to combine all the counts into one table.

### Extract required columns

For first .tsv file: extract columns 1 (IDs) and 4 (counts)

```bash
cut -f1,4 -d$'\t' abundance.tsv > Lb1.tsv
```

Then for other .tsv files, just extract column 4 (counts):

```bash
cut -f4 -d$'\t' abundance.tsv > condx_repx.tsv
```


### Paste together

```bash
paste newfile1 newfile2 > counts.tsv
```

### Remove existing header

```bash
tail -n +2 counts.tsv > headless.tsv
```

### Add new header

```bash
echo -e "gene\tcondition1\tcondition2" | cat - headless.tsv > newtable.tsv
```

### Examine the file

```
less newtable.tsv
```

Look at columns 2 and 3 of the counts of reads against particular transcripts. Do any look much higher in a particular condition?

## Test for DGE

Degust:

upload counts file.

Name the file.
Info columns = gene
Analyze server side - tick
Min gene read coutn - 10
Add condition - tick Lb1, Lb2, Lb3; condition will be named Lb
Add condition - tick Mg1, Mg2, Mg3; condition will be named Mg
Save changes; view

## Links

Kallisto paper
http://www.nature.com.ezp.lib.unimelb.edu.au/nbt/journal/v34/n5/full/nbt.3519.html

Kallisto + sleuth paper
http://biorxiv.org/content/biorxiv/early/2016/06/10/058164.full.pdf
