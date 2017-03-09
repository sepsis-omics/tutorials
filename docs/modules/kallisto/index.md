#Kallisto


## Background

Two conditions
RNA-seq from each (as a proxy for which genes are expressed)
Are genes expressed more in one condition?

In this tutorial, we have
- one bacterial strain.
- its reference transcriptome (we can get the from the reference genome: extract out transcriptomic regions)
- the RNA-seq reads in condition 1 and condition 2 (with x replicates in each)

Process:
- Take one set of RNA-seq reads and count how many per gene (or region): use Kallisto
- Repeat for every set of reads; Collect into a table of counts
- Test whether counts are different between conditions (statistically): use Degust (or, sleuth)

## Get data

e.g. E_coli set
a reference transcriptome, e.g. the .ffn from Prokka

** Need to get these names as setA, ptsG etc? **



## Run Kallisto to get counts

### Index the transcripts file
```
kallisto index -i transcripts.idx [this is the output file] transcripts.ffn
```

### Run kallisto for every read set
```
kallisto quant -i transcripts.idx -o output_LB1 -b 100 --single -l 100 -s 10 Ecoli1.fastq.gz
```

** estimated length (l) and s (SD) . need to check this

** Rename the col header here for each


## Combine the counts

We need to combine all the counts into one table.

### Extract required columns

Condition 1/ rep1: extract columns 1 (IDs) and 4 (counts)
```
cut -f1,4 -d$'\t' abundance.tsv > Lb1.tsv
```

Then for other conditions/reps, just extract column 4 (counts):

```
cut -f4 -d$'\t' abundance.tsv > condx_repx.tsv
```


## Paste together

```
paste newfile1 newfile2 > counts.tsv
```

### Remove existing header

```
tail -n +2 counts.tsv > headless.tsv
```

### Add new header

```
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





## Links

Kallisto paper
http://www.nature.com.ezp.lib.unimelb.edu.au/nbt/journal/v34/n5/full/nbt.3519.html

Kallisto + sleuth paper
http://biorxiv.org/content/biorxiv/early/2016/06/10/058164.full.pdf
