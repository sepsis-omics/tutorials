# Salmon

Available on the ABRPI mGVL command line.

Not yet available in mGVL Galaxy.

What does Salmon do?

* quantifies transcripts from RNA-seq data
* uses "quasi-mapping" yo map reads to reference - more detail [here](http://bioinformatics.oxfordjournals.org/content/32/12/i192.full)




##Start

Open a terminal.
Log in to the mGVL. [is Salmon on the general mGVL?]
Make a folder called RNA_seq
```text
mkdir salmon
cd salmon
```



##Input

* RNA-seq reads in FASTA/FASTQ format
* Reference in FASTA format -- need to get the ecoli transcriptome - containing genes and RNAs

how to get?


```genbank2fasta.pl < file.gbk > file.ffn
```
 

get these files from [Nectar container] - how to scp these in using terminal


##Run

Index the reference.
[Note: the reads need to be in random order with respect to the reference. Are they?]

```text
salmon index -t ecoli.fasta -i ecoli_index
```

**-t** specifies the input transcript fasta file
**-i** specifies the name of the output index (folder)

Quantify the reads
[Salmon says "quantify, don't count" ? ]

```text
salmon quant -i ecoli_index -l A -r rna_seqs.fasta.gz -p 8 --numAuxModelSamples 50000 -o quants
```

**-i** is the name of the index for the reference file
**-l A** means that salmon will determine the library type for the input sequencing reads (whether they are stranded, anything else? page just says etc.)
**-p 8** use 8 threads (possible on this machine?)
**numAuxModelSamples** default was 5 million but there are only ~2million reads per file; set this to 50,000 (not sure if ok)
**-o quants** output directory for results is called quants

Warming:
Detected a *potential* strand bias > 1% in an unstranded protocol check the file: quants/lib_format_counts.json for details

(is this a problem?)

##Output

```text
less quant.sf
```

shows the output:

TPM - Transcripts per Million
NumReads - estimated number of reads for this reference
