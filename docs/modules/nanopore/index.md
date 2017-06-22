# Nanopore data - bacterial genome assembly

Base calling the reads, then canu, then Polish, or use racon


## Get data

MinION data => Metrichor server => calls bases => fast5 files

- check the output from Metrichor
- we could assume we have fastq files, but we need original (fast5?) files to use in Nanopolish


Enterobacter kobei

- from DOI: 10.1099/mgen.0.000085
-  ERR1341575 (MinION pass reads) - fastq


Alternative data:
E coli K-12 MG1655 from Loman
http://www.nature.com.ezp.lib.unimelb.edu.au/nmeth/journal/v12/n8/full/nmeth.3444.html



## Assess data: poretools

https://poretools.readthedocs.io/en/latest/

### Convert fast5 to fastq

```
poretools fastq fast5/
```

(is fast5 a dir?)

### plot histogram of read sizes

```
poretools hist fast5/
```

- how to viz?



## Assemble

- canu


```
canu -p Ekobei -d Ekobei -genomeSize =?m -nanopore-raw fastq


```
=> 14 contigs, max 1.5m
need fewer contigs.





## Polish

- nanopolish
- uses info from raw signal data (before they have been called as bases)
-



http://simpsonlab.github.io/2015/03/30/optimizing-hmm/

https://github.com/jts/nanopolish

- this can polish the assembly



## Evaluate
