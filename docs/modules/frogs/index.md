#Metagenomics

-microbial only?

Data: 16S rRNA
Tool:FROGS

Find Rapidly OTU with Galaxy Solution

##Metagenomics background

- uses environmental samples (not cultured)
- uses 16S rRNA if classifying bacteria in sample
- uses WGS if trying to use DNA to identify the functions of genes (and to see what pathways are happening)


- if just using 16S, uses "amplicons"


##notes

swarm clustering is hierarchical

why frogs is a good way to do this

look up frogs at galaxy last years conf

## Overview of FROGS

There are five steps:

- Get data
- Pre-process
- Cluster the sequences
- Remove chimera
- Affiliation (assign identity)

##Get data



##Pre-process

What it does:
- if data is not contiged, will overlap read1 and read2 (allows some mismatch in overlapping region)
- filters out contigs that are too big or small
- if using illumina standard protocol: looks for those primers (primers provided?) - filters out contigs without primers. cuts primer seqs.
- filter out seqs that are too small or poor qual.
- de-replication - remove duplicates but keep a copy of the number of counts for each duplicate

Input: R1 and R2 reads (multiple samples).

Sequencer - Illumina
Input type - Files by samples
Reads alread contiged? - No
Samples - Name - sample1
Reads1 - 01R1.fq
Reads2 - 01R2.fq
Click plus sign for Insert Samples
Reads1 - 02R1.fq
Reads2 - 02R2.fq
Reads 1 size - 250
Reads 2 size - 250
Expected amplicon size - 420
Minimum amplicon size - 380
Maximum amplicon size - 460
Sequencing protocol - Illumina standard
5' primer - GGCGVACGGGTGAGTAA
3' primer - GTGCCAGCNGCNGCGG //note needs to be in 5' to 3' orientation
Execute

Output: 3 files
deprelicated.fasta
count.tsv
report.html

1. report.html

- shows how samples were filtered. eg. the number of reads kept at each filtering stage.
- (what would be a problem to see here? lots of reads removed at a particular filtering stage?) -- e.g. if many (>20%) are lost at the overlapping stage the data could be poor quality.

2. counts.tsv

- number of (the same seq) in sample 1 and sample 2

3. dereplicated.fasta
this is the sequences (just one copy if there are dups)


##Clustering swarm

Sequences are clustered into groups using Swarm.

Input: the fasta and counts file from pre-processing.

What it does:
Sorts reads by abundance.
Clusters the reads into pre-clusters using Swarm and distance parameter of 1.
Sorts these pre-clusters by abundance.
Cluster the pre-clusters using Swarm and a user-specified distance.

settings:

Aggregation distance: 3
Perform deionising clustering step?: Yes
Execute

Outputs:

- abundance file in biom format
- seed_sequences.fasta - the cluster (OTU) representative sequence
- swarms.composition.tsv - what is in each cluster.



##Remove chimera




##Affliation OTU
OTUs- Operation Taxonomic Unit
- this is a cluster of sequences. at 97%


## Comparative metagenomics

compare microbial communities from different environments

do this by comparing the metagenome.
e.g. GC-content, k-mers, size, how many spp, functions of genes



eg for function, can compare to a database (COG or KEGG) and see whether metagenomes differ (statistically test) -



## Other software

QIIME
UPARSE
MOTHUR
MG-RAST



##Links
Slides by LeBras Yvan

https://f1000research.com/slides/5-1832
