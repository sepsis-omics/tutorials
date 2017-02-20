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

Closely-related sequences may form chimeras (mixed sequences) during PCR (libray prep). This step removes these sequences.


How it works:
splits input data into samples
uses vsearch to find chimeras in each sample
removes these


Sequences file: seed_sequences.fasta
Abundance type: BIOM file
Abundance file: abundance.biom
Execute



Outputs:

A filtered file containing no chimeras non_chimera.fasta
A filtered abundance file containing no chimeras non_chimera.biom
Summary  report.html

seems to be a lot removed?

## Filters

The OTUs have now been clustered.
In this step, we will filter out some of the OTUs. (? we have set 2 samples per OTU - i.e. OTU must be in both sample? OTU must contain at least 0.005 percent of all the seqs?)

Sequences file: non_chimera.fasta
Abundance file:  non_chimera_abundance.biom
*** THE FILTERS ON OTUS IN SAMPLES, OTUS SIZE and SEQUENCE PERCENTAGE: Apply filters
Minimum number of samples: 2
Minimum proportion/number of sequences to keep OTU: 0.00005
N biggest OTU: leave blank?
*** THE FILTERS ON RDP: No filters
*** THE FILTERS ON BLAST: No filters
*** THE FILTERS ON CONTAMINATIONS: No filters
Execute

Outputs:
a filtered abundance and fasta file
an excluded.txt
a report.html

only kept %20 of OTUs?



##Affliations OTU
OTUs- Operation Taxonomic Unit
- this is a cluster of sequences.
- this step adds the taxonomy to the abundance file. (why?)
- uses the SILVA database for rRNA- that has been filtered for only (16S or all 30S?), non-redundant, keeping certain taxonomic levels, and split into pro and eukaryotes.
- this step uses blastn+ to align each OTU against seqs in the dbase, keeping the best.
- It can return multi-affiliation - see notes below tool panel.

Using reference database: silva123
Also perform RDP assignation: No
OTU seed sequence: output.fasta from step 4
Abundance file: output.biom from step 4
Execute

Outputs:
abundance file with affiliation   affiliation.biom - note this biom file is not human-readable. You can convert it with the frogs biom to tsv tool.

report.html


##Affiliations stat
- Computes some statistics
- generates a report of the OTUs/taxonomy found

Abundance file: affiliation.biom from step 5
Rarefaction ranks: Class Order Family Genus Species
Affiliation processed: FROGS blast
Execute

Outputs:
summary.html

- click on Display global distribution
- Pie chart thing - how to read - start in the centre. major groups each have a segment by colour. then more detail as you go outwards.
- click cross to exit
- but does this show sample 1 and sample 2 combines? 

- tick the boxes next to the samples
- then with selection order - click Display rarefaction
- Rarefaction curve - what does it mean




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


The SILVA database
https://www.arb-silva.de/

SILVA paper 10.1093/nar/gks1219

More about 16S rRNA:
The 16S sequence is RNA that makes up part of the bacterial ribosome (together with other rRNA and protein). The operon is 16S, 23S, 5S - and there may be many copies of the operon (7 in E coli).
