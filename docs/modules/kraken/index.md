# Species identification with Kraken

Kraken: a tool to taxonomically identify the sample from sequence reads, or a mixed set of sequence reads.

This can be used to confirm the species identification, or for metagenomics.

- e.g. reads from one sample &rarr; kraken &rarr; *Staphylococcus aureus*

- e.g. mixed reads &rarr; kraken &rarr; 50% from *Staphylococcus aureus*, 40% from *Campylobacter concisus*, 10% from Genus *Rothia*, unidentified species.

is this a proportion of reads?
is proportion of reads a proxy for relative abundance?
(e.g. something might not sequence well - be truly abundant but rare in the sample)

We will use Kraken to confirm the identify of reads from one bacterial isolate.

## Get data

R1 and R2  fastq.


## Run Kraken

Galaxy / Kraken


### How it works

Generally: compare sequence to database of known sequence identities.

In detail:

- uses k-mers (default k = 31)
- user to specify a library of genomes
- query database for each k-mer in a sequence
- summarize the hits for all the k-mers from a seq to give the ID
- if seq has no kmers in dbase, it is left unclassified.






## View results

Pie chart


## Next



## Links

[Kraken paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46)
