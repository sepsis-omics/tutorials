# Snippy-core


Run Snippy
  - ref genome + sets of reads
  - one output is a whole-genome alignment of snps (not just the core-genome snps)
  - => gubbins



## Get data

Upload all these to galaxy:

the Staph aureus BPH2986.gbk (from prokka annotation) =25749
(from sepsis done)

  - manually change datatype to "genbank"

the read files from samples 25745-25748 (paired end)

  - change datatype from fastq to fastqsanger

### Make dataset collections

In Galaxy - see the tick box at the top of the history pane.
Click on this - select all the files for the collection

change the search bar to say R1 instead of -1
check the pairs
auto pair
create list
name list
save



https://galaxyproject.org/tutorials/collections/





## Run Snippy

Reference type: Genbank. choose genbank file

Single or Paired-end reads: Paired collection
