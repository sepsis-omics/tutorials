# Spades - commandline

This tutorial follows on from "PacBio assembly with commandline tools".


## Short-read assembly: a comparison

So far, we have assembled the long PacBio reads into one contig (the chromosome) and found an additional plasmid in the Illumina short reads.

If we only had Illumina reads, we could also assemble these using the tool Spades.

You can try this here, or try it later on your own data.

## Get data

We will use the same Illumina data as we used above:

- <fn>illumina_R1.fastq.gz</fn>: the Illumina forward reads
- <fn>illumina_R2.fastq.gz</fn>: the Illumina reverse reads

This is from Sample 25747.

## Assemble

Run Spades:

```text
spades.py -1 illumina_R1.fastq.gz -2 illumina_R2.fastq.gz --careful --cov-cutoff auto -o spades_assembly_all_illumina
```

- `-1` is input file of forward reads
- `-2` is input file of reverse reads
- `--careful` minimizes mismatches and short indels
- `--cov-cutoff auto` computes the coverage threshold (rather than the default setting, "off")
- `-o` is the output directory

## Results

Move into the output directory and look at the contigs:

```text
infoseq contigs.fasta
```

<!-- ### Questions

How many contigs were found by Spades?

- many

How does this compare to the number of contigs found by assembling the long read data with Canu?

- many more.

Does it matter that an assembly is in many contigs?

- Yes

  - broken genes => missing/incorrect annotations
  - less information about structure: e.g. number of plasmids

- No

  - Many or all genes may still be annotated
  - Gene location is useful (e.g. chromosome, plasmid1) but not always essential (e.g. presence/absence of particular resistance genes)

How can we get more information about the assembly from Spades?

- Look at the assembly graph <fn>assembly_graph.fastg</fn>, e.g. in the program Bandage. This shows how contigs are related, albeit with ambiguity in some places.
-->

## Next

Run "Prokka" to annotate the contigs.
