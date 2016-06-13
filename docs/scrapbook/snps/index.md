# Variant calling

## Background

Variant calling is the process of identifying differences between to genome samples.
Usually differences are limited to single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels). Larger structural variation such as inversions, duplications and large deletions are not typically covered by "variant calling".

## Learning Objectives

At the end of this tutorial you should be able to:

1. Find variants between a reference genome and a set of reads
2. Determine the effect of those variants on genomic features
3. Understand if the SNP is potentially affecting the phenotype

## Experiment

FIXME: You are working on a bug and you notice one of them is forming smaller colonies than normal. You want to find out why this msall colony vairant (SCV) is doing at the DNA level.

## Prepare reference

FIXME: need FASTA or prefer an *annotated* genome eg. GBK or GFF3 of the original strain you used

!!! note
    Please make sure your reference genome includes all chromosomes and plasmids

## Align reads

FIXME: BWA MEM align the SCV reads

## Call variants

FIXME: freebayes?  varscan2 ?
FIXME: talk about multimapping reads?

## Filter variants

FIXME:  vcffilter? something else?  mindepth, homozygous?

## Annotate consequencs

FIXME: snpEff - but it is hard to add a genome

!!! hint
    Just use Snippy and all this will happen magically?

## What next?

* SNPs can be used to build [phylogentic trees](/trees/index.md).
