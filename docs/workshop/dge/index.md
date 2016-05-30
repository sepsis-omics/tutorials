# Differential Gene Expression

## Background

Differential Gene Expression (DGE) is the process of determining whether any genes were expressed at a different level between two conditions. For example, the conditions could be wildtype versus mutant, or two growth conditions. Usually multiple biological replicates are done for each condition - these are needed to separate variation within the condition from that between the conditions.

## Learning Objectives

At the end of this tutorial you should be able to:

1. Align RNA-Seq data to a reference genome  
2. Count transcripts for each sample
3. Perform statistical analysis to obtain a list of differentially expressed genes
4. Interpret the DGE list
5. Visualize the results in Degust

## Input data

A typical experiment will have 2 conditions each with 3 replicates, for a total of 6 samples. Each sample will be [RNA-Seq data](/rna/data.md), either as one file per sample (single-end reads / SE) or two files (paired-end reads / PE).

|             | Condition 1 | Condition 2 |
|-------------|-------------|-------------|
| Replicate 1 |     1       |      4      |
| Replicate 2 |     2       |      5      |
| Replicate 3 |     3       |      6      |

FIXME:  load reads

## Prepare reference

FIXME: tricky need it in correct format for htseq-count GFF2 ?

## Align reads

FIXME: use BWA MEM with defaults? for each sample, should be able to use the "6 at once" feature of Galaxy?  

## Count reads

FIXME: htseq-count, use the CDS and RNA features to align to, get count table

## DGE Analysis

### Within Galaxy

FIXME: Need to use Voom/Limma here

### Within Degust

FIXME: need to combine each of the results in count-reads section into a single table (using galaxy table tools?) but need to munge in the annotation as well, so i or simon will need to add new tools to toolshed to do this

## What next?

FIXME
