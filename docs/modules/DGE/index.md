# Differential Gene Expression

for bacteria, using Galaxy

FIXME: Degust to be installed on mGVL Galaxy

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

Sample data set:

FIXME: can this have a link to swift or a galaxy history URL
http://45.113.234.197/galaxy/history/export_archive?id=6d4c8d216c52b289


E coli
glycolysis growth phenotype
dataset from simon
2 conditions, 3 reps for each condition







## Prepare reference

FIXME: tricky need it in correct format for htseq-count GFF2 ?
(No? seems like GTF is ok? - Anna)

- a reference genome
E coli in fna and gtf formats

## Align reads

FIXME: use BWA MEM with defaults? for each sample, should be able to use the "6 at once" feature of Galaxy?  


- We need to map the transcripts to a reference genome.
- Bacteria don't need splice-aware mapping (don't have introns)
- Galaxy: tools: <ss>NGS Analysis: NGS Mapping: Map with BWA-MEM</ss>
- <ss>Will you select a reference genome from your history or use a built-in index?: Use a genome from history and build index</ss>
- <ss>Use the following dataset as the reference sequence:</ss> E_coli_ref_genome
- <ss>Single or Paired-end reads: single </ss>
- <ss>Select FASTQ dataset</ss>:
    - click on the multiple files button (image) in centre
    - make sure all 6 FASTQ files are in there  
    - hold down shift to select them all (they turn blue)
    - this will map each set of transcripts to the ref genome, so there will be 6 output files in BAM format



    <!-- ## Visualize the mapped reads
    - The mapped reads are now as .bam files which can't be viewed by just clicking on them.
    - could use JBrowse -->

## Count reads

FIXME: htseq-count, use the CDS and RNA features to align to, get count table


HTSeq-count is a tool to count how the reads overlap with genes.
the gene info is from the annotated ref genome.



NGS: RNA Analysis &rarr; htseq-count

- For <ss>Aligned SMA/BAM File</ss>, first select the multiple files icon, then select the 6 bam files.

- For <ss>GFF File</ss>, select the GTF file.

- For <ss>Stranded</ss> select reverse (or if not, yes, or if not, no)

- Click <ss>Execute</ss>

Output:
feature (eg the exons) with the counts for each replicate (so X 6)




## DGE Analysis

### Within Galaxy

FIXME: Need to use Voom/Limma here

<ss>NGS: RNA Analysis &rarr; Differential Count models</ss>
- seems to be options for using any/all of edgeR, DESeq, Voom ?
- so tick option for Voom?




### Within Degust

FIXME: need to combine each of the results in count-reads section into a single table (using galaxy table tools?) but need to munge in the annotation as well, so i or simon will need to add new tools to toolshed to do this


get EC numbers in so Kegg pathway



General degust options:
- Input: read counts (htseq output)
- Configuration: settings
- Name - for data
- Info columns - info for each gene - how to choose
- Add condition ( for all conditions)
- Set min read count to 10 (< 10 reads = gene omitted)
- Save, view
- Execute


## Output

Compare the expression of genes between conditions 1 and 2

## What next?

FIXME

-------------
