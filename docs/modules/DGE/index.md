# Differential Gene Expression

- for bacteria, using Galaxy tools and Degust (web)
- FIXME: Degust to be installed on mGVL Galaxy

## Background

Differential Gene Expression (DGE) is the process of determining whether any genes were expressed at a different level between two conditions. For example, the conditions could be wildtype versus mutant, or two growth conditions. Usually multiple biological replicates are done for each condition - these are needed to separate variation within the condition from that between the conditions.

## Learning Objectives

At the end of this tutorial you should be able to:

1. Align RNA-Seq data to a reference genome  
2. Count transcripts for each sample
3. Perform statistical analysis to obtain a list of differentially expressed genes
4. Interpret the DGE list
5. Visualize the results in Degust

## Input data - rna-seq reads

A typical experiment will have 2 conditions each with 3 replicates, for a total of 6 samples. Each sample will be [RNA-Seq data](/rna/data.md), either as one file per sample (single-end reads / SE) or two files (paired-end reads / PE).

|             | Condition 1 | Condition 2 |
|-------------|-------------|-------------|
| Replicate 1 |     1       |      4      |
| Replicate 2 |     2       |      5      |
| Replicate 3 |     3       |      6      |

Sample data set:

- E coli, dataset from simon, 2 conditions, 3 reps for each condition
- FIXME: link to swift URL?
- make sure they are all fastqsanger - (change datatype if required)

## Prepare reference

We have two reference genomes for E. coli.

- in FASTA format - nucleotide sequence - for mapping RNA-Seq reads against
- in GTF format (FIXME says GTF but filetype is GFF -- change name to GFF?) which contains gene annotations (counting the mapped reads against to match to gene names)
- FIXME: how to get EC numbers in here for later use in KEGG pathways? they should be in the GTF file shouldn't they (ideally?)

## Align reads

FIXME: use BWA MEM with defaults? for each sample, should be able to use the "6 at once" feature of Galaxy?  [I'm not sure - Anna]

We need to map the RNA-seq reads to a reference genome.

Bacteria don't need splice-aware mapping (don't have introns)

In Galaxy:

- Go to <ss>Tools: NGS Analysis: NGS Mapping: Map with BWA-MEM</ss>
- Under <ss>Will you select a reference genome from your history or use a built-in index?</ss>: *Use a genome from history and build index*
- <ss>Use the following dataset as the reference sequence</ss>: <fn>E_coli fasta file </fn>
- <ss>Single or Paired-end reads</ss>: single
- <ss>Select FASTQ dataset</ss>:
    - click on the multiple files button (image) in centre
    - make sure all 6 FASTQ files are in there  
    - hold down shift to select them all (they turn blue)
    - this will map each set of transcripts to the ref genome, so there will be 6 output files in BAM format
- Click <ss>Execute</ss>

FIXME: takes ages

Note: can't view the BAM files directly. If you want to view, need to make a file using the JBrowse tool.

## Count reads using htseq-count

HTSeq-count is a tool to count how the reads overlap with genes.
The gene info is from the annotated ref genome.

In Galaxy,

- Go to <ss>NGS: RNA Analysis &rarr; htseq-count</ss>
- For <ss>Aligned SMA/BAM File</ss>, first select the multiple files icon, then select the 6 bam files.
- For <ss>GFF File</ss>, select the GTF file.
- For <ss>Stranded</ss> select FIXME ?? [still can't work out what it is]
- Click <ss>Execute</ss>

Output:

- feature (eg the exons) with the counts for each replicate (so X 6)

Get these counts into one file:

- combine into one table: one column of gene IDs, then columns of counts per sample
- best way to do this? (rather than excel)

FIXME: need to combine each of the results in count-reads section into a single table (using galaxy table tools?) but need to munge in the annotation as well, so i or simon will need to add new tools to toolshed to do this)

## DGE Analysis Within Degust

Degust is a tool to test for DGE, or just to view the results of previous testing.
[Link to Degust on github](https://github.com/Victorian-Bioinformatics-Consortium/degust#degust-formerly-known-as-dge-vis)

We will use Degust here on the web.

We will use the table of HTSeq counts prepared above.

### Upload file

1. Go to .....
2. upload file of counts
2. give it a name
3. click on add condition - e.g. control
4. in the adjacent box, use the drop-down menu to select the columns that correspond to this condition
5. click on add condition again - e.g. treatment
6. repeat 4
7. save changes
8. view - this brings up the degust viewing window.

### Test for DGE

Degust statistically tests the counts (choose edgeR or voom/limma- then apply) [is this correct]

Make a note here: the counts alone are not enough because:

- we want to see if some genes are expressed more - in condition x - *compared to other genes* (not how a whole genome is expressed more/less due to the changed condition / randomly / technically) [is this correct?]

### View output

Things to look at:

- We want to see whether the gene expression is higher in different conditions.
Ideally, replicates within each condition should be similar.

Viewing options overview:

- Top panel at right: Configure - click to change the input settings / toggles with view
- Conditions (left hand, top): These are the diff exp conditions
- Centre - top: There are three plots (with options at the right for each)
- Under this is a heat map of gene expression
- Under this is a list of the genes (choose to sort by FDR or by the ?most highly expressed genes?  from each condition)
- Show R Code to see [stuff]

**Plots: Expression: Parallel Coordinates**

- data has been scaled.
- Each condition has a vertical, parallel axis
- we are trying to see the relationship between conditions (axes)
- if lines between axes are parallel, = pos relationship (high A = highB)
- if lines cross between axes, = neg relationship (high A = low B)
change order of axes by dragging the name at the top of each axis. as there is no natural order (as there would be in a time series) we can try a few ways. no absolute way to find "best" order.

if you want to look at a particular section (e.g. low expression in condition A) drag the rectangle that sits on that axis to cover the area of interest. then only genes that connect to that area will be shown. to turn off, drag the top of the rectangle right down past the base of the axis.
what does it mean if lines converge at 0 on one graph - is that the point all other condtions are compared to?
good discussion of PC.  https://eagereyes.org/techniques/parallel-coordinates

- options.....

**Plots: Expression: MA plot**

- what does it stand for
- each gene is a dot on the graph
- y axis: fold.
- x axis: expression (averaged over samples -- does this mean all reps in all conditions?)
fold change vs expression.
called MA bc y M (y) is the ratio (folds) and A (x) is the avg expression
what are we looking for? high fold changes at low and high expressions?

"The MA-plot represents each gene with a dot. The x axis is the average expression over all samples, the y axis the log2 fold change of normalized counts (i.e the average of counts normalized by size factor) between treatment and control. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call."

FC = fold change: doubling = fold change of 1, one-fold increase. 100% increase.
20 -> 80 = 3 fold
80 -> 20 = neg ¾ fold
bc 80->40 is neg ½ fold

difference/ orig         x100
so 20 -> 80     60/20 = 3   = 300% increase
so 80 -> 20     -60/80 = -0.75  =  75% decrease

but note there is no single fold change of x with biological importance

- options...

**Plots: MDS plot**
- stuff
- options...

**Kegg Pathway**

- in options.

A pathway is a drawn network to show interaction between molecules, including some or all of genes, proteins, RNAs, chemical reactions. E.g.
The Kegg pathway database: (Kyoto Encyclopedia of Genes and Genomes)
numbers are EC numbers - enzyme commission ->  enzyme/s that catalyze a reaction (might be >1)
Click on a pathway in the drop down menu
highlighted numbers (relate to lines in the parallel coords graph - but what is it?)


FIXME: to do:
get EC numbers in so Kegg pathway
can these be obtained from the gene IDs in the annotation?

**Heatmap**
no heading, but this is the patchwork under the plots
conditions on the left
don't know what the blocks are - one per gene?
colours: fold change - e.g. blue, lower expression, red, increased expression.

**Genes**
list of genes.

### Within Galaxy

FIXME: Need to use Voom/Limma here

<ss>NGS: RNA Analysis &rarr; Differential Count models</ss>
- seems to be options for using any/all of edgeR, DESeq, Voom

## What next?

FIXME

## More information

http://rnaseq.uoregon.edu/#exp-design-differential-gene-expression



voom paper
 https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29
