# Galaxy - Velvet

Modified from LSCC tutorial by [Simon Gladman](mailto:simon.gladman@unimelb.edu.au) - VLSCI

## Tutorial Overview

In this tutorial we cover the concepts of Microbial de novo assembly using a very small synthetic dataset from a well studied organism.

## Background

**Where is the data in this tutorial from?**

The data for this tutorial is from a whole genome sequencing experiment of a multi-drug resistant strain of the bacterium *Staphylococcus aureus*. The DNA was sequenced using an Illumina GAII sequencing machine. The data we are going to use consists of about 4 million x 75 base-pair, paired end reads (two FASTQ read files, one for each end of a DNA fragment.) The data was downloaded from the [NCBI Short Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra/) (http://www.ncbi.nlm.nih.gov/sra/). The specific sample is a public dataset published in April 2012 with SRA accession number ERR048396.

We will also use a FASTA file containing the sequences of the Illumina adapters used in the sequencing process. It is desirable to remove these as they are artificial sequences and not part of the bacterium that was sequenced.

We will use software called [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) (Zerbino et al 2008) for the main de novo assembly, as well as some other peripheral software for pre- and post-processing of the data. Details of these can be found in the background document linked above.

**The protocol:**

We are performing a de novo assembly of the read data into contigs and then into scaffolds (appropriately positioned contigs loosely linked together). We firstly need to check the quality of the input data as this will help us choose the most appropriate range of input parameters for the assembly and will guide us on an appropriate quality trimming/cleanup strategy. We will then use an iterative method to assemble the reads using the [Velvet Optimiser](http://www.vicbioinformatics.com/software.velvetoptimiser.shtml) (a program that performs lots of Velvet assemblies searching for an optimum outcome.) Once this is complete we will obtain summary statistics on the final results (contigs) of the assembly.

More information about this protocol at the end of this tutorial.

**The protocol in a nutshell:**

*Input:* Raw reads from sequencer run on microbial DNA sample.

*Output:* File of assembled scaffolds/contigs and associated information.


## Input data

  1. On the Galaxy tools panel, click on **Get data -> Upload File**.
  2. Click on the **Paste/Fetch Data** button.
  3. Paste the URL: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_1.fastq.gz* into the text box. Change the type to *fastqsanger* (Not *fastqcsanger*).
  4. Click on the **Paste/Fetch Data** button again.
  5. Paste the URL: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_2.fastq.gz* into the text box and change it's type to *fastqsanger* as well.
  6. Repeat the process for the last URL: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/illumina_adapters.fna* , but make it's type *fasta*
  7. Click on the **Start** button. Once all of the uploads are at 100%, click on the **Close** button.
  8. When the files have finished uploading, rename them to ‘ERR048396_1.fastq’, ‘ERR048396_2.fastq’ and ‘illumina_adapters.fna’ respectively by clicking on the <img src="../media/Galaxy-edit.png" width=20 /> icon to the top right of the file name in the right hand Galaxy panel (the history panel)

You should now have the following files in your Galaxy history:

* *ERR048396_1.fastq* - forward reads in fastq format
* *ERR048396_2.fastq* - reverse reads in fastq format
* *illumina_adapters.fa* - Illumina adapter sequences in fasta format

Click on the <img src="../media/Galaxy-view.png" width=20 /> icon to the top right of each fastq file to view the first part of the file

## Section 1: Quality control

The basic process here is to collect statistics about the quality of the reads in the sample FASTQ readsets. We will then evaluate their quality and choose an appropriate regime for quality filtering using Trimmomatic (a FASTQ read quality trimmer.)

### Run FastQC on both input read files
1. From the tools menu in the left hand panel of Galaxy, select **NGS QC and manipulation > FastQC: Comprehensive QC** (down the bottom of this category) and run with these parameters:
    * "FASTQ reads": *ERR048396_1.fastq*
    * Use default for other fields
2. Click **Execute**
3. Now repeat the above process on the second read file: *ERR048396_2.fastq*

It is important to do both read files as the quality can be very different between them.

##### Figure 1: Screenshot of FastQC interface in Galaxy

<img src="../media/screenshot-fastqc-interface.png" width=80% />

### Examine the FastQC output

You should have two output objects from the first step:

* *FastQC_ERR048396_1.fastqc.html*
* *FastQC_ERR048396_2.fastqc.html*

These are a html outputs which show the results of all of the tests FastQC performed on the read files.

1. Click on the <img src="../media/Galaxy-view.png" width=20 /> icon of each of these objects in turn to see the FastQC output.

  The main parts of the output to evaluate are:

  * Basic statistics. This section tells us that the ASCII quality encoding format used was Sanger/Illumina 1.9 and the reads are length 75 and the percent GC content of the entire file is 35%.
  * Per base sequence quality. In the plot you should see that most of the early bases are up around the '32' mark and then increase to 38-40, which is very high quality; The spread of quality values for the last few bases increases and some of the outliers have quality scores of less than 30. This is a very good quality dataset. 20 is often used as a cutoff for reliable quality.

##### Figure 2: Screenshot of FastQC output in Galaxy

  <img src="../media/screenshot-fastqc-output.png" width=80% />

### Quality trim the reads using Trimmomatic.

1. From the tools menu in the left hand panel of Galaxy, select **NGS QC and manipulation > Trimmomatic** and run with these parameters (only the non-default selections are listed here):

    * "Input FASTQ file (R1/first of pair)": *ERR048396_1.fastq*
    * "Input FASTQ file (R2/second of pair)": *ERR048396_2.fastq*
    * "Perform initial ILLUMINACLIP step?": *Yes*
    * "Adapter sequences to use": *TruSeq3 (additional seqs) (paired end, for MiSeq and HiSeq)*
    * "How accurate ... read alignment": *40*
    * "How accurate ... against a read": *15*
    * We will use the default settings for the SLIDING_WINDOW operation but we need to add a few more Trimmomatic operations.
    * Click **Insert Trimmomatic Operation**
        * Add *Cut bases ... (LEADING)*
        * "Minimum quality required to keep a base": *15*
    * Repeat the **Insert Trimmomatic Operation** for:
        * Trim trailing bases, minimum quality: *15*
        * Minimum length read: *35*

2. Click **Execute**
      <br></br>

##### Figure 3: Screenshot of Trimmomatic inputs in Galaxy

<img src="../media/screenshot-trimmomatic-inputs.png" width=80% />

### Examine the Trimmomatic output FastQ files.

  You should have 4 new objects in your history from the output of Trimmomatic:

  * *Trimmomatic on data 2 and data 1 (R1 Paired)*
  * *Trimmomatic on data 2 and data 1 (R1 Unpaired)*
  * *Trimmomatic on data 2 and data 1 (R2 Paired)*
  * *Trimmomatic on data 2 and data 1 (R2 Unpaired)*

  Click on the <img src="../media/Galaxy-view.png" width=20 /> on one of the objects to look at its contents. You’ll notice that not all of the reads are the same length now, as they have had the illumina adapters cut out of them and they’ve been quality trimmed.

----------------------------------------------

## Section 2: Assemble reads into contigs with Velvet and the Velvet Optimiser

The aim here is to assemble the trimmed reads into contigs/scaffolds using Velvet and the Velvet Optimiser.

We will use a single tool, Velvet Optimiser, which takes the trimmed reads from Trimmomatic and performs numerous Velvet assemblies to find the best one. We need to add the reads in two separate libraries. One for the still paired reads and the other for the singleton reads orphaned from their pairs by the trimming process.

[Click here for a more detailed explanation of Velvet assemblies and the Velvet Optimiser](assembly-background.md)

### De novo assembly of the reads into contigs

1. From the tools menu in the left hand panel of Galaxy, select **NGS: Assembly -> Velvet Optimiser** and run with these parameters (only the non-default selections are listed here):
    * "Start k-mer value": *55*
    * "End k-mer value": *69*
    * In the input files section:
        * "Select first set of reads": *Trimmomatic on data 2 and data 1 (R1 paired)*
        * "Select second set of reads": *Trimmomatic on data 2 and data 1 (R2 paired)*
    * Click the **Insert Input Files** button and add the following:
        * "Single or paired end reads": *Single*
        * "Select the reads": *Trimmomatic on data 2 and data 1 (R1 unpaired)*
    * Repeat the above process to add the other unpaired read set *Trimmomatic on data 2 and data 1 (R2 unpaired)* as well.
2. Click **Execute**.

##### Figure 4: Screenshot of Velvet Optimiser inputs in Galaxy

<img src="../media/screenshot-velvetoptimiser-inputs.png" width=80% />


### Examine assembly output
Once step 1 is complete, you should now have 2 new objects in your history:
* *VelvetOptimiser on data 9, data 7, and others: Contigs*
* *VelvetOptimiser on data 9, data 7, and others: Contig Stats*

Click on the <img src="../media/Galaxy-view.png" width=20 /> icon of the various objects.

* Contigs: You’ll see the first MB of the file. Note that the contigs are named NODE_XX_length_XXXX_cov_XXX.XXX. This information tells you how long (in k-mer length) each contig is and what it’s average k-mer coverage is. (See detailed explanation of Velvet and Velvet Optimiser for explanation of k-mer coverage and k-mer length.)

* Contig stats: This shows a table of the contigs and their k-mer coverages and which read library contributed to the coverage. It is interesting to note that some of them have much higher coverage than the average. These are most likely to be repeated contigs. (Things like ribosomal RNA and IS elements.)

##### Figure 5: Screenshot of assembled contigs (a) and contig stats (b)

##### a

<img src="../media/screenshot-contigs.png" width=80% />

##### b

<img src="../media/screenshot-contig-stats.png" width=80% />


### Calculate some statistics on the assembled contigs

1. From the tools menu in the left hand panel of Galaxy, select **FASTA Manipulation -> Fasta Statistics** and run with these parameters:
    * "Fasta or multifasta file": *Velvet Optimiser ... Contigs*
2. Click **Execute**
3.  Examine the Fasta Stats output

You should now have one more object in your history: *Fasta Statistics on data 10: Fasta summary stats*

  Click on the <img src="../media/Galaxy-view.png" width=20 /> icon next to this object and have a look at the output. You’ll see a statistical summary of the contigs including various length stats, the % GC content, the n50 as well as the number of contigs and the number of N bases contained in them.

----------------------------------------------

## Section 3: Extension.

Examine the contig coverage depth and blast a high coverage contig against a protein database.

### Examine the contig coverage depth.

Look at the Contig Stats data (Velvet Optimiser vlsci on data 8, data 9, and data 7: Contig stats) by clicking on the <img src="../media/Galaxy-view.png" width=20 /> icon. Note that column 2 contig length (lgth), shows a number of very short contigs (some are length 1).

  * We can easily filter out these short contigs from this information list by using the **Filter and Sort -> Filter tool.**
  * Set the following:
    * "Filter": *Velvet Optimiser on data 8, data 7 and others: Contig stats*
    * "With the following condition": *c2 > 100*
  * Click **Execute**


The new data object in the history is called: *Filter on data 11*.

Click on its <img src="../media/Galaxy-view.png" width=20 /> icon to view it. Look through the list taking note of the coverages. Note that the average of the coverages (column 6) seems to be somewhere between 16 and 32.  There are a lot of contigs with coverage 16. We could say that these contigs only appear once in the genome of the bacteria. Therefore, contigs with double this coverage would appear twice. Note that some of the coverages are >400! These contigs will appear in the genome more than 20 times!

Lets have a look at one of these contigs and see if we can find out what it is.

### Extract a single sequence from the contigs file.
Note the contig number (column 1 in the Contig stats file) of a contig with a coverage of over 300. There should be a few of them. We need to extract the fasta sequence of this contig from the contigs multifasta so we can see it more easily.

To do this we will use the tool:

  * **Fasta manipulation -> Fasta Extract Sequence**
  * Set the following:
    * "Fasta or multifasta file": *Velvet Optimiser ... : Contigs*
    * "Sequence ID (or partial): *NODE\_1_...* (for example)
  * Click **Execute**

The new data object in the history is called: *Fasta Extract Sequence on data 10: Fasta*.

Click on its <img src="../media/Galaxy-view.png" width=20 /> icon to view it. It is a single sequence in fasta format.

### Blast sequence to determine what it contains.

We want to find out what this contig is or what kind of coding sequence (if any) it contains. So we will blast the sequence using the NCBI blast website. (External to Galaxy). To do this:

  * Bring up the sequence of the contig into the main window of the browser by clicking on the <img src="../media/Galaxy-view.png" width=20 /> icon if it isn’t already.
  * Select the entire sequence by clicking and dragging with the mouse or by pressing ctrl-a in the browser.
  * Copy the selected sequence to the clipboard.
  * Open a new tab of your browser and point it to: http://blast.ncbi.nlm.nih.gov/Blast.cgi
  * Under the BASIC BLAST section, click “blastx”.
  * Paste the sequence into the large text box labelled: Enter Accession number(s), gi(s) or FASTA sequence(s).
  * Change the Genetic code to: Bacteria and Archaea (11)
  * Click the button labelled: BLAST

After a while the website will present a report of the blast run. Note that the sequence we blasted (if you chose NODE_1) is identical to part of a transposase gene (IS256) from a similar Staphylococcus aureus bacteria. These transposases occur frequently as repeats in bacterial genomes and so we shouldn’t be surprised at its very high coverage.

##### Figure 6: Screenshot of the output from the NCBI Blast website

<img src="../media/screenshot-blast.png" width=80% />



----------------------------------------------

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does it's job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)





## De novo assembly with Velvet and the Velvet Optimiser.

### Velvet

Velvet is software to perform dna assembly from short reads by manipulating de Bruijn graphs. It is capable of forming long contigs (n50 of in excess of 150kb) from paired end short reads. It has several input parameters for controlling the structure of the de Bruijn graph and these must be set optimally to get the best assembly possible. Velvet can read Fasta, FastQ, sam or bam files. However, it ignores any quality scores and simply relies on sequencing depth to resolve errors. The Velvet Optimiser software performs many Velvet assemblies with various parameter sets and searches for the optimal assembly automatically.

### de Bruijn graphs

A de Bruijn graph is a directed graph which represents overlaps between sequences of symbols. The size of the sequence contained in the nodes of the graph is called the word-length or k-mer size. In Figure 2, the word length is 3. The two symbols are 1 and 0. Each node in the graph has the last two symbols of the previous node and 1 new symbol. Sequences of symbols can be produced by traversing the graph and adding the “new” symbol to the growing sequence.

**Figure 2: A de Bruijn graph of word length 3 for the symbols 1 and 0.**

<img src="../media/debruijn.jpg" width=80% />

*From: https://cameroncounts.wordpress.com/2015/02/28/1247/*

Velvet constructs a de Bruijn graph of the reads. It has 4 symbols (A, C, G and T - N’s are converted to A’s) The word length (or k-mer size) is one of Velvet’s prime parameters.

Velvet is not the only assembly software that works in this manner. Euler, Edena and SOAP de novo are examples of others.

### The Velvet algorithm

#### *Step 1: Hashing the reads.*

* Velvet breaks up each read into k-mers of length k.
  * A k-mer is a k length subsequence of the read.
  * A 36 base pair long read would have 6 different 31-mers.
* The k-mers and their reverse complements are added to a hash table to categorize them.
* Each k-mer is stored once but the number of times it appears is also recorded.
* This step is performed by “velveth” - one of the programs in the Velvet suite.

#### *Step 2: Constructing the de Bruijn graph.*

* Velvet adds the k-mers one-by-one to the graph.
* Adjacent k-mers overlap by k-1 nucleotides.
* A k-mer which has no k-1 overlaps with any k-mer already on the graph starts a new node.
* Each node stores the average number of times its k-mers appear in the hash table.
* Figure 3 shows a section of a de Bruijn graph constructed by Velvet for k=5.
* Different sequences can be read off the graph by following a different path through it. (Figure 3)

**Figure 3: Section of a simple de Bruijn graph of reads with k-mer size 5. Coloured sequences are constructed by following the appropriately coloured line through the graph.** *(Base figure Zerbino et al 2008.)*

<img src="../media/velvet-graph.png" width=80% />

#### *Step 3: Simplification of the graph.*

* Chain merging: When there are two connected nodes in the graph without a divergence, merge the two nodes.
* Tip clipping: Tips are short (typically) chains of nodes that are disconnected on one end. They will be clipped if their length is < 2 x k or their average k-mer depth is much less than the continuing path.
* Bubble removal: Bubbles are redundant paths that start and end at the same nodes (Figure 4.) They are created by sequencing errors, biological variants or slightly varying repeat sequences.
  * Velvet compares the paths using dynamic programming.
  * If they are highly similar, the paths are merged.
* Error removal: Erroneous connections are removed by using a “coverage cutoff”. Genuine short nodes which cannot be simplified should have a high coverage. An attempt is made to resolve repeats using the “expected coverage” of the graph nodes.
* Paired end read information: Velvet uses algorithms called “Pebble” and “Rock Band” (Zerbino et al 2009) to order the nodes with respect to one another in order to scaffold them into longer contigs.

**Figure 4: Representation of “bubbles” in a Velvet de Bruijn graph.**  *(Base figure Zerbino et al 2008.)*

<img src="../media/velvet-bubbles.png" width=80% />

#### *Step 4: Read off the contigs.*

* Follow the chains of nodes through the graph and “read off” the bases to create the contigs.
* Where there is an ambiguous divergence/convergence, stop the current contig and start a new one.

### K-mer size and coverage cutoff values

The size of the k-mers that construct the graph is very important and has a large effect on the outcome of the assembly. Generally, small k-mers create a graph with increased connectivity, more ambiguity (more divergences) and less clear “paths” through the graph. Large k-mers produce graphs with less connectivity but higher specificity. The paths through the graph are clearer but they are less connected and prone to breaking down.

The coverage cutoff c used during the error correction step of Velvet also has a significant effect on the output of the assembly process. If c is too low, the assembly will contain nodes of the graph that are the product of sequencing errors and misconnections. If c is too high, it can create mis-assemblies in the contigs and destroys lots of useful data.

Each dataset has its own optimum values for the k-mer size and the coverage cutoff used in the error removal step. Choosing them appropriately is one of the challenges faced by new users of the Velvet software.

### Velvet Optimiser

The Velvet Optimiser chooses the optimal values for k and c automatically by performing many runs of Velvet (partially in parallel) and interrogating the subsequent assemblies.  It uses different optimisation functions for k and c and these can be user controlled.

It requires the user to input a range of k values to search (to cut down on running time).

## References
http://en.wikipedia.org/wiki/Sequence_assembly

Zerbino DR, Birney E, Velvet: algorithms for de novo short read assembly using de Bruijn graphs, Genome Research, 2008, 18:821-829

Zerbino DR, McEwen GK, Margulies EH, Birney E, Pebble and rock band: heuristic resolution of repeats and scaffolding in the velvet short-read de novo assembler. PLoS One. 2009; 4(12):e8407.

Gladman SL, Seemann T, Velvet Optimiser, http://www.vicbioinformatics.com/software.shtml 2009.
