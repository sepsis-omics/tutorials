# Assembly with Spades in Galaxy

FIXME: This tutorial includes the Workshop 2a "Assembly with Spades" but also some extra info:

- section: Pre-requisites
- section: How does Spades work  
- more detail on output files
- section: Filter output
- section: questions  


## Background
Spades is one of a number of *de novo* assemblers that use short read sets as input (e.g. Illumina Reads), and the assembly method is based on de Bruijn graphs. For information about Spades see this [link](http://bioinf.spbau.ru/spades). A protocol for assembling with Velvet (another *de novo* assembler) is available [here](https://docs.google.com/document/d/1xs-TI5MejQARqo0pcocGlymsXldwJbJII890gnmjI0o/pub).

In this activity, we will perform a *de novo* assembly of a short read set (from an Illumina sequencer) using the Spades assembler. The output from Spades that we are interested in is a multifasta file that contains the draft genome sequence.

The read set for today is from an imaginary *Staphylococcus aureus* bacterium with a miniature genome.

We have a closed, annotated genome sequence for a closely related *wildtype* strain.

The whole genome shotgun method used to sequence our mutant strain read set was produced on an Illumina DNA sequencing instrument.

-   The reads are paired-end
-   Each read is 150 bases (before trimming)
-   The number of bases sequenced is equivalent to 19x the genome sequence of the wildtype strain. (Read coverage 19x - rather low!).

## Learning objectives
At the end of this tutorial you should be able to:

1. import data into Galaxy  
2. view the files
3. evaluate the read quality
4. assemble the reads using Spades, and
5. examine the output assembly.

##Pre-requisites

- Galaxy
- de novo assembly
- QC
- Trimming

## Login to Galaxy
-  Go to this Galaxy address: <http://43.240.98.1/galaxy>  (FIXME: or alternative)
- [Remind me how to logon.](https://docs.google.com/document/d/1LAQvhIG8s-vv6T14bb8lGRkmoNha7E3bHf9kAgUwMs0/pub)
FIXME: note this contains the same galaxy address as above - change?

## Import data
- Click on the <ss>Analyze Data</ss> menu at the top of the page.    
- Click on the <ss>History options</ss> button the ![history button](./images/image02.png) on the top right of the history pane.
- Click <ss>Import from File</ss> (at the bottom of the list).  
- A new page will appear with a text box for the URL of the history to import.  
- Copy the following URL into the text box: <http://43.240.98.1/public/dieter/Galaxy-History-Colombiaworkshopstart.tar.gz>  
- Click <ss>Submit</ss>.  
- Galaxy will download the data files from the internet and will be available as an additional history (takes about one minute).  
- To view this new history, click the <ss>View all histories</ss> button![Histories button](./images/image01.png) (top right of the history pane).  
  - If the history has finished downloading it will appear as <fn>imported from archive: Colombia_workshop_start</fn>
- Click on the ![Switch button](./images/image06.png) button above the <fn>imported from archive:Colombia_workshop_start</fn> then the ![Done button](./images/image05.png) button.
- You should now have four files in the history pane as follows:

![Files in history](./images/image07.png)

## View files
All the files are text files.

- <fn>mutant_R1.fastq</fn> and <fn>mutant_R2.fastq</fn>: a paired-end read set  
- <fn>wildtype.fna</fn>: a file that contains the genome sequence of the wildtype strain in fasta format (a header line, then the nucleotide sequence of the genome)
- <fn>wildtype.gff</fn>: a file that contains the genome sequence of the wildtype strain in general feature format. (a list of features - one feature per line, then the nucleotide sequence of the genome)

Look at the contents of these files

- Click on the View Data button (the ![Eye icon](./images/image04.png)) next to each of the files in turn.
- Brief Discussion about the GFF format (FIXME: add?)

![GFF format](./images/image08.png)

## Evaluate the input reads

Questions you might ask about your input reads include:

- How good is my read set?
- Do I need to ask for a new sequencing run?  
- Is it suitable for the analysis I need to do?

We will evaluate the input reads using the FastQC tool.

- This runs a standard series of tests on your read set and returns a relatively easy to interpret report.
- We will use the FASTQC tool in Galaxy to evaluate the quality of one of our fastq files.
- Go to <ss>Tools &rarr; NGS:Analysis &rarr; NGS: QC and Manipulation &rarr; FastQC</ss>
- Select <fn>mutant_R1.fastq</fn>
- <ss>Execute</ss>
- Once finished, examine the output called <fn>FastQC on data1:webpage</fn> (Hint:![Eye icon](./images/image04.png)). It has a summary at the top of the page and a number of graphs.

Some of the important outputs of FastQC for our purposes are:

-   <ss>Basic Statistics: Sequence length</ss>: will be important in setting maximum k-mer size value for assembly
-   <ss>Basic Statistics: Encoding</ss>: Quality encoding type: important for quality trimming software
-   <ss>Basic Statistics: % GC</ss>: high GC organisms don’t tend to assemble well and may have an uneven read coverage distribution.
-   <ss>Basic Statistics: Total sequences</ss>: Total number of reads: gives you an idea of coverage.
-   <ss>Per base sequence quality</ss>: Dips in quality near the beginning, middle or end of the reads: determines possible trimming/cleanup methods and parameters and may indicate technical problems with the sequencing process/machine run.
-   <ss>Per base N content</ss>: Presence of large numbers of Ns in reads: may point to poor quality sequencing run. You would need to trim these reads to remove Ns.
-   <ss>Kmer content</ss>: Presence of highly recurring k-mers: may point to contamination of reads with barcodes, adapter sequences etc.

Although we have warnings for two outputs (per base sequence content; Kmer content), we can ignore these for now. For a fuller discussion of FastQC outputs and warnings, see the [FastQC website link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), including the section on each of the output [reports](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/), and examples of ["good"](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and ["bad"](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) Illumina data. We won’t be doing anything to these data to clean it up as there isn’t much need. Therefore we will get on with the assembly!

##How does Spades work?

1. As with several other de novo assembly programs (e.g. Velvet) Spades uses an algorithm based on [de Bruijn graphs.](http://www.nature.com/nbt/journal/v29/n11/full/nbt.2023.html) Such graphs use sub-lengths of sequence reads to build an overall genome assembly. The span of the sub-length is called a k-mer, where "k" is the number of nucleotides (e.g. k=21). The user chooses three values of k and Spades makes three assemblies based on these.

2. For the first value of k, each read is broken into as many fragments as possible. For example, if the input read is 22 nucleotides long, and the chosen value of k is 21, then there are two possible fragments (positions 1-21 and 2-22).

3. One randomly-chosen fragment becomes the first node on the de Bruijn graph.

4. A second fragment is connected to this node if it overlaps.

5. Repeat until all fragments are connected. Output &rarr; de Bruijn graph.

6. Find a connected pathway through this graph. Output &rarr; a pathway (sequence) known as a contig. Because of poor or incorrect sequencing, not all the fragments can be joined together. There will be several de Bruijn graphs and so several contigs, usually of different sizes.

7. Repeat these steps for a further two values of k (e.g. k = 33, k = 55). Output &rarr; Three (sets of) contigs.

8. Merge the three (sets of) contigs to get one. Output &rarr; one set of contigs.

9. For paired-end reads (as in this tutorial), the two reads are sequenced from each end of a longer DNA fragment. The middle part of the fragment is not sequenced, but information about the distance between the reads can be used by Spades to join contigs into larger sequences, called scaffolds. Output &rarr; one set of scaffolds.

10. To fix any errors map the original sequence reads onto the scaffolds with the program BWA. Output &rarr; assembled genome.

## Assemble reads with Spades

- We will perform a *de novo* assembly of the mutant fastq reads into long contiguous sequences (in fasta format.)

- Spades produces both contigs and scaffolds. Ask your demonstrator if you would like to know the difference between contigs and scaffolds.

- Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: Assembly &rarr; spades</ss>
- Set the following parameters:

    - <ss>Run only Assembly</ss>: *Yes*  
    - <ss>Kmers to use separated by commas:</ss> *33,55,91*  no spaces  
    - <ss>Coverage cutoff:</ss> *auto*  
    - <ss>Files &rarr; Forward reads:</ss> <fn>mutant_R1.fastq</fn>  
    - <ss>Files &rarr; Reverse reads:</ss> <fn>mutant_R2.fastq</fn>  

- Your tool interface should look like this:

![Spades interface](./images/image03.png)

-  Click <ss>Execute</ss>

## Examine the output

- Galaxy is now running Spades on the reads for you.
- When it is finished, you will have five new files in your history.  
- <fn>contig stats</fn>: There are x contigs. Look at the variation in length and coverage. A short contig with high coverage could be a result of contamination, a collapsed repeat, or a plasmid.
- <ss>contigs</ss>: Each contig is listed, followed by its sequence in fasta format.
- <ss>scaffold stats</ss>: There are x scaffolds.
- <ss>scaffolds</ss>: Each scaffold is listed, followed by its sequence in fasta format.
- <ss>log</ss>: The specific actions performed in the analysis.
- Click on the View Data button ![Eye icon](./images/image04.png) on each of the files.
- Note that the short reads have been assembled into much longer contigs.
- (However, in this case, the contigs have not been assembled into larger scaffolds.)
- The stats files will give you the length of each of the contigs.

## Filter output

Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: Assembly &rarr; Filter SPAdes output</ss>.

This is a quick way to discard contigs that are too short (e.g., they might be contamination) or contigs that do not have enough coverage (e.g., they might be too unreliable).

- Under <ss>Sequences</ss>, choose the contigs fasta file.
- Under <ss>Contig stats</ss> choose the contigs stats file. Change the cut-off values for length and coverage or leave them as they are.
- For <ss>Save filtered-out sequences?</ss> click <ss>Yes</ss>.
- Click <ss>Execute</ss>. A new fasta file with only the filtered sequences will be saved in the right-side history pane.

## Questions

<details> <summary>
How does SPAdes differ from other genome assembly programs?</summary>
It uses multiple values of k in de Bruijn graphs. Larger fragment sizes will more accurately position sections of duplicated DNA (repeats), but these larger fragments will only overlap well in densely-sequenced (high-coverage) areas of the genome. Because bacterial genomes may have low-coverage regions, using smaller fragments can increase the potential for overlaps (joins) in these low-coverage regions. Using a range of fragment sizes will therefore get the benefit from both approaches. More information [here](http://thegenomefactory.blogspot.com.au/2013/08/how-spades-differs-from-velvet.html>).</details>

<details> <summary>
How do I choose values of k?  </summary>
The k values need to be odd numbers, and shorter than the read lengths.  A good strategy could be to choose some that are small, medium and large. e.g. if the read is 150 nucleotides, k values could be 33, 55, 91. There is no absolute rule; rather, the aim is to get a good spread of k values across the read length. </details>

<details> <summary>
What can I do with my assembled genome?</summary>
This tutorial used a subset of a real dataset, so is not a complete genome (is it?). You could re-try it with short reads from a whole genome, at NCBI SRA. You can [annotate] (describe) the genomic features such as genes or [compare] it to other genomes to see variation in structure.  </details>

## What Next?
Annotate the genome, e.g. with [Prokka](../anno/prokka.md).
