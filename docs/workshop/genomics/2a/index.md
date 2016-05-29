# Assembly using Spades

## Background
Spades is one of a number of *de novo* assemblers that use short read sets (Illumina Reads) as an input and the method for assembly is based on de Bruijn graphs. For information about Spades see this [link](http://bioinf.spbau.ru/spades). A protocol for assembling with Velvet (another *de novo* assembler) is available [here](https://docs.google.com/document/d/1xs-TI5MejQARqo0pcocGlymsXldwJbJII890gnmjI0o/pub).

In this activity, we will perform a *de novo* assembly of a short read set (from an Illumina sequencer) using the Spades assembler. The output from Spades that we are interested in is a multifasta file that contains the draft genome sequence.

The read set for today is from an imaginary *Staphylococcus aureus* bacterium with a miniature genome.

We have a closed, annotated genome sequence for a closely related *wildtype* strain.

The whole genome shotgun method used to sequence our mutant strain The read set was produced on an Illumina DNA sequencing instrument.

-   The reads are paired-end
-   Each read is 150 bases (before trimming)
-   The number of bases sequenced is equivalent to 19x the genome sequence of the wildtype strain. (Read coverage 19x - rather low!)

## Learning objectives
At the end of this tutorial you should be able to:

1. import data into Galaxy  
2. view the files
3. evaluate the read quality
4. assemble the reads using Spades, and
5. examine the output assembly.

## Login to Galaxy
-  **Go to the Galaxy Page**: Web address: <http://43.240.98.1/galaxy>  (FIXME: this one?)
- [Remind me how to logon](https://docs.google.com/document/d/1LAQvhIG8s-vv6T14bb8lGRkmoNha7E3bHf9kAgUwMs0/pub)

## Import data
- Click on the **Analyze Data** menu at the top of the page.    
- Click on the History menu button the ![history button](./images/image02.png) on the top right of the history pane.
- Click **Import from File** (at the bottom of the list).  
- A new page will appear with a text box for the URL of the history to import.  
- Copy the following URL into the text box: <http://43.240.98.1/public/dieter/Galaxy-History-Colombiaworkshopstart.tar.gz>  
- Click **Submit**.  
- Galaxy will download the data files from the internet and will be available as an additional history (takes about one minute).  
- **To make the newly imported history appear as the current history**    
  - Click on the View all Histories button ![Histories button](./images/image01.png) (the on the top right of the history pane).  
  - If the history has finished downloading it will appear as **imported from archive: Colombia_workshop_start**
- Click on the ![Switch button](./images/image06.png) button above the **imported from archive:Colombia_workshop_start** then the ![Done button](./images/image05.png)button.
- You should now have four files in the history pane as follows:

![Files in history](./images/image07.png)

## View files
All the files are are text files.

- <fn>mutant_R1.fastq</fn> and <fn>mutant_R2.fastq</fn>: a paired-end read set  
- <fn>wildtype.fna</fn>: a file that contains the genome sequence of the wildtype strain in fasta format (a header line, then the nucleotide sequence of the genome)
- <fn>wildtype.gff</fn>: a file that contains the genome sequence of the wildtype strain in general feature format. (a list of features - one feature per line, then the nucleotide sequence of the genome)

Look at the contents of these files

- Click on the View Data button (the ![Eye icon](./images/image04.png)) next to each of the files in turn.
- Brief Discussion about the GFF format (FIXME: add)

![GFF format](./images/image08.png)

## Evaluate the input reads

Questions you might ask about your input reads include:

- How good is my read set?
- Do I need to ask for a new sequencing run?  
- Is it suitable for the analysis I need to do?

We will evaluate the input reads using the FastQC tool. [FastQC website link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- This runs a standard series of tests on your read set and returns a relatively easy to interpret report.
- We will use the FASTQC tool in Galaxy to evaluate the quality of one of our fastq files.
- Open the FASTQC tool interface (on the tool pane - **NGS: QC and Manipulation -> FastQC: Comprehensive QC)**
- Select mutant_R1.fastq
- **Execute**
- Once finished, examine the output (Hint:![Eye icon](./images/image04.png)). It has a summary at the top of the page and a number of graphs.

Some of the important outputs of FastQC for our purposes are:

-   Read length: will be important in setting maximum k-mer size value for assembly
-   Quality encoding type: important for quality trimming software
-   % GC: high GC organisms don’t tend to assemble well and may have an uneven read coverage distribution.
-   Total number of reads: gives you an idea of coverage.
-   Dips in quality near the beginning, middle or end of the reads: determines possible trimming/cleanup methods and parameters and may indicate technical problems with the sequencing process/machine run.
-   Presence of highly recurring k-mers: may point to contamination of reads with barcodes, adapter sequences etc.
-   Presence of large numbers of Ns in reads: may point to poor quality sequencing run. You need to trim these reads to remove Ns.

We won’t be doing anything to these data to clean it up as there isn’t much need. Therefore we will get on with the assembly!

## Assemble reads with Spades

- We will perform a *de novo* assembly of the mutant fastq reads into long contiguous sequences (in fasta format.)

- Spades produces both contigs and scaffolds. Ask your demonstrator if you would like to know the difference between contigs and scaffolds.

- Open the Spades assembler tool interface (on the tool pane - **NGS: Assembly -> spades**)
- Set the following parameters:

    - **Run only Assembly**: *Yes*  
    - **Kmers to use separated by commas:** *33,55,91*  no spaces  
    - **Coverage cutoff:** *auto*  
    - **Forward reads:** <fn>mutant_R1.fastq</fn>  
    - **Reverse reads:** <fn>mutant_R2.fastq</fn>  

- Your tool interface should look like this:

![Spades interface](./images/image03.png)

-  Click **Execute**

## Examine the output

- Galaxy is now running Spades on the reads for you.
- When it is finished, you will have 5 new files in your history.
- Fasta files of the resulting contigs and scaffolds, some statistics on each and the SPAdes logfile.
- Click on the View Data button ![Eye icon](./images/image04.png) on each of the files.
- Note that the short reads have been assembled into much longer contigs.
- The stats files will give you the length of each of the contigs.
