# FastQC

After sequencing, the reads should be checked for their quality. This tutorial demonstrates how to use the tool called FastQC to examine bacterial paired-end sequence reads from Illumina.

[FastQC website link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Pre-requisites
- mGVL instance
- knowledge: de novo assembly
- knowledge: galaxy or cmdline

****

# FastQC on Galaxy

FIXME: include screenshots

## Start
- open your mGVL galaxy instance

## Input
- required input data: R1 and R2 reads
- Go to: history: import from file [public URL of a history file]
- [FIXME: how to transfer a history to that spot ?]
- [FIXME: At present I am using Pasteurella multocida but would be good to have a data set with contaminants, adapters - how to find?]

## Run
- Tools pane: <ss>NGS Analysis: NGS QC and manipulation: FastQC</ss>.
- Under <ss>Short read data from your current history</ss>: click the multiple files button: select R1 and R2 files. [FIXME: is this correct? you would want to view both?]
- Leave the other settings as they are.
- <ss>Execute</ss>.

## Output
- For each of the inputs (R1 and R2), there will two output files at the top of the history pane (right).
- We will look at the "webpage" output (which is the output displayed graphically).
Click on the eye icon next to the file called <ss>FastQC on data1: Webpage</ss>.
- Look at the following information:  
    - Basic Statistics: total sequences - to understand coverage.
    - Basic Statistics: sequence length - to know how to set k values later.
    - Basic Statistics: %GC - if high, assembly may be more difficult.  
    - Per base sequence quality: might be low at the start and the end, which may support some trimming of reads
    - Per base N content: if high, may indicate poor quality and a need for trimming the Ns.
    - Kmer content: should be no spikes
- Look at the same information for the R2 file.
- [Example of good illumina data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
- [Example of bad illumina data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

****

# FastQC on cmdline

## Start
- on your local machine, make sure XQuartz is installed (but it doesn't have to be open - it will open automatically later).
- in terminal, ssh to your virtual machine with -X and -Y, e.g. `ssh -X -Y ubuntu@111.111.111.111` (the -X and -Y means it will use your local XQuartz to display some files).
- `module load fastqc_dist_0_10_1` [FIXME: will this have a different name? eg just fastqc?]
- navigate to where you want to make a FastQC analysis folder.
- Make a folder: `mkdir fastqc_analyses`
- Move to that folder: `cd fastqc_analyses`

## Input
- [FIXME: get the data into this folder]

## Run
- `fastqc R1reads.fastq` [runs fastqc]
- type `fastqc --help` to see settings that you can change, and defaults
- FIXME: any to change
- FIXME: repeat for R2reads?

## Output
- R1reads_fastqc: folder containing the output e.g. fastqc_report.html
- to view this, type: firefox fastqc_report.html - firefox should open and display the report
(you may get an error message in terminal but ignore this)

****

## More information

FIXME: include these?

- link to a fastqc protocol:
http://vlsci.github.io/lscc_docs/tutorials/assembly/assembly-protocol/#section-1-read-quality-control

- more detailed information:
https://docs.google.com/document/pub?id=16GwPmwYW7o_r-ZUgCu8-oSBBY1gC97TfTTinGDk98Ws

## Next
- Trim reads: Trimmomatic (link)
