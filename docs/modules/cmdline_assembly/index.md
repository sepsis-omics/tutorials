#PacBio reads: assembly with command line tools

This tutorial demonstrates how to use long sequence reads to assemble a bacterial genome.

## Learning objectives

At the end of this tutorial, be able to use command line tools to:

1. assemble long (PacBio) sequence reads
2. trim overhang
3. correct the draft assembly with short (Illumina) reads
4. circularise the corrected assembly
5. &rarr; produce a finished assembly.

##Input files

- Open the mGVL command line
- Navigate to the directory in which you want to work.
- Make a directory called "staph"

```
mkdir staph
cd staph
```
- Make a subdirectory with the sample name "25745"

```
mkdir 25745
cd 25745
```

**Find the PacBio files for this sample**

- Obtain the input files. e.g. from the BPA portal.
- PacBio files are often stored in the format:
Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq.gz
- We will use the <fn>fastq.gz</fn> files.
- The reads are usually split into three separate files because they are so large.
- Right click on the first <fn>subread.fastq.gz</fn> file and "copy link address".

- In the cmd line,

```
wget --user username --password password [paste link URL for file]
```
- repeat for the other two <fn>subread.fastq.gz</fn> files.

**Join PacBio fastq files**

- If the files are gzipped, type:

```
cat filepath/filep0.*.subreads.fastq.gz > subreads.fastq.gz
```

- If the files are not gzipped, type:

```
cat filepath/filep0.*.subreads.fastq | gzip > subreads.fastq.gz
```

- We now have a file called <fn>subreads.fastq.gz</fn>.

**Find the Illumina files for this sample**

- We will also use 2 x Illumina (Miseq) fastq.gz files.
- These are the <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> files.
- Right click on the file name and "copy link address".
- In the cmd line,

```
wget --user username --password password [paste link URL for file]
```
- Repeat for the other read.fastq.gz file.
- Shorten the name of each of these files:

```
mv longfilename_R1.fastq.gz R1.fastq.gz
mv longfilename_R2.fastq.gz R2.fastq.gz
```

**View files**

- Type "ls" to display the folder contents.

```
ls
```

- The 3 files we will use in this analysis are:
    - <fn>subreads.fastq.g</fn> (the PacBio reads)
    - <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> (the Illumina reads)

##Assemble with Canu

- Run canu with these commands:

```Unix
canu -p staph -d staph_outdir corMhapSensitivity=high corMinCoverage=0 genomeSize=2.8m -pacbio-raw subreads.fastq
```

- **staph** is the prefix given to output files
- **staph_outdir** is the output directory
- **corMhapSensitivity** and **corMinCoverage** is a recommended option for sensitivity
- **genomeSize** only has to be approximate. e.g., for staph, 2.8M.


- Note: it may say "Finished..."  but it is probably still running. Type:

```
squeue
```

- This will show you what is running.

- Canu will correct, trim and assemble the reads.
- This will take ~ x minutes.

**Check the output**

```
cd staph_outdir
```

- The <fn>contigs.fasta</fn> are the assembled sequences.
- The <fn>unassembled.fasta</fn> are the reads that could not be assembled.
- The <fn>file.gfa</fn> is the graph of the assembly.

**Look at the contigs file**

- Display summary information about the contigs:

```
fa -f staph.contigs.fasta
```

![contigs_info_screenshot](images/contigs_info.png)

- There are two contigs. Note down their lengths: 2,725,231 and 43,991
- We will split these into two separate fasta files.

**Separate the contigs into single files**

- Index the contigs file:

```
samtools faidx staph.contigs.fasta
```

- this makes an indexed file with the suffix -fai

- send each contig to a new file:

```
samtools faidx staph.contigs.fasta tig00000000 > contig1.fa
samtools faidx staph.contigs.fasta tig00000001 > contig2.fa
```

**What is in the unassembled reads?**

- Are they contaminants?
- Blast against NCBI
    - Use Cyberduck to transfer the file to your local computer
    - Go to NCBI and upload the file
    - (doesn't work with this file - too big)
    - (or: explain how to do so on command line?)

## Trim overhang in the large contig

The assembly graph will be split at an ambiguous/unstable node. However, this area of the graph likely overlaps in the bacterial chromosome, but has not aligned with itself completely. This is called an overhang. We need to identify these overhangs and trim them, for the chromosome and any plamsids.

Contig 1.

Take the first 30,000 bases:

```
head -n 501 contig1.fa > contig1.fa.head
```

- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:

```
formatdb -i contig1.fa -p F
```

- blast the start of the assembly against all of the assembly:

```
blastall -p blastn -i contig1.fa.head -d contig1.fa -e 1e-10 -F F -o contig1.bls
```

- look at contig1.bls to see hits:

```
less contig1.bls
```

The first hit is against the start of the chromosome, as expected.

![screeshot of blast](images/contig1_bls.png)


- This position, near the end of the contig, is where the overhang starts: position 2725159
- check that this match goes to the end of the contig.2725231 - yes.
- we will now trim the contig
- first, index the contig1 fasta file

```
samtools faidx contig1.fa
```

- this makes contigs.fasta.fai
(samtools will find this automatically though in the next step)


```
samtools faidx contig1.fa tig00000000:1-2725158 > contig1.fa.trimmed
```

- open the contig1.fa.trimmed file

```
nano contig1.fa.trimmed
```

delete the header info except contig name (e.g. tig00000000)
exit.






## Trim overhang in the plamsid

contig 2.





make single FASTA file
chromosome
plasmid

##Correct with Pilon

canu doesn't do any polishing, so there will be thousands of errors, almost all insertions.
input: draft pacbio assembly (overhang trimmed) and illumina reads (in bam format)
output: corrected assembly


**Align Illumina reads => bam**

Illumina reads are shorter but have fewer errors (and/or there are just more of them? so consensus is more accurate)

These can be aligned to the PacBio assembly to correct any errors.

- trimmed contigs + Illumina R1 & R2 => align with BWA-MEM => hybrid contigs => samtools sort => bam
- samtools index bam

```Unix
bwa index contigs.fasta
bwa mem -t 8 contigs.fasta r1s r2s | samtools sort > aln.bam
samtools index aln.bam
samtools faidx contigs.fasta
```

-t is the number of cores (e.g. 8)
to find out how many you have,
grep -c processor /proc/cpuinfo


**Look at the assembly graph**

- Copy the <fn>file.gfa</fn> to your local computer
    - e.g. using a file transfer program such as Cyberduck
    - Install Cyberduck
    - Open Cyberduck
    - click on "open connection"
    - choose SFTP from drop down menu
    - server = your virtual machine IP address (e.g. abrpi.genome.edu.au)
    - username = your username
    - password = your password
    - You should now see a window showing the folders and files on your virtual machine.
    - You can drag and drop files into the preferred folder.

- Examine the assembly in the program Bandage.
    - Install Bandage.
    - File: Load graph: <fn>file.gfa</fn>
    - In the left hand panel, click "Draw graph"
    - How many contigs? Sizes?


**Correct**

contigs & bam => pilon => corrected contigs

pilon

--genome contigs.fasta
--frags aln.bam [that we just created]
--prefix corrected
--changes
--threads 72
--verbose
[does have some options for mindepth - can change these if we want]
--minmq 60  [min mapping aln quality - default 0 is silly- bwa-mem uses 60, bowtie uses 42; 60 should remove most multi mapping stuff]

don't run on things with lots (10+) of contigs, but 1 or 2 is ok

```bash
pilon --genome contigs.fasta --frags aln.bam --output corrected --changes --threads 8 --verbose --minmq 60
```

output:
lots of files. look at .pilon.fasta
so cf to input file, this one could be longer
can look at the changes file

less corrected.changes


a lot are single, but there are some bigger ones.

=> corrected assembly


(Correct again?)

##Circularise
corrected contigs => circlator => circular genome starting at dnaA

```Unix
circlator fixstart trimmed.fa outprefix
```


orients at DNAa

=> circular, corrected assembly


for plasmid:
circlator fixstart --genes_fa filename contig outprefix

need to give it the file (DNA?) of repA
need to download


## Next

Annotate - in Galaxy, obtain this file

## Links

[Canu manual](http://canu.readthedocs.io/en/stable/quick-start.html)

[Canu code](https://github.com/marbl/canu)

[Pilon article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963)

[Pilon on github](https://github.com/broadinstitute/pilon/wiki)

Circlator
http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0
http://sanger-pathogens.github.io/circlator/

Finishing

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies

Files

reference to bas.h5 files (details) https://s3.amazonaws.com/files.pacb.com/software/instrument/2.0.0/bas.h5+Reference+Guide.pdf
