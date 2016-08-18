#PacBio reads: assembly with command line tools

This tutorial demonstrates how to use long sequence reads to assemble a bacterial genome.

## Learning objectives

At the end of this tutorial, be able to use command line tools to:

1. assemble long (PacBio) sequence reads
2. trim overhang
3. correct the draft assembly with short (Illumina) reads
4. circularise the corrected assembly.

=> finished assembly

##Input files

- Open the mGVL command line
- Navigate to the directory in which you want to work.

You will also need to know the location of the sequences (the file path).
TBA: filepath of E. coli tutorial data

PacBio files are often stored in the format:

Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq

The PacBio sequences are usually stored in three fastq files.

This command takes the three FASTQ files, pipes them to gzip, and stores the gzipped concatenated sequences in the file subreads.fastq.gz

Type in:

```Unix
cat FILEPATH/*.fastq | gzip > subreads.fastq.gz
```

##Assemble with Canu

Type in:

```Unix
canu -p ecoli -d outdir corMhapSensitivity=high corMinCoverage=0 genomeSize=2.8m -pacbio-raw subreads.fastq.gz
```

- **ecoli** is the prefix given to output files (check)
- **outdir** is the output directory
- **corMhapSensitivity** and **corMinCoverage** is a recommended option for (to be explained)
- **genomeSize** only has to be approximate

Canu will correct, trim and assemble the reads.

**Check output**

Look at the outdir.

- contigs.fasta : assembled sequences
- unassembled.fasta
- file.gfa : the graph of the assembly (could examine in Bandage. should show...?)

**Look at the contigs file**


```Unix
fa -f contigs.fasta
```

(Is "fa" available during tutorial? )

e.g. 2 contigs (1 chromosome, 1 plasmid)

what else to look for?

**What is in the unassembled reads**

- are they contaminants? blast against ncbi?
(explain how to do so on command line)

output => assembled contigs

## Trim overhang
The assembly graph will be split at an ambiguous/unstable node. However, this area of the graph likely overlaps in the bacterial chromosome, but has not aligned with itself completely. This is called an overhang.

To find out whether the ends match, blast the first bit of the assembly with the last bit. e.g. the first 60kb.

Find where these should properly overlap, and trim. (explain how)

- take the first 501 lines (there are 60 bases per line) and send them to a new file
- this is the start of the assembly. want to see if it matches the end (overhang)

```Unix
head -n 501 ass.contigs.fasta > start.fasta
```

- format the assembly file for blast
```Unix
formatdb -i ass.contigs.fasta -p F
```

- blast the start against the whole assembly, send to file called start.bls

```Unix
blastall -p blastn -i start.fasta -d ass.contigs.fasta -e 1e-100 -F F -o start.bls
```

- look at start.bls to see hits
- note where the overlap starts

```Unix
less start.bls
```

- use the position at the start of the overhang; trim

```Unix
fa-subseq.pl -end 1825546 ass.contigs.fasta > trimmed.fa
```

- open the trimmed.fa file and remove the header info except contig name

```Unix
nano trimmed.fa
```
output => trimmed contigs


##Correct with Pilon

canu doesn't do any polishing, so there will be thousands of errors, almost all insertions.


**Align Illumina reads**

TBA: get the Illumina seqs into this folder

Illumina reads are shorter but have fewer errors (and/or there are just more of them? so consensus is more accurate)

These can be aligned to the PacBio assembly to correct any errors.

- trimmed contigs + Illumina R1 & R2 => align with BWA-MEM => hybrid contigs => samtools sort => bam
- samtools index bam

```Unix
bwa index contigs.fasta
bwa mem -t 72 contigs.fasta r1s r2s | samtools sort > aln.bam
samtools index aln.bam
samtools faidx contigs.fasta
```

**Correct**

hybrid contigs & bam => pilon => corrected contigs



Pilon
- draft assembly (contigs) + read alignments (bam file)
- if evidence that draft is wrong, change
- how: looks at each base position and all the aligned reads. ignores ends of reads as less reliable. measures evidence for each position using count of reads + base quality from seq instrument + alignment mapping quality.
- the either keeps base in assembly, or changes it, or notes that it is     ambiguous, or notes there is not enough info.
    - ambigous: due to diploid; or difficult-to-seq region meaning lots of errors; or collapsed repeat (i.e. two repeats being considered as one).
- pilon also has ways to identify local problem areas and uses trusted alignments either side to try to draw in better seqs to fix.
- output => fasta file of improved assembly

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

output:
lots of files. look at .pilon.fasta
so cf to input file, this one could be longer
can look at the changes file
a lot are single, but there are some bigger ones.

=> corrected assembly

##Circularise
corrected contigs => circlator => circular genome starting at dnaA

circlator fixstart trimmed.fa dnaa

orients at DNAa

=> circular, corrected assembly

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
