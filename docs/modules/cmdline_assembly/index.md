<br>
# Pacbio reads: assembly with command line tools

Keywords: de novo assembly, PacBio, PacificBiosciences, Illumina, command line, Canu, Circlator, BWA, Spades, Pilon, Microbial Genomics Virtual Laboratory

This tutorial demonstrates how to use long Pacbio sequence reads to assemble a bacterial genome and plasmids, including correcting the assembly with short Illumina reads.

## Resources

Tools (and versions) used in this tutorial include:

- canu 1.5 [recently updated]
- infoseq and sizeseq (part of EMBOSS) 6.6.0.0
- circlator 1.5.1 [recently updated]
- bwa 0.7.15
- samtools 1.3.1
- makeblastdb and blastn (part of blast) 2.4.0+
- pilon 1.20

<!-- check these are the correct version - updates? -->

## Learning objectives

At the end of this tutorial, be able to:

1. Assemble, circularise and polish PacBio data for a prokaryote
2. Retrieve small plasmids missed by long read sequencing, using Illumina data
3. Explore the effect of polishing assembled sequences with a different data set.

## Overview

Simplified version of workflow:

![workflow](images/flowchart.png)

## Get data

- Open the command line. <!-- own machine, mGVL or BPA VM -->
- Navigate to or create the directory in which you want to work.
- If you already have the files ready, skip forward to next section, [Assemble](#assemble).

### Find the PacBio files for this sample

- If the files are already on your server, you can symlink by using

```text
ln -s real_file_path [e.g. data/sample_name/pacbio1.fastq.gz] chosen_symlink_name [e.g. pacbio1.fastq.gz]
```

- Alternatively, obtain the input files from elsewhere, e.g. from the BPA portal. (You will need a password.)

- Pacbio files are often stored in the format:
    - <fn>Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq.gz</fn>

- We will use the <fn>longfilename.subreads.fastq.gz</fn> files.

- The reads are usually split into three separate files because they are so large.

- Right click on the first <fn>subreads.fastq.gz</fn> file and "copy link address".

- In the command line, type:

```text
wget --user username --password password [paste link URL for file]
```
- Repeat for the other two <fn>subreads.fastq.gz</fn> files.

### Join Pacbio fastq files

- If the files are gzipped, type:

```text
cat pacbio*.fastq.gz > subreads.fastq.gz
```

- If the files are not gzipped, type:
```text
cat pacbio*.fastq | gzip > subreads.fastq.gz
```

- We now have a file called <fn>subreads.fastq.gz</fn>.

### Find the Illumina files for this sample

- We will also use 2 x Illumina (Miseq) fastq.gz files.
- These are the <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> files.
- Symlink or "wget" these files as described above for PacBio files.
- Shorten the name of each of these files:

```text
mv longfilename_R1.fastq.gz R1.fastq.gz
mv longfilename_R2.fastq.gz R2.fastq.gz
```

### View files

- Type "ls" to display the folder contents.
```text
ls
```
- The 3 files we will use in this analysis are:
    - <fn>subreads.fastq.gz</fn> (the Pacbio reads)
    - <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> (the Illumina reads)

<!--
Find information about read lengths: fq subreads.fastq.gz
- Look at the average and maximum lengths.
-->

- In this tutorial we will use *Staphylococcus aureus* sample 25745.

or: sample 25747 (add details).

### Sample information

Sample 25745:
The sample used in this tutorial is a gram-positive bacteria called *Staphylococcus aureus*. Some strains of *S. aureus* are resistant to antibiotics. For example, methicillin-resistant *S. aureus* (or MRSA) cannot be treated with the antibiotic methicillin (a type of penicillin). The sample used here however is a methicillin-sensitive (MSSA). It was isolated from (human) blood and caused bacteraemia, an infection of the blood.

## Assemble<a name="assemble"></a>

- We will use the assembly software called [Canu](http://canu.readthedocs.io/en/stable/).
- Run Canu with these commands:

```text
canu -p canu -d canu_output genomeSize=2.8m -pacbio-raw subreads.fastq.gz
```

- the first **canu** tells the program to run
- **-p canu**: names prefix for output files ("canu")
- **-d canu_output**: names output directory ("canu_output")
- **genomeSize** only has to be approximate.
    - e.g. *Staphylococcus aureus*, 2.8m
    - e.g. *Streptococcus pyogenes*, 1.8m
- the **reads** can be unzipped or .gz

- Canu will correct, trim and assemble the reads.
- Various output will be displayed on the screen.
- This will take ~ 30 minutes.

### Check the output

Move into the output directory and "ls" to see the output files.

- The <fn>canu.contigs.fasta</fn> are the assembled sequences.
- The <fn>canu.unassembled.fasta</fn> are the reads that could not be assembled.
- The <fn>canu.correctedReads.fasta.gz</fn> are the corrected Pacbio reads that were used in the assembly.
- The <fn>canu.file.gfa</fn> is the graph of the assembly.
- Display summary information about the contigs:
```
infoseq canu.contigs.fasta
```

- (infoseq is a tool from [EMBOSS](http://emboss.sourceforge.net/index.html))

- This will show the number of contigs, e.g.

```text
    - tig00000031	49397
    - tig00000056	2748030
```

This looks like a chromosome of 2.7 million bases and a plasmid of 49 k bases.

[explain some of these features: %GC content, len=49397 reads=61 covStat=3.95 gappedBases=no class=contig suggestRepeat=no suggestCircular=no]

Or: Sample 25747:  tig00000001 2,851,895



### Change Canu parameters if required

- If the assembly is poor with many contigs, re-run Canu with extra sensitivity parameters; e.g.
```text
canu -p prefix -d outdir corMhapSensitivity=high corMinCoverage=0 genomeSize=2.8m -pacbio-raw subreads.fastq
```

## Questions

Q: How do long- and short-read assembly methods differ? A: short reads: De Bruijn graphs; long reads: a move back towards simpler overlap-layout-consensus methods.

Q: Where can we find out the what the approximate genome size should be for the species being assembled? A: NCBI Genomes - enter species name - click on Genome Assembly and Annotation report - sort table by clicking on the column header Size (Mb) - look at range of sizes in this column.

Q: In the assembly output, what are the unassembled reads? Why are they there?

Q: What are the corrected reads? How did canu correct the reads?

Q: Where could you view the output .gfa and what would it show?

## Trim and circularise

### Run Circlator
Circlator identifies and trims overhangs (on chromosomes and plasmids) and orients the start position at an appropriate gene (e.g. dnaA). It takes in the assembled contigs from Canu, as well as the corrected reads prepared by Canu.

Move back into your main folder.

To run:

```text
circlator all --threads 8 --verbose canu_output/canu.contigs.fasta canu_output/canu.correctedReads.fasta.gz circlator_output
```
- **&#45;&#45;threads** is the number of cores
- **&#45;&#45;verbose** prints progress information to the screen

- **<fn>canu_output/canu.contigs.fasta</fn>** is the file path to the input multi-fasta assembly
- **<fn>canu_output/canu.correctedReads.fasta.gz</fn>** is the file path to the corrected Pacbio reads - note, fasta not fastq
- **circlator_output** is the name of the output directory.

Some output will print to screen. When finished, it should say "Circularized x of x contig(s)".

### Check the output

Move into the circlator_output directory and "ls" to list files.

*Were the contigs circularised?* :

```text
less 04.merge.circularise.log
```

- Yes, both were circularised (last column).
- Type "q" to exit.

*Where were the contigs oriented (which gene)?* :

```text
less 06.fixstart.log
```
- Look in the "gene_name" column.
- Contig tig00000056 (the chromosome) has been oriented at tr|A0A090N2A8|A0A090N2A8_STAAU, which is another name for dnaA. <!-- (search swissprot - uniprot.org) -->
- Contig tig00000031 (the plasmid) has been oriented at a gene predicted by prodigal.

*What are the trimmed contig sizes?* :

```text
infoseq 06.fixstart.fasta
```

- tig00000031 25012  (24835 bases trimmed - almost half)[explain why]

- tig00000056 2725223 (22807 bases trimmed)

or: sample 25747: 2,823,331 (28564 bases trimmed)

*Re-name the contigs file*:

- The trimmed contigs are in the file called <fn>06.fixstart.fasta</fn>.
- Re-name it <fn>contig_1_2.fa</fn>:

or: sample 25747: contig1.fa

```text
mv 06.fixstart.fasta contig_1_2.fa
```

to do: rename contig "chromosome.fasta" and in fasta file header here with nano. e.g. chromosome

### Options:

[explain this option better]

**&#45;&#45;b2r_length_cutoff** using approximately 2X average read length (could be omitted at first; if all contigs don't circularise, include this option to see if any improvement)

e.g. setting as: 20000


## Questions

Q: Were all the contigs circularised? Why/why not?

Q: Circlator can set the start of the sequence at a particular gene. Which gene does it use? Is this appropriate for all contigs? A: Uses dnaA for the chromosomal contig. For other contigs, uses a centrally-located gene. However, ideally, plasmids would be oriented on a gene such as repA. It is possible to provide a file to Circlator to do this.


## Find smaller plasmids
Pacbio reads are long, and may have been longer than small plasmids. We will look for any small plasmids using the Illumina reads.

This section involves several steps:

1. Use the multifasta canu-circlator output of trimmed assembly contigs.
2. Map all the Illumina reads against these Pacbio assembled contigs.
3. Extract any reads that *didn't* map and assemble them together: this could be a plasmid, or part of a plasmid.
5. Look for overhang: if found, trim.


**Optional**:

If no overhang found in assembled contigs from the unmapped reads:

6. Search Genbank for any matching proteins: a replication protein found.  
7. Assemble all the Illumina reads and produce an assembly graph.
8. Search the graph for a match to the replication protein and its adjoining regions.
9. Extract this longer sequence from the Illumina assembly: this is the small plasmid.
10. Check for overhang in this plasmid and trim.




### Align Illumina reads to PacBio contigs

- Copy the trimmed, circularised contigs from circlator back to your main folder.
- Align illumina reads to these contigs

- First, index the contigs file. [explain what this does]

```text
bwa index contig_1_2.fa
```

- then, align using bwa mem:

```text
bwa mem -t 8 contig_1_2.fa R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
```
- **bwa mem** is the alignment tool
- **-t 8** is the number of cores
- **contig_1_2.fa** is the input assembly file
- **R1.fastq.gz R2.fastq.gz** are the Illumina reads
- ** | samtools sort** pipes the output to samtools to sort
- **> aln.bam** sends the alignment to the file <fn>aln.bam</fn>


### Extract unmapped Illumina reads

- Index the alignment file
```text
samtools index aln.bam
```
- extract the fastq files from the bam alignment - those reads that were unmapped to the Pacbio alignment - and save them in various "unmapped" files:
```text
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam
```
- **fastq** is a command that coverts a <fn>.bam</fn> file into fastq format
- **-f 4** : only output unmapped reads
- **-1** : put R1 reads into a file called <fn>unmapped.R1.fastq</fn>
- **-2** : put R2 reads into a file called <fn>unmapped.R2.fastq</fn>
- **-s** : put singleton reads into a file called <fn>unmapped.RS.fastq</fn>
- **aln.bam** : input alignment file

We now have three files of the unampped reads:

- <fn> unmapped.R1.fastq</fn>
- <fn> unmapped.R2.fastq</fn>
- <fn> unmapped.RS.fastq</fn>

### Assemble the unmapped reads

- assemble with spades
```text
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq --careful --cov-cutoff auto -o spades_assembly
```

- **-1** is input file forward
- **-2** is input file reverse
- **-s** is unpaired
- **&#45;&#45; careful** : minimizes mismatches and short indels
- **&#45;&#45; cov-cutoff auto** : computes the coverage threshold (rather than the default setting, "off")
- **-o** is the output directory

Look at the output:

```text
cd spades_assembly
infoseq contigs.fasta
```
- shows 11 contigs were assembled, with the max length of 2229.     

<!--
- sort fasta by size of seqs:
```text
sizeseq
input sequence set: contigs.fasta
return longest sequence first [N]: Y
output sequence(s) [contigs.fasta]: sorted_contigs.fasta
```
-->

Print the first row of each seq to see coverage:
```text
grep cov contigs.fasta  
```
- result: NODE_1_length_2229_cov_610.298
    - longest contig is 2229 and high coverage
- all other contigs are short so we will disregard.

<!--
- see if any other nodes have high coverage
    - e.g. NODE_135_length_78_cov_579
- look at the sequence of this contig:
```text
tail sorted_contigs.fasta
```
- This is a homopolymer, so disregard.
-->

or sample 25747: NODE_1_length_2550_cov_496.613
all other nodes are < 650kb

- We will extract the first sequence (NODE_1):
```text
samtools faidx contigs.fasta
samtools faidx contigs.fasta NODE_1_length_2229_cov_610.298 > contig3.fa
```
- this is now saved as <fn>contig3.fa</fn>
- open this file in nano, make the header ">contig3", save (Ctrl-X).

or sample 25747: extracted as plasmid.fa

### Investigate the small plasmid (contig3)

- Blast the start of contig3 against itself
- Take the start of the contig:
```text
head -n 10 contig3.fa > contig3.fa.head
```
- We want to see if it matches the end (overhang)
- Format the assembly file for blast:
```text
makeblastdb -in contig3.fa -dbtype nucl
```
- blast the start of the assembly (.head file) against all of the assembly:
```text
blastn -query contig3.fa.head -db contig3.fa -evalue 1e-3 -dust no -out contig3.bls
```
- look at <fn>contig3.bls</fn> to see hits:
```text
less contig3.bls
```
- the first hit is against itself, as expected
- there are no few further hits, so we assume there is no overhang that needs trimming.
- however, the sequence is likely then to be longer than this.


or: sample 25747:
first hit at start
second hit is at 2474 all the way to the end - 2550.
so trim to 2473

- Index the plasmid.fa file:
```text
samtools faidx plasmid.fa
```
- Trim:
```text
samtools faidx plasmid.fa plasmid:1-2473 > plasmid.fa.trimmed
```
- Open this file in nano and change the header to ">plasmid_trimmed", save.
- We now have a trimmed plasmid.

to do: end up with naming this plasmid.fasta

to do: orient this at repA.



-----------------------------
or continuing on from sample 25745:

```text
less contig3.fa
```
- Copy the sequence
- Go to NCBI: <tt>https://blast.ncbi.nlm.nih.gov/Blast.cgi</tt>; choose blastx
- Paste the sequence from contig3.fa
- Choose genetic code = 11
- Blast
![ncbi blast](images/ncbi_blast_contig3.png)
- This hits a replication (plasmid) protein. Hypothesise that	this is a small plasmid; search for the entire sequence within the assembly of all the Illumina reads (next step).

###Assemble *all* the illumina reads
- Assemble all the Illumina reads with spades (not just those reads that did not map to the Pacbio assembly).

```text
spades.py -1 R1.fastq -2 R2.fastq --careful --cov-cutoff auto -o spades_assembly_all_illumina
```
<!--
Alternatively, if you have spades-fast, you can run with these options:

```text
spades-fast --R1 R1.fastq.gz --R2 R2.fastq.gz --gsize 2.8M --outdir spades_assembly_all_illumina --cpus 8
```
-->
Navigate to the output:

```text
cd spades_assembly_all_illumina
```
- in here is the <fn>assembly_graph.fastg</fn>
- Transfer this file to your local computer (e.g. using the file transfer program [Cyberduck](https://cyberduck.io/?l=en)).
- Examine the assembly in the program [Bandage](https://rrwick.github.io/Bandage/).
    - File: Load graph: <fn>assembly_graph.fastg</fn>
    - In the left hand panel, click "Draw graph"
    - Your assembly graph may look like this:
![bandage pic](images/illumina_assembly_bandage.png)
- Blast the small plasmid sequence in this assembly
    - In the left hand panel: Blast: create/view BLAST search
    - Build blast database
    - Paste in the sequence of contig3
    - Blast
    - The main hit is around node 10.

- Go to the main Bandage window
    - "find nodes" in right hand panel - 10
    - This node is slightly longer: 2373: this could be the plasmid
    - Extract this node in fasta format: In the top panel, go to Output: Save selected node sequences; save as <fn>contig3b.fa</fn>
    - Open this file in nano and change the header to ">contig3b", save.

### Trim small plasmid
- Take the start of the sequence and see if it matches the end:
```text
head -n 10 contig3b.fa > contig3b.fa.head
makeblastdb -in contig3b.fa -dbtype nucl
blastn -query contig3b.fa.head -db contig3b.fa -evalue 1e-3 -dust no -out contig3b.bls
less contig3b.bls
```
- The first hit is against the start of the chromosome, as expected.
- The last hit starts at position 2253; we will trim the plasmid to position 2252
- Index the contig3b.fa file:
```text
samtools faidx contig3b.fa
```
- Trim:
```text
samtools faidx contig3b.fa contig3b:1-2252 > contig3b.fa.trimmed
```
- Open this file in nano and change the header to ">contig3b", save.
- We now have a trimmed contig3b.

### Collect all contigs in one file
```text
cat contig_1_2.fa contig3b.fa.trimmed > all_contigs.fa
```
- See the three contigs and sizes:
```text
infoseq all_contigs.fa
```

or with sample 25747: cat contig1 plasmid1 > contigs.fasta

## Questions

Q: Why is this section so complicated? A: Finding small plasmids is difficult for many reasons! This paper has a nice summary: On the (im)possibility to reconstruct plasmids from whole genome short-read sequencing data. doi: https://doi.org/10.1101/086744

Q: Why can PacBio sequencing miss small plasmids? A: Library prep size selection

Q: We extract unmapped Illumina reads and assemble these to find small plasmids. What could they be missing? A: Repeats that have mapped to the PacBio assembly.

Q: How do you find a plasmid in a Bandage graph? A: It is probably circular, matches the size of a known plasmid, has a rep gene...

Q: Are there easier ways to find plasmids? A: Possibly. One option is the program called Unicycler which may automate many of these steps. https://github.com/rrwick/Unicycler

## Re-orient the plasmids at repA?


##Correct

We will correct the Pacbio assembly with Illumina reads.

<!--
First, we will change some of the mGVL settings so that we can use a worker node with more CPUs.

Type in

```text
sinteractive --cpus=8 --mem=10g
```

(To later exit out of this worker node, if you want to, type in "exit").
--->

Align the Illumina reads (R1 and R2) to the draft PacBio assembly, e.g. <fn>contigs.fasta</fn>:

```text
bwa index contigs.fasta
bwa mem -t 32 contigs.fasta R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
samtools index aln.bam
samtools faidx contigs.fasta
```

- **-t** is the number of cores (e.g. 8)
    - to find out how many you have, grep -c processor /proc/cpuinfo
- now we have an alignment file to use in pilon: <fn>aln.bam</fn>

Look at how the illumina reads are aligned:

```text
samtools tview -p contig1 aln.bam contigs.fasta
```
note: **contig1** is the name of the contig to view; e.g. tig00000000.

[explain settings in tview - how to scroll, find]

Run pilon:

```text
pilon --genome contigs.fa --frags aln.bam --output pilon1 --fix all --mindepth 0.5 --changes --verbose
```
- **&#45;&#45;genome** is the name of the input assembly to be corrected
- **&#45;&#45;frags** is the alignment of the reads against the assembly
- **&#45;&#45;output** is the name of the output prefix
- **&#45;&#45;fix** is an option for types of corrections
- **&#45;&#45;mindepth** gives a minimum read depth to use
- **&#45;&#45;changes** produces an output file of the changes made
- **&#45;&#45;verbose** prints information to the screen during the run
- if you are using pilon on a different machine and you want to specify the number of CPUs, type in **&#45;&#45;threads** number (e.g. 32).

Look at the changes file:

```text
less pilon1.changes
```

<!-- fix this to show this file
- This shows the corrections made by Pilon:

![Pilon](images/pilon.png)
 -->

Look at the fasta file:

```text
less pilon1.fasta
```

Look at the details of the fasta file:

```text
infoseq pilon1.fasta
```

If there are more than 2 changes, run Pilon again, using the pilon1.fasta file as the input assembly, and the Illumina reads to correct.



eg sample 25747: ~ 10 changes. (compared to snippy  ~ 10 changes)
re-ran pilon. only 2 changes made.
saved pilon2.fasta as staph_assembly.fasta


note: annotated with prokka. plasmid only has two proteins. ermC, and a hypothetical protein. protein blast genbank: matches a staph replication and maintenance protein. (rep?)


**Final output:**

- the corrected genome assembly of *Staphylococcus aureus* in .fasta format, containing three contigs: chromosome, large plasmid and small plasmid.  


## Questions

Q: Why don't we correct earlier in the assembly process? A: We need to circularise the contigs and trim overhangs first.

Q: Why can we use some reads (Illumina) to correct other reads (PacBio) ? A: Illumina reads have higher accuracy

Q: Could we just use PacBio reads to assemble the genome? A: Yes, if accuracy adequate.



##Next

**Further analyses:**

- Annotate with Prokka.
- Comparative genomics, e.g. with Roary.

**Links:**

- [Details of bas.h5 files](https://s3.amazonaws.com/files.pacb.com/software/instrument/2.0.0/bas.h5+Reference+Guide.pdf)
- Canu [manual](http://canu.readthedocs.io/en/stable/quick-start.html) and [gitub repository](https://github.com/marbl/canu)
- Circlator [article](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0) and [github repository](http://sanger-pathogens.github.io/circlator/)
- Pilon [article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963) and [github repository](https://github.com/broadinstitute/pilon/wiki)
- Notes on [finishing](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes) and [evaluating](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies) assemblies.

<!--
## Commands summary <a name="summary"></a>

```text
#get the sequencing reads
#joined Pacbio reads: subreads.fastq
#illumina reads: R1.fastq.gz, R2.fastq.gz

#run canu to assemble
canu -p staph -d output genomeSize=2.8m -pacbio-raw subreads.fastq
#option: use sensitivity settings
#output: staph.contigs.fasta

#run circlator to trim and orient
circlator all --verbose ../staph.contigs.fasta ../staph.corrected.reads.fastq.gz circlator_all_output
#option: use branch lengths option
#output: 06.fixstart.fasta; rename as contig_1_2.fa

#find smaller plasmids

#align illumina to Pacbio contigs
bwa index contig_1_2.fa
bwa mem -t 8 contig_1_2.fa R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
#output is aln.bam

#extract unmapped reads
samtools index aln.bam
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam

#assemble these reads with spades
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq -o spades_assembly
#option: use careful and cutoff options
#save any long contigs; e.g. as contig3.fa

#blast this contig to find overhang
head -n 10 contig3.fa > contig3.fa.head
makeblastdb -in contig3.fa -dbtype nucl
blastn -query contig3.fa.head -db contig3.fa -evalue 1e-3 -dust no -out contig3.bls
less contig3.bls
#no overhang found

#assemble all illumina reads
spades-fast --R1 R1.fastq.gz --R2 R2.fastq.gz --gsize 2.8M --outdir spades_fast --cpus 32
#blast the assembly graph with contig3.fa; extract out node and save as contig3b.fa

#trim overhang
head -n 10 contig3b.fa > contig3b.fa.head
makeblastdb -in contig3b.fa -dbtype nucl
blastn -query contig3b.fa.head -db contig3b.fa -evalue 1e-3 -dust no -out contig3b.bls
less contig3b.bls
samtools faidx contig3b.fa
samtools faidx contig3b.fa contig3b:1-2252 > contig3b.fa.trimmed
#open contig3b.fa in nano and shorten header name

#join all contigs: all_contigs.fa

#correct with pilon

#1. correct using trimmed, corrected Pacbio reads
bwa index all_contigs.fa
bwa mem -t 32 all_contigs.fa canu.trimmedReads.fasta.gz | samtools sort > aln.bam
samtools index aln.bam
samtools faidx all_contigs.fa
pilon --genome all_contigs.fa --unpaired aln.bam --output corrected --fix bases --mindepth 0.5 --changes --threads 32 --verbose

#2. correct using illumina reads - use output from 1 as the contigs file
# use --frags not --unpaired

#3. repeat correction using Illumina reads - use output from 2 as the contigs file

#output: corrected genome assembly of Staphylococcus aureus in .fasta format, containing three contigs: chromosome, large plasmid and small plasmid.
```
-->
