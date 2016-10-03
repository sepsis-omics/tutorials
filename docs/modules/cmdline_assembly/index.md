#PacBio reads: assembly with command line tools
This tutorial demonstrates how to use long PacBio sequence reads to assemble a bacterial genome and plasmids, including correcting the assembly with short Illumina reads. The last section contains a [summary](#summary) of the main commands used.

## Learning objectives
At the end of this tutorial, be able to use command line tools to produce a bacterial genome assembly using the following workflow:

1. get data
2. assemble long (PacBio) reads
3. trim overhangs and circularise
4. search for smaller plasmids
5. correct with short (Illumina) reads

##Get data
- Open the mGVL command line
- Navigate to or create the directory in which you want to work.
e.g.
```
mkdir staph
cd staph
```
###Find the PacBio files for this sample
- If the files are already on your server, you can symlink by using
```
ln -s real_file_path chosen_symlink_name
```
- Alternatively, obtain the input files from elsewhere, e.g. from the BPA portal. (You will need a password.)
- PacBio files are often stored in the format: Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq.gz
- We will use the <fn>longfilename.subreads.fastq.gz</fn> files.
- The reads are usually split into three separate files because they are so large.
- Right click on the first <fn>subreads.fastq.gz</fn> file and "copy link address".
- In the command line, type:
```bash
wget --user username --password password [paste link URL for file]
```
- Repeat for the other two <fn>subreads.fastq.gz</fn> files.

###Join PacBio fastq files
- If the files are gzipped, type:
```
cat filepath/filep0.*.subreads.fastq.gz > subreads.fastq.gz
```
- If the files are not gzipped, type:
```
cat filepath/filep0.*.subreads.fastq | gzip > subreads.fastq.gz
```
- We now have a file called <fn>subreads.fastq.gz</fn>.
###Find the Illumina files for this sample
- We will also use 2 x Illumina (Miseq) fastq.gz files.
- These are the <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> files.
- Right click on the file name and "copy link address".
- In the command line, type:
```
wget --user username --password password [paste link URL for file]
```
- Repeat for the other read.fastq.gz file.
- Shorten the name of each of these files:
```
mv longfilename_R1.fastq.gz R1.fastq.gz
mv longfilename_R2.fastq.gz R2.fastq.gz
```
###View files
- Type "ls" to display the folder contents.
```
ls
```
- The 3 files we will use in this analysis are:
    - <fn>subreads.fastq.gz</fn> (the PacBio reads)
    - <fn>R1.fastq.gz</fn> and <fn>R2.fastq.gz</fn> (the Illumina reads)

- In this tutorial we will use *Staphylococcus aureus* sample 25745.

##Assemble
- We will use the assembly software called [Canu](http://canu.readthedocs.io/en/stable/).
- Run Canu with these commands:
```
canu -p staph -d output genomeSize=2.8m -pacbio-raw subreads.fastq
```
- **staph** is the prefix given to output files
- **output** is the name of the output directory
- **genomeSize** only has to be approximate.
    - e.g. *Staphylococcus aureus*, 2.8m
    - e.g. *Streptococcus pyogenes*, 1.8m
- the **reads** can be unzipped or .gz

- Canu will correct, trim and assemble the reads.
- This will take ~ 30 minutes.
- Type
```
squeue
```
to check that the job is running.

###Check the output
```
cd output
```
- The <fn>staph.contigs.fasta</fn> are the assembled sequences.
- The <fn>staph.unassembled.fasta</fn> are the reads that could not be assembled.
- The <fn>staph.correctedReads.fasta.gz</fn> are the corrected pacbio reads that were used in the assembly.
- The <fn>staph.file.gfa</fn> is the graph of the assembly.
- Display summary information about the contigs:
```
fa -f staph.contigs.fasta
```

- This will show the number of contigs, e.g.
    - tig00000000	dna	2746242
    - tig00000001	dna	48500

###Change Canu parameters if required
- If the assembly is poor with many contigs, re-run Canu with extra sensitivity parameters; e.g.
```
canu -p prefix -d outdir corMhapSensitivity=high corMinCoverage=0 genomeSize=2.8m -pacbio-raw subreads.fastq
```

## Trim and circularise

### Run Circlator
Circlator identifies and trims overhangs (on chromosomes and plasmids) and orients the start position at an appropriate gene (e.g. dnaA). It takes in the assembled contigs from Canu, as well as the corrected reads prepared by Canu.

To run:
```
circlator all --verbose ../staph.contigs.fasta ../staph.corrected.reads.fastq.gz circlator_all_output
```
**circlator_all_output** is the name of the output directory.

Check the output: contig sizes, whether they were circularised, trimmed and the start position chosen.
```
less 04.merge.circularise.log
less 06.fixstart.log
```
The trimmed contigs are in the file called <fn>06.fixstart.fasta</fn>. Re-name it <fn>contig_1_2.fa</fn>:
```
mv 06.fixstart.fasta contig_1_2.fa
```

<!--
###Separate the contigs into single files
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
- change contig names:
```
nano filename.fa
```
- change header to >contig1 or >contig2
- save
- We now have two files:
    - <fn>contig1.fa</fn> : the chromosome
    - <fn>contig2.fa</fn> : a plasmid

-->

<!-- ###What is in the unassembled reads?

- Are they contaminants?
- Blast against NCBI
    - Use Cyberduck to transfer the file to your local computer
    - Go to NCBI and upload the file
    - (doesn't work with this file - too big)
    - (or: explain how to do so on command line?)
-->


<!--## Alternatives for assembly
- canu without sensitivity settings
- HGAP2

##Use Quiver/Arrow here?

to correct assembly with the raw reads.
-->

<!---
## Trim - manuallly

The bacterial chromosome and plasmids are circular.

The canu assembly graph will be split at an ambiguous/unstable node. However, this area of the graph likely overlaps in the bacterial chromosome, but has not aligned with itself completely. This is called an overhang. We need to identify these overhangs and trim them, for the chromosome and any plamsids.
### Chromosome: identify overhang
- Take the first 30,000 bases of contig1.
```
head -n 501 contig1.fa > contig1.fa.head
```
- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:
```
formatdb -i contig1.fa -p F

change to:

makeblastdb -in contig1.fa -dbtype nucl
```
- blast the start of the assembly (.head file) against all of the assembly:
```
blastall -p blastn -i contig1.fa.head -d contig1.fa -e 1e-10 -F F -o contig1.bls

change to:

blastn -query contig1.fa.head -db contig1.fa -evalue 1e-10 -dust no -out contig1.bls
```
- look at <fn>contig1.bls</fn> to see hits:
```
less contig1.bls
```
- The first hit is against the start of the chromosome, as expected.
- To find the next hit, type:
```
/ Score
```
then
```
n
```
to scroll through the remaining hits.

- The next hit is from position 2,725,223 until the end of the chromosome (position 2,746,242).

![screeshot of blast](images/blast_chrm.png)

- Therefore, the overhang starts at position 2,725,233 and we can trim the contig to this position minus 1.
### Chromosome: trim
- First, index the <fn>contig1.fa</fn> file
```
samtools faidx contig1.fa
```
- this makes an index file called <fn>contig1.fa.fai</fn>
- next, extract all the sequence except for the overhang. (We don't have to specify the name of the *index* file; it will be found automatically):
```
samtools faidx contig1.fa contig1:1-2725222 > contig1.fa.trimmed
```
- open the <fn>contig1.fa.trimmed</fn> file
```
nano contig1.fa.trimmed
```
- delete the header info except contig name (e.g. contig1)
- save, exit.
- we now have a trimmed contig1.
### Plasmid: examine
- Contig 2 is 48,500 bases.
- This seems long for a plasmid in this species.
- Do a dot plot to examine the sequence.
###Transfer file to local computer
- Use Cyberduck to copy <fn>contig2.fa</fn> to your local computer
    - Install Cyberduck
    - Open Cyberduck
    - click on "open connection"
    - choose SFTP from drop down menu
    - server = your virtual machine IP address (e.g. abrpi.genome.edu.au)
    - username = your username
    - password = your password
    - You should now see a window showing the folders and files on your virtual machine.
    - You can drag and drop files into the preferred folder.
###View dot plot
- Install Gepard
- Open Gepard
- Sequences - sequence 1 - select file - <fn>contig2.fa</fn>
- Sequences - sequence 2 - select file - <fn>contig2.fa</fn>
- Create dotplot
![Gepard dotplot screenshot](images/gepard_contig2.png)
- the sequence starts to repeat at position 24880
- there are probably two plasmids combined into one contig in tandem.
- the whole length is 48500
- we will cut 10,000 from each end

###Trim plasmid
- use samtools to extract the region of the contig except for 10k on each end
```
samtools faidx contig2.fa
samtools faidx contig2.fa contig2:10000-38500 > contig2.fa.half
```
- we now have a ~ 28500 bases plasmid
```
nano contig2.fa.half
```
- delete the rest of the name after >contig2
###Trim overhang in the plasmid

- Take the first x bases:
```
head -n 10 contig2.fa.half > contig2.fa.half.head
```
- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:
```
formatdb -i contig2.fa.half -p F
```
- blast the start of the assembly against all of the assembly:
```
blastall -p blastn -i contig2.fa.half.head -d contig2.fa.half -e 1e-3 -F F -o contig2.bls
```
- look at contig2.bls to see hits:
```
less contig2.bls
```
- The first hit is against the start of the chromosome, as expected.
- The last hit starts at position 24886.
- We will trim the plasmid to position 24885
- first, index the contig2.fa.half file
```
samtools faidx contig2.fa.half
```
- trim
```
samtools faidx contig2.fa.half contig2:1-24885 > contig2.fa.half.trimmed
```
- open the contig2.fa.half.trimmed file
```
nano contig2.fa.half.trimmed
```
- delete the header info except contig name (e.g. contig2)
exit.
- we now have a trimmed contig2.
--->

## Find smaller plasmids
Pacbio reads are long, and may have been longer than small plasmids. We will look for any small plasmids using the Illumina reads.

This section involves several steps:

1. Use the multifasta canu-circlator output of trimmed assembly contigs.
2. Map all the Illumina reads against these pacbio assembled contigs.
3. Extract any reads that *didn't* map and assemble them together: this could be a plasmid, or part of a plasmid.
5. Look for overhang: if found, trim. If not, continue:
6. Search Genbank for any matching proteins: a replication protein found.  
7. Assemble all the Illumina reads and produce an assembly graph.
8. Search the graph for a match to the replication protein and its adjoining regions.
9. Extract this longer sequence from the Illumina assembly: this is the small plasmid.
10. Check for overhang in this plasmid and trim.

<!-- combine contigs 1 and 2
take chr and plasmid => one fasta file
```
cat contig1.fa.trimmed contig2.fa.half.trimmed > contig_1_2.fa
```
```
fa -f contig_1_2.fa
```
- Make directory and move in combined contigs file and the illumina reads
```
cp contig_1_2.fa ../
cd ..
mkdir find_contig_3
cp contig_1_2.fa R1.fastq.gz R2.fastq.gz find_contig_3/
cd find_contig_3
```
-->

###Align Illumina with BWA
- Align illumina reads to these contigs
- First, index the contigs file
```
bwa index contig_1_2.fa
```
- then, align using bwa mem
```
bwa mem -t 8 contig_1_2.fa R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
```
- the output alignment is <fn>aln.bam</fn>

<!--
###Or align with bowtie2
-->

###Extract unmapped illlumina reads
- Index the alignment file
```
samtools index aln.bam
```
- extract the fastq files from the bam alignment - those reads that were unmapped to the pacbio alignment - and save them in various "unmapped" files:
```
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam
```

<!--
- check R1 and R2 are approx same size:
```
ls -l
```
- check first reads have same name in R1 and R2:
```
head unmapped.R1.fastq unmapped.R2.fastq
```
-->

- we now have three files of the unampped reads:
- <fn> unmapped.R1.fastq</fn>
- <fn> unmapped.R2.fastq</fn>
- <fn> unmapped.RS.fastq</fn>
###Assemble the unmapped reads
- assemble with spades
```
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq -o spades_assembly
```

- **-1** is input file forward
- **-2** is input file reverse
- **-s** is unpaired
- **-o** is the output directory
```
cd spades_assembly
fa contigs.fasta
```
- shows how many assembled:
    - e.g. no=135
    - max = 2229
- sort fasta by size of seqs:
```
sizeseq
```
- Input sequence set: contigs.fasta
- Return longest sequence first [N]: Y
- output sequence(s) [contigs.fasta]: sorted_contigs.fasta

Print the first row of each seq to see coverage:
```
grep cov sorted_contigs.fasta  
```
- result: NODE_1_length_2229_cov_610.583
    - longest contig is 2229 and high coverage
- all the nodes are listed
- see if any other nodes have high coverage
    - e.g. NODE_135_length_78_cov_579
- look at the sequence of this contig:
```
tail sorted_contigs.fasta
```
- This is a homopolymer, so disregard.
- We will extract the first sequence (NODE_1):
```
samtools faidx sorted_contigs.fasta
samtools faidx sorted_contigs.fasta NODE_1_length_2229_cov_610.583 > contig3.fa
```
- this is now saved as <fn>contig3.fa</fn>
- open this file in nano, make the header ">contig3", save.
###Investigate the small plasmid (contig3)
- Blast the start of contig3 against itself
- Take the start of the contig:
```
head -n 10 contig3.fa > contig3.fa.head
```
- We want to see if it matches the end (overhang)
- Format the assembly file for blast:
```
makeblastdb -in contig3.fa -dbtype nucl
```
- blast the start of the assembly (.head file) against all of the assembly:
```
blastn -query contig3.fa.head -db contig3.fa -evalue 1e-3 -dust no -out contig3.bls
```
- look at <fn>contig3.bls</fn> to see hits:
```
less contig3.bls
```
- the first hit is against itself, as expected
- there are no few further hits, so we assume there is no overhang that needs trimming.
- however, the sequence is likely then to be longer than this.

```
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
- Assemble all the Illumina reads with spades (not just those reads that did not map to the pacbio assembly).

```
spades-fast --R1 R1.fastq.gz --R2 R2.fastq.gz --gsize 2.8M --outdir spades_fast --cpus 32
```
```
cd spades_fast
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

<!--- - there are two hits
this is probably due to overhang? eg this is where the plasmid overhang occurs and this is where it would later be trimmed?
--->

- Go to the main Bandage window
    - "find nodes" in right hand panel - 10
    - This node is slightly longer: 2373: this could be the plasmid
    - Extract this node in fasta format: In the top panel, go to Output: Save selected node sequences; save as <fn>contig3b.fa</fn>
    - Open this file in nano and change the header to ">contig3b", save.

### Trim small plasmid
- Take the start of the sequence and see if it matches the end:
```
head -n 10 contig3b.fa > contig3b.fa.head
makeblastdb -in contig3b.fa -dbtype nucl
blastn -query contig3b.fa.head -db contig3b.fa -evalue 1e-3 -dust no -out contig3b.bls
less contig3b.bls
```
- The first hit is against the start of the chromosome, as expected.
- The last hit starts at position 2253; we will trim the plasmid to position 2252
- Index the contig3b.fa file:
```
samtools faidx contig3b.fa
```
- Trim:
```
samtools faidx contig3b.fa contig3b:1-2252 > contig3b.fa.trimmed
```
- Open this file in nano and change the header to ">contig3b", save.
- We now have a trimmed contig3b.

### Collect all contigs in one file
```
cat contig_1_2.fa contig3b.fa.trimmed > all_contigs.fa
```
- See the three contigs and sizes:
```
fa -f all_contigs.fa
```

##Correct
<!-- check: Canu doesn't do any polishing, so there will be thousands of errors, almost all insertions. -->
We will correct the pacbio assembly, first with pacbio corrected reads (from Canu) and then with Illumina reads.

- inputs:
    - Draft pacbio assembly (overhang trimmed from each of the three replicons)
    - corrected pacbio reads from Canu.
    - illumina reads (aligned to pacbio assembly: in bam format)
- output: corrected assembly

### 1. Correct with pacbio corrected reads

Align reads to assembly:

```bash
bwa index all_contigs.fa
bwa mem -t 32 all_contigs.fa canu.correctedReads.fasta.gz | samtools sort > aln.bam
samtools index aln.bam
samtools faidx all_contigs.fa
```

- **-t** is the number of cores (e.g. 8); to find out how many you have, grep -c processor /proc/cpuinfo
- now we have an alignment file to use in pilon: <fn>aln.bam</fn>

Run pilon:

```bash
pilon --genome all_contigs.fa --frags aln.bam --output corrected --fix all --mindepth 0.5 --changes --threads 32 --verbose
```
- **--genome** is the name of the input assembly to be corrected
- **--frags** is the alignment of the reads against the assembly
- **--output** is the name of the output directory
- **--fix** is an option for types of corrections
- **--mindepth** gives a minimum read depth to use
- **--changes** produces an output file of the changes made
- **--threads** is the number of cores
- **--verbose** prints information to the screen during the run


Look at the .changes file:
```bash
less corrected.changes | column -t
```
- This shows the corrections made by Pilon:

![Pilon](images/pilon.png)

- look at the .fasta file.


###Align Illumina reads => bam

used the corrected assembly from the previous step:


```
bwa index all_contigs.fa
bwa mem -t 32 all_contigs.fa R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
samtools index aln.bam
samtools faidx all_contigs.fa
```


- look at how the illumina reads are aligned:
```
samtools tview -p contig1 aln.bam all_contigs.fa
```


<!--
contigs & bam => pilon => corrected contigs
don't run on things with lots (10+) of contigs, but 1 or 2 is ok
-->
- run pilon:

```bash
pilon --genome all_contigs.fa --frags aln.bam --output corrected --fix bases --mindepth 0.5 --changes --threads 32 --verbose
```
- look at the changes file
```bash
less corrected.changes
```
- look at the .fasta file.

- option: re run pilon to correct again.
  

**Final output:**

- the corrected genome assembly of *Staphylococcus aureus* in .fasta format, containing three contigs: chromosome, large plasmid and small plasmid.  

<!--
-but there are some bigger ones.
    - e.g. at position 2,463,699 there is a 23-bp seq that is deleted.

###View the deletion in a pileup

tig00000000:2463600-2463622 tig00000000_pilon:2463699 GTTAAAGGTTATTTGAATGATCA .

- view the illumina reads mapped against the original assembly:
```
samtools tview -p tig00000000:2463600 aln.bam all_contigs.fa
```

- -p gives the chromosome position that we want to view
- then input file; then reference file

view:
- reference sequence along the top
- then consensus from the aligned illumina reads
- a dot is a match on F
- a comma is a match on R
- capital letter is a correction on F
- small letter is a correction on R
- asterisk
- underlining

view in IGV:
- transfer aln.bam, aln.bam.bai, contigs.fa, contigs.fa.fai to local computer
- open IGV
- File - genomes - load genome from file: contigs
- File - load from file - aln.bam file
- select the main chromosome eg tig00000000
- go to coordinate 2463600
- there is a dip in coverage here.
- pilon has identified this as a local misassembly and has deleted this section in the corrected assembly

![IGV screenshot](images/IGV_first_pilon_correction.png)

- output:corrected assembly

- re run pilon to correct again.
    - but this time, the reference is the first pilon correction

```
bwa index corrected.fasta
bwa mem -t 8 corrected.fasta R1.fastq.gz R2.fastq.gz | samtools sort > aln_to_corrected.bam
samtools index all_to_corrected.bam
samtools faidx corrected.fasta
```

- separate into three fasta files

-->
<!--
###or, use Bowtie for alignment

index the fasta file

```
bowtie2-build all_contigs.fa bowtie
```

input file
output ref name

align:

```
bowtie2 --end-to-end --threads 72 -x bowtie -1 25745_1_PE_700bp_SEP_UNSW_ARE4E_TCGACGTC_AAGGAGTA_S5_L001_R1.fastq.gz -2 25745_1_PE_700bp_SEP_UNSW_ARE4E_TCGACGTC_AAGGAGTA_S5_L001_R2.fastq.gz | samtools sort > bowtie2.bam

```

- "bowtie" means use the index files that were generated above that we called bowtie
- use 72 threads on mdu but <8 on sepsis
- then index the bam file
```
samtools index bowtie2.bam
samtools faidx contigs.fasta
```

- then view in tview
```
samtools tview -p tig00000000:2463600 bowtie2.bam all_contigs.fa
```
<!--
###or use BLASR for alignment: pacbio reads to pacbio assembly

(in MDU marvin)
```
blasr subreads.fastq all_contigs.fa -nproc 72 -sam | samtools sort > blasr_aln.bam
```

blasr subreads.fastq all_contigs.fa -nproc 72 -sam | samtools sort > blasr_aln.bam
samtools index blasr_aln.bam
samtools faidx all_contigs.fa
samtools tview -p tig00000000:2463600 blasr_aln.bam all_contigs.fa

###or align pacbio reads to pacbio assembly with bwa mem

canu: pacbio corrected reads
bwa mem
align to pacbio assembly


bam index all_contigs.fasta

bwa mem -t 72 all_contigs.fa staph.correctedReads.fasta | samtools sort > pacbio_aln.bam


(leave out -x pacbio option as these are corrected reads)
--->


<!--

##Examine 23bp

align corrected pacbio reads to pacbio assembly

make a 23bp seq

```
nano 23bp.fa
```

paste in the seq

-->

<!--

### Split into three final contigs

- Examine the all_contigs.fa file
```
fa -f
```
- contig1	dna	2725158
- contig2	dna	24884
- contig3b	dna	2252
- split these
```
samtools faidx all_contigs.fa contig1: > chromosome.fa
samtools faidx all_contigs.fa contig2: > plasmid.fa
samtools faidx all_contigs.fa contig3b: > small_plasmid.fa
```

##Circularise
corrected contigs => circlator => circular genome (e.g. starting at dnaA)

###Chromosome

```
circlator fixstart contig.fa outprefix
```

circlator fixstart chromosome.fa chrm



default: orients at DNAa

=> circular, corrected assembly

###Plasmid 1

blast against ncbi
find hit
what have they done
download the start of their assembly (may be a gene or may be tandem repeats etc. )

```
circlator fixstart --genes_fa filename contig outprefix
```


 -or, whatever has been done with the plasmid it most closely blast matches to
eg
  https://www.ncbi.nlm.nih.gov/nucleotide/260066114?report=genbank&log$=nuclalign&blast_rank=1&RID=WTK2RJSZ014


need to download

###Plasmid 2

-->
<!--
##Quickstart overview
A summary of the commands used in this example (not all directory changes and minor steps listed):

(Note, filenames in this section and the next section are not yet finalised).

example with Staph aureus. Sample 25745
```bash
# join files and analyse with canu
cat pacbio* > subreads.fastq.gz
# gunzip
# find info about read lengths
fq subreads.fastq
# average length is 10k; mode is 20k
mkdir 01.canu
cd 01.canu
canu -p canu -d canu.default genomeSize=2.8m -pacbio-raw ../subreads.fastq.gz
# examine how many contigs found and their lengths
cd canu.default/
fa -f canu.contigs.fasta
# 1 chrm; 1 plasmid
# if too many contigs found, repeat canu with sensitive options
# use circlator to trim overhangs and orient start of assemblies
mkdir 02.circlator
cd circlator
# using a length cutoff of approx 2X mode read length
circlator all --threads 72 --verbose --b2r_length_cutoff 40000 ../01.canu/canu.default/canu.contigs.fasta ../01.canu/canu.default/canu.correctedReads.fasta.gz circlator_all_output
# check the output to see if contigs were circularised
cd circlator_all_output
less 04.merge.circularise.log
# check where the contigs were oriented
less 06.fixstart.log
# look at new contig trimmed lengths
fa -f 06.fixstart.fasta
# rename file and copy up
mv 06.fixstart.fasta canu.circlator.contigs.fa
# use illumina reads to find smaller plasmids
mkdir 03.find_smaller_plasmids
bwa index ../canu.circlator.contig.fa
bwa mem -t 32 ../canu.circlator.contig.fa ../R1 ../R2 | samtools sort > aln.bam
samtools index aln.bam
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq -o spades_assembly
# use sizeseq to find longest contig
# blast start agains all; if overhang found; trim

    # if no overhang found: ncbi blastx hits a rep protein suggesting this is a small plasmid
    # assemble all illumina reads to find flanking regions
    # mkdir 04.all_illumina_spades
    spades-fast --R1 ../R1 --R2 ../R2 --gsize 2.8M --outdir spades_fast --cpus 32
    # view assembly_graph.fastg in bandagae and blast using the contig
    # it matches a node of length 2373: extract this node from the all-illumina contigs
    samtools faidx contigs.fasta
    samtools faidx contigs.fasta NODE_11_length_2373_cov_417.492 > contig3b.fa
    # blast start against all to find overhang
    head -n 10 contig3b.fa > contig3b.fa.head
    makeblastdb -in contig3b.fa -dbtype nucl
    blastn -query contig3b.fa.head -db contig3b.fa -evalue 1e-3 -dust no -out contig3b.bls
    samtools faidx contig3b.fa
    samtools faidx contig3b.fa contig3b:1-2252 > contig3b_trimmed.fa
    # ncbi blast x: find this matches rep protein; save as nucleotides as rep_protein.fa

# fixstart
circlator fixstart --genes_fa rep_protein.fa contig3b_trimmed.fa fixstart
# join all contigs and plasmids
cat canu.circlator.contigs.fa small_plasmid.fa > pre_pilon_staph.fa
#run pilon to correct using illumina reads
mkdir 05.pilon/
cd 05.pilon/
bwa index ../pre_pilon_staph.fa
bwa mem -t 32 ../pre_pilon_staph.fa ../R1 ../R2 | samtools sort > aln.bam
samtools index aln.bam
samtools faidx ../pre_pilon_staph.fa
pilon --genome ../pre_pilon_staph.fa --frags aln.bam --output pilon --fix all --mindepth 0.5 --changes --threads 32 --verbose
#check the output: one large deletion - this is a tandem repeat in the pacbio assembly; unsupported by the illumina reads.
# final assembly pilon.fasta changed name to staph_25745.fasta

options:

#re-orient larger plasmid, rather that use circlator's choice
samtools faidx staph_25745.fasta
samtools faidx staph_25745.fasta tig00000001_pilon > large.plasmid.fa
samtools faidx staph_25745.fasta tig00000000_pilon > chrm.fa
samtools faidx staph_25745.fasta contig3b:1-2252_pilon > small.plasmid.fa
# blastn on ncbi : matches staph aureus plasmid CP002115.1 by Tim Stinear's group
# this is oriented at start of repA; saved this as a fasta file
circlator fixstart --genes_fa repA_CP002115.fa ../large.plasmid.fa fixstart
# rejoin all three contigs into a multifasta file
```

-->
##Next

**Further analyses:**

- Annotate with Prokka.
- Comparative genomics, e.g. with Roary.

**Links:**

Canu

- [Canu manual](http://canu.readthedocs.io/en/stable/quick-start.html)
- [Canu code](https://github.com/marbl/canu)

Circlator

- http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0
- http://sanger-pathogens.github.io/circlator/

Pilon

- [Pilon article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963)
- [Pilon on github](https://github.com/broadinstitute/pilon/wiki)

Finishing

- https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes
- https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies

Details of bas.h5 files

- https://s3.amazonaws.com/files.pacb.com/software/instrument/2.0.0/bas.h5+Reference+Guide.pdf


## Commands summary <a name="summary"></a>

```bash
//get the sequencing reads
//joined pacbio reads: subreads.fastq
//illumina reads: R1.fastq.gz, R2.fastq.gz

//run canu to assemble
canu -p staph -d output genomeSize=2.8m -pacbio-raw subreads.fastq
//output: staph.contigs.fasta

//run circlator to trim and orient
circlator all --verbose ../staph.contigs.fasta ../staph.corrected.reads.fastq.gz circlator_all_output
//output: 06.fixstart.fasta; rename as contig_1_2.fa

//find smaller plasmids

//align illumina to pacbio contigs
bwa index contig_1_2.fa
bwa mem -t 8 contig_1_2.fa R1.fastq.gz R2.fastq.gz | samtools sort > aln.bam
//output is aln.bam

//extract unmapped reads
samtools index aln.bam
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam

//assemble these reads with spades
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq -o spades_assembly
//save any long contigs; e.g. as contig3.fa

//blast this contig to find overhang
head -n 10 contig3.fa > contig3.fa.head
makeblastdb -in contig3.fa -dbtype nucl
blastn -query contig3.fa.head -db contig3.fa -evalue 1e-3 -dust no -out contig3.bls
less contig3.bls
//no overhang found

//assemble all illumina reads
spades-fast --R1 R1.fastq.gz --R2 R2.fastq.gz --gsize 2.8M --outdir spades_fast --cpus 32
//blast the assembly graph with contig3.fa; extract out node and save as contig3b.fa

//trim overhang
head -n 10 contig3b.fa > contig3b.fa.head
makeblastdb -in contig3b.fa -dbtype nucl
blastn -query contig3b.fa.head -db contig3b.fa -evalue 1e-3 -dust no -out contig3b.bls
less contig3b.bls
samtools faidx contig3b.fa
samtools faidx contig3b.fa contig3b:1-2252 > contig3b.fa.trimmed
//open contig3b.fa in nano and shorten header name

//join all contigs: all_contigs.fa

//correct with pilon

//1. correct using corrected pacbio reads
bwa index all_contigs.fa
bwa mem -t 32 all_contigs.fa canu.correctedReads.fasta.gz | samtools sort > aln.bam
samtools index aln.bam
samtools faidx all_contigs.fa
pilon --genome all_contigs.fa --frags aln.bam --output corrected --fix bases --mindepth 0.5 --changes --threads 32 --verbose

//2. correct using illumina reads - use output from 1 as the contigs file

//3. repeat correction using Illumina reads - use output from 2 as the contigs file

//output: corrected genome assembly of Staphylococcus aureus in .fasta format, containing three contigs: chromosome, large plasmid and small plasmid.
```
