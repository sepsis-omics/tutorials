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

###Find the PacBio files for this sample

- Obtain the input files. e.g. from the BPA portal.
- PacBio files are often stored in the format:
Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq.gz
- We will use the <fn>fastq.gz</fn> files.
- The reads are usually split into three separate files because they are so large.
- Right click on the first <fn>subread.fastq.gz</fn> file and "copy link address".

- In the command line, type:

```
wget --user username --password password [paste link URL for file]
```
- repeat for the other two <fn>subread.fastq.gz</fn> files.

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
- This will take ~ 30 minutes.

###Check the output

```
cd staph_outdir
```

- The <fn>contigs.fasta</fn> are the assembled sequences.
- The <fn>unassembled.fasta</fn> are the reads that could not be assembled.
- The <fn>file.gfa</fn> is the graph of the assembly.
- Display summary information about the contigs:

```
fa -f staph.contigs.fasta
```

![contigs_info_screenshot](images/contigs_info.png)

- There are two contigs. Note down their lengths: 2,725,231 and 43,991
- We will split these into two separate fasta files.

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
- change header to contig1 or contig2
- save
- We now have two files:
    - <fn>contig1.fa</fn>
    - <fn>contig2.fa</fn>

<!-- ###What is in the unassembled reads?

- Are they contaminants?
- Blast against NCBI
    - Use Cyberduck to transfer the file to your local computer
    - Go to NCBI and upload the file
    - (doesn't work with this file - too big)
    - (or: explain how to do so on command line?)
-->

## Alternatives for assembly
- canu without sensitivity settings
- HGAP2

##Use Quiver/Arrow here?

to correct assembly with the raw reads.

## Trim overhang in chromosome 1

Chromosome 1 is in the <fn>contig1.fa</fn> file.

The canu assembly graph will be split at an ambiguous/unstable node. However, this area of the graph likely overlaps in the bacterial chromosome, but has not aligned with itself completely. This is called an overhang. We need to identify these overhangs and trim them, for the chromosome and any plamsids.

- Take the first 30,000 bases of contig1.

```
head -n 501 contig1.fa > contig1.fa.head
```

- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:

```
formatdb -i contig1.fa -p F
```

- blast the start of the assembly (.head file) against all of the assembly:

```
blastall -p blastn -i contig1.fa.head -d contig1.fa -e 1e-10 -F F -o contig1.bls
```

- look at <fn>contig1.bls</fn> to see hits:

```
less contig1.bls
```

- The first hit is against the start of the chromosome, as expected.
- This position, near the end of the contig, is where the overhang starts: position 2725159

![screeshot of blast](images/contig1_bls.png)


- check that this match goes to the end of the contig: 2725231 - yes.
- we will now trim the contig
- first, index the <fn>contig1.fa</fn> file

```
samtools faidx contig1.fa
```

- this makes an index file called <fn>contig1.fa.fai</fn>
- next, extract all the sequence except for the overhang. (We don't have to specify the name of the *index* file; it will be found automatically):
```
samtools faidx contig1.fa tig00000000:1-2725158 > contig1.fa.trimmed
```

- open the <fn>contig1.fa.trimmed</fn> file

```
nano contig1.fa.trimmed
```

- delete the header info except contig name (e.g. tig00000000)
exit.
- we now have a trimmed contig1.

##Investigate the plasmid

- Contig 2 is 43,991 bases.
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

![Gepard dotplot screenshot](images/gepard_plasmid.png)

- the sequence starts to repeat at position ~ 24850
- there are probably two plasmids combined into one contig in tandem.
- the whole length is 43,991
- we will cut 9000 from each end
- use samtools to extract the region of the contig except for 9k on each end

```
samtools faidx contig2.fa
samtools faidx contig2.fa tig00000001:9000-35000 > contig2.fa.half
```

- we now have a ~ 26000 bases plasmid
```
nano contig2.fa.half
```
- delete the rest of the name after tig00000001

## Trim overhang in the plasmid

contig 2.

Take the first x bases:

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
- The last hit starts at position 24885.
- We will trim the plasmid to position 24884

- first, index the contig2.fa.half file

```
samtools faidx contig2.fa.half
```

- trim

```
samtools faidx contig2.fa.half tig00000001:1-24884 > contig2.fa.half.trimmed
```

- open the contig2.fa.half.trimmed file

```
nano contig2.fa.half.trimmed
```

- delete the header info except contig name (e.g. tig00000000)
exit.

- we now have a trimmed contig2.

## Combine contigs 1 and 2

take chr and plasmid => one fasta file
```
cat contig1.fa.trimmed contig2.fa.half.trimmed > contig_1_2.fa
```

```
fa -f contig_1_2.fa
```
## Look for smaller plasmids

- Make directory and move in combined contigs file and the illumina reads

```
cp contig_1_2.fa ../

cd ..

mkdir find_contig_3

cp contig_1_2.fa R1.fastq.gz R2.fastq.gz find_contig_3/

 cd find_contig_3
```

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

###Or align with bowtie2

Add info here


###Extract unmapped illlumina reads

- Index the alignment file


```
samtools index aln.bam

```

- extract the fastq files from the bam alignment - those reads that were unmapped to the pacbio alignment.

```
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam

ls -l # check R1 and R2 are approx same size

head unmapped.R* #check first reads have same name in R1 and R2

head unmapped.R1.fastq unmapped.R2.fastq

```
- we now have three files of the unampped reads
###Assemble the unmapped reads

- assemble with spades

```
spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq -o spades_assembly
```
-1 is input file forward
-2 is input file reverse
-s is unpaired
-o is output_dir

```
cd spades_assembly

fa scaffolds.fasta
```

- shows how many assembled:
    - e.g. no=136
    - max = 2229

- sort fasta by size of seqs:

```
sizeseq
Input sequence set: scaffolds.fasta
Return longest sequence first [N]: Y
output sequence(s) [scaffolds.fasta]: sorted_scaffolds.fasta
grep cov sorted_scaffolds.fasta  #just print the first row of each seq to see coverage
```

- result: NODE_1_length_2229_cov_610.033
    - longest contig is 2229 and high coverage

- all the nodes are listed
- see if any other ones have high coverage
    - e.g. >NODE_136_length_78_cov_579
- look at the sequence of this contig:

```
tail sorted_scaffolds.fasta
```

- seq looks bad for node136 so disregard (seq is only Cs)

- we will extract the first seq

```
samtools faidx sorted_scaffolds.fasta
samtools faidx sorted_scaffolds.fasta NODE_1_length_2229_cov_610.033 > contig3.fa
```
- this is now saved as <fn>contig3.fa</fn>

###Investigate the small plasmid (contig3)

- blast start of contig3 against itself
- Take the first x bases:

```
head -n 10 contig3.fa > contig3.fa.head
```

- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:

```
formatdb -i contig3.fa -p F
```

- blast the start of the assembly against all of the assembly:

```
blastall -p blastn -i contig3.fa.head -d contig3.fa -e 1e-3 -F F -o contig3.bls
```

- look at contig3.bls to see hits:

```
less contig3.bls
```
- the first hit is against itself, as expected
- there are a few further hits to small sections, but none to the end, so we assume there is no overhang that needs trimming.
- however, the sequence is likely then to be longer than this.

- copy the sequence:
```
less contig3.fa
```

- copy the sequence
- go to ncbi blast
- nuc to protein
- paste seq
- choose genetic code = 11
- blast
- hits: replication protein
- so this is a protein that has not been found in the pacbio assembly.
- hypothesise that	this might be a true small plasmid but the rest of its seq is in common with other parts of the staph genome, so they haven't been assembled with the rep protein

##Assemble *all* the illumina reads

- assemble all the illumina reads with spades (not just those reads unmapped to pacbio assembly)
- resulting assembly: search for small plasmid sequence
- use spades-fast on MDU server
- mkdir
- wget illumina reads
- run spades-fast:

```bash
spades-fast --R1 file --R2 file --gsize 2.8M --outdir spades_fast --cpus 32
```
```
cd spades_fast
```

- in here is the <fn>assembly_graph.fastg</fn>
- use cyberduck to transfer this file to local computer
- Examine the assembly in the program Bandage.
    - Install Bandage.
    - File: Load graph: <fn>file.gfa</fn>
    - In the left hand panel, click "Draw graph"

![bandage pic](images/illumina_assembly_bandage.png)

- see: main chromosome, plasmid, small plasmid (annotate pic)

- blast the small plasmid sequence in this assembly
- left hand panel: Blast: create/view BLAST search
- build blast database
- paste in nuc seq of contig3
- blast

- the main hit is around node 10
<!--- - there are two hits
this is probably due to overhang? eg this is where the plasmid overhang occurs and this is where it would later be trimmed?
--->

- go to bandage window
- "find nodes" in right hand panel - 10
- this node is slightly longer: 2373
- this could be the plasmid
- find the sequence: NODE_10
- extract this seq
- this is called NODE_10 in bandage but in the assembly_graph.fastg file it is called EDGE_10
- go to the folder for spades_fast

```
less contigs.fasta
/2373
```
- search for length of this node
- see that this is actually called NODE_11_length_2373_cov_417.492
- extract out node 11

```
samtools faidx contigs.fasta
samtools faidx contigs.fasta NODE_11_length_2373_cov_417.492 > contig3b.fa
```

- we will call it contig 3 "b" because it is larger than our original contig3.
- open it - copy the seq
- move this file back to sepsis abrpi machine:
    - open a nano file: nano contig3a.fa
    - paste the seq  
    - save

## Trim small plasmid
- next: check for overhang
- Take the first x bases:

```
head -n 10 contig3b.fa > contig3b.fa.head
```

- this is the start of the assembly
- we want to see if it matches the end (overhang)
- format the assembly file for blast:

```
formatdb -i contig3b.fa -p F
```

- blast the start of the assembly against all of the assembly:

```
blastall -p blastn -i contig3b.fa.head -d contig3b.fa -e 1e-3 -F F -o contig3b.bls
```
- look at contig3b.bls to see hits:
```
less contig3b.bls
```
- The first hit is against the start of the chromosome, as expected.
- The last hit starts at position 2253
- We will trim the plasmid to position 2252
- first, index the contig3b.fa file
```
samtools faidx contig3b.fa
```
- trim
```
samtools faidx contig3b.fa contig3b:1-2252 > contig3b.fa.trimmed
```
- open the contig3b.fa.trimmed
```
nano contig3b.fa.trimmed
```
- delete the header info except contig name (e.g. contig3b)
exit.
- we now have a trimmed contig3b.

## Collect all contigs in one file

- move up to main analysis folder
- mkdir pilon
- copy the trimmed contigs 1, 2, 3b into this folder
- cat them
```
cat contig_1_2.fa contig3b.fa.trimmed > all_contigs.fa

fa -f all_contigs.fa
```

- see the three contigs and sizes
- rename all the contigs
- contig1
- contig2
- contig3b

##Correct with Pilon

Canu doesn't do any polishing, so there will be thousands of errors, almost all insertions.

We will use the illumina reads to correct the pacbio assembly.

- inputs:
    - draft pacbio assembly (overhang trimmed from each of the three replicons)
    - illumina reads (aligned to pacbio assembly: in bam format)

- output: corrected assembly

###Align Illumina reads => bam
- might have already been done above
```Unix
bwa index contigs.fasta
bwa mem -t 8 contigs.fasta r1s r2s | samtools sort > aln.bam
samtools index aln.bam
samtools faidx contigs.fasta
```

- -t is the number of cores (e.g. 8)
- to find out how many you have, grep -c processor /proc/cpuinfo

###Correct

<!--
contigs & bam => pilon => corrected contigs
don't run on things with lots (10+) of contigs, but 1 or 2 is ok
-->

- run pilon:
```bash
pilon --genome all_contigs.fa --frags all_aln.bam --output corrected --fix bases --changes --threads 8 --verbose
```

- look at the output file called .pilon.fasta
    - compared to input file, this one could be longer
- look at the changes file
```
less corrected.changes
```

- a lot are single, but there are some bigger ones.
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
###or, use Bowtie for alignment

index the fasta file

```
bowtie2-build all_contigs.fa bowtie
```

input file
output ref name

align:

```bash
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
###or use BLASR for alignment: pacbio reads to pacbio assembly

(in MDU marvin)
```
blasr subreads.fastq all_contigs.fa -nproc 72 -sam | samtools sort > blasr_aln.bam
```

##Circularise
corrected contigs => circlator => circular genome (e.g. starting at dnaA)

###Chromosome

```
circlator fixstart contig.fa outprefix
```

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




## Next

Annotate with prokka


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
