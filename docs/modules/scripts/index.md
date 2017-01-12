<br>
#Bash scripts

It is possible to put some commands into a script, to run automatically.

Here, we will put commands for several tools into a bash script and use it to assemble a bacterial genome.

Bash is a programming language, and we can run a bash script in the bash shell (e.g. the terminal).

##Prepare the bash script

- Open the terminal.
- Open a text editor, e.g., nano.

```bash
nano myscript.sh
```

- We now have a script called myscript.sh
- In the script, write

```bash
#!/bin/bash  
echo This is my assembly script
```
- Ctrl-X to exit nano.

- Make the script executable:

```bash
chmod 755 myscript.sh
```

- Run the script to test.

```bash
./myscript.sh
```

- The text should print to the screen:
``` bash
This is my assembly script
```



##Edit the bash script

- Now we will add to the bash script.
- This will be a simple script that is designed to assemble, trim, and correct.
- It will likely only work in simple cases.
- Open nano

```bash
nano myscript.sh
```

- and type in the following
(#comments and echo statements are optional but probably useful):


```bash
#!/bin/bash
#a script for assembling a bacterial genome

#set a variable for the number of CPUS
#this can then be easily changed if needed
CPUS=16

#this script requires 5 arguments
#we will test whether 5 were entered
if [ ! $# -eq 5 ] #if there are not 5 arguments
then
    echo "ERROR: needs 5 arguments"
    echo "usage: $0 outdir genomeSize (e.g. 2.1m) pacbio.fq R1.fq R2.fq"
    #the usage shows what arguments are required
    exit 1 #exits
fi

#report back the arguments entered:
#each argument is assigned a variable name from 1 to 5
#these variables are referred to with the $ sign
echo "your directory with output is called $1"
echo "your genomeSize is $2"
echo "your pacbio reads are in file $3"
echo "your R1 reads are in file $4"
echo "your R2 reads are in file $5"
echo #echo on its own puts in a blank line on the screen

#run canu to assemble the reads
mkdir -p "$1"  #the p creates intermediate directories where required
echo "now running canu"
echo "your canu output is in directory "$1/canu""
nice canu -p canu -d "$1/canu" genomeSize="$2" -pacbio-raw "$3"
#commands are preceeded by 'nice' in this script for particular servers

#test canu output: number of contigs
#make a variable called num_contigs
num_contigs=$(grep -c '>' "$1/canu/canu.contigs.fasta")
echo "number of contigs found by canu is: $num_contigs"
echo
if [ $num_contigs -gt 10 ]
then
    echo "there are more than 10 contigs"
    echo "the analysis is stopping"
    exit 1
#in more complex scripts, you can add the option to repeat the analysis with different parameters, etc.
fi

#run circlator to trim the contigs
echo "now running circlator to trim and orient contigs"
echo "your circlator output is in directory "$1/circlator" "
nice circlator all --threads $CPUS --verbose "$1/canu/canu.contigs.fasta" "$1/canu/canu.correctedReads.fasta.gz" "$1/circlator"
#in more complex scripts, you can add information about which contigs were merged; where contigs were oriented; how contigs were trimmed, etc.
#rename the contigs file
cp "$1/circlator/06.fixstart.fasta" "$1/contigs.fa"
echo "your contigs from circlator are now in the file called contigs.fa"
#note: in this case, these contigs may not have been trimmed.

#find smaller plasmids
echo "now looking in illumina reads for small plasmids"
echo "indexing the contigs.fa file"
nice bwa index "$1/contigs.fa"
echo "aligning illumina reads to the contigs.fa file"
nice bwa mem -t $CPUS "$1/contigs.fa" "$4" "$5" | samtools sort > "$1/aln.bam"
echo "indexing the aln.bam file"
nice samtools index "$1/aln.bam"
echo "extracting the unmapped illumina reads"
samtools fastq -f 4 -1 "$1/unmapped.R1.fastq" -2 "$1/unmapped.R2.fastq" -s "$1/unmapped.RS.fastq" "$1/aln.bam"
echo "assembling the unmapped reads with spades"
echo
nice spades.py -t $CPUS -1 "$1/unmapped.R1.fastq" -2 "$1/unmapped.R2.fastq" -s "$1/unmapped.RS.fastq" -o "$1/spades_assembly"
echo
echo "your output is in the spades_assembly directory"
echo

#check if any contigs are > 1500bp
seqtk seq -L 1500 "$1/spades_assembly/contigs.fasta" > "$1/illumina_contigs.fasta"
echo "extracting contigs > 1500bp into illumina_contigs.fasta"
echo
#in more complex scripts, you can trim and orient these contigs, or do further investigation into whether they are true plasmids

#join the canu pacbio contigs and the small plasmids
cat "$1/contigs.fa" "$1/illumina_contigs.fasta" > "$1/allcontigs.fasta"

echo
echo "joining small contigs to the pacbio contigs as allcontigs.fasta"
echo

#correct this assembly with Pilon
echo "preparing to correct contigs with illumina reads, using Pilon"
echo "indexing allcontigs.fasta"
nice bwa index "$1/allcontigs.fasta"
echo "aligning illumina reads to allcontigs.fasta"
nice bwa mem -t $CPUS "$1/allcontigs.fasta" "$4" "$5" | samtools sort > "$1/aln.bam"
echo "indexing aln.bam"
nice samtools index "$1/aln.bam"
echo "indexing allcontigs.fasta"
nice samtools faidx "$1/allcontigs.fasta"
echo "correcting with Pilon"
nice pilon --genome "$1/allcontigs.fasta" --frags "$1/aln.bam" --output "$1/pilon1" --fix all --mindepth 0.5 --changes --verbose --threads $CPUS
echo
echo "Pilon finished"
echo

#count how many changes were made
pilon_changes=$(wc -l < "$1/pilon1.changes")
echo "there were $pilon_changes corrections made"
echo
echo "Assembly complete"
echo "The sequence is saved as "$1/pilon1.fasta""

exit 1
```

##Run

To run, type in

```bash

./myscript.sh outdir 1.9m subreads.fastq R1.fastq R2.fastq

```
- **outdir** is the name of the output directory
- **1.9m** is the approximate genome size for the species
- **subreads.fastq** is the name of the Pacbio reads file
- **R1.fastq** is the name of the R1 Illumina reads file
- **R2.fastq** is the name of the R2 Illumina reads file

- If the Pacbio and Illumina reads files are not in the current directory, put in the whole path to those files

- The script should run to the end or may exit early if necessary (e.g. the canu assembly found too many contigs).

- There are many refinements and additions that could be made to this script, and some of these are noted as comments above.
