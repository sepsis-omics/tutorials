# Galaxy - Spades

This tutorial will demonstrate how to assemble a bacterial genome de novo using the program Spades in the mGVL.

The Spades webpage is [here](http://bioinf.spbau.ru/spades).

##Pre-requisites

- Galaxy
- de novo assembly
- QC
- Trimming

##Start

- open your mGVL galaxy

##Input data

FIXME: get different data that shows how contigs are joined into larger scaffolds.

Bacterial DNA sequences:

- Illumina paired-end reads of approximately 150 base pairs, from *Staphylococcus aureus*. Subset of a real dataset.

Get data files into Galaxy:

1. In the history pane, click on the **History Options** button at the top right corner.

2. A list appears: at the bottom of the list, click on **Import from File**.

3. In the centre pane, a text box appears called "Archived History URL:". Into this box, paste the web address of the stored data that we will use. This is: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Galaxy-History-adelaideworkshopstart.tar.gz  (FIXME: rename something simpler and only keep the fastq files) and then click **Submit**. Galaxy will download these files.

4. In the history pane, click on the **View all Histories** button at the top right corner. The imported files will be a new history called "imported from archive: adelaide_workshop_start" (rename). Above this, click on the **Switch to** button to make this the current history. Then click **Done** in the top left corner.

5. In the history pane, we now have two files: aaaR1.fastq and aaaR2.fastq. (FIXME: to rename)

6. Click on the eye icon in the history panel to view the contents of each file. This sequencing produced paired-end reads, which means each input DNA fragment was sequenced at both the start and the end (and not sequenced in the middle). All the start sequences are grouped into one file (aaaR1.fastq), as are all the end sequences (aaaR2.fastq). Underneath each sequence is a string of letters/symbols, which correspond to estimates of the quality of each nucleotide call. For more information on fastq file format, click [here](https://en.wikipedia.org/wiki/FASTQ_format).

##How does Spades work?

1. As with several other de novo assembly programs (e.g.Velvet) Spades uses an algorithm based on [de Bruijn graphs.](http://www.nature.com/nbt/journal/v29/n11/full/nbt.2023.html) Such graphs use sub-lengths of sequence reads to build an overall genome assembly. The span of the sub-length is called a k-mer, where "k" is the number of nucleotides (e.g. k=21). The user chooses three values of k and Spades makes three assemblies based on these.

2. For the first value of k, each read is broken into as many fragments as possible. For example, if the input read is 22 nucleotides long, and the chosen value of k is 21, then there are two possible fragments (positions 1-21 and 2-22).

3. One randomly-chosen fragment becomes the first node on the de Bruijn graph.

4. A second fragment is connected to this node if it overlaps.

5. Repeat until all fragments are connected. Output: de Bruijn graph.

6. Find a connected pathway through this graph. Output: a pathway (sequence) known as a contig. Because of poor or incorrect sequencing, not all the fragments can be joined together. There will be several de Bruijn graphs and so several contigs, usually of different sizes.

7. Repeat these steps for a further two values of k (e.g. k = 33, k = 55). Output: Three (sets of) contigs.

8. Merge the three (sets of) contigs to get one. Output: one set of contigs.

9. For paired-end reads (as in this tutorial), the two reads are sequenced from each end of a longer DNA fragment. The middle part of the fragment is not sequenced, but information about the distance between the reads can be used by Spades to join contigs into larger sequences, called scaffolds. Output: one set of scaffolds.

10. To fix any errors map the original sequence reads onto the scaffolds with the program BWA. Output: assembled genome.

## Run Spades

1. In the tools pane, click on **NGS Analysis -> NGS: Assembly -> spades**. Leave all settings as they are except these:

2. Under **Run only assembly?** click **Yes**.

3. Under **K-mers to use, separated by commas**: enter in 33,55,91. Don't put any spaces after the commas.

4. For **Coverage cutoff**, choose **Auto**.

5. For **Libraries -> Files -> 1:Files -> Select file format -> Forward reads** : R1.fastq.

6. For **Libraries -> Files -> 1:Files -> Select file format -> Reverse reads** : R2.fastq.

7. Click **Execute**.

FIXME: include screenshot

## Output

There are five output files that will appear in the history pane. Click on the eye icon to view each file. These are the assembled genome sequence, and various statistics/information on how the program ran.

1. **contig stats**: There are x contigs. Look at the variation in length and coverage. A short contig with high coverage could be a result of contamination, a collapsed repeat, or a plasmid.

2. **contigs**: Each contig is listed, followed by its sequence in fasta format.

3. **scaffold stats**: There are x scaffolds.

4. **scaffolds**: Each scaffold is listed, followed by its sequence in fasta format.

5. **log**: The specific actions performed in the analysis.

## Filter output

1. In the tools pane, click on **NGS Analysis -> NGS: Asssembly -> Filter SPAdes output**. This is a quick way to discard contigs that are too short (e.g., they might be contamination) or contigs that do not have enough coverage (e.g., they might be too unreliable).

2. Under **Sequences**, choose the contigs fasta file.

3. Under **Contig stats** choose the contigs stats file. Change the cut-off values for length and coverage or leave them as they are.

4. For **Save filtered-out sequences?** click **Yes**.

5. Click **Execute**. A new fasta file with only the filtered sequences will be saved in the right-side history pane.

##Questions

<details> <summary>
How does SPAdes differ from other genome assembly programs?</summary>
It uses multiple values of k in de Bruijn graphs. Larger fragment sizes will more accurately position sections of duplicated DNA (repeats), but these larger fragments will only overlap well in densely-sequenced (high-coverage) areas of the genome. Because bacterial genomes may have low-coverage regions, using smaller fragments can increase the potential for overlaps (joins) in these low-coverage regions. Using a range of fragment sizes will therefore get the benefit from both approaches.</details>

<details> <summary>
How do I choose values of k?  </summary>
The k values need to be odd numbers, and shorter than the read lengths.  A good strategy could be to choose some that are small, medium and large. e.g. if the read is 150 nucleotides, k values could be 33, 55, 91. There is no absolute rule; rather, the aim is to get a good spread of k values across the read length. </details>

<details> <summary>
What can I do with my assembled genome?</summary>
This tutorial used a subset of a real dataset, so is not a complete genome (is it?). You could re-try it with short reads from a whole genome, at NCBI SRA. You can [annotate](../annotation_viz/annotation_viz.md) (describe) the genomic features such as genes or [compare](../mauve/mauve.md) it to other genomes to see variation in structure.  </details>

## More information

http://thegenomefactory.blogspot.com.au/2013/08/how-spades-differs-from-velvet.html

## Next
- Annotation
