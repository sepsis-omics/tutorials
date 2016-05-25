<span class="c3 c21"></span>

<span class="c3 c23">Assembling using SPAdes</span> {#h.9cr0qpqpazk2 .c12 .c25}
===================================================

<span></span>

<span class="c3 c17">Background</span> {#h.da5zfog6qur3 .c25 .c12}
--------------------------------------

<span class="c3 c22"></span>

<span>SPAdes is one of a number of </span><span class="c19">de
novo</span><span> assemblers that use short read sets (Illumina Reads)
as an input and the method for assembly is based on de Bruijn graphs.
For information about SPAdes see this </span><span
class="c8">[link](https://www.google.com/url?q=http://bioinf.spbau.ru/spades&sa=D&ust=1464150946842000&usg=AFQjCNFTKOTw_OCIWN7RU8IKiwvrSDNH4Q){.c9}</span><span>.
A protocol for assembling with Velvet (another </span><span
class="c19">de novo </span><span>assembler) is available </span><span
class="c8">[here](https://www.google.com/url?q=https://docs.google.com/document/d/1xs-TI5MejQARqo0pcocGlymsXldwJbJII890gnmjI0o/pub&sa=D&ust=1464150946843000&usg=AFQjCNFk1fsIALzFm3Lwamirt0OfdErc6w){.c9}</span><span>.</span>

<span></span>

<span>In this activity, w</span><span>e will perform a </span><span
class="c19">de novo</span><span> assembly of a short read set (from an
Illumina sequencer) using the SPAdes assembler. The output from SPAdes
that we are interested in is a multifasta files that contains the draft
genome sequence.</span>

<span></span>

<span>The read set for today is from an imaginary </span><span
class="c19">Staphylococcus aureus</span><span> bacterium with a
miniature genome. </span>

<span></span>

<span>We have a closed, annotated genome sequence for a closely related
</span><span class="c19">wildtype</span><span> strain.</span>

<span></span>

<span>The whole genome shotgun method used to sequence our mutant strain
 The read set was produced on an Illumina DNA sequencing
instrument.</span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 109.20px; height: 59.50px;">![](images/image00.png)</span>

-   <span>The reads are paired-end</span>
-   <span>Each read is 150 bases (before trimming)</span>
-   <span>The number of bases sequenced is equivalent to 19x the genome
    sequence of the wildtype strain. (Read coverage 19x - rather low!)
    </span>

------------------------------------------------------------------------

<span class="c17 c3"></span> {#h.29ghte8svr23 .c12 .c20}
----------------------------

<span class="c17 c3">Section 1: Import data into Galaxy</span> {#h.2av488nv3p6q .c25 .c12}
--------------------------------------------------------------

<span></span>

<span>Summary of Section 1 activity:</span>

1.  <span>Log in to your Galaxy server</span>
2.  <span>Import files required for the activity</span>
3.  <span>View imported files </span>

<span></span>

1.  <span class="c3">Go to the Galaxy Page</span>

<!-- -->

1.  <span class="c18">Web address: 43.240.98.1/galaxy</span>
2.  <span class="c8 c3">[Remind me how to
    logon](https://www.google.com/url?q=https://docs.google.com/document/d/1LAQvhIG8s-vv6T14bb8lGRkmoNha7E3bHf9kAgUwMs0/pub&sa=D&ust=1464150946852000&usg=AFQjCNE64WW4yg7wALRjj8DmbYeiixh02w){.c9}</span>

<span class="c3">        </span>

<span class="c3"></span>

2.  <span class="c3">Import Files to Galaxy</span>

<!-- -->

1.  <span>Click on the </span><span class="c3">Analyze
    Data</span><span> menu at the top of the page.</span>
2.  <span>Click on the History menu button (the </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 22.00px; height: 20.00px;">![Screen
    Shot 2015-09-17 at 1.34.01
    pm.png](images/image02.png)</span><span>on the top right of the
    history pane)</span>
3.  <span>Click </span><span class="c3">Import from File
    </span><span>(at the bottom of the list)</span>
4.  <span>A new page will appear with a text box for the URL of the
    history to import.</span>

<!-- -->

1.  <span>Copy the following URL into the text box</span>

<span></span>

<span
class="c8">[http://43.240.98.1/public/dieter/Galaxy-History-Colombiaworkshopstart.tar.gz](https://www.google.com/url?q=http://43.240.98.1/public/dieter/Galaxy-History-Colombiaworkshopstart.tar.gz&sa=D&ust=1464150946856000&usg=AFQjCNEjS_bqNpCv6YF6F48n657lkRIwqg){.c9}</span>

<span></span>

<span></span>

2.  <span>Click </span><span class="c3">Submit</span>

<span>Galaxy will download the data files from the internet and will be
available as an additional history (takes about one minute).</span>

<span class="c3"></span>

3.  <span class="c3">To make the newly imported history appear as the
    current history</span>

<!-- -->

1.  <span>Click on the View all Histories button (the </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 26.00px; height: 22.00px;">![Screen
    Shot 2015-09-17 at 1.39.45
    pm.png](images/image01.png)</span><span> on the top right of the
    history pane.)</span>

<span>If the history has finished downloading it will appear as
“</span><span class="c3 c5 c24">imported from archive:
Colombia\_workshop\_start</span><span>“</span>

2.  <span>Click on the </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 73.00px; height: 21.00px;">![](images/image06.png)</span><span
    class="c3"> </span><span> button above the “</span><span
    class="c19">imported from archive:
    Colombia\_workshop\_start</span><span>” then the </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 50.00px; height: 20.00px;">![](images/image05.png)</span><span> button.</span>

<span></span>

<span>You should now have 4 files in the history pane as follows:</span>

<span></span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 259.00px; height: 156.00px;">![Screen
Shot 2015-09-17 at 1.41.07 pm.png](images/image07.png)</span>

<span></span>

<span class="c3">4.        View imported Files</span>

<span>All the files are are text files </span>

<span class="c6 c3 c5">mutant\_R1.fastq</span><span
class="c3 c28"> </span><span>and</span><span
class="c3 c28 c29"> </span><span
class="c6 c3 c5">mutant\_R2.fastq</span><span> a pair-end read
set</span>

<span class="c3 c5 c6">Wildtype.fna</span><span> a file that contains
the genome sequence of the wildtype strain in fasta format (a header
line, then the nucleotide sequence of the genome)</span>

<span class="c6 c3 c5">Wildtype.gff</span><span> a file that contains
the genome sequence of the wildtype strain in general feature format. (a
list of features - one feature per line, then the nucleotide sequence of
the genome)</span>

<span></span>

<span class="c3">Look at the contents of these file</span><span> </span>

1.  <span>Click on the View Data button (the </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 23.00px; height: 20.00px;">![](images/image04.png)</span><span>)
    next to each of the files in turn.</span>

<span></span>

<span>Brief Discussion about the GFF format</span>

<span></span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 601.70px; height: 336.00px;">![](images/image08.png)</span>

<span></span>

<span>5.        </span><span class="c3">Evaluate the quality of the read
set</span>

<span>How good is my read set? Do I need to ask for a new sequencing
run? Is it suitable for the analysis I need to do?</span>

<span></span>

<span
class="c3 c11">[FastQC](https://www.google.com/url?q=http://130.56.248.135/galaxy/tool_runner?tool_id%3Dtoolshed.g2.bx.psu.edu%252Frepos%252Fiuc%252Ffastqc%252Ffastqc%252F0.53&sa=D&ust=1464150946869000&usg=AFQjCNFM7cN7deOwzeheFHnk3tTRpurQvQ){.c9}</span><span> is
a tool that runs a standard series of tests on your read set and returns
a relatively easy to interpret report.</span>

<span></span>

<span>We will use the FASTQC tool in Galaxy to evaluate the quality of
one of our fastq files.</span>

1.  <span>Open the FASTQC tool interface (on the tool pane -
    </span><span class="c3">NGS: QC and Manipulation -&gt; </span><span
    class="c3 c11">[FastQC: Comprehensive
    QC](https://www.google.com/url?q=http://130.56.248.135/galaxy/tool_runner?tool_id%3Dtoolshed.g2.bx.psu.edu%252Frepos%252Fiuc%252Ffastqc%252Ffastqc%252F0.53&sa=D&ust=1464150946871000&usg=AFQjCNEwAdeMDZHDTA9RTT_a8g6q2Dn4hw){.c9}</span><span>)</span>
2.  <span>Select </span><span
    class="c6 c3 c5">mutant\_R1.fastq</span><span> and click
    </span><span class="c3">Execute</span>
3.  <span>Once finished, examine the output (Hint: </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 23.00px; height: 20.00px;">![](images/image04.png)</span><span>).
    It has a summary at the top of the page and a number
    of graphs.</span>

<span></span>

<span>Brief Discussion about the FastQC output</span>

<span></span>

<span>Some of the important outputs of FastQC for our purposes
are:</span>

-   <span>Read length - Will be important in setting maximum k-mer size
    value for assembly</span>
-   <span>Quality encoding type - Important for quality trimming
    software</span>
-   <span>% GC - High GC organisms don’t tend to assemble well and may
    have an uneven read coverage distribution.</span>
-   <span>Total number of reads - Gives you an idea of coverage..</span>
-   <span>Dips in quality near the beginning, middle or end of the
    reads - Determines possible trimming/cleanup methods and parameters
    and may indicate technical problems with the sequencing
    process/machine run.</span>
-   <span>Presence of highly recurring k-mers - May point to
    contamination of reads with barcodes, adapter sequences etc.</span>
-   <span>Presence of large numbers of N’s in reads - May point to poor
    quality sequencing run. You need to trim these reads to
    remove N’s.</span>

<span></span>

<span>We won’t be doing anything to these data to clean it up as there
isn’t much need. Therefore we will get on with the assembly!</span>

<span></span>

<span class="c17 c3">Section 2: Fastq reads assembled with SPAdes</span> {#h.3xqhbjrbggh8 .c25 .c12}
------------------------------------------------------------------------

<span></span>

<span>In this section we will perform a </span><span class="c19">de
novo</span><span> assembly of the mutant fastq reads into long
contiguous sequences (in fasta format.) SPAdes produces both contigs and
scaffolds. Ask your demonstrator if you would like to know the
difference between contigs and scaffolds.</span>

<span></span>

1.  <span class="c3">Perform an assembly</span>

<!-- -->

1.  <span>Open the SPAdes assembler tool interface (on the tool pane -
    </span><span class="c3">NGS: Assembly -&gt;
    spades</span><span>)</span>
2.  <span>Set the following parameters:</span>

<!-- -->

1.  <span> </span><span class="c3">Run only Assembly</span><span>:
    </span><span class="c19">Yes</span>
2.  <span> </span><span class="c3">Kmers to use separated by
    commas:</span><span> </span><span class="c19">33,55,91</span>
3.  <span> </span><span class="c3">Coverage
    cutoff:</span><span> </span><span class="c19">auto</span>
4.  <span> </span><span class="c3">Forward
    reads:</span><span> </span><span
    class="c6 c3 c5">mutant\_R1.fastq</span>
5.  <span> </span><span class="c3">Reverse
    reads:</span><span> </span><span
    class="c6 c3 c5">mutant\_R2.fastq</span>

<span></span>

<span>Your tool interface should look like this:</span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 474.04px; height: 474.70px;">![Screen
Shot 2015-09-17 at 2.01.38 pm.png](images/image03.png)</span>

3.  <span>Click </span><span class="c3">Execute</span>

<span class="c3"></span>

2.  <span class="c3">Examine the SPAdes output files.</span>

<span>Galaxy is now running SPAdes on the reads for you. When it is
finished, you will have 5 new files in your history. Fasta files of the
resulting contigs and scaffolds, some statistics on each and the SPAdes
logfile.</span>

1.  <span>Click on the View Data button (</span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 23.00px; height: 20.00px;">![](images/image04.png)</span><span>)
    on each of the files. Note that the short reads have been assembled
    into much longer contigs. The stats files will give you the length
    of each of the contigs.</span>

<span></span>
