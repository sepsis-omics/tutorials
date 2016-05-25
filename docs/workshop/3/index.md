<span class="c10 c6"></span>

<span class="c13 c6 c26">Comparative Genomics</span> {#h.h35wz63prixt .c4 .c27}
====================================================

<span class="c17 c13 c6">Background</span> {#h.6q2o04s2nbfz .c4 .c19}
------------------------------------------

<span>In this activity we will use the BWA short read mapper and
Freebayes variant caller to identify ‘micro’ differences between genome
sequences</span><span class="c6">.</span>

<span></span>

<span>The genome sequences being compared are those of the ‘wildtype’
and ‘mutant’ strains.</span>

<span>\*\*The relevant files should already be available on
Galaxy\*\*</span>

<span></span>

<span>Just for a recap:</span>

<span></span>

<span>We have a closed, annotated genome sequence for the
</span><span>wildtype</span><span> strain. </span>

<span>Files: wildtype.fna and wildtype.gff</span>

<span>For the mutant strain we have whole genome shotgun short sequence
reads from an Illumina DNA sequencing instrument</span>

<span>        Files: mutant\_R1.fastq and mutant\_R2.fastq (fastq
format)</span>

<span></span>

-   <span>The reads are paired-end</span>
-   <span>Each read is 150 bases</span>
-   <span>The reads coverage depth is estimated at 19x.</span>

<span></span>

<span>After investigating the ‘micro’ differences (SNPs / INDELS), we
will attempt to detect larger ‘macro’ differences using Mauve.</span>

<span class="c6"></span>

<span class="c6"></span>

<span class="c17 c13 c6">Section 1: Mapping Reads to Reference</span> {#h.ia6szkwynb7f .c4 .c19}
---------------------------------------------------------------------

1.  <span class="c6">Map the reads on to the reference sequence</span>

<span>Several programs could be used for this but we will use
BWA-MEM.</span>

1.  <span>Open the BWA-MEM tool interface (from the tool panel
    </span><span class="c6">NGS: Mapping -&gt; Map with
    BWA-MEM</span><span>.)</span>
2.  <span>Set the following parameters:</span>

<!-- -->

1.  <span> </span><span class="c13 c6">Will you select a reference
    genome from your history or use a built-in index?: </span><span
    class="c1">Use a genome from history and build index</span>
2.  <span class="c13"> </span><span class="c6 c13">Use the following
    dataset as the reference sequence: </span><span
    class="c1">wildtype.fna</span>
3.  <span class="c13"> </span><span class="c13 c6">Select first set of
    reads: </span><span class="c1">mutant\_R1.fastq</span>
4.  <span class="c13"> </span><span class="c13 c6">Select second set of
    reads:</span><span class="c13"> </span><span
    class="c1">mutant\_R2.fastq</span>

<!-- -->

3.  <span class="c13">Click </span><span class="c13 c6">Execute</span>

<span></span>

<span class="c6"></span>

2.  <span class="c6">Examine the contents of the BAM file</span>

<span>NOTE: The BAM file is a Binary Compressed Datafile and cannot be
viewed directly. If you attempt to view it using the view data button it
will be downloaded to your local computer. Instead we must convert it to
a non-compressed text format (SAM) first.</span>

1.  <span>First we have to convert the BAM file to a SAM file.</span>
2.  <span>Open the BAM-to-SAM tool interface: (on the tool panel
    </span><span class="c6">NGS: SAM tools -&gt;
    BAM-to-SAM</span><span>.)</span>
3.  <span>View the resultant SAM file by clicking on the View
    Data button.</span>

<!-- -->

1.  <span>Have a look at the fields in the file. </span>

<span></span>

<span>The demonstrator will now point out what all the fields
are.</span>

<span></span>

<span class="c17 c13 c6">Section 2: Viewing the BAM file using Artemis</span> {#h.nc855y9q41uk .c4 .c19}
-----------------------------------------------------------------------------

<span>In this section we will view the BAM file produced we produced
above in Artemis.</span>

<span></span>

<span></span>

1.  <span class="c6">Download the BAM file to your local computer</span>

<!-- -->

1.  <span>Click on the name of the BAM file that you created in
    Section 1.</span>

<!-- -->

1.  <span>Click on the download button </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 27.00px; height: 26.00px;">![Screen
    Shot 2015-09-18 at 3.50.18 pm.png](images/image02.png)</span><span>,
    you need to download both the BAM file and the bam\_index.</span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 249.00px; height: 166.00px;">![](images/image03.png)</span>

<span></span>

<span></span>

2.  <span class="c6">Start Artemis and load the wildtype.gff</span>

<!-- -->

1.  <span>From the Artemis menu, Click </span><span class="c6">File
    -&gt; Open ...</span>

<!-- -->

1.  <span>Load </span><span class="c20">wildtype.gff</span>

<span></span>

<span>You should now have the wildtype’s annotated sequence loaded into
the Artemis genome browser.</span>

<span></span>

3.  <span class="c6">Load the BAM file into Artemis</span>

<!-- -->

1.  <span>Click </span><span class="c6">File -&gt; Read BAM / VCF</span>
2.  <span>Select: </span><span class="c20">Galaxy … .bam</span>
3.  <span>Click </span><span class="c6">Ok</span>

<span></span>

<span></span>

<span>You should see something like this:</span>

<span></span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 361.50px; height: 437.23px;">![](images/image00.png)</span>

<span></span>

<span>Can you find a SNP? </span>

<span>Demonstration of the ways that the view of the BAM file can be
enhanced!</span>

<span></span>

<span></span>

<span></span>

<span>Imagine finding each SNP manually -  Luckily this can be automated
using a tool available on Galaxy.</span>

<span></span>

<span class="c13 c6 c17">Section 3: Variant Calling</span>

<span></span>

<span>We will now search through our alignment file (BAM) for
statistically valid SNPs using the Freebayes variant calling
tool.</span>

<span></span>

1.  <span class="c6">Run Freebayes</span>

<!-- -->

1.  <span>Load the Freebayes tool interface (on the tool panel:
    </span><span class="c6">NGS: Variant Detection -&gt;
    Freebayes</span><span>.)</span>
2.  <span>Set the following parameters:</span>

<!-- -->

1.  <span> </span><span class="c6">Load reference genome
    from:</span><span> </span><span class="c20">History</span>
2.  <span> </span><span class="c6">Sample BAM file:</span><span
    class="c20"> Map with BWA-MEM on data … BAM format)</span>
3.  <span> </span><span class="c6">Use the following dataset as the
    reference sequence: </span><span class="c20">wildtype.fna</span>

<!-- -->

3.  <span>Click </span><span class="c6">Execute</span>

<span class="c6"></span>

2.  <span class="c6">Examine the Freebayes output</span>

<!-- -->

1.  <span>Freebayes will create a VCF file. This stands for Variant
    Calling Format.</span>
2.  <span>Click on its View Data button and have a look at the file.
    There is a lot of header information, the variants appear
    lower down.</span>
3.  <span>Can you spot a SNP?</span>
4.  <span>What about an insertion? A deletion?</span>

<span></span>

<span></span>

<span></span>

<span></span>

<span class="c17 c13 c6">Section 4: Investigation of Variants</span>

<span>What is the impact of the differences we have observed?</span>

<span></span>

<span>In this section we will use some simple investigation strategies
to predict the impact of the difference on the function of the gene and
perhaps even the strain itself.</span>

<span></span>

<span class="c6">Artemis</span><span> - the annotated draft genome
sequence of the mutant strain - what is the impact the protein coding
region? what is the predicted function? </span>

<span class="c6">blastp</span><span> -
(http://blast.ncbi.nlm.nih.gov/Blast.cgi)  the protein domain display -
are any major protein domains truncated by the difference?</span>

<span class="c6">LipoP/SignalP/TmHMM</span><span> - (</span><span
class="c22">[http://www.cbs.dtu.dk/services/](https://www.google.com/url?q=http://www.cbs.dtu.dk/services/&sa=D&ust=1464150964648000&usg=AFQjCNG3wnXby1GVaooVIBrx-flME1kq2Q){.c5}</span><span>)
membrane location prediction - has the change had an impact on the
membrane location of the protein?</span>

<span class="c6">Literature?</span>

<span></span>

<span>Can you suggest a type of nucleotide sequence that might have no
impact on the function of the encoded protein?</span>

<span></span>

<span>In this section we will investigate a few variants together as a
demonstration </span>

<span>\*\*perhaps a few individually too??\*\*</span>

<span class="c21 c6"></span>

<span class="c6 c8">Section 5: </span><span class="c6 c23">Detection of
‘macro’ INDELS and rearrangement using Mauve</span>

<span></span>

<span>We will now examine our earlier assembly and compare it with the
reference on a genome wide basis using Mauve.</span>

<span></span>

<span>Download and install Mauve. More information on Mauve and its use
can be found </span><span
class="c22">[here](https://www.google.com/url?q=http://darlinglab.org/mauve/mauve.html&sa=D&ust=1464150964651000&usg=AFQjCNEpl5L36SMHuxVG1WQuPyZYBnVrVA){.c5}</span><span>.</span>

<span></span>

<span>You will then need to load both the reference </span><span
class="c20">wildtype.gff</span><span> file and the mutant gff file that
you downloaded earlier.</span>

<span></span>

<span
style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 601.70px; height: 184.00px;">![](images/image01.png)</span>

<span></span>
