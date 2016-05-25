<span class="c9 c13"></span>

<span class="c1">Genome Annotation using Prokka</span> {#h.9cr0qpqpazk2 .c2}
======================================================

<span></span>

<span class="c9 c19">Background</span> {#h.twbpq7rc0ic0 .c2}
--------------------------------------

<span class="c10">In this section we will use a software tool called Prokka to annotate the draft genome sequence produced in Activity 2a. Prokka is a “wrapper”, it collects together several pieces of software (from various authors) - avoids “re-inventing the wheel”. </span> {#h.qww5912pyrnk .c2}
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

<span class="c10">Prokka finds and annotates features (both protein coding regions and RNA genes i.e. tRNA, rRNA) present on on a sequence. Note: Prokka uses a two step process for the annotation of protein coding regions, first protein coding regions on the genome are identified using </span><span class="c16 c10">[Prodigal](https://www.google.com/url?q=http://prodigal.ornl.gov/&sa=D&ust=1464150956262000&usg=AFQjCNHsb-lcoZ7BADy23CpkSlxwsUDhuQ){.c6}</span><span class="c4 c14"> and the function of the encoded protein is predicted by similarity to proteins in one of many protein or protein domain databases. Prokka is a software tool that can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats. More information about Prokka can be found </span><span class="c4 c16">[here](https://www.google.com/url?q=https://github.com/tseemann/prokka&sa=D&ust=1464150956263000&usg=AFQjCNF70QjNHbyPrRo4eQzMVN7UKJezOw){.c6}</span><span class="c4 c14">.</span> {#h.xrw96u8qe6qz .c2}
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

<span class="c10"></span>

1.  <span class="c10 c9">Run Prokka</span>

<!-- -->

1.  <span class="c10">Load the Prokka interface from the tool panel
    (</span><span class="c10 c9">NGS: Annotation -&gt;
    Prokka</span><span class="c10">)</span>
2.  <span class="c10">Set the following parameters:</span>

<!-- -->

1.  <span class="c10 c9">Contigs to annotate: </span><span
    class="c10 c12">SPAdes contigs</span>
2.  <span class="c4 c9">Force GenBank/ENA/DDJB compliance (--compliant):
    </span><span class="c3">Yes</span>
3.  <span class="c4 c9">Genus Name: </span><span
    class="c3">Staphylococcus</span>
4.  <span class="c4 c9">Strain Name:</span><span
    class="c4"> </span><span class="c3">aureus</span>
5.  <span class="c4 c9">Use genus-specific BLAST database: </span><span
    class="c3">No</span>
6.  <span class="c4 c9">Locus tag prefix: </span><span
    class="c3">P</span>
7.  <span class="c4 c9">Sequencing Centre: </span><span
    class="c3">V</span>

<!-- -->

3.  <span class="c4">Click </span><span class="c4 c9">Execute</span>

<span class="c4"></span>

2.  <span class="c4 c9">Examine the output files</span>

<span class="c4">Once Prokka has finished, examine each of its output
files.</span>

1.  <span class="c4">The gff and gbk files contains all of the
    information about all of the features annotated (in
    different formats.)</span>
2.  <span class="c4">The txt file contains a summary of the number of
    features annotated.</span>
3.  <span class="c4">The faa file contains the protein sequences of the
    genes annotated.</span>
4.  <span class="c4">The ffn file contains the nucleotide sequences of
    the genes annotated.</span>

<span class="c4"></span>

3.  <span class="c4 c9">Download the gff file to your local
    computer</span>

<span class="c4">Now that we have annotated the draft genome sequence,
we would like to view the sequence in the Artemis genome viewer.</span>

5.  <span class="c4">The file is downloaded from the history panel.
    Identify the file with the .gff extension, expand file header to
    reveal the expanded file header and download the file using the
    </span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 25.00px; height: 20.00px;">![](images/image00.png)</span><span
    class="c4"> icon.</span><span
    style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 138.83px; height: 178.50px;">![](images/image01.png)</span>
6.  <span class="c4">Open Artemis and load the downloaded .gff file.
    </span>

<span class="c4"></span>

<span class="c4">Discussion - a closer look at the annotated features
</span>
