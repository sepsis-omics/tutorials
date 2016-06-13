# JBrowse

## Background

JBrowse is a genome browser available through the Galaxy interface or as a stand-alone tool.

http://jbrowse.org/

http://genome.cshlp.org/content/19/9/1630.full

Today we will use JBrowse within Galaxy.

A variety of file types can be viewed using JBrowse, including GFF3, BED, FASTA, BAM, VCF.

structure:
- top: a reference bar showing chromosome position
- feature tracks based on genome annotations:
   - particular annotations: e.g. tRNA, non coding RNA ?
   - all annotations in a row?

navigation:
- dragging
- use navigation buttons (left, right, zoom)
- type in coordinates or feature name

makefile:
- transform datatypes


- go to <ss>Graph/Display Data &rarr; JBrowse</ss>
- choose the ref sequence? only allows the fasta file? will track options then appear?
- chose wildtype.fna
but: 	python: can't open file 'jbrowse.py': [Errno 2] No such file or directory


## Learning objectives

At the end of this tutorial, you should be able to:

1. download a bacterial genome sequence, and
2. open the genome sequence in JBrowse
3. and look at the sequence features.

## Download a bacterial genome

We will download the sequence of *Leptospira borgpetersenii* serovar Hardjobovis Chromosome II from the NCBI website in GenBank format.  

- Go to <http://www.ncbi.nlm.nih.gov>
- Select the “Taxonomy” database from the list on the left hand side.
- Enter “Leptospira borgpetersenii” in the search box and click <ss>Search</ss>.

![NCBI search box](./images/image00.png)

- Click on the species name:

![Species name link](./images/image08.png)

&nbsp;

- Click on the species name again (at the top of the list):

![Species name link](./images/image06.png)

&nbsp;

- Next to <ss>Genome</ss>, click on the far right number <ss>1</ss>.

<!---FIXME: add arrow to point to it.  
--->

![link to genome](./images/image02.png)

&nbsp;

- Click on <ss>Genome Assembly and Annotation report</ss>:

&nbsp;

![link to report](./images/image04.png)

&nbsp;

- Then, for the first species listed, scroll to the far right of the table, see chromosome 2, and click the link to the annotated assembly number <ss>NC_008509.1</ss> (not CP0003):

![link to report](./images/image01.png)

&nbsp;

- Then, click <ss>Send</ss>
- choose <ss>Complete Record</ss>, <ss>Destination &rarr; File</ss>, <ss>Format &rarr; GenBank (full)</ss>
- click <ss>Create File</ss>.
- Note the download location (e.g. Downloads folder).

![Send button](./images/image05.png)

&nbsp;

## Open the GenBank file in JBrowse




Overview:  

- 6-frame annotation
- annotated genome features are highlighted
- black lines are stop codons

DNA view:

- 6-frame translation
- DNA sequence in the middle
- amino acid translations above and below

Text summary:

- text summary of features

Navigation:

- Go to the overview pane
- click on one annotated feature (highlighted in blue)
- it will be summarized on the top line ("selected feature")
- the corresponding sequence will be highlighted in the DNA view pane
- the corresponding feature will be higlighted in the text summary pane
- now double click on the same annotated feature
- all three panes will be centred for this feature
- to move left or right, use the horizontal scroll bars under each pane
- to zoom, use the vertical scroll bars on the right

## What next  
- [Assemble a bacterial genome using Spades.](../spades/index.md)
