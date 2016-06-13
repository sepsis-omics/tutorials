# HGAP

this means SMRT Suite

Assembly with PacBio data

PacBio is a company that uses a technology called single molecule real time (SMRT) sequencing. This produces very long sequencing reads (up to xx?). These raw sequencing reads can be assembled into genomes by using the software from PacBio, called the SMRT Portal.

link: http://www.pacb.com/smrt-science/
link: http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/


info:
https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Large-Genome-Assembly-with-PacBio-Long-Reads


Background:
assembly
galaxy
cmdline?

Start
Open your GVL dashboard.

Go to Admin. There is a list of packages. Find SMRT Analysis. On the right, click on "Install".  (Note: you will need 16 cores in your GVL).

open up the link and log on (create a username and password I think)
FIXME

Sepsis Data
e.g. https://downloads-qcif.bioplatforms.com/bpa/sepsis/pacbio/

e.g. one sample
it will be split into three files ?
for each sample: need the bas file? and the xml

but not the subreads



bas.h5 Reference Guide (PDF):
Describes the main output files produced by the primary analysis pipeline: bas.h5,.1.bax.h5, .2.bax.h5, and .3.bax.h5. The bax.h5 files contain base call information from the sequencing run. The bas.h5 file is essentially a pointer to the three bax.h5 files.
Metadata Output Guide (PDF): Describes the file metadata.xml, which contains top-level information about the data, including what sequencing enzyme and chemistry were used, sample name, and other metadata.


which files to put over into sepsis gvl





Input data
from here [link]

Import and Manage
Import SMRT cells: SMRT Cells
shows file paths that are searched
add
(so should have made folder, put smrt cells data in there, and then tell it this path using "add")








How it works
The PacBio SMRT analysis software: de novo assembly.

link: http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/analysis-applications/de-novo-assembly/
Several options, but for example HGAP + BridgeMapper:
RS_HGAP Assembly.3
pre-assembly
de novo assembly with AssembleUnitig
finalise assembly with Quiver

RS_Bridgemapper
assesses assembly quality by comparing to a reference genome.






Run

give job name and comments
choose a protocol
click on the data you want (how to get in) and then arrow to transfer it to right pane
run
click on the monitor tab
as it runs, new items will appear on the left hand side under reports, e.g. starting with filtering
graphs will appear for some items, e.g. mapped subread length, what does it mean
new data sets will appear under on the left hand side under data eg polished assembly fastq.
click on top right corner - log (to check for anything?)

Output
when finished, go to view data tab; click on job name, open (what is SMRT view? another option)

shows all reports and data files
e.g. assembled genome

download any or leave in here?
main assembly parameters /reports - what to check
how do you know assembly is good



https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes

-- good detail in here


Next

Links to more information
