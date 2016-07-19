#Assembly with PacBio data and SMRT Portal

This tutorial will show you how to assemble and finish a bacterial genome de novo, using the PacBio SMRT Portal on the mGVL.

## Start
- Open your mGVL dashboard.
- Go to Admin. There is a list of packages. Find SMRT Analysis. On the right, click on <ss>Install</ss>.
- You should see SMRT Portal as one of the instance services on your GVL dashboard.
- Open up the SMRT portal web link (to the right) and register/log on.



##PacBio files

FIXME: check/ add more detail here.

- what are SMRT cells?

Types of files:

- filename.bax.h5  x3
- subreads.fasta/fastq  x3
- filename.xml x1

What we need for analysis:

## Input
- Choose your data. (FIXME: e.g. on GenomeSpace?)
- Get the data you want to use onto your mGVL. (FIXME: explain how)
- In the SMRT Portal, go to <ss>Design Job</ss>, the top left tab.
- Go to <ss>Import and Manage: Import SMRT cells: SMRT Cells</ss>. Work out where you put the data on your GVL, and make sure the file path is showing. If not, click <ss>Add</ss> and list the file path to the data.
- Click on the file path and then <ss>Scan</ss> to check for new data.

## Run

The Hierarchical Genome Assembly Process (HGAP) is one protocol used to assemble bacterial genomes from PacBio data (long reads).

- Go back to the top tab <ss>Design Job</ss>.
- Go to <ss>Create New</ss>.
- An <ss>Analysis</ss> window should appear. Check the box next to <ss>De novo assembly</ss>, then <ss>Next</ss>.
- Under <ss>Job Name</ss> enter a name.
- Under <ss>Protocols</ss> choose <ss>RS_HGAP_Assembly.3</ss>.

![smrt portal screenshot](/media/screenshots/smrt1.png)

- There is an ellipsis underneath <ss>Protocols</ss> - click on the ellipsis.

This brings up the settings.

- This protocol is a collection of other protocols for filtering, assembly, mapping and consensus. (These are the .xml files).

- Filtering:
    - take a long circular polymerase PacBio read
    - adapters will be removed automatically I think?
    - a "subread" is the read left between two adapters
    - set the <ss>Minimum Polymerase Read Length</ss>: *default*
        - note: this is listed last but will be filtered first?
        - default is shorter than subreads so disregarded anyway?
    - set the <ss>Minimum Polymerase Read Quality</ss>: *default* (0.80)
    - set the <ss>Minimum Subread Length</ss>: *default*

- Assembly:
    - Pre-assembly:
        - maps subreads to long subreads of chosen length
        - For <ss>Compute Minimum Seed Read Length</ss>: ensure box is ticked -- it will compute automatically?
        - Does this split the reads into longer seed reads used for making consensus preassembled reads; and shorter reads that will be used in the next Mapping stage? Or is the seed read actually quite long and is used for subreads to just map to? (in which case - what length are reads being split on for sending to BLASR mapping?)
        - For <ss>Number of Seed Read Chunks</ss>: enter *12*
        - uses BLASR to map subreads to each other to produce consensus "preassembled reads"

    - Assembly:
        - is this stage = Celera?
        - Celera uses Overlap-Layout-Consensus?
        - Change the <ss>Genome Size</ss> to an approximately correct size for the species.
        - For <ss>Target Coverage</ss>: enter *10*
        - For <ss>Overlapper Error Rate</ss>: enter *0.04*
        - Click <ss>Ok</ss>.  

FIXME: add new screenshot showing options selected/changed

- Mapping:
    - uses BLASR to map the shorter subreads (that were not used for assembly) back to the assembly
    - leave defaults
    - ok

- Consensus:
    - uses Quiver to call a consensus from the previous mapping step.
    - uses original quality values from the sequencing reads.
    - ok


- In the <ss>SMRT Cells Available</ss> window, select the file to be used. Click on the arrow to transfer these files to the SMRT Cells in Job window.

![smrt portal screenshot](/media/screenshots/smrt3.png)

- Click <ss>Save</ss>.
- Next to <ss>Save</ss>, click <ss>Start</ss>.
- The <ss>Monitor Jobs</ss> window should open.
    - As each step proceeds, new items will appear under the <ss>Reports</ss> and <ss>Data</ss> tabs on the left.
    - Click on each of these items to see the details and graphs available, which will appear in the main pane.
    - The default display in the main pane is "Overview".

![smrt portal screenshot](/media/screenshots/smrt6.png)

##Output

- The current running jobs will be under the <ss>Monitor Jobs</ss> tab. Click on the job to see the reports and data.
- The finished jobs will be under the <ss>View Data</ss> tab.


General: Filtering:

- we removed short polymerase reads - graph shows read lengths post-filter
- we removed low quality reads - graph starts at 0.8 level
- anything to check here? probably not?
- look at general stats - e.g. number of reads left post-filter. add screenshot


![smrt portal screenshot](/media/screenshots/smrt5.png)

General: subread filtering:

- we removed short subreads - graph shows subread length post-filter (I assume)
- stats: see number of subreads left; average length.

Assembly: Pre-Assembly:

- stats: length cutoff (this split the reads: long => assembly; short => mapping)
- pre-assembled reads - the number of reads

Assembly: Polished assembly:

- graph: confidence vs depth (for each contig):
-  look to see if you can see chromosomes and plasmids. e.g. chromosomal contigs should cluster together; plasmid ones may have higher coverage if they are multi copy (eg further along the x axis).

Assembly: Corrections:

- expect a lot of corrections but should be random
- if all in one spot, could be a problem

Assembly: Top corrections:

- if > one page = a lot. Need to do another correction.

##RS_Resequencing
= BLASR + Quiver

##RS_Bridgemapper
= BLASR + Quiver + Bridgemapper  
=> view using SMRTview



## Links to more information:


- [A full ist of reports and terminology is here](http://files.pacb.com/software/smrtanalysis/2.3.0/doc/smrtportal/help/Webhelp/SMRT_Portal.htm)



- [Finishing bacterial genomes](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes)


[Link to PacBio analysis software](http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)

Link to info:
https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Large-Genome-Assembly-with-PacBio-Long-Reads
