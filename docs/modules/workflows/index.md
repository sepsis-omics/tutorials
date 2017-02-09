<br>
# Galaxy workflows

This tutorial assumes you have used Galaxy before.

Although we can use tools in Galaxy to analyse data and create a history, there is also a way to create a workflow of files, tools, settings and outputs.

You can then input different datasets and run the workflow.

## Start

Go to your Galaxy instance and Register/Login.

Import a history of data files:

- Click on the <ss>History</ss> cog ![cog icon](images/image02.png)
- Select <ss>Import from File</ss>
- In the box called <ss>Archived History URL</ss>, paste in this link address to the Galaxy history of input files:

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/Galaxy_history_input_files.tar.gz</tt>

- Click <ss>Submit</ss>
- Wait a few seconds.
- Click on the "view all histories" button ![histories icon](images/histories.png)
- See if the Galaxy history has been imported: it will be called <fn>imported from archive: Data</fn>
- Above that pane, click on the <ss>Switch to</ss> button.
- Then click <ss>Done</ss> (in the top left corner).
- You should now have a list of five files in your current history.

- Re-name this history "Workflows".

## Build a workflow

We will write a workflow for assembling a bacterial genome with Spades and annotating it with Prokka.

In the top menu bar in Galaxy, click on "Workflow".

workflow_menu.png

Click on <ss>Create new workflow</ss>

create_new.png

Under "Workflow Name:" put in e.g. Reads to Annotation.

Click <ss>Create</ss>

This will bring up the "Workflow Canvas", a grid where you can arrange the workflow.

### Add inputs

In the Tools panel, click <ss>Inputs: Input datset</ss> (at the very top of the list).

A box will appear: drag it to the left and there will be another box underneath it. Drag this also to the left.

inputs.png

Click on the first box. Look in the right hand panel (now called "Details") and change the name of the Input dataset to R1.fastq. Press Enter for the change to be saved.
R1fastq.png

Repeat for the second input dataset box, naming that one R2.fastq.

### Add the tool "spades"

In the tools panel, click on NGS Analysis: NGS Assembly: spades
This brings up the spades tool box:

spades.png

Click on that box and look in the Details pane on the right. This shows all the options in spades. Choose:

- <ss>Run only Assembly</ss>: *Yes* [the *Yes* button should be darker grey]
- <ss>Kmers to use separated by commas:</ss> *33,55,91*  [note: no spaces]  
- <ss>Coverage cutoff:</ss> *auto*  

### Join inputs to the tool

Now tell spades which input files to use. Look at the input dataset box called R1.fastq and find the small arrow: >

Click on this and drag the arrow over to the spades box input arrow > next to "Libraries 1 > Files 1 > Forward reads".

Repeat for R2.fastq, joining to "Reverse reads":

join.png

### Save it and run

Click on the cog at the top right of the workflow canvas and "Save".

save.png

Click the cog again and choose "Run".

Under <ss>Step1: Input dataset</ss> choose <fn>mutant_R1.fastq</fn>.
Under <ss>Step2: Input dataset</ss> choose <fn>mutant_R2.fastq</fn>.

Click <ss>Run workflow</ss>.

This will run the workflow (spades) and save the output to the top of your current history in the right hand panel.

View some of the files with the eye icon to check that the workflow (in this case, just spades) ran correctly.

## Add annotation to the worfklow

Go to the top panel and click "Workflow".

Your workflow "Reads to Annotation" should be in the list. Click on the drop-down arrow next to this workflow and choose "Edit". This will bring up the Workflow Canvas where we can add more inputs and tools.

Click on prokka.

join spades output: contigs to prokka input: contigs.

prokka.png

Click on the prokka box and change some of the settings in the right hand "Details" panel:

- Set the following parameters (leave everything else unchanged):
    - <ss>Locus tag prefix (--locustag)</ss>: P
    - <ss>Force GenBank/ENA/DDJB compliance (--compliant)</ss>: *No*
    - <ss>Sequencing Centre ID (--centre)</ss>: V
    - <ss>Use genus-specific BLAST database</ss> *No*  


Click on the cog to save.

Click on the cog to run. Again choose the R1 and R2 input files and Run workflow.

## Add visualization to the workflow

Add JBrowse.
