#Assembly with PacBio data and SMRT Portal

This tutorial will show you how to assemble a bacterial genome *de novo*, using the PacBio SMRT Portal on the mGVL. We will use an analysis pipeline called HGAP, the Hierarchical Genome Assembly Process.

<!-- FIXMEs:

- will SMRT portal be available on all training GVLs
- will students use existing SMRT portal registrations, or will each set up their own
- PacBio data to be loaded into SMRT portal(s) -->


## Start
- Open your mGVL dashboard.
- You should see SMRT Portal as one of the instance services on your GVL dashboard.
- Open up the SMRT portal web link (to the right) and register/log on.

## Data
We will use a publicly available data-set of PacBio reads from the bacteria *E. coli* (reference link below). This has already been added to SMRT Portal for this tutorial.

## Input
If using this tutorial during a training session, proceed to the next step (Assembly).

Otherwise:

- Load the PacBio data (your own, or the training *E. coli* dataset) onto your GVL.
- In the SMRT Portal, go to <ss>Design Job</ss>, the top left tab.
- Go to <ss>Import and Manage</ss>.
![smrt portal screenshot](images/image03.png)
- Click <ss>Import SMRT cells</ss>.
![smrt portal screenshot](images/image04.png)
- Work out where you put the data on your GVL, and make sure the file path is showing.
    - If not, click <ss>Add</ss> and enter the file path to the data.
- Click on the file path and then <ss>Scan</ss> to check for new data.

## Assembly
- In the SMRT Portal, go to the top left tab, <ss>Design Job</ss>.
- Go to <ss>Create New</ss>.
- An <ss>Analysis</ss> window should appear. Check the box next to <ss>De novo assembly</ss>, then <ss>Next</ss>.
- Under <ss>Job Name</ss> enter a name.
- Under <ss>Protocols</ss> choose <ss>RS_HGAP_Assembly.3</ss>.
- There is an ellipsis underneath <ss>Protocols</ss> - click on the ellipsis.

![smrt portal screenshot](images/smrt1.png)

This brings up the settings. Click on <ss>Assembly</ss>.

- For <ss>Compute Minimum Seed Read Length</ss>: ensure box is ticked
- For <ss>Number of Seed Read Chunks</ss>: enter *12*
- Change the <ss>Genome Size</ss> to an approximately correct size for the species.
- For <ss>Target Coverage</ss>: enter *10*
- For <ss>Overlapper Error Rate</ss>: enter *0.04*
- Leave all other settings as they are.
- Click <ss>Apply</ss>

 Your protocol window should look like this:

![smrt portal screenshot](images/image02.png)

- Click <ss>Ok</ss>.  

- In the <ss>SMRT Cells Available</ss> window, select the file to be used. Click on the arrow to transfer these files to the SMRT Cells in Job window.

![smrt portal screenshot](images/smrt3.png)

- Click <ss>Save</ss> (bottom right hand side).
- Next to <ss>Save</ss>, click <ss>Start</ss>.
- The <ss>Monitor Jobs</ss> window should open.
    - As each step proceeds, new items will appear under the <ss>Reports</ss> and <ss>Data</ss> tabs on the left.

![smrt portal screenshot](images/smrt6.png)

## Output

- Click on the top right tab, <ss>View Data</ss>.
    - Double click on the job name to open its reports.
- Click on different <ss>Reports</ss> in the left hand panel.
- Click on <ss>Assembly &rarr; Pre-Assembly</ss>
    - see *Length Cutoff* for the calcuated sequence read length that was used: shorter reads were mapped against reads of at least this length to make pre-assembled reads.
    - see number of *Pre-Assembled Reads*
    - see average *Pre-Assembled Reads Length*    
- Click on <ss>Assembly &rarr; Polished Assembly</ss>.
    - See how many *Polished Contigs* were found.
    - See the *Max Contig Length*.
    - Are there multiple contigs?
      - Are some of these plasmids?

<!--   
- Poor assembly?
      - Does the sample require different assembly parameters?
      - Does the sample require new sequencing? -->

##Further Polishing

During polishing, raw reads are used to correct the assembly.
During HGAP, the assembly was polished once but may need further corrections.

- From the previous step, Go to <ss>Data &rarr; Assembly &rarr; Polished Assembly</ss> and download the FASTA file by clicking on it.
    - Unzip the .gz file
- Go to <ss>Design Job &rarr; Import and Manage</ss> and click <ss>New</ss> on the bottom right hand side. Then, select that FASTA assembly file to upload.
    - creates a new reference.
- Go to <ss>Design Job &rarr; Create New</ss>
    - choose reference-based
    - Select protocol: RS_Resequencing.1
    - Leave all settings.
    - Select your reference from the drop down menu.
    - Click <ss>Save</ss> and <ss>Start</ss>.
- Examine the output assembly and repeat if necessary.

<!-- - When complete, see Reports.

    - Variants: how many found? if less than 2, does not need any more polishing.
    - If 2+ variants found, repeat the polishing step (including adding a new reference).
-->

<!-- ## Next
Correct with Illumina reads <link to tutorial> -->

## Links to more information

[PacBio *E. coli* data set](https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly)


[HGAP overview](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP)

[A full ist of reports and terminology](http://files.pacb.com/software/smrtanalysis/2.3.0/doc/smrtportal/help/Webhelp/SMRT_Portal.htm)

[Video overview of HGAP on SMRT portal](http://www.pacb.com/training/BacterialAssemblyandEpigeneticAnalysis/story.html)
