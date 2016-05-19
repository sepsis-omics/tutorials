# PacBio SMRT Portal

**Bacterial de novo assembly with PacBio data**

PacBio is a company that uses a technology called single molecule real time (SMRT) sequencing. This produces very long sequencing reads. These raw sequencing reads can be assembled into genomes by using the software in the PacBio SMRT Portal.

[PacBio link](http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)

## Background:
- de novo assembly
- mGVL
- Galaxy

## Start
- Open your GVL dashboard.
- Go to Admin. There is a list of packages. Find SMRT Analysis. On the right, click on <ss>Instal<ss>. (Note: you will need 16 cores in your GVL: FIXME: so will people be able to do this? Is this the plan?).
- You should see SMRT Portal as one of the instance services on your GVL dashboard.
- Open up the SMRT portal web link (to the right) and register/log on.

## Input data
- Choose your data. (FIXME: e.g. on GenomeSpace?)
- Get the data you want to use onto your mGVL. (FIXME: explain how?)
- In the SMRT Portal, go to <ss>Design Job</ss>, the top left tab.
- Go to <ss>Import and Manage: Import SMRT cells: SMRT Cells</ss>. Work out where you put the data on your GVL, and make sure the file path is showing. If not, click <ss>Add</ss> and list the file path to the data.
- Click on the file path and then <ss>Scan</ss> to check for new data.

## Run
- Go back to the top tab <ss>Design Job</ss>.
- Go to <ss>Create New</ss>.
- An <ss>Analysis</ss> window should appear. Check the box next to <ss>De novo assembly</ss>, then <ss>Next</ss>.
- Under <ss>Job Name</ss> enter a name.
- Under <ss>Protocols</ss> choose <ss>RS_HGAP_Assembly.3</ss>.
![smrt portal screenshot](/media/screenshots/smrt1.png)
- There is an ellipsis underneath <ss>Protocols</ss> - click on the ellipsis. This brings up the settings. Leave everything as is, except for: Click on <ss>Assembly</ss>. Change the <ss>Genome Size</ss> to an approximately correct size for the sample. Click <ss>Ok</ss>.  
![smrt portal screenshot](/media/screenshots/smrt2.png)
- In the <ss>SMRT Cells Available</ss> window, select the file to be used. Click on the arrow to transfer these files to the SMRT Celles in Job** window.
![smrt portal screenshot](/media/screenshots/smrt3.png)
- Click <ss>Save</ss>.
- Next to <ss>Save</ss>, click <ss>Start</ss>.
- The <ss>Monitor Jobs</ss> window should open. As each step proceeds, new items will appear under the <ss>Reports</ss> and <ss>Data</ss> tabs on the left.
![smrt portal screenshot](/media/screenshots/smrt4.png)
- FIXME: how long will it take for this example data.
- FIXME: are any reports important to look at /check during the run.
- [information about report files and what they mean - FIXME: maybe expand on this](http://files.pacb.com/software/smrtanalysis/2.3.0/doc/smrtportal/help/Webhelp/SMRT_Portal.htm)

## Outputs
- The current running jobs will be under the <ss>Monitor Jobs</ss> tab. Click on the job to see the reports and data.
- The finished jobs will be under the <ss>View Data</ss> tab.
- FIXME: which reports are important to look at.
- FIXME: what are the files under <ss>Data</ss> for - further analyses later? where would these be saved if we want to use later.
- When the assembly finishes, look at the <ss>View Data</ss> tab for all the reports and ?data.

## Links to more information:
- [Finishing bacterial genomes](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes)
