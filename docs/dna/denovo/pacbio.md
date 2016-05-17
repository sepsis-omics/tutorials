# PacBio SMRT Portal

**Bacterial de novo assembly with PacBio data**

PacBio is a company that uses a technology called single molecule real time (SMRT) sequencing. This produces very long sequencing reads. These raw sequencing reads can be assembled into genomes by using the software in the PacBio SMRT Portal.

link: http://www.pacb.com/smrt-science/

link: http://www.pacb.com/products-and-services/analytical-software/smrt-analysis/

info: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Large-Genome-Assembly-with-PacBio-Long-Reads


## Background:
- de novo assembly
- mGVL
- Galaxy

## Start
- Open your GVL dashboard.

- Go to Admin. There is a list of packages. Find SMRT Analysis. On the right, click on "Install". (Note: you will need 16 cores in your GVL: FIXME: so will people be able to do this? Is this the plan?).

- You should see SMRT Portal as one of the instance services on your GVL dashboard.

- Open up the SMRT portal web link (to the right) and register/log on.

## Input data

- Choose your data. (FIXME: e.g. on GenomeSpace?)

- Get the data you want to use onto your mGVL. (FIXME: explain how?)

- In the SMRT Portal, go to **Design Job**, the top left tab.

- Go to **Import and Manage: Import SMRT cells: SMRT Cells**. Work out where you put the data on your GVL, and make sure the file path is showing. If not, click **Add** and list the file path to the data.

- Click on the file path and then **Scan** to check for new data.

## Run

- Go back to the top tab **Design Job**.

- Go to **Create New**.

- An **Analysis** window should appear. Check the box next to **De novo assembly**, then **Next**.

- Under **Job Name** enter a name.

- Under **Protocols** choose **RS_HGAP_Assembly.3**.

![smrt portal screenshot](/media/screenshots/smrt1.png)

- There is an ellipsis underneath **Protocols** - click on the ellipsis. This brings up the settings. Leave everything as is, except for: Click on **Assembly**. Change the **Genome Size** to an approximately correct size for the sample. Click **Ok**.  

![smrt portal screenshot](/media/screenshots/smrt2.png)

- In the **SMRT Cells Available** window, select the file to be used. Click on the arrow to transfer these files to the **SMRT Celles in Job** window.

![smrt portal screenshot](/media/screenshots/smrt3.png)

- Click **Save**.

- Next to **Save**, click **Start**.

- The **Monitor Jobs** window should open. As each step proceeds, new items will appear under the **Reports** and **Data** tabs on the left.

![smrt portal screenshot](/media/screenshots/smrt4.png)

- FIXME: how long will it take for this example data.

- FIXME: are any reports important to look at /check during the run.

- [information about report files and what they mean - FIXME: maybe expand on this](http://files.pacb.com/software/smrtanalysis/2.3.0/doc/smrtportal/help/Webhelp/SMRT_Portal.htm)

## Outputs

- The current running jobs will be under the **Monitor Jobs** tab. Click on the job to see the reports and data.

- The finished jobs will be under the **View Data** tab.

- FIXME: which reports are important to look at.

- FIXME: what are the files under **Data** for - further analyses later? where would these be saved if we want to use later.

- When the assembly finishes, look at the **View Data** tab for all the reports and ?data.

## Links to more information:

https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes
