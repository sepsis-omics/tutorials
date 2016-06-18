# Dataset

## Galaxy histories

* [Galaxy history of input files](Data.tar.gz)
* [Galaxy history: FastQC](FastQChistory.tar.gz)
* [Galaxy history: Spades](Spadeshistory.tar.gz)
* [Galaxy history: Prokka](Prokkahistory.tar.gz)
* [Galaxy history: Snippy](Snippyhistory.tar.gz)

## Individual input files
###Wildtype reference

* [wildtype.fna](wildtype.fna)
* [wildtype.gbk](wildtype.gbk)
* [wildtype.gff](wildtype.gff)

### Mutant Illumina sequence

* [mutant_R1.fastq.gz](mutant_R1.fastq.gz)
* [mutant_R2.fastq.gz](mutant_R2.fastq.gz)

## Get tutorial data into Galaxy

To get the saved tutorial history (a set of files) into Galaxy:

- Right-click on <fn>Galaxy history of input files</fn> above and copy link address.
- Go to your Galaxy instance. Make sure you are registered and logged in. Refresh the page.
- Click on the <ss>History</ss> cog ![cog icon](images/image02.png)
- Select <ss>Import from File</ss>

![history options](images/image03.png)

- In the box called <ss>Archived History URL</ss>, paste in the link address to the Galaxy history.
- Click <ss>Submit</ss>
- Wait a few seconds.
- Click on the "view all histories" button ![histories icon](images/image11.png)
- See if the Galaxy history has been imported: it will be called <fn>imported from archive: Data</fn>
- Above that pane, click on the <ss>Switch to</ss> button.
- Then click <ss>Done</ss> (in the top left corner).
- You should now have a list of five files in your current history.

![files in galaxy history](images/datafiles.png)

## What next?

- If you are are working through the tutorials in a different order or want to see the completed history for another section, additional Galaxy histories are available above.
- Next: [Learn about quality control](../fastqc/index.md).
