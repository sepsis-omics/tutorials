<br>
# Dataset

## Upload a Galaxy history

- Copy this link:
    - **<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/Galaxy_history_input_files.tar.gz</tt>**

- Go to your Galaxy instance. Make sure you are registered and logged in. Refresh the page.
- Click on the <ss>History</ss> cog ![cog icon](images/image02.png)
- Select <ss>Import from File</ss>

![history options](images/image03.png)

- In the box called <ss>Archived History URL</ss>, paste in the link address to the Galaxy history (that you copied above).
- Click <ss>Submit</ss>
- Wait a few seconds.
- Click on the <ss>view all histories</ss> button ![histories icon](images/image11.png)
- See if the Galaxy history has been imported: it will be called <fn>imported from archive: Data</fn>
- Above that pane, click on the <ss>Switch to</ss> button.
- Then click <ss>Done</ss> (in the top left corner).
- You should now have a list of five files in your current history. We will use these for the Genomics Workshop; or see below for additional files.

![files in galaxy history](images/datafiles.png)

## Additional Galaxy histories

If you are using only part of the Genomics Workshop, you can upload any required histories listed here. Follow the instructions above.

### Galaxy history: FastQC

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/FastQChistory.tar.gz</tt>

### Galaxy history: Spades

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/Spadeshistory.tar.gz</tt>

### Galaxy history: Prokka

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/Prokkahistory.tar.gz</tt>

### Galaxy history: Snippy

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/Snippyhistory.tar.gz</tt>

## Additional files

If you need individual files, you can upload any of the files listed here. The instructions are listed below.

###Wildtype reference

* wildtype.fna

<tt> https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/wildtype.fna</tt>

* wildtype.gbk

<tt> https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/wildtype.gbk</tt>

* wildtype.gff

<tt> https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/wildtype.gff</tt>

### Mutant Illumina sequence

* mutant_R1.fastq.gz

<tt> https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/mutant_R1.fastq.gz</tt>

* mutant_R2.fastq.gz

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/mutant_R2.fastq.gz</tt>

###Assembled contigs

* SPAdes_contigs.fasta

<tt>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Microbial_tutorials/SPAdes_contigs.fasta</tt>

### Upload to Galaxy

There are two ways to upload these files to Galaxy. You can either download to your local computer and upload to Galaxy, or you can tell Galaxy to directly upload the file from an external source.

Download and upload:

-  Download required file(s) to your computer.
-  From the Galaxy tool panel, click on <ss>Get Data &rarr; Upload File</ss>  
-  Click the <ss>Choose local file</ss> button  
-  Find and select the <fn>file</fn> you downloaded and click <ss>Open</ss>  
-  Set the <ss>Type</ss> correctly.  
-  Click the <ss>Start</ss> button.  
-  Once the progress bar reaches 100%, click the <ss>Close</ss> button  
- The file will now upload to your current history.

Or, tell Galaxy to find the file from an external source:

-  From the Galaxy tool panel, click on <ss>Get Data &rarr; Upload File</ss>  
-  Click the <ss>Paste/Fetch data</ss> button  
-  Paste the URL into the box.
-  Click the <ss>Start</ss> button.  
-  Once the progress bar reaches 100%, click the <ss>Close</ss> button  
- The file will now upload to your current history.
