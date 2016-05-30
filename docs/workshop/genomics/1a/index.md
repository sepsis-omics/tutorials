# Introduction to Galaxy

## Background

Galaxy is a web-based analysis and workflow platform designed for biologists to analyse their own data. It can be used to run a variety of bioinformatics tools. The selection of bioinformatics tools installed on the Galaxy instance we are using today caters for the analysis of bacterial genomics data sets.

Bioinformatics tools can be added from the Galaxy ‘tool shed’ or removed as necessary from your Galaxy instance.

Galaxy is an open, web-based platform. Details about the project can be found [here](https://galaxyproject.org/).

The Galaxy interface is separated into three parts. The <ss>Tools</ss> list on the left, the <ss>Viewing</ss> panel in the middle and the analysis and data <ss>History</ss> on the right. We will be looking at all three parts in this tutorial.

![galaxy overview screenshot](images/image05.png)

This activity will familiarise you with the Galaxy interface. It will cover the following operations:

-   Logging in to the server
-   Putting data onto Galaxy
-   Using some common tools

## Learning Objectives

At the end of this tutorial you should be able to:

1.  Register and login to a Galaxy server.
2.  Upload data to a Galaxy server from:
    -   A file on your local computer.
    -   A file on a remote datastore with an accessible URL.  
3.  Use tools in Galaxy by:
    -   Accessing the tool via the tool menu.
    -   Using the tool interface to run the particular tool.
    -   Viewing/accessing the tool output.

## Login to Galaxy

- Open a new tab or window on your web browser.
- Use Firefox or Chrome - Please don’t use Internet Explorer or Safari.
- Type in the following address: 43.240.98.1/galaxy (FIXME or your personal mGVL instance / sepsis mGVL / TBA - if this changes, change screenshot too)

&nbsp;

![Galaxy URL](images/image01.png)

&nbsp;

Click on <ss>User</ss> button on the right and either register or login.

&nbsp;

![Register or Login screenshot](images/image04.png)

&nbsp;

If you haven't yet registered, <ss>Register:</ss>

- Select: <ss>User &rarr; Register</ss>
- Enter your email, choose a password, repeat it and add a (all lower case - FIXME:why?) one word name
- Click <ss>Submit</ss>

If you have already registered, <ss>Login:</ss>

- Select: <ss>User &rarr; Login</ss>
- Enter your username & password
- Click <ss>Submit</ss>

## Put data onto Galaxy

There are two main ways to put your data onto Galaxy; this section will run through both ways. First, we need to make a new history.

### Make a new history

Note: Make a new folder to store the work we are about to perform.

- Click on the history menu button ![history icon](images/image02.png) at the top of the <ss>History</ss> panel.
- Select <ss>Create New</ss>
- Click on <ss>Unnamed history</ss> to rename.

### Datatypes

What sort of file is being uploaded?

We need to tell Galaxy what sort of file is being uploaded. Some common datatypes (file formats) are: text, fasta, fastq, vcf, GFF, GenBank, tabular. (FIXME: determine correct format (e.g. capitals) and link to file formats page)

### Upload a file from your own computer

With this method you can get most of the files on your own computer into Galaxy.

#### Download the following file to your computer:

- Copy this URL and paste it into the address bar in your web browser: <https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/galaxy101/Contig_stats.txt.gz>  
This is a file stored on the cloud. Press <ss>Enter</ss>, and the file will download (FIXME: on Mac but others?). Note the download location.

&nbsp;

![Browser address bar](./images/image06.png)

&nbsp;

#### Upload the file to Galaxy

-  From the Galaxy tool panel, click on <ss>Get Data &rarr; Upload File</ss>  
-  Click the <ss>Choose local file</ss> button  
-  Find and select the <fn>Contig_stats.txt.g</fn> file you downloaded and click <ss>Open</ss>  
-  Set the <ss>Type</ss> to *tabular*  
-  Click the <ss>Start</ss> button  
-  Once the progress bar reaches 100%, click the <ss>Close</ss> button  
- The file will now upload to your current history.

### Upload a file from a URL

If a file exists on a web resource somewhere and you know its URL (Unique Resource Location - a web address) you can directly load it into Galaxy.

- From the tool panel, click on <ss>Get Data &rarr; Upload File</ss>
- Click on the <ss>Paste/Fetch Data</ss> button
- Copy and paste the following web address into the URL/Text box:
<https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq.gz>
- Set the <ss>Type</ss> to *fastqsanger* (CAREFUL: not fastqCsanger)
- Click <ss>Start</ss>
- Once the progress bar has reached 100%, click <ss>Close</ss>.
- Note that Galaxy is smart enough to recognize that this is a compressed file and so it will uncompress it as it loads it.

### Upload another file from a URL

Now we are going to upload another file from the remote data source.

- Repeat the above for: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/MRSA0252.fna
- Note: this file <ss>Type</ss> is *fasta*.
- The DNA sequence of *Staphylococcus aureus MRSA252* will be loaded into your history as a fasta file.
- Your <ss>History</ss> should now look like this:

![History screenshot](./images/image03.png)

### The data

A brief description of each of the three files uploaded to Galaxy.  

<fn>Contigs_stats.txt</fn>  

- this file contains a table of summary data from a *de novo* genome assembly (the process of attempting to recover the full genome of an organism from the short read sequences produced by most DNA sequencing machines).
- The columns contain a lot of information but the ones we will be using indicate the amount of data (or coverage) that went into making up each piece of the final assembly.

<fn>bacterial_std_err_1.fastq.gz</fn>  

- This file contains sequence reads, in the format produced by Illumina sequencing machines. Read more about the
[fastq](https://en.wikipedia.org/wiki/FASTQ_format) format at Wikipedia.

<fn>MRSA0252.fna</fn>

  - This file contains the genome sequence of *Staphylococcus aureus MRSA252*. Read more about the [fasta](https://en.wikipedia.org/wiki/FASTA_format) format at Wikipedia.

## Galaxy tools

The purpose of this section is to help you become familiar with the way
tools are run on Galaxy.

We will see how to:

- rename files
- summarize assembly statistics
- convert file formats, and
- find features in a DNA sequence.

### Rename files

Two of the files in the <ss>History</ss> have very long and confusing names. File names can be changed by taking the following steps:

- Click on the edit icon ![edit icon](./images/image07.png) next to the file in the <ss>History</ss> called: <fn>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq</fn>
- In the <ss>Name</ss> text box, give it a new name. Rename it to: <fn>typical.fastq</fn>
- Click the <ss>Save</ss> button.

Repeat the process with another file:

- Find the file called: <fn>https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/MRSA0252.fna</fn>  
- Rename it to <fn>MRSA252.fna</fn>

Much better. There is other functionality hidden behind that edit icon ![edit icon](./images/image07.png)

You can change a file’s datatype, convert its format and many other things. Feel free to play around with them at a later date.

### Summarize assembly statistics

We are going to produce a histogram of contig read-coverage depths and calculate the summary statistics from the <fn>Contig_stats.txt</fn> file.

To do this we need to make some changes to the <fn>Contig_stats.txt</fn> file:

- cut out a couple of columns from the file
- remove a line from the file
- produce a histogram

#### Cut out and keep two columns

- Click on the eye icon ![eye icon](./images/image00.png) of the <fn>Contig_stats.txt</fn> file to have a look at it.  
- Note that there are 18 columns in this file. We want column 1 and column 6.
- Go to <ss>Tools &rarr; Text Manipulation &rarr; Cut</ss> and set the following:
- Set <ss>Cut Columns</ss> to: *c1,c6*
- <ss>Delimited by</ss>: *Tab*
- <ss>From</ss>: <fn>Contig_stats.txt</fn>
- Click <ss>Execute</ss>
- Examine the new file by clicking on its eye icon ![eye icon](./images/image00.png). We now have 2 columns instead of the 18 in the original file.

#### Remove the Header lines of the new file

- Go to <ss>Tools &rarr; Text Manipulation &rarr; Remove beginning</ss> and set the following:
- <ss>Remove First</ss>: *1*
- <ss>from</ss>: <fn>Cut on data1</fn>
- click <ss>Execute</ss>
- Note the the new file is the same as the previous one without the header line.

#### Make a histogram

- Go to <ss>Tools &rarr; Graph/Display Data &rarr; Histogram</ss> and set the following:
- <ss>Dataset</ss> <fn>Remove beginning on Data 4</fn>
- <ss>Numerical column for X axis</ss> Column: 2
- <ss>Number of breaks</ss> *25*
- <ss>Plot title</ss> *Histogram of Contig Coverage*
- <ss>Label for X axis</ss> *Coverage depth*
- Click <ss>Execute</ss>
- Click on the eye icon ![eye icon](./images/image00.png) of the histogram to have a look at it. Note there are a few peaks. Maybe these correspond to single, double and triple copy number of these contigs.

#### Calculate summary statistics for contig coverage depth

- Go to <ss>Tools &rarr; Statistics and Visualisation &rarr; Statistics &rarr; Summary Statisitics</ss> and set the following:  
- <ss>Summary statistics on</ss> <fn>Remove beginning on Data 4</fn>
- <ss>Column or expression</ss> *c2*
- Click <ss>Execute</ss>
- You’ll note that the summary statistics tool failed (red background in the <ss>History</ss>). There was an error!
- If you click on the filename, and then the bug symbol ![bug icon](./images/image08.png), it will tell you what went wrong. (There is a missing python library).
- At this point, you would normally contact your Galaxy server administrator.

### Convert file formats

This shows how to convert a fastq file to a fasta file. The tool creates a new file with the converted data.

- Go to <ss>Tools &rarr; Basic Tools &rarr; Convert Formats &rarr; FASTQ to FASTA</ss> and set the following:
- <ss>FASTQ file to convert</ss>: <fn>typical.fastq</fn>
- Click <ss>Execute</ss>
- The output is a new Fasta file called <fn>FASTQ to FASTA on data 2</fn>.

### Find features

This example shows how to use a tool called “barrnap” to search for rRNAs in a DNA sequence.

#### Find all of the ribosomal RNAs in a sequence

- Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: Annotation &rarr; barrnap</ss> and set the following:
- <ss>Fasta file</ss>: <fn>MRSA252.fna</fn>
- Click <ss>Execute</ss>
- The output is <fn>barrnap on data 3</fn> It is a gff3 format file. (general feature format version 3). Each line in the file describes a feature in the DNA sequence.

#### Filter the annotations to get the 23S RNAs

- Make a file with only the 23S rRNA features
- Go to <ss>Tools &rarr; Basic Tools &rarr; Filter and Sort &rarr; Select</ss> and set the following:
- <ss>Select lines from</ss>: (whatever you called the barrnap gff3 output)
- <ss>the pattern</ss>: *23S* (this will look for all the lines in the file that contain “23S”)
- Click <ss>Execute</ss>
- Now you have a gff3 file with just the 23S annotations!

## What now?

Remember how we started a new <ss>History</ss> at the beginning? If you want to see any of your old histories, click on the History menu button ![history button](./images/image02.png) at the top of the <ss>History</ss> panel and then select “Saved Histories.” This will give you a list of all the histories you have worked on in this Galaxy server.
