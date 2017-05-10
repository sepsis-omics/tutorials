<br>
# Metabolomics with XCMS

## Overview

Metabolomics is the study of metabolites: small molecules (smaller than proteins) produced by organisms. One technique to identify and quantify metabolites is LC/MS (liquid chromatography-mass spectrometry).

In liquid chromatography, metabolites are sent through a column and are separated by various chemical properties.

- At each time point, the abundance of the particular group of metabolites is measured (the intensity).
- This is called a chromatogram:

<img src="../images/chrom.png" width="500" height ="300" />

*image: C Wenger, Wikipedia*

At each time point in this chromatogram, the group of metabolites is then ionized (charged) and fired through a mass spectrometer.

- These are separated based on their mass-to-charge ratio.
- This adds another dimension to the graph: an axis with the mass-to-charge ratio of the metabolites found at that time point, and their intensities.

![peaks](images/peaks.png)

*image: Daniel Norena-Caro, Wikipedia.*


Each of these peaks is a "feature": an ion with a unique m/z and retention time:

Mass spectra data is often simplified into a graph of peaks only (local maxima):

<img src="../images/centroid.png" width="600" height ="250" />

These masses can then be matched to a database to identify the metabolite.

To summarize:

one chromatogram - from the LC
multiple specta - from the MS


A common aim is to compare metabolites between two samples (e.g. two different bacterial strains), and from there, to understand which biological pathways may be up- or down-regulated.

## XCMS

In this tutorial, we will use XCMS online to analyse metabolite data.

- Input: raw data from mass spectrometry
- Output: identifies the metabolites, and compares their abundances between samples.
- Go to: <https://xcmsonline.scripps.edu> and sign up.

![startpage](images/start_page.png)

## Get data

The data we will use today is from two bacterial strains of *Klebsiella pneumoniae*.

- strain AJ218 (check AMR) - x6 replicates
- strain KPC2 (antibiotic resistance) - x6 replicates

[Note: there may not be any metabolomic diffs between these. If not, look at Staph vs Klebs]

The raw data output from the mass spectrometer are points in a 3D graph: intensity, m/z, time (retention time).

Download these to your local computer from swift.[if data release policy allows. upload to swift and put link here].

Data format:
- We will use .mzML format. [this format is going to be uploaded to the Sepsis Data]
- The machine used to produced this data originally produced .d files. These have been converted into .mzML format using the Proteowizard MSConvert program.

To use Proteowizard: (Windows only)

-  Download the program and open the MSConvert program.
- Browse. Add files.
- Change output format to .mzML
- Leave default settings.
- Click Start in the bottom right hand corner.

An alternative is the commerially-available Qualitative Analysis Software.

## Upload data

In the top panel, go to <ss>Stored Datasets</ss>. We will upload some data here.

![datasets](images/stored_data.png)

In the top right corner, click on <ss>Add Dataset(s)</ss>.
<img src="../images/add_data.png" width="400" height ="200" />

- Drag the <fn>.mzML</fn> files into the <ss>Drop Here</ss> box. - 6 reps for each strain.

- Wait until all files have a green tick (scroll down to check all).
- Name the datset (e.g. Sample AJ218 or Sample KPC2) and click <ss>Save</ss>.
- Click <ss>Save Dataset & Proceed</ss>.


![upload](images/file_upload.png)


Repeat with the second strain.


## Set up job

In the top panel, click on <ss>Create Job</ss> and select <ss>Pairwise Job</ss>.

![pairwise](images/pairwise.png)

On the right hand side, under Job Summary, Job Name: enter job name. e.g. Klebsiella metabolomics.

Under <ss>Dataset 1</ss> click on <ss>Select Dataset</ss>.

- Choose Sample AJ218.

Under <ss>Dataset 2</ss> click on <ss>Select Dataset</ss>.

- Choose Sample KPC2.

We now need to set parameters that correspond with the machine on which the data was generated.

In a typical analysis, we would look at the raw output files and examine the chromatograms and mass spectra to inform some of the settings. Here we have chosen appropirate settings for this data set.

e.g. use SeeMS (part of proteowizard)

- internal standard - a known molecule added
- e.g. a labelled standard (lablled with C13 etc)
- use the expected and actual mass of labelled standard to calculate the appropriate error setting (ppm - setting )


Under <ss>Parameters</ss> select HPLC / UHD Q-TOF (HILIC, neg. mode).

- Click <ss>View/Edit</ss>.
- This brings up a window to change some settings.
- First, click <ss>Create New</ss> in the bottom right hand corner. [add image]
- Give it a name. e.g. Agilent 6545
- (Don't click Save Current yet).
- See the tabs along the top: we will change some of these settings.

- Note: even if you use an existing setup, some values will need to be re-set - e.g. Statistical test, View pairs (I set up as 1 matches 1 etc?); also biosource will need to be re-set.

![parameters](images/params1.png)

**General**

- Retention time format: seconds
- Polarity: negative

**Feature Detection**

Method: centWave

ppm: 50
minimum peak width: 10
maximum peak width: 50

**Retention Time Correction**

This is to correct shift as run progresses? yes.
averages the RT for each feature. incl. min and max RTs

Method: obiwarp
profStep: 0.1




**Alignment**

same as ret time correction

mzwid: 0.5
minfrac: 0.5
bw: 20


**Statistics**

Statistical test: Paired non-parametric (Wilcoxon signed rank test)
[this didn't keep previously set stats - are these ok? didn't change any thresholds etc]

**Annotation**

Search for: isotopes + adducts
m/z absolute error: 0.05
ppm: 50



can have diff adducts (eg NH44)

helps with ID

the adducts depend on the solvent used in LC





**Identification**

ppm: 50
adducts: [M-H]-
sample biosource: K. pneumo [then SELECT button on left]







also: biosource - choose a K pneumo? yes. but would maybe jsut be for downstream pathway stuff.



ppm - same as feature detection ppm eg 7



**Visualization**

which settings
200 ok. jsut width of window for chromatogram.


**Miscellaneous**

Bypass file sanity check: tick

Then:

- <ss>Save Current</ss>
- <ss>Submit Job</ss>

This will now bring up the "View Results" page.

- The current job will be listed as "Processing" with a % completion bar.


## View results

Click on <ss>View</ss>.

First check samples are different enough? PCA and MDS? iPCA?



**Graphs: Total Ion Chromatograms**

All the ions detected. Their intensity vs retention time. Original and corrected. Also, a correction curve graph.



**Results Table**

A table of features (ion with unique m/z and retention time).

sort by heading or click on a feature row.

to filter, click on the magnifying glass at the top left.
e.g. filter by p-value or fold change. (or multiple things)

The headings mean:

- **fold**: fold change (log?). ratio of mean intensities.
- **p-value**: for whatever test was performed?
- no q value (or is this an option to toggle p value?)
- **updown**: up or down regulated (compared to the other sample)
- **mzmed**: median value for m/z
- **rtmed**: median value for retention time
- **maxint**: max intensity of the feature
- **dataset1_mean**: average intensity of this feature in dataset1
- **isotopes**: isotopes found (why not all have at least one? )
- **adducts**: ? other things that formed?
- **peakgroup**: abritrary number for different groups?  where from?
- no Metlin column?  only have as a graph to the right?

Graphs to the right: (only have EIC ?)

- Extracted Ion Chromatogram (EIC): intensity vs. retention time
- Mass spectrum: intensivy vs. m/z
- Metlin IDs

(how can you match to Metlin with only MS and not MS/MS data?)


The [METLIN database](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage) contains data on metabolites, their mass, their known and predicted fragment masses.

**Metabolomic Cloud Plot**

Cloud plot: m/z vs retention time
What are the circles - abundance of metabolites?
Can adjust:
m/z
retention time
intensity
-- why can u adjust? is it for zooming in?  


**Interactive Heatmap**



**iPCA**


**Conclusion**

We found the metabolites A, B, C were upregulated in strain KPC2, and they are part of the XYZ pathway involved in <something>.

## Next

More complex analyses: eg compare wildtype with 5 mutants. Metabolites in common to the 5 mutants identified; then investigate biochemical pathway and function.


## Links

[XCMS Online](https://xcmsonline.scripps.edu)

["XCMS Online" documentation](https://xcmsonline.scripps.edu/landing_page.php?pgcontent=documentation)

Smith, R. et al. (2014) Proteomics, lipidomics, metabolomics: a mass spectrometry tutorial from a computer scientist's point of view. *BMC Bioinformatics*. DOI: 10.1186/1471-2105-15-S7-S9.

- See Figure 2 for an excellent explanation of the various graphs produced from MS.
