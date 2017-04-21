<br>
# Metabolomics with XCMS

## Overview

Metabolomics is the study of metabolites: small molecules (smaller than proteins) produced by organisms. One technique to identify and quantify metabolites is mass spectrometry. Metabolites are ionized, which breaks them into charged fragments. Next, the abundance of each fragment is measured. A graph shows abundance peaks at various mass-to-charge ratios. This spectrum of peaks (the mass spectrum) is matched to a database of known spectra, to identify the metabolite from which it was produced.


LC separates the metabolites and then these are sent to MS.

Molecular ions can have +1 or more charged, but we mainly look at the M+1

This ion (what value - its mass?) is entered into Metlin.

"each peak is associated with a mass-to-charge ratio, retention time,
fold change, p-value, and relative intensity "

identify features (ions with unique m/z and retention time)
and then filter eg based on:
- those that differ bn treatment and control with eg certain fold change (eg >3) and certain p value (eg < 0.0001)
- then put these in a cloud plot. cloud plots:
useful bc they display info about the biochemistry of the molecules. eg they show the m/z value and retention time (eg this differs for fat-soluble and water-soluble molecules).



## XCMS

inputs: raw data from mass spec

output: data converted to identified metabolites (via METLIN db?)



## Get data

Data from Agilent (is this correct) is in .D format.

Upload these to xcms, or convert these to .mzdata or .mzXML files first. (e.g. with Proteowizard (http://proteowizard.sourceforge.net/downloads.shtml)




data files:
eg sample from ABR project (if not under embargo)
eg public data (all v big though? make smaller?)


ABR project data:
acquired on a HPLC/6545 Q-TOF




## Steps

XCMS online.
- Work out on Galaxy

Compare metabolites in two groups (e.g. 2 different bacterial strains).

This is a pairwise comparison job.

Select data files.  - What format are ABR files.
    -  Glu_1_1.mzData
    -  
Select Instrument. - the instrument it was run on or closest available.

Select machine parameters -

Stats:
- are 2 samples different? (is that right)

Identification:
- is this converting pattern to metabolite?


The [METLIN database](https://metlin.scripps.edu/landing_page.php?pgcontent=mainPage) contains data on metabolites, their mass, their known and predicted fragment masses.

Viz:
- what is EIC width

## Results
- The main result is a matrix of metabolite values. (amts?)


Things to know:<!-- correct? -->
Retention time: the gaseous, charged fragment will be retained in the detector chamber for longer if it is heavier. Thus, retention time is a proxy for the fragment's mass.  
Peak: In a graph of mass vs ?, the peak is a detection of a fragment. Compare height of peaks between 2 groups. Just of the metabolite, prior to fragmentation?


mass and the fragmentation data (MS/MS spectra) for each metabolite peak. ??



GRAPHS
Pairwise Jobs will have an initial results screen that provides a summary with the following graphs:

Total Ion Chromatograms (TIC) Before retention time alignment

Retention Time Deviation vs Retention Time

LTotal Ion Chromatograms (TIC) post retention time correction


Cloud plot: m/z vs retention time
What are the circles - abundance of metabolites?
Can adjust:
m/z
retention time
intensity
-- why can u adjust? is it for zooming in?  



Multidimensional Scaling (MDS)



Principal Component Analysis (PCA)
- can re-scale? why?


Optionally a MS/MS scan location plot if MS/MS data was included in file upload


TABLE
The View Results Table is different than the View Results page and can be found by clicking the “VIEW” button of a job. On the left side is a button named "View Results Table."




## Next

More complex analyses: eg compare wildtype with 5 mutants. Metabolites in common to the 5 mutants identified; then investigate biochemical pathway and function.


## Links

[XCMS Online](https://xcmsonline.scripps.edu)

["XCMS Online" documentation](https://xcmsonline.scripps.edu/landing_page.php?pgcontent=documentation)
