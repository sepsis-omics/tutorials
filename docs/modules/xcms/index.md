<br>
# Metabolomics with XCMS


## Overview

Metabolomics is the study of metabolites: small molecules (smaller than proteins) produced by organisms. One technique for detecting metabolites is to use mass spectrometry. The metabolites are ionized, which breaks them into charged fragments. <!-- not always fragmented? eg only in MS/MS? --> These fragments are detected in certain patterns, which are compared with a database of patterns from known metabolites. The metabolites are then identified and quantified.

## XCMS

inputs: raw data from mass spec

output: data converted to identified metabolites (via METLIN db?)



## Get data

data files:
eg sample from ABR project (if not under embargo)
eg public data (all v big though? make smaller?)

## Steps

XCMS online.
- Work out on Galaxy

Compare metabolites in two groups (e.g. 2 different bacterial strains).

This is a pairwise comparison job.

Select data files.  - What format are ABR files.
    -  Glu_1_1.mzData
    -  
Select Instrument. - the instrument it was run on?
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
