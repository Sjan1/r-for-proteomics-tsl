Plan for Brussell meeting
=========================

Here is the summary of our progress, my thoughts on what we could achieve in the upcoming meeting, and what Laurent could help with. In general, we made a good  progress and now is the right time to review it and make further plans. The hope is that the meeting will help us to shift the focus on the presentation for the Open Plant we plan for March.
   
In the last couple of months we have worked on PDLP1 dataset. For the simplicity's sake, I suggest to stick with it for now; it should be fairly easy to apply it to other datasets.


#### What we have done is here:  
https://github.com/Sjan1/r-for-proteomics-tsl/tree/clean/rfp4

----- 

### 1--- S06-tab.R
Here is our table that describes PDLP1 experiment: "experiment-sample-factor(s)-run-search" data. We create unique strings to label each sample with, that can be parsed to merge fractions or find out sample replicates to calculate statistics. We still have table to do so. 	
TO DO:
This script should ideally merge with design.Rmd


### 2--- S07-msnset-build.R
This generates MSnSet with peptides on rows for all the samples that are in the table with exception of the peptides that are commented out in the column "not used". Then it combines peptides to proteins. Expression values are spectral counts based on unique peptide sequence.
It seems to me, for what follows after peptides are loaded, QC plots and differentiation on peptide level, it would be better to split it and move combining the proteins to a new script.
TO DO:
More different counts would be useful, esp. emPAI, protein count, unique peptide count. 
We do not use meta data yet, at the bottom of MSnSet! What is it good for?


### 3--- S08-msnset-QC.R
Here we read MSnSet expression values to convert data into long format for QC plots with ggplot on the peptide level. 
TO DO:
This could be further improved - more plots using feature data and raw data. (one of the Open Plant objectives)
We need QC plots on three levels: raw data, peptide, protein. I reckon we keep these levels separately, then combine for some final report at the end.

### 4--- The next step should be proper differential analysis and visualization.
TO DO:
Scatter plot, MA plot or Volcano plot with coloured proteins differentially changing between samples. Possibly heatmap and sample clustering would be nice too. Please mind there is a deliberate catch in PDLP1 dataset that makes differentiation and visualization more difficult – there are four factors describing the experiment: C,P,H,I. The scientific questions in this experiment were: which proteins are enriched in P, while C is a negative control. Is the result different when both C and P samples are infected? CH- healthy control, CI- infected control, PH - healthy PDLP1 phenotype , PI - infected PDLP1 phenotype. Also not the different number of sample in each replicate groups.

### 5--- When we have the differential, 
we need to mine raw data to print the spectra of selected proteins we found, label fragment ions and asses PSMs visually.

----- 

##### Other useful things:
(I think to do the above well is plenty of work, but just in case...)

- Management of the project. How to make scripts we generate easily reusable. How to link them and document the usage. Anything that helps to present the project and preserve it for future. 

- XIC - extracted ion chromatograms of selected features, zoom in, label position of the PSM in several samples.

- Visualize experimental factors with data.tree(s).

- iRTs and retention time alignment.

- Contingency table of spectral counts, peptides, PTMs and samples.

- vi Check the integrity of the individual measurements (LC-MS runs) to find out retention time drift, signal (sensitivity) drop, MS1 or MS2 level problem, distribution changes (charge states, precursors, MS2 fragments, infusion times, RTs, etc.) 
 
  
 ##### more notes in Brussels
- read PD output
- how to MSnSet->ggplot using pData
-feature data for QC, print spectra with scores, precursor mass, precursor error 
 
 -----