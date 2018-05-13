
# Notes from 4th session (R-project: RfP3)
------------------------------------------------


## General thoughts
Our OP project is starting to take shape and it looks like that when six-month period finishes, we should be able to have three (see below), well documented  work-flows we are going to test and document.

We have to start asking ourselves the following questions:

- How are we going to finish OpenPlant?
- How to get maximum out of the OP project?
- How are we going to secure what we have been shown?
- How shall we benefit in future? (In the high school language, how to move from "beginner" to "secure" or "confident" level.)

To find the answers, I suggest we start meeting and talking more often here in Norwich. Weekly may not be realistic for our other engagements, but at least every fortnight in between sessions with Laurent should be possible. I hope a short, one to two hour seminar should be sufficient. 

What could we do?

- Have a little workshop
- Try what we have learned
- Discuss and troubleshoot
- Teach each other

This could, for the sake of OP project, help us finish six months period. (If we like it, we can carry on after OP in future.) The main goal now is to find the way to complete our work-flows, document them, preserve for the future use, and finally present what we will have achieved.



## Datasets to work with
We have three datasets to play with (I believe it is enough, but let me know is you think otherwise):

- TASK1: Quality control (QC); HeLa QC runs; time course visualization
- TASK2: Spectral counting (SPC): PDLP1 dataset, Differential abundance analysis
- TASK3: Label based quantitative proteomics (TMT): Gerhard's data; differential abundance analysis



## Problems encountered during April session

### I. 'mzid' problem
We found a problem with Mascot export files in the form of MzIdentML. They are not readable by RfP, function readMzIdData. 
1. Would more recent version of Mascot produce a compatible mzid format? => No, the most recent version behave the same way.
2. Scaffold can also produce mzIdentML, is this format compatible?
3. csv exported search is limited search output, a last resort possibility. 
4. Is there a converter 'dat -> mzid' available?
5. Would Marielle and Govind like to look at Mascot dat parser?
6. We will set up MSGFPlus search engine on nbi cluster and search the data there in the mean time.
	
### II. We need an Experiment-Sample table
We need this table to start processing many files.
File name in the first column, any number of associated factors (including notes) in the other columns; file name independent system.


### III. Differential abundance analysis of SPC by statistical and non statistical approach
1. We will use msmsTest and MSnID packages to calculate emPAi and to differentiate between sample groups. (statistical approach)
2. Sample groups will be provided in "Experiment-Sample Table". Ideally a shiny application could help to dynamically change sample grouping (samples, fractions, and replicates).
3. Specificity and SPC based filter for groups of interactors - implementation my Sum-Count based filtering.  
4. Visualization of hits of interest in e.g. MA plot.
5. Any SP counting differential expression analysis make sense on protein level only.
6. The next step in our analysis is list peptides for every protein, and then PSMs - raw spectra with fragment ion annotation. We rely on manual analysis for selected proteins to avoid overly optimistic PSMs. 
7. PTMs; another reason for the passion for well annotated spectra is our interest in PTMs. It is a hot topic in our lab, with several project depending on confident identification and quantitation of e.g. phosphorylation of certain residue. We need to have a great confidence in PSM with PTM, for then we confirm the structure and quantify using synthetic peptide on QQQ. That is quite laborious and costly, and only probability of PSM is not enough. We check the spectral quality manually; fragment ion should be high above noise, most of them assigned, error distribution of both precursors and the fragments smooth, the neutral losses should make sense considering the sequence assigned, etc. 


### IV. QC pipeline - which plots and statistics we wish to produce.
1. Visualize SPC, UPC, PC (spectal counts, unique peptide counts, protein counts) in as a time course.
2. Add notes from the instrument history file -> Experiment-Sample table.
3. Score distributions for the time points.
4. XIC (extracted ion chromatogram) of selected peptides and corresponding MS2 spectra for every QC file.
5. XIC areas in time: peak area and FWHM (full width half maximum).
6. RT (retention time) stability in time.
7. Correlation of the last run to a short/long-term average (last week, recent best performance, best performance ever).
8. Sequence coverage plot would be extremely nice visualization, preferably with peptide coverage (i) and PTMs (ii), miscleavages (iii), heat map of spectral count intensity plotted over the protein sequence (iv). The final goal is to compare many samples is such a plot. 

### V. TMT pipeline - we should start working on it as soon as a progress is made on TASK1&2. 
	
## Other points discussed in April
1. Plot theoretical fragment masses with search intervals over MS2 spectra. That would require a shiny spectrum zooming functionality.
2. Charge deconvoluted MS2 spectra plots (with annotation) are a convenient tool for manual spectrum validation.
3. We found 1Da error in precursor masses in MS.GF+ result. That can only mean that Asp/Glu deamidation took place or precursor masses were no recognized correctly. Mascot has a function to corrects the latter and deamidation could be added as a variable modification in new search. This would be interesting to compare, I guess when all the other things are done. Do I remember well that in another Laurent's data file this problem was not observed?
4. Compare the results of two search engines. To develop confidence in PSMs when using MS.GF+ we are not used to, the best is to compare with Mascot first. It could give us a few extra PSMs and increase confidence. In the long run, running a second search engine on HPC might be of interest.
5. I should show you Scaffold and what it provides. I mentioned well annotated spectra already, but there is still more: E.g protein similarity plot, quantitative functions that we should also compare with our R-calculated result.


--------------------------------------------------------------------------

## Priorities for the next session (discussed on 8th May, Norwich)
1. Generate MSGFPlus searches for data in all the tasks need to be done before the next meeting.
2. Use sample-experiment input table in QC and PDLP1 Tasks.
3. Use the packages for spectral counting to quantify differences between the samples. (QC - time course, PDLP1 - four categories)
4. Start with TMT quantitation. (Task3)
