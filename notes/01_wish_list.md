# RfP training


## Top priorities

As soon as possible we need to learn how to parse the search results,
combine it with the study factors from experiment documentation and
read the raw data. Then we should focus on these following goals:

### QC workflow; spectral counting (SPC) quantitation

 * **Data:** HeLa standard runs from January 2018
 * Parse (Mascot) search engine output and visualize
 * Start independent work-flow using either Mascot or MSGF+ search
   engines. The latter should be able to achieve higher throughput.
 * Package for spectral counting calculating e.g. emPAI
 * Visualize SPC, UPC, PC in as a time course
 * Add factors from the instrument history file
 * Visualize score distributions for the time points
 * XIC of selected peptides, MS2 of selected peptides
 * Reproducibility of selected peptides
 * RT stability
 * Correlation of the last run to a short/long-term average (last
   week, recent best performance, best performance ever)
 * Chromatographic peak widths (FWHM)
 * **To discuss:** How realistic are the latter two points without PRM
   data? I guess we need to run Skyline to process PRM runs? Could we
   calculate LF precursor area?

### Find differentially changing proteins; SPC quantitation

 * **Data:** Importin dataset -- two sample types (Orbitrap XL)
 * **Data:** PDLP1 -- several sample types (Orbitrap XL)
 * **Data:** Erin/Yasin data for many sample types (Orbitrap Fusion)
 * **Data:** Alternatively Will's vesicles. (Orbitrap XL)
 * Differential changes in SPC in multiple sample types.
 * How to get proteins enriched in certain sample type(s)?
 * Heatmap of the ratios over negative control, missing data replaced
   by very small number. Is there any other way? How to get proteins
   from the heatmap: most abundant, occuring/not occuring in selected
   bio-samples.
 * How to extract patterns when more samples are present? We tested
   `library("geneplotter")`.
 * Is there a way to normalize data from samples prepared by
   immuno-affinity enrichment?
 * **To discuss:** To process the latter two datasets in our current
   Mascot-Scaffold work-flow takes about a week due to long processing
   time and many interruptions. Any improvement that would increase
   throughput welcome.

### PTM identification; SPC quantitation

 * We use error tolerant search to hunt for PTMs. Can we identify
   modified peptides in search output, print the spectra with search
   engine annotation and compare spectral count in sample categories
   to identify how are peptides differentially changing? Currently I
   use Scaffold exported csv then pivot table in Excel then I go to
   Scaffold for the spectra, copy paste to a document.
 * **Data:** HeLa run in error tolerant mode
 * **Data:** Zane's Cell paper
 * **Data:** Possibly Hailong

### More search engines

 * When we have a work-flow automated, could we add another search
   engine(s) to maximize the result?

### HPC computing

 * Do we need it?
 * Do we want to compare with our current work-flows?
 * We have slow throughput with challenging datasets, such as
   Erin/Yasin example or many error tolerant searches.

### How to collect experimental meta information most efficiently?

 * **To discuss:** There is an underlying data management question,
   important in all work-flows that that we are going to generate. We
   discussed with Marielle we could populate a list of input files
   with Mascot search results. This should be linked with the study
   factors that are provided with every experiment. How should be the
   files containing annotations formatted? We should also collect
   search metadata and how many searches we have carried out with
   every sample. The data management can get quite complicated
   easily. We are seeking the lightest set-up.
