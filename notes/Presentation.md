# R for Proteomimcs
## Open Plant project (2018-2019)
## Final presentation

Jan Sklenar, Lauent Gatto, Marielle Vigouroux, Govind Chandra

---

## Outline

- Jan
	+ Introduction of OP project
		- conception
		- goals
 	+ Using example of PDLP1 experiment show the work-flow
	+ Run script(s) or show results ... finish with S-curve
 	+ Explain step by step important bits:
		+ inputs, pre-requisites
		+ exp-table
		+ MSnSet object
		+ filters, decoys, sanity check
		+ rtslprot library
		+ QC plots, SPC, UPC
		+ FC analysis
		+ msmsTests

- Laurent
	+ History and philosophy of RfP
	+ Documentation, vignettes, glossy pictures 
	+ What is actually under the hood, hidden details, what is missing in our example above (do not be afraid to be critical)
	+ MSnSet layout
	+ Can rtslprot be of any use for the RfP in general?
 	+ Future direction of RfP
 	+ What more we can do now, how to expand. (us here in Norwich)
	+ etc. 

- Everybody - natural start of the discussion here
	+ Personal experience
	+ Learning curve
	+ Like/don't like
	+ We learned something (by-products)
		- R, R-packages
		- git
		- symbolic links in Windows
		- installed MSGF+ on our cluster
	+ etc. 

---

## Details
---

## RfP pipeline description

[github.com/Sjan1/r-for-proteomics-tsl/](https://github.com/Sjan1/r-for-proteomics-tsl/rfp6)

**The pipeline:**

* S01_code.R
* S03_code_prot.R
* S05_FC_analysis.R

**Inputs:**

* mzid search files in e.g. *mascot/*
* project meta data e.g. *SampleExperimentTable.csv*

**Description:**

1.  S01 strats working on experimental design, removes unnecessary factors, and defines sample names (from the experimet factors), then it takes all the peptide data in one MSnSet object => e. Functions used here we believe to be reusable are in our library(“rtslprot”)
2.  S03 groups peptides to form proteins => eprot
3.  S05 calculated fold changes for one selected pair of project factors

**Saved results:**
e_mascot_fdr1pc.rds
eprot_mascot_fdr1pc.rds

**Markdown scripts:**

* Design.Rmd
* QC_plots.Rmd
 
The first assists with experimental design and data integrity (names of files and factors).
The second produces QC plots from data processed in our pipeline and provides a bigger picture of the results.

While the scripts run they can generate plots, pdf or html documents, or presentations as required.

---

## Inputs in detail => run design.Rmd
---
## QC plots => run QC_plots.Rmd 
---
## Run the scripts
---

