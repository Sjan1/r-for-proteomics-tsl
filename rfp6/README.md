## Project description

### Experimental design:
- biorep: biological replicates
- techrep: technical replicates
- dates: dates experiments where run

### Directories:
- raw: contains converted raw files
- mzml: contains mzML files or it is symlink to the data
- mascot: contains mzid files or it is symlink to the data
- msgf: contains mzid files or it is symlink to the data

### Files:
- README.md
- SampleExperimentTable.csv
- S01-2.R
	Latest version made during Laurent's visit 9th July 2019
- S01-2_psm.R
	To get psm MSnSet
- S01-3.R	
	Tests, checks, etc.