# ResCap Repository

This is the ResCap repository for software, raw data tables and data bases.


Repository structure:
* **db**	Folder with the databases used for the analysis (bowtie2 index format).
	- _argannot_ (ARG-Annot: Antibiotic Resistance Database)
	- _bacmet_ (BACMET: Biocide & Metal Resistance Database)
	- _relax_trim_ (ConjDB: Relaxases database)
* **data**	Folder with the raw data tables and other results
	- _AlleleNetwork.cys_ (Allele Network in Cytoscape format)
	- _AssembliesData.tsv_ (Tabular data table with the assemblies summary)
	- _BacMetFunctionallData.tsv_ (Tabular data table with the functional annotation of MGC)
	- _BlastnResultsAgainstResCapDB.tsv.gz_ (Blast results of the ORF of metagenomic assemblies)
	- _BlastpResultsAgainstUniProtKB.tsv.gz_ (Blast results of the ORF of metagenomic assemblies)
	- _DiffAnalysisMSS.tsv_ (Differential analysis from DESeq2 for MSS samples)
	- _DiffAnalysisRescap.tsv_ (Differential analysis from DESeq2 for ResCao samples)
	- _FullDataTable.tsv_ (All samples data table in tidy format)

