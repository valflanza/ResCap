# ResCap Repository

This is the ResCap repository for software, raw data tables and data bases.


Repository structure:
* **db**	Folder with the databases used for the analysis (bowtie2 index format).
	- argannot (ARG-Annot: Antibiotic Resistance Database)
	- bacmet (BACMET: Biocide & Metal Resistance Database)
	- relax_trim (ConjDB: Relaxases database)
* **data**	Folder with the raw data tables and other results
	- AlleleNetwork.cys (Allele Network in Cytoscape format)
	- AssembliesData.tsv (Tabular data table with the assemblies summary)
	- BacMetFunctionallData.tsv (Tabular data table with the functional annotation of MGC)
	- BlastnResultsAgainstResCapDB.tsv.gz (Blast results of the ORF of metagenomic assemblies)
	- BlastpResultsAgainstUniProtKB.tsv.gz (Blast results of the ORF of metagenomic assemblies)
	- DiffAnalysisMSS.tsv (Differential analysis from DESeq2 for MSS samples)
	- DiffAnalysisRescap.tsv (Differential analysis from DESeq2 for ResCao samples)
	- FullDataTable.tsv (All samples data table in tidy format)

