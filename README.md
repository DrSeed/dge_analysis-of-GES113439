Dioffrentail gene expression analyis of GSEE113439 (Microarray Data)
 The project involves data preprocessing, quality control, differential expression analysis, and pathway enrichment analysis to identify potential biomarkers for PAH.

Getting Started
Prerequisites
Before you begin, ensure you have R installed on your system. Additionally, the Bioconductor package manager, BiocManager, is required to install the necessary bioinformatics packages.

Data Preparation
The analysis begins by setting the working directory and creating necessary folders for storing raw data and results. Array data can be imported using the ArrayExpress package, and sample metadata is read and processed.

Quality Control and Preprocessing
Quality control metrics are assessed using arrayQualityMetrics. Preprocessing steps including background correction, normalization, and summarization are performed with the oligo package.

Differential Expression Analysis
limma is used to identify differentially expressed genes between PAH and control samples. The results are further explored to identify significant genes based on adjusted p-values and log fold changes.

Pathway Enrichment Analysis
Enrichment analysis for the identified differentially expressed genes is conducted using topGO, ReactomePA, and clusterProfiler to understand the biological implications of the findings.

Visualization
Various visualization techniques such as PCA plots, boxplots, and heatmaps are utilized to illustrate the data quality, expression patterns, and results of the differential expression analysis.


Output
The scripts generate a variety of outputs, including preprocessed datasets, lists of differentially expressed genes, enriched pathways, and various plots to visualize the data and analysis results. 
