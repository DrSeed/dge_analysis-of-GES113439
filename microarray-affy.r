###Required packages 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") ##only required once

#General Bioconductor packages
BiocManager::install("Biobase") #core package that contains multiple genomic data structure classes
BiocManager::install("oligoClasses") #package contains class definitions, validity checks, and initialization methods for classes used by the oligo and crlmm packages.

#Annotation and data import packages
BiocManager::install("ArrayExpress") #package to import data from ArrayExpression 
BiocManager::install("pd.hugene.1.0.st.v1") #human gene annotation database package
BiocManager::install("pd.hugene.1.0.st.v2") 
BiocManager::install("hugene10sttranscriptcluster.db") #human gene annotation database package with probe information

#Quality control and pre-processing packages
BiocManager::install("oligo") #helps in generating expressionset
BiocManager::install("arrayQualityMetrics") #package that helps in automated quality control analysis of expression dataset.

#Analysis and statistics packages
BiocManager::install("limma") #package with multiple linear models to apply on the expression data to find differentially expressed genes
BiocManager::install("topGO") #package for topGO gene ontology annotation 
BiocManager::install("ReactomePA") #package for reactome pathway annotation 
BiocManager::install("clusterProfiler")

#Plotting and color options packages
BiocManager::install("gplots") #package for data visualization
BiocManager::install("ggplot2")  #package for data visualization
BiocManager::install("geneplotter") #package for microarray data visualization
BiocManager::install("RColorBrewer") #package for color brewer, color generation 
BiocManager::install("pheatmap") #package to generate heat map

#Formatting/documentation packages
#BiocManager::install("rmarkdown")
#BiocManager::install("BiocStyle")
BiocManager::install("dplyr") #package for data manipulation and preprocessing 
BiocManager::install("tidyr") #package for data manipulation and preprocessing 

#Helpers:
BiocManager::install("stringr") #package with functions to help in working with string datasets
BiocManager::install("matrixStats") #package with functions that Apply to Rows and Columns of Matrices ("and to Vectors")
BiocManager::install("genefilter") #package with methods for filtering genes from high-throughput experiments
BiocManager::install("openxlsx") #package for reading/writing of excel files
#BiocManager::install("devtools")


###########################
# Calling Libraries
###########################

#General Bioconductor packages
library(Biobase) #core package that contains multiple genomic data structure classes
library(oligoClasses) #package contains class definitions, validity checks, and initialization methods for classes used by the oligo and crlmm packages.

#Annotation and data import packages
library(ArrayExpress) #package to import data from ArrayExpression 
library(pd.hugene.1.0.st.v1) #human gene annotation database package
library(pd.huex.1.0.st.v2)
library(hugene10sttranscriptcluster.db) #human gene annotation database package with probe information

#Quality control and pre-processing packages
library(oligo) #helps in generating expressionset
library(arrayQualityMetrics) #package that helps in automated quality control analysis of expression dataset.

#Analysis and statistics packages
library(limma) #package with multiple linear models to apply on the expression data to find differentially expressed genes
library(topGO) #package for topGO gene ontology annotation 
library(ReactomePA) #package for reactome pathway annotation 
library(clusterProfiler)

#Plotting and color options packages
library(gplots) #package for data visualization
library(ggplot2)  #package for data visualization
library(geneplotter) #package for microarray data visualization
library(RColorBrewer) #package for color brewer, color generation 
library(pheatmap) #package to generate heat map

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr) #package for data manipulation and preprocessing 
library(tidyr) #package for data manipulation and preprocessing 

#Helpers:
library(stringr) #package with functions to help in working with string datasets
library(matrixStats) #package with functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(genefilter) #package with methods for filtering genes from high-throughput experiments
library(openxlsx) #package for reading/writing of excel files
#library(devtools)


#Actual Folder for Analysis
raw_data_dir <- "GSE113439" #defining directory 

#Directory Creation
if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir) #if directory doesn't exist, then it'll be created
}

# Fetching ArrayExpression/Dataset
#downloaded_AE <- getAE("E-GEOD-38974", path = raw_data_dir, type = "raw")

##Reading SDRF File
sdrf <- file.path(raw_data_dir, "samples.txt")
sdrf <- read.csv(sdrf)
#Row names of SDRF from Array.Data.File column of SDRF
rownames(sdrf) <- sdrf$File
#AnnotatedDataFrame is being formed for SDRF 
sdrf <- AnnotatedDataFrame(sdrf)

#Reading the raw sample files according to the phenotypic data (SDRF)
raw_data <- read.celfiles(filenames = file.path(raw_data_dir, 
                                                       sdrf$File),
                                 verbose = FALSE, phenoData = sdrf)
stopifnot(validObject(raw_data))

#Retaining only specific columns of pData (Samples)
pData(raw_data) <- pData(raw_data)[, c("File",
                                       "Phenotype")]
#Logarithmic transformation of the expression data
expressionData <- log2(exprs(raw_data))

###Basic Quality Control
#Transposition and PCA construction of the log2 expression data.
PCA <- prcomp(t(expressionData), scale. = FALSE)

#Finding the percentage/ratio of the first 2 PCs.
percentage <- round(100*PCA$sdev^2/sum(PCA$sdev^2), 1)
standard_deviation_ratio <- sqrt(percentage[2]/percentage[1])

PCAframe <- data.frame(PC1 = PCA$x[,1], PCA$x[,2],
                       Phenotype = pData(raw_data)$Phenotype)
colnames(PCAframe)[2] <- "PC2"
PCAframe$Phenotype[1:15] <- "PAH"
PCAframe$Phenotype[16:26] <- "Control"

#PCA Plot Visualization
pca_plot <- ggplot(data = PCAframe) + geom_point(mapping = aes(x= PC1, y = PC2, color = Phenotype))
pca_plot + labs(title = "PCA Plot of Log2 Transformed Expression Data") + xlab (paste0("PC1, Exp: ", percentage[1], "%")) + ylab (paste0("PC2, Exp: ", percentage[2], "%")) + scale_shape_manual(values = c(4,15)) + scale_color_manual(values = c("red", "darkgreen"))

ggsave(
  "pca_raw.png",
  plot = last_plot(),
  dpi = 900,
)

#Box plot construction of Log2 Intensities of Raw Data
oligo::boxplot(raw_data, target = "core", main = "Log2-Intensities of Raw Data")

#Generating automated quality controls
#arrayQualityMetrics(expressionset = raw_data, outdir = getwd(), force = TRUE, do.logtransform = TRUE, intgroup = c("Phenotype", "Factor.Value.phenotype."))


#Summarization and Background Correction Without Normalziation
rma_Data <-rma(raw_data, normalize = FALSE)
row_medians <- Biobase::rowMedians(as.matrix(exprs(rma_Data)))
rle <- sweep(Biobase::exprs(rma_Data), 1, row_medians)
rle <- as.data.frame(rle)
rle_gathered <- gather(rle, sample_name, log2_expression_deviation)
rle_plot <- ggplot(rle_gathered, aes(sample_name,
                                          log2_expression_deviation))
rle_plot + geom_boxplot(outlier.shape = NA) + ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = "aquamarine4",
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))


##Summarization and Background Correction With Normalziation
rma_Data_norm <-rma(raw_data, normalize = TRUE)
rma_norm_expression <- Biobase::exprs(rma_Data_norm)
PCA_rma_norm <- prcomp(t(rma_norm_expression), scale. = FALSE)

oligo::boxplot(rma_Data_norm, target = "core", main = "Log2-Intensities of Normalized Data")

percentage_rma_norm <- round(100*PCA_rma_norm$sdev^2/sum(PCA_rma_norm$sdev^2), 1)
standard_deviation_ratio_rma_norm <- sqrt(percentage[2]/percentage[1])

PCAframe_rma_norm <- data.frame(PC1 = PCA_rma_norm$x[,1], PC2 = PCA_rma_norm$x[,2],
                                Phenotype = pData(rma_Data_norm)$Phenotype
                                )


pca_plot_norm <- ggplot(data = PCAframe_rma_norm) + geom_point(mapping = aes(x= PC1, y = PC2, color = Phenotype))
pca_plot_norm + labs(title = "PCA Plot of Summarized, Calibrated and Normalized Data") + xlab (paste0("PC1, Exp: ", percentage_rma_norm[1], "%")) + ylab (paste0("PC2, Exp: ", percentage_rma_norm[2], "%")) + coord_fixed(ratio = standard_deviation_ratio_rma_norm) + scale_shape_manual(values = c(4,15)) + scale_color_manual(values = c("red", "darkgreen"))

#Data Annotation for Heatmap Visualization
phenotype_names <- ifelse(str_detect(pData(rma_Data_norm)$Phenotype,"PAH"), "PAH", "Control")

data_annotation <- data.frame(Phenotype = phenotype_names)
row.names(data_annotation) <- row.names(pData(rma_Data_norm))

distances <- as.matrix(dist(t(rma_norm_expression), method = "manhattan"))
rownames(distances) <- row.names(pData(rma_Data_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal("YlOrRd"))(255))
colnames(distances) <- NULL
diag(distances) <- NA #don't need any color in diagonal, because it's self-self
ann_colors <- list(Phenotype = c(Control = "green", PAH = "red"))

map <- pheatmap(distances, annotation_row = data_annotation, annotation_colors = NA,
legend = TRUE,
treeheight_row = 0,
legend_breaks = c(min(distances, na.rm = TRUE),
max(distances, na.rm = TRUE)),
legend_labels = (c("Small distance", "Large distance")),
main = "Heatmap for Calibrated Samples (Normalized)")


#Intesnity-based Filtering of Low-intensity Transcripts
threshold_set <- 4
transcript_medians <- rowMedians(rma_norm_expression)
histo_transcript_medians <- hist(transcript_medians, 100, col = "lightblue", freq = FALSE,
                 main = "Histogram of the Median Intensities Per Gene With Threshold",
                 border = "black",
                 xlab = "Median Intensities")
abline(v = threshold_set, col = "red", lwd = 2)

#Filtering out the Genes that are Above Threshold
no_of_samples <- table(pData(rma_Data_norm)$Phenotype)
samples_cutoff <- min(no_of_samples)
filtered_genes_data <- apply(Biobase::exprs(rma_Data_norm), 1,
                             function(x){
                               sum(x > threshold_set) >= samples_cutoff})
filtered_genes_table <- table(filtered_genes_data)
filtered_genes <- subset(rma_Data_norm, filtered_genes_data)

#Annotating the data with identifiers and names
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)

annotated_data <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                        keys = (featureNames(filtered_genes)),
                                        columns = c("SYMBOL", "GENENAME"),
                                        keytype = "PROBEID")

annotated_data <- subset(annotated_data, !is.na(SYMBOL))

#Removal of the PROBEIDS that match to multiple genes. 
grouped_genes <- dplyr::group_by(annotated_data, PROBEID)
grouped_genes_summarized <- dplyr::summarise(grouped_genes, total = n_distinct(SYMBOL))
filtered_annotated_genes <- filter(grouped_genes_summarized, total > 1)
probe_statistics <- filtered_annotated_genes
dim(probe_statistics)

#Generation of ExpressionSet without the probe IDs that have multiple mappings
ids_to_exclude <- (featureNames(filtered_genes) %in% probe_statistics$PROBEID)
table(ids_to_exclude)
filtered_genes_set_final <- subset(filtered_genes, !ids_to_exclude)
validObject(filtered_genes_set_final)

#Generation of filtered data
fData(filtered_genes_set_final)$PROBEID <- rownames(fData(filtered_genes_set_final))
fData(filtered_genes_set_final) <- left_join(fData(filtered_genes_set_final), annotated_data)
rownames(fData(filtered_genes_set_final)) <- fData(filtered_genes_set_final)$PROBEID

write.csv(filtered_genes_set_final,"filtered_genes_set_final.csv")

###Linear Model Data Preparation 

individual <- as.character(row.names(pData(filtered_genes_set_final)))

disease <- ifelse(str_detect(Biobase::pData(filtered_genes_set_final)$Phenotype, "PAH"), "PAH", "N")

##Factors Preparation
#Disease 
design_palmieri <- model.matrix(~ 0 + disease)
colnames(design_palmieri)[1:2] <- c("N","PAH")
rownames(design_palmieri) <- individual

design <- select(as.data.frame(design_palmieri), PAH, N)
design_sorted <- arrange(design, desc(PAH))

write.csv(design_sorted,"design_matrix.csv", row.names = TRUE)

design_palmieri <- as.matrix(design_sorted)

##Application of limma on all genes 

contrast_matrix <- makeContrasts(PAH-N, levels = design_palmieri)
palmieri_fit <- eBayes(contrasts.fit(lmFit(filtered_genes_set_final, design_palmieri), contrast_matrix))


#Extraction of Data from the model
table_CD <- topTable(palmieri_fit, number = Inf)
head(table_CD)
#Verification visualization
hist(table_CD$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "DEGs Frequency Per P-Value Bin PC vs. NPC", xlab = "adjusted p-values")

write.csv(table_CD, "results_without_threshold.csv")


#Threshold for genes that are significant
degs_CD <- subset(table_CD, P.Value < 0.05)
length(degs_CD$SYMBOL)

##DEGs by setting cutoff
DEG_0_5 <- subset(degs_CD, logFC < -0.5 | logFC > 0.5)
DEG_1 <- subset(degs_CD, logFC < -1 | logFC > 1)

write.csv(DEG_0_5, 'DEGs_with_0.5.csv')
write.csv(DEG_1, 'DEGs_with_1.csv')

DEG_overexpressed_0_5 <- subset(DEG_0_5, logFC > 0.5)
DEG_overexpressed_1 <- subset(DEG_1, logFC > 1)

write.csv(DEG_overexpressed_0_5, 'Overexpressed_0.5.csv')
write.csv(DEG_overexpressed_1, 'Overexpressed_1.csv')

DEG_underexpressed_0_5 <- subset(DEG_0_5, logFC < -0.5)
DEG_underexpressed_1 <- subset(DEG_1, logFC < -1)

write.csv(DEG_underexpressed_0_5, 'underexpressed_0.5.csv')
write.csv(DEG_underexpressed_1, 'underexpressed_1.csv')

#Visualization of the DEGs with FC >= 1 and only 100 lowest P-value genes
library(EnhancedVolcano)

table_CD <- subset(table_CD, !is.na(table_CD$SYMBOL))
EnhancedVolcano(table_CD ,
                
                lab = table_CD$SYMBOL,
                
                x = "logFC",
                
                y = "P.Value",
                ylim = c(0, -log10(10e-12)),
                
                pCutoff = 0.05,
                
                FCcutoff = 1.0,
                
                title = "PAH vs. Control")


library(ggplot2)

analyze_enrichment <- function(data_file, dbs) {
  data <- read.csv(data_file, header = TRUE)
  #all <- read.csv("final_", header = TRUE)
  
  up.idx <- which(data$P.Value < 0.05 & data$logFC > 1)
  dn.idx <- which(data$P.Value < 0.05 & data$logFC < -1)
  
  all.genes <- data$SYMBOL
  up.genes <- data[up.idx,]$SYMBOL
  dn.genes <- data[dn.idx,]$SYMBOL
  
  library(enrichR)
  
  for (db in dbs) {
    
    enriched_pw <- enrichr(genes = up.genes, databases = db)
    
    plotEnrich(enriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", 
               title = paste0(db, " [Up]"))
    ggsave(paste0(db, "_up.png"), dpi = 300, width = 8, height = 6)
    
    write.csv(enriched_pw[[1]], paste0(db, "_up.csv"))
    
    enriched_pw <- enrichr(genes = dn.genes, databases = db)
    
    plotEnrich(enriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", 
               title = paste0(db, " [Down]"))
    ggsave(paste0(db, "_dn.png"), dpi = 300, width = 8, height = 6)
    
    write.csv(enriched_pw[[1]], paste0(db, "_dn.csv"))
  }
}


analyze_enrichment("results_without_threshold.csv", c("KEGG_2019_Human", "GO_Biological_Process_2021", 
                                               "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"))
