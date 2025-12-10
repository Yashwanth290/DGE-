Differential Gene Expression (DGE) Analysis of SARSCoV2 Infection


 Analysis Overview
1. Data Source
•	Host transcriptomic data from SARSCoV2 infected vs control samples
 2. Preprocessing & Quality Control
•	Sample filtering (low read counts, library size)
•	Gene filtering (low expression)
     QC metrics:
•	PCA / clustering
•	Library size distributions
•	Expression density plots
 3. Normalization
•	Method: DESeq2 ( edgeR/limma )
•	Variance stabilizing transformation (VST) for PCA & visualization.
 4. Differential Expression
•	Using DESeq2:
•	Model: ~ condition
•	Contrast: infected vs control
  Outputs:
•	log2 fold change (LFC)
•	pvalues
•	Benjamini–Hochberg adjusted pvalues (FDR)
 5. Visualization
•	Volcano plots for significant genes
•	Heatmaps of top DE genes
•	PCA of samples
•	MA plots
 6. Pathway & Functional Enrichment
•	GO: Biological Processes
•	KEGG: Viral infection, immune response
•	Reactome pathways: Interferon signaling, cytokine responses


Requirements
•	R Packages
•	DESeq2
•	ggplot2
•	pheatmap
•	clusterProfiler
•	EnhancedVolcano
•	tidyverse

License
This project is licensed under the MIT License.
