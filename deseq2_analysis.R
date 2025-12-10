#!/usr/bin/env Rscript
###############################################################################
# DESeq2 Differential Gene Expression Analysis
# Modified for GSE164073_raw_counts_GRCh38.p13_NCBI.tsv
# Metadata: metadata.csv
# Working directory:  setwd("//wsl.localhost/Ubuntu-24.04/home/yashu/24MSBI116/question2")
###############################################################################

# =============================================================================
# 1. LOAD REQUIRED PACKAGES
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

required_pkgs <- c("DESeq2", "apeglm", "ggplot2", "dplyr", "data.table",
                   "pheatmap", "RColorBrewer", "tidyverse")

for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE)
    }
    library(pkg, character.only = TRUE)
}

# =============================================================================
# 2. DEFINE FILE PATHS AND OUTPUTS
# =============================================================================

setwd("//wsl.localhost/Ubuntu-24.04/home/yashu/24MSBI116/question2")

counts_file <- "GSE164073_raw_counts_GRCh38.p13_NCBI.tsv"
metadata_file <- "metadata.csv"
output_dir <- "output/deseq2_results"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# 3. READ AND CLEAN COUNT MATRIX
# =============================================================================

cat("Reading count matrix...\n")
counts_raw <- fread(counts_file, header = TRUE, sep = "\t")

# Display first few columns to understand structure
cat("Original column names:\n")
print(colnames(counts_raw))

# Assuming first column is gene IDs (adjust if needed)
gene_ids <- counts_raw[[1]]
counts_data <- counts_raw[, -1]
rownames(counts_data) <- gene_ids

# Convert to numeric and handle NAs
counts_data <- as.data.frame(sapply(counts_data, function(x) as.integer(as.character(x))))
rownames(counts_data) <- gene_ids

# Replace NAs with zeros
counts_data[is.na(counts_data)] <- 0

# Remove genes with all zeros
counts_data <- counts_data[rowSums(counts_data) > 0, ]

# Verify cleanup
cat("NA values after cleanup:", sum(is.na(counts_data)), "\n")
cat("Genes remaining after removing all-zero rows:", nrow(counts_data), "\n")
cat("Count matrix dimensions:", dim(counts_data), "\n")
cat("Samples detected:\n")
print(colnames(counts_data))

# =============================================================================
# 4. READ AND ALIGN METADATA
# =============================================================================

cat("\nReading metadata...\n")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Display metadata structure
cat("Metadata columns:\n")
print(colnames(metadata))
cat("Metadata preview:\n")
print(head(metadata))

# Adjust column names based on your metadata structure
# Assuming columns are: SampleID, Group, Replicate
colnames(metadata) <- c("SampleID", "Group", "Replicate")

# Keep only samples present in count matrix
metadata <- metadata[metadata$SampleID %in% colnames(counts_data), ]

# Match order of samples between metadata and count matrix
metadata <- metadata[match(colnames(counts_data), metadata$SampleID), ]

if (!all(metadata$SampleID == colnames(counts_data))) {
    stop("ERROR: Sample IDs in metadata do not match count matrix columns!")
}

rownames(metadata) <- metadata$SampleID
metadata$Group <- factor(metadata$Group, levels = unique(metadata$Group))

cat("Metadata loaded successfully:\n")
print(metadata)

# =============================================================================
# 5. CREATE DESEQ2 DATASET
# =============================================================================

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                               colData = metadata,
                               design = ~ Group)

# Pre-filter genes with low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("\nDESeq2 dataset ready.\n")
cat("Genes:", nrow(dds), "  Samples:", ncol(dds), "\n")

# =============================================================================
# 6. VARIANCE STABILIZATION AND QC
# =============================================================================

cat("\nPerforming variance stabilization and PCA...\n")
vsd <- vst(dds, blind = FALSE)

# PCA Plot
pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = -1.2, size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_minimal(base_size = 14) +
    ggtitle("PCA Plot - Sample Clustering by Group")

yes

# Sample distance heatmap
cat("Generating sample distance heatmap...\n")
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(vsd$Group, rownames(sample_dist_matrix), sep = "_")

pdf(file.path(output_dir, "sample_distance_heatmap.pdf"), width = 9, height = 7)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distance")
dev.off()

# =============================================================================
# 7. RUN DESEQ2 ANALYSIS
# =============================================================================

cat("\nRunning DESeq2 differential expression analysis...\n")
dds <- DESeq(dds)

# Get group levels for contrast
groups <- levels(metadata$Group)
cat("Comparing:", groups[2], "vs", groups[1], "\n")

res <- results(dds, contrast = c("Group", groups[2], groups[1]))
res_lfc <- lfcShrink(dds, coef = 2, type = "apeglm")

summary(res)

# Count significant genes
cat("\nSignificant genes summary:\n")
cat("padj < 0.05:", sum(res$padj < 0.09, na.rm = TRUE), "\n")
cat("padj < 0.05 & |log2FC| > 1:", sum(res$padj < 0.09 & abs(res$log2FoldChange) > 1, na.rm = TRUE), "\n")

# =============================================================================
# 8. VISUALIZATIONS
# =============================================================================

# MA plots
pdf(file.path(output_dir, "MA_plot.pdf"))
plotMA(res, ylim = c(-5, 5), main = "MA Plot (Unshrunken LFC)")
dev.off()

pdf(file.path(output_dir, "MA_plot_shrunk.pdf"))
plotMA(res_lfc, ylim = c(-5, 5), main = "MA Plot (Shrunken LFC)")
dev.off()

# Volcano plot
cat("Generating volcano plot...\n")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.09 & res_df$log2FoldChange > 1] <- "Up-regulated"
res_df$Significance[res_df$padj < 0.09 & res_df$log2FoldChange < -1] <- "Down-regulated"

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal(base_size = 14) +
    ggtitle(paste("Volcano Plot:", groups[2], "vs", groups[1])) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")

ggsave(file.path(output_dir, "volcano_plot.pdf"), volcano, width = 9, height = 7)
ggsave(file.path(output_dir, "volcano_plot.png"), volcano, width = 9, height = 7, dpi = 300)

# Dispersion plot
pdf(file.path(output_dir, "dispersion_plot.pdf"), width = 8, height = 6)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("\nSaving DESeq2 results...\n")

res_ordered <- res[order(res$padj), ]

# Write all results
write.table(as.data.frame(res_ordered),
            file = file.path(output_dir, paste0(Sys.Date(), "_all_genes_DESeq2_results.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Significant subsets
res_sig <- subset(res_ordered, padj < 0.09)
res_sig_fc <- subset(res_ordered, padj < 0.09 & abs(log2FoldChange) > 1)

write.table(as.data.frame(res_sig),
            file = file.path(output_dir, paste0(Sys.Date(), "_significant_genes_padj0.05.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)

write.table(as.data.frame(res_sig_fc),
            file = file.path(output_dir, paste0(Sys.Date(), "_significant_genes_padj0.05_FC2.txt")),
            sep = "\t", quote = FALSE, row.names = TRUE)

# Save up/down lists
up_genes <- rownames(subset(res_ordered, padj < 0.09 & log2FoldChange > 1))
down_genes <- rownames(subset(res_ordered, padj < 0.09 & log2FoldChange < -1))

write.table(up_genes, file = file.path(output_dir, paste0(Sys.Date(), "_upregulated_genes.txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(down_genes, file = file.path(output_dir, paste0(Sys.Date(), "_downregulated_genes.txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Up-regulated:", length(up_genes), " Down-regulated:", length(down_genes), "\n")

# =============================================================================
# 10. HEATMAP OF TOP 50 DEGs
# =============================================================================

cat("\nGenerating heatmap for top 50 DEGs...\n")
top_genes <- head(rownames(res_ordered), 50)
norm_counts <- counts(dds, normalized = TRUE)
top_counts <- log2(norm_counts[top_genes, ] + 1)

annotation_col <- data.frame(Group = dds$Group)
rownames(annotation_col) <- colnames(dds)

pdf(file.path(output_dir, "top50_DEGs_heatmap.pdf"), width = 10, height = 12)
pheatmap(top_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_col = annotation_col,
         main = "Top 50 Differentially Expressed Genes",
         fontsize_row = 6)
dev.off()

# =============================================================================
# 11. SAVE R OBJECTS AND SESSION INFO
# =============================================================================

saveRDS(dds, file = file.path(output_dir, "dds_object.rds"))
saveRDS(res, file = file.path(output_dir, "results_object.rds"))
saveRDS(res_lfc, file = file.path(output_dir, "results_lfc_shrunk.rds"))

cat("\nAnalysis complete! Results saved in:", output_dir, "\n")

cat("\n=== Session Information ===\n")
sessionInfo()
