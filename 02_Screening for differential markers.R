##############################################################
## Differential Expression Analysis Using Fold Change & p-value
## Author: Honghao Zhang
## Date: 2025-03-15
##############################################################

################################################################
## Step1: Data Loading - Reuse Previous Data
################################################################
# Note: Uses same data loading as PCA analysis
# Define input file paths
inputfiles <- "./inputdata.csv"  # Path to expression matrix (features x samples)
groups <- "./sample.csv"        # Path to sample metadata with group info

# Read expression matrix (features as rows, samples as columns)
data <- read.csv(file = inputfiles, row.names = 1, header = TRUE)

# Read experimental group information
groups <- read.csv(file = groups, row.names = 1, header = TRUE)

################################################################
## Step2: Differential Expression Analysis with DESeq2
################################################################
library(DESeq2)

# Create DESeq2 dataset object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = data,        # Raw count data matrix
  colData = groups,        # Sample group information
  design = ~ group         # Experimental design formula
)

# Perform differential expression analysis
dds_res <- DESeq2::DESeq(
  object = dds, 
  sfType = "poscounts"     # Estimation method for size factors
)

# Extract results for specific contrast
res <- DESeq2::results(
  object = dds_res,
  tidy = TRUE,             # Return tibble format
  format = "DataFrame",    # Output as DataFrame
  contrast = c("group", "TreatA", "TreatB")  # Comparison groups
)

# Rename first column for consistency
colnames(res)[1] <- "ID"  # Unique feature identifier

################################################################
## Step3: Initial Feature Filtering
################################################################
# Filter significantly upregulated features (log2FC > 2)
res_up <- res[res$log2FoldChange > 2 & !is.na(res$padj), ]

# Filter significantly downregulated features (log2FC < -2)
res_down <- res[res$log2FoldChange < -2 & !is.na(res$padj), ]

################################################################
## Step4: Comprehensive Feature Annotation
################################################################
# Set significance thresholds
cut_off_stat <- 0.05   # Adjusted p-value cutoff (FDR)
cut_off_logFC <- 1     # Minimum absolute log2 fold change

# Prepare results dataframe
diff_res <- as.data.frame(res)

# Handle missing values
diff_res[is.na(diff_res)] <- 0  # Impute NA values for visualization

# Create differential expression labels
diff_res$label <- ifelse(
  diff_res$padj < cut_off_stat & abs(diff_res$log2FoldChange) >= cut_off_logFC,
  ifelse(diff_res$log2FoldChange > cut_off_logFC, "Up", "Down"),
  "No-significant"
)

# Create feature labels for top significant features
diff_res$name <- ifelse(
  diff_res$pvalue < cut_off_stat & abs(diff_res$log2FoldChange) >= 8,
  as.character(diff_res$ID), 
  ""
)

################################################################
## Step5: Volcano Plot Visualization
################################################################
library(ggplot2)
library(ggrepel)

ggplot(data = diff_res, aes(x = log2FoldChange, y = -log10(pvalue), color = label)) +
  # Core visualization elements
  geom_point(size = 3.5, alpha = 0.4) +  # Semi-transparent points
  # Reference lines for thresholds
  geom_vline(xintercept = c(-cut_off_logFC, cut_off_logFC), 
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(cut_off_stat), 
             linetype = "dashed", color = "black") +
  # Color mapping for expression states
  scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) +
  # Axis labels with scientific notation
  labs(x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("p-value"))) +
  # Theme customization
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank()) +
  # Feature labeling system
  geom_text_repel(
    aes(label = name),
    size = 3,
    box.padding = unit(0.8, "lines"),   # Padding around labels
    point.padding = unit(0.8, "lines"),
    max.overlaps = Inf,                # Show all labels
    show.legend = FALSE
  )
