########################################################
## PCA Analysis Script for Experimental Data Validation
## Author: Honghao zhang
## Date: 2025-03-15
########################################################

################################################################
## Step1: Data Loading
################################################################
# Define input file paths
inputfiles <- "./inputdata.csv"  # Path to main data matrix (features x samples)
groups <- "./sample.csv"        # Path to sample metadata file

# Read data matrix (samples in columns, features in rows)
data <- read.csv(file = inputfiles, row.names = 1, header = T)  # Use first column as row names

# Read sample group information
groups <- read.csv(file = groups, row.names = 1, header = T)    # Use first column as row names

################################################################
## Step2: Data Preprocessing
################################################################
# Filter low-abundance features (rows with total counts > 10)
data_clean <- data[rowSums(data) > 10, ]

# Apply log2 transformation (add 1 to avoid log(0))
data_log <- log2(data_clean + 1)  # Pseudocount addition for zero values

# Z-score normalization (scale rows to mean=0, variance=1)
data_scale <- scale(data_log)     # Standardization for PCA

################################################################
## Step3: Principal Component Analysis
################################################################
# Perform PCA on transposed data (samples in rows, features in columns)
data_pca <- prcomp(t(data_scale),   # Transpose matrix for PCA (samples x features)
                   center = TRUE,   # Center already scaled data (redundant but safe)
                   scale. = FALSE)  # No scaling needed after standardization

# Create PCA results matrix with grouping information
data_pca_mat <- data.frame(
  pc1 = data_pca$x[, 1],   # First principal component scores
  pc2 = data_pca$x[, 2],   # Second principal component scores
  group = groups           # Experimental groups from metadata
)

################################################################
## Step4: Variance Explanation Calculation
################################################################
# Calculate variance explained by each principal component
var_percent <- data_pca$sdev^2 / sum(data_pca$sdev^2) * 100  # sdev = SD of principal components

# Extract variance explained for visualization labels
pc1_var <- round(var_percent[1], 2)  # Variance explained by PC1
pc2_var <- round(var_percent[2], 2)  # Variance explained by PC2

# Display variance explanation results
cat("PC1 Variance Contribution:", pc1_var, "%\n")
cat("PC2 Variance Contribution:", pc2_var, "%\n")

################################################################
## Step5: Visualization with ggplot2
################################################################
library(ggplot2)
library(ggforce)
library(ggrepel)

# Create PCA plot with enhanced visualization elements
ggplot(data_pca_mat, aes(x = pc1, y = pc2, color = group)) +
  # Add confidence ellipses for group clusters
  ggforce::geom_mark_ellipse(
    aes(fill = group),
    expand = unit(4, "mm"),   # Ellipse expansion size
    alpha = 0.15,             # Ellipse transparency
    n = 100                   # Smoothness of ellipse
  ) +
  # Add sample points
  geom_point(size = 5) +      # Point size for samples
  # Add non-overlapping sample labels
  geom_text_repel(
    aes(label = rownames(data_pca_mat)),
    size = 4,                 # Label text size
    max.overlaps = Inf        # Allow unlimited label overlaps
  ) +
  # Add reference lines
  geom_hline(
    yintercept = 0,           # Horizontal reference line
    linetype = "dotdash",
    color = "grey80"
  ) +
  geom_vline(
    xintercept = 0,           # Vertical reference line
    linetype = "dotdash",
    color = "grey80"
  ) +
  # Set coordinate limits
  scale_y_continuous(limits = c(-25, 25)) +
  # Configure axis labels with variance percentages
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),  # X-axis label with variance
    y = paste0("PC2 (", pc2_var, "%)")   # Y-axis label with variance
  ) +
  # Set plot theme
  theme_bw() +                # Black-and-white theme
  # Custom color scheme
  scale_color_manual(values = c("#8dcdd5", "#e6846d")) +  # Point colors
  scale_fill_manual(values = c("#8dcdd5", "#e6846d")) +   # Ellipse fill colors
  # Theme customization
  theme(
    plot.title = element_text(hjust = 0.5),  # Center plot title
    legend.position = "right"                # Legend positioning
  )