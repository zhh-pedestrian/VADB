##############################################################
## Single-cell RNA-seq Data Processing Pipeline
## Author: Honghao Zhang
## Date: 2025-03-15
##############################################################

################################################################
## Step1: Data Loading & Initial Processing
################################################################

# Load raw count matrices from 10X Genomics format
test1 <- Seurat::Read10X("./scRNAseq/test1/")  # Sample1 raw count matrix
test2 <- Seurat::Read10X("./scRNAseq/test2/")  # Sample2 raw count matrix

# Create Seurat objects with basic QC thresholds
data_test1 <- Seurat::CreateSeuratObject(
  counts = test1, 
  project = "test1",
  min.cells = 3,     # Retain genes detected in ≥3 cells
  min.features = 200 # Retain cells with ≥200 genes detected
)

data_test2 <- Seurat::CreateSeuratObject(
  counts = test2, 
  project = "test2",
  min.cells = 3,
  min.features = 200
)

# Add sample metadata
data_test1@meta.data$Sample <- "test1"  # Add sample identifier
data_test2@meta.data$Sample <- "test2"

################################################################
## Step2: Dataset Integration & QC Metrics Calculation
################################################################

# Create merged Seurat object
data_list <- list(test1 = data_test1, test2 = data_test2)
df <- Seurat::merge(
  x = data_test1,
  y = list(data_test2),
  add.cell.ids = names(data_list)  # Preserve sample origin in cell IDs
)

# Calculate quality control metrics
mTOR_genes <- c("MTOR", "RAPTOR", "RICTOR", "ULK1")  # mTOR pathway markers
df[["percent.mTOR"]] <- Seurat::PercentageFeatureSet(df, features = mTOR_genes)

# Mouse mitochondrial genes (use "^MT-" for human data)
df[["percent.mt"]] <- Seurat::PercentageFeatureSet(df, pattern = "^mt-")  

# Hemoglobin genes (erythrocyte marker)
df[["percent.hb"]] <- Seurat::PercentageFeatureSet(df, pattern = "^Hbb-")  

################################################################
## Step3: Quality Control Visualization
################################################################

# Violin plots for key metrics
Seurat::VlnPlot(
  df, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mTOR"),
  ncol = 3, 
  pt.size = 0  # Remove individual data points
)

# Scatter plot for feature relationships
Seurat::FeatureScatter(
  df, 
  feature1 = "nCount_RNA",  # Total UMI count per cell
  feature2 = "nFeature_RNA", # Number of detected genes
  pt.size = 0.5
)

################################################################
## Step4: Cell Filtering Strategy
################################################################

# Calculate statistical thresholds
feature_stats <- quantile(df$nFeature_RNA, c(0.05, 0.95))  # 5-95% percentile
count_stats <- quantile(df$nCount_RNA, c(0.05, 0.95))

# Apply multi-criteria filtering
df_filter <- subset(df, 
                    subset = 
                      nFeature_RNA > feature_stats[1] &    # Remove low complexity cells
                      nFeature_RNA < feature_stats[2] &    # Remove potential doublets
                      nCount_RNA > count_stats[1] &        # Remove low-quality cells
                      nCount_RNA < count_stats[2] &        # Remove potential multiplets
                      percent.mt < 20 &                    # Control mitochondrial content
                      percent.hb < 5                       # Control erythrocyte contamination
)

# Verify filtering effect
Seurat::VlnPlot(
  df_filter, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, 
  pt.size = 0
)

# Display dimension changes
message("Original dimensions: ", paste(dim(df), collapse = " x "))
message("Filtered dimensions: ", paste(dim(df_filter), collapse = " x "))