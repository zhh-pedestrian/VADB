##############################################################
## Comprehensive Enrichment Analysis Pipeline (KEGG GO AND GSEA)
## Author: Honghao Zhang
## Date: 2025-03-15
##############################################################

################################################################
## Section 1: Data Preparation
################################################################

# Load differential expression results
inputfile <- "./volcano_data"
volcano_data <- read.csv(file = inputfile, header = TRUE, row.names = 1)

# Extract significant genes (Up/Down regulated)
gene <- data.frame(
  Gene = volcano_data$GeneID[volcano_data$group %in% c("Up", "Down")]
)

# Display available ID types in human annotation database
keytypes(org.Hs.eg.db)

# Convert ENSEMBL IDs to ENTREZ IDs for enrichment analysis
gene_vector <- as.character(gene$Gene)
gene_entrez <- clusterProfiler::bitr(
  geneID = gene_vector,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

################################################################
## Section 2: KEGG Pathway Enrichment Analysis
################################################################

# Perform KEGG enrichment analysis
KEGG <- clusterProfiler::enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",               # Human genome
  pvalueCutoff = 0.05,           # Significance threshold
  qvalueCutoff = 0.05,           # Adjusted p-value threshold
  pAdjustMethod = "BH"           # Benjamini-Hochberg correction
)

# Convert ENTREZ IDs to gene symbols for interpretation
KEGG <- DOSE::setReadable(
  KEGG, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID"
)

# Calculate enrichment metrics
KEGG_frame <- as.data.frame(KEGG) %>%
  dplyr::mutate(
    RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
    FoldEnrichment = (as.numeric(sub("/.*", "", GeneRatio)) * 
                        as.numeric(sub("/.*", "", BgRatio))) /
      (as.numeric(sub(".*/", "", GeneRatio)) * 
         as.numeric(sub(".*/", "", BgRatio)))
  )

# Visualize top pathways
dotplot(
  KEGG, 
  x = "GeneRatio", 
  color = "p.adjust",
  showCategory = 20,
  font.size = 12,
  title = "Top 20 KEGG Pathways",
  label_format = 30
)

################################################################
## Section 3: Gene Ontology (GO) Enrichment Analysis
################################################################

# Perform comprehensive GO analysis
go_enrich <- clusterProfiler::enrichGO(
  gene = gene_entrez$ENTREZID,
  ont = "all",                   # Analyze all ontologies
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Simplify redundant GO terms
go_enrich_simple <- clusterProfiler::simplify(
  go_enrich,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# Calculate enrichment metrics
go_result <- go_enrich_simple@result %>%
  dplyr::mutate(
    RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
    FoldEnrichment = (as.numeric(sub("/.*", "", GeneRatio)) * 
                        as.numeric(sub("/.*", "", BgRatio))) /
      (as.numeric(sub(".*/", "", GeneRatio)) * 
         as.numeric(sub(".*/", "", BgRatio)))
  )

# Basic visualization
dotplot(go_enrich_simple,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        split = 'ONTOLOGY') +
  facet_grid(ONTOLOGY ~ ., scales = "free_y")

################################################################
## Section 4: Advanced GO Visualization
################################################################

# Load and prepare data for advanced visualization
data <- read.csv("./output/1_GO_enrichment.csv")

# Select top terms for each ontology
top_terms <- function(ontology) {
  data %>%
    dplyr::filter(ONTOLOGY == ontology) %>%
    dplyr::arrange(-Count) %>%
    head(10)
}

df <- purrr::map_dfr(c("MF", "CC", "BP"), ~{
  top_terms(.) %>%
    dplyr::mutate(Category = .)
})

# Create multi-layered visualization
ggplot() +
  ggnewscale::new_scale_fill() +
  # Molecular Function layer
  geom_point(
    data = dplyr::filter(df, ONTOLOGY == "MF"),
    aes(x = Count, y = interaction(forcats::fct_reorder(Description, Count), ONTOLOGY), 
        fill = p.adjust, size = Count), 
    shape = 21
  ) +
  scale_fill_gradient(low = "#a1d99b", high = "#238b45", name = "MF p.adjust") +
  # Cellular Component layer
  ggnewscale::new_scale_fill() +
  geom_point(
    data = dplyr::filter(df, ONTOLOGY == "CC"),
    aes(x = Count, y = interaction(forcats::fct_reorder(Description, Count), ONTOLOGY), 
        fill = p.adjust, size = Count), 
    shape = 21
  ) +
  scale_fill_gradient(low = "#fdd49e", high = "#d7301f", name = "CC p.adjust") +
  # Biological Process layer
  ggnewscale::new_scale_fill() +
  geom_point(
    data = dplyr::filter(df, ONTOLOGY == "BP"),
    aes(x = Count, y = interaction(forcats::fct_reorder(Description, Count), ONTOLOGY), 
        fill = p.adjust, size = Count), 
    shape = 21
  ) +
  scale_fill_gradient(low = "#8c96c6", high = "#8c6bb1", name = "BP p.adjust") +
  # Theme and layout configuration
  theme_bw() +
  labs(title = "GO Annotation Landscape", x = "Count", y = "Description") +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

################################################################
## Section 5: Gene Set Enrichment Analysis (GSEA)
################################################################

# Prepare ranked gene list
data_GSEA <- volcano_data %>%
  dplyr::select(GeneID, log2FoldChange) %>%
  dplyr::inner_join(
    clusterProfiler::bitr(
      geneID = volcano_data$GeneID,
      fromType = "ENSEMBL",
      toType = "SYMBOL",
      OrgDb = 'org.Hs.eg.db'
    ),
    by = c("GeneID" = "ENSEMBL")
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# Create ranked list
genelist <- data_GSEA$log2FoldChange
names(genelist) <- data_GSEA$SYMBOL
genelist <- sort(genelist, decreasing = TRUE)

# Perform GSEA analysis
gsea <- clusterProfiler::GSEA(
  genelist,
  TERM2GENE = gene_go2[c('GOALL', 'SYMBOL')],  # Gene set definitions
  TERM2NAME = gene_go2[c('GOALL', 'Term')],    # Gene set names
  pvalueCutoff = 1,            # Include all results for exploration
  pAdjustMethod = 'BH',        # Multiple testing correction
  eps = 0,                     # Boundary for calculating p-values
  seed = TRUE                  # Ensure reproducibility
)
