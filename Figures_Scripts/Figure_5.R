# =============================================================================
# FIGURE 5: CHRNA1+ scRNA Clusters and Myogenesis Signature (GSE115978)
# Panels:
#   A – UMAP of CHRNA1+ malignant cells (clusters 0 and 1) with batch correction
#   B – Violin plots: CHRNA1, CHRNB1, CHRNG expression per cluster
#   C – Hallmark GSVA heatmap per cluster (myogenesis highlighted)
#   D – UMAP of cell cycle phases + bar plot of phase percentages
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "Seurat", "tidyverse", "ggplot2", "grid", "pheatmap",
  "clusterProfiler", "org.Hs.eg.db", "escape", "GSVA"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# 1. LOAD GSE115978 COUNT MATRIX AND ANNOTATION
# =============================================================================

GSE115978 <- read.csv("GSE115978_counts.csv", header = TRUE, sep = ",", row.names = 1)
GSE115978.annotation <- read.csv("GSE115978_annotations.csv", header = TRUE, sep = ",")

# =============================================================================
# 2. SUBSET CHRNA1+ MALIGNANT CELLS AND BUILD SEURAT OBJECT – PANEL 5A
# =============================================================================

# Cells expressing CHRNA1
GSE115978.OnlyA1 <- GSE115978[rownames(GSE115978) == "CHRNA1", ]
GSE115978.OnlyA1 <- GSE115978.OnlyA1[, colSums(GSE115978.OnlyA1) > 0]

GSE115978.OnlyA1.count <- GSE115978[, colnames(GSE115978.OnlyA1)]

# Annotation for CHRNA1+ cells
GSE115978.Annot.A1 <- GSE115978.annotation[
  GSE115978.annotation$cells %in% colnames(GSE115978.OnlyA1.count),
]
GSE115978.Annot.A1 <- GSE115978.Annot.A1[
  GSE115978.Annot.A1$cell.types == "Malignant",
]

# Adjusted patient grouping (batch correction)
GSE115978.Annot.A1$Adjusted.Patient <- ifelse(
  GSE115978.Annot.A1$samples == "Mel110",
  "Mel110",
  "Combined.Samples"
)

# Reassign counts to malignant CHRNA1+ cells only
GSE115978.OnlyA1.count <- GSE115978.OnlyA1.count[
  , colnames(GSE115978.OnlyA1.count) %in% GSE115978.Annot.A1$cells
]

# Create Seurat object
GSE115978.A1 <- CreateSeuratObject(
  counts  = GSE115978.OnlyA1.count,
  project = "GSE115978_CHRNA1"
)

GSE115978.A1$Condition       <- GSE115978.Annot.A1$type
GSE115978.A1$Cell.Type       <- GSE115978.Annot.A1$cell.types
GSE115978.A1$Patient         <- GSE115978.Annot.A1$samples
GSE115978.A1$Adjusted.Patient <- GSE115978.Annot.A1$Adjusted.Patient
GSE115978.A1$DataSet         <- "GSE115978"

# =============================================================================
# 3. INTEGRATION / BATCH CORRECTION AND CLUSTERING – PANEL 5A
# =============================================================================

split.list <- SplitObject(GSE115978.A1, split.by = "Adjusted.Patient")

split.list <- lapply(
  X   = split.list,
  FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x
  }
)

Features <- SelectIntegrationFeatures(object.list = split.list)
Find.Anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = Features)

GSE115978.A1.combined <- IntegrateData(anchorset = Find.Anchors, k.weight = 50)

DefaultAssay(GSE115978.A1.combined) <- "integrated"

GSE115978.A1.combined <- ScaleData(GSE115978.A1.combined, verbose = FALSE)
GSE115978.A1.combined <- RunPCA(GSE115978.A1.combined, npcs = 10, verbose = FALSE)
GSE115978.A1.combined <- RunUMAP(GSE115978.A1.combined, reduction = "pca", dims = 1:10)
GSE115978.A1.combined <- FindNeighbors(GSE115978.A1.combined, reduction = "pca", dims = 1:10)
GSE115978.A1.combined <- FindClusters(GSE115978.A1.combined, resolution = 0.4)

# UMAP of clusters (Panel 5A)
Idents(GSE115978.A1.combined) <- "seurat_clusters"

p5A <- DimPlot(
  GSE115978.A1.combined,
  label.size = 7,
  repel      = TRUE,
  label      = TRUE,
  reduction  = "umap",
  pt.size    = 5
) +
  theme(
    plot.title  = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 12, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
    legend.text = element_text(size = 13, colour = "black", face = "bold")
  )

ggsave("Fig5A_CHRNA1_Clusters_UMAP.png", p5A, width = 8, height = 5.5, dpi = 1200)

# =============================================================================
# 4. VIOLIN PLOTS: CHRNA1, CHRNB1, CHRNG – PANEL 5B
# =============================================================================

DefaultAssay(GSE115978.A1.combined) <- "RNA"
GSE115978.A1.combined <- ScaleData(GSE115978.A1.combined, verbose = FALSE)

p5B <- VlnPlot(
  GSE115978.A1.combined,
  features = c("CHRNA1", "CHRNB1", "CHRNG"),
  pt.size  = 0.4
)

ggsave("Fig5B_CHRNs_Violin.png", p5B, width = 7, height = 5, dpi = 1200)

# Optional EMT TF DotPlot (MITF, SOX10, ZEB1) – not strictly part of panel B
dp_emt <- DotPlot(
  GSE115978.A1.combined,
  cols      = c("blue", "red"),
  dot.scale = 9,
  features  = c("MITF", "SOX10", "ZEB1")
) +
  theme_classic() +
  theme(
    axis.title   = element_text(size = 15, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 30, colour = "black", face = "bold"),
    axis.text.y  = element_text(size = 30, colour = "black", face = "bold.italic"),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold")
  ) +
  ylab("Cluster") + xlab("") + coord_flip()

ggsave("Fig5B_MITF_SOX10_ZEB1_DotPlot.png", dp_emt, width = 8, height = 5, dpi = 600)

# =============================================================================
# 5. CLUSTER MARKERS AND GO ENRICHMENT (USED FOR PANEL 5C CONTEXT)
# =============================================================================

DEG.Markers.Batch.Malig <- FindAllMarkers(
  GSE115978.A1.combined,
  only.pos        = TRUE,
  logfc.threshold = 0.25,
  min.pct         = 0.25,
  p.val.cutoff    = 0.05
)

write.csv(DEG.Markers.Batch.Malig, "Fig5_cluster_markers_CHRNA1_scRNA.csv", row.names = FALSE)

top10.markers <- DEG.Markers.Batch.Malig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  as.data.frame()

p_heat <- DoHeatmap(
  GSE115978.A1.combined,
  features = top10.markers$gene,
  angle    = 0,
  label    = FALSE
) +
  theme(
    axis.text.y   = element_text(size = 18, colour = "black", face = "bold.italic"),
    legend.key.size = unit(0.5, "cm"),
    legend.text  = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15, face = "bold")
  ) +
  scale_color_discrete(name = "Cluster", na.translate = FALSE)

ggsave("Fig5_cluster_markers_heatmap.png", p_heat, width = 8, height = 6, dpi = 1200)

# GO Enrichment (Biological Process) – cluster 0 vs 1
cluster0.genes <- DEG.Markers.Batch.Malig$gene[DEG.Markers.Batch.Malig$cluster == "0"]
cluster1.genes <- DEG.Markers.Batch.Malig$gene[DEG.Markers.Batch.Malig$cluster == "1"]

list.for.GO <- list("0" = cluster0.genes, "1" = cluster1.genes)

CHRNA1.GO <- compareCluster(
  geneCluster = list.for.GO,
  fun         = enrichGO,
  OrgDb       = "org.Hs.eg.db",
  keyType     = "SYMBOL",
  ont         = "BP"
)

write.csv(as.data.frame(CHRNA1.GO), "Fig5_GO_Enrichment_CHRNA1_clusters.csv", row.names = FALSE)

png("Fig5_GO_BP_CHRNA1_clusters.png", width = 11, height = 9, units = "in", res = 500)
dotplot(CHRNA1.GO, title = "GO Enrichment: Biological Process", showCategory = 6) +
  theme(
    plot.title   = element_text(size = 22, colour = "black", face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 19, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 19, colour = "black", face = "bold"),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 17, face = "bold")
  ) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))
dev.off()

# =============================================================================
# 6. HALLMARK GSVA PER CLUSTER – PANEL 5C
# =============================================================================

# Use integrated assay averages per cluster for GSVA
DefaultAssay(GSE115978.A1.combined) <- "integrated"
Average <- AverageExpression(GSE115978.A1.combined, return.seurat = TRUE)

gene.sets.H <- getGeneSets(library = "H", species = "Homo sapiens")

GSVA.Enrich <- gsva(
  Average@assays$integrated@data,
  gene.sets.H,
  min.sz = 5,
  max.sz = 500,
  kcdf  = "Gaussian"
)
colnames(GSVA.Enrich) <- paste0("Cluster ", colnames(GSVA.Enrich))

ROWnames <- lapply(
  substr(gsub("_", " ", gsub("", "", rownames(GSVA.Enrich))), 1, 500),
  function(x) bquote(bold(.(x)))
)
COLnames <- lapply(
  colnames(GSVA.Enrich),
  function(x) bquote(bold(.(x)))
)

png("Fig5C_CHRNA1_Hallmarks_GSVA.png", width = 8, height = 10, units = "in", res = 1200)
pheatmap::pheatmap(
  GSVA.Enrich,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  angle_col    = 90,
  scale        = "column",
  clustering_method = "ward.D2",
  fontsize     = 14,
  labels_col   = as.expression(COLnames),
  labels_row   = as.expression(ROWnames)
)
dev.off()

# =============================================================================
# 7. CELL-CYCLE SCORING, UMAP, AND BAR PLOT – PANEL 5D
# =============================================================================

s.genes  <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cell.cycle <- GSE115978.A1.combined
cell.cycle <- CellCycleScoring(
  cell.cycle,
  s.features  = s.genes,
  g2m.features = g2m.genes,
  set.ident   = TRUE
)

cell.cycle <- RunPCA(cell.cycle, features = c(s.genes, g2m.genes))

p5D_umap <- DimPlot(
  cell.cycle,
  cols       = c("red", "navy", "green"),
  label.size = 8,
  pt.size    = 4
) +
  theme(
    plot.title  = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 12, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
    legend.text = element_text(size = 20, colour = "black", face = "bold")
  )

ggsave("Fig5D_CHRNA1_Cell_Cycle_UMAP.png", p5D_umap, width = 7, height = 5, dpi = 600)

# Cell cycle percentages
cell.cycle.count.table <- table(cell.cycle@meta.data$old.ident, cell.cycle@meta.data$Phase)
Total <- rowSums(cell.cycle.count.table)

# Manual counts from your script
Cluster.name      <- c("Cluster 0", "Cluster 0", "Cluster 0", "Cluster 1", "Cluster 1", "Cluster 1")
Cell.Cycle.State  <- c("G1", "S", "G2M", "G1", "S", "G2M")
Cell.Numbers      <- c(5, 67, 27, 5, 12, 22)
Cell.Percentage   <- c(5, 68, 27, 13, 31, 56)

Cell.Cycle.dataframe <- data.frame(
  Cluster.name,
  Cell.Cycle.State,
  Cell.Numbers,
  Cell.Percentage
)

Cell.Cycle.dataframe$Cluster.name <- factor(
  Cell.Cycle.dataframe$Cluster.name,
  ordered = TRUE,
  levels  = c("Cluster 0", "Cluster 1")
)
Cell.Cycle.dataframe$Cell.Cycle.State <- factor(
  Cell.Cycle.dataframe$Cell.Cycle.State,
  ordered = TRUE,
  levels  = c("S", "G2M", "G1")
)

p5D_bar <- ggplot(Cell.Cycle.dataframe, aes(x = Cluster.name, y = Cell.Percentage, fill = Cell.Cycle.State)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("red", "navy", "green")) +
  ylab("Cell Cycle Phase %") +
  labs(fill = "Cell Cycle Phase") +
  theme_classic() +
  geom_text(
    aes(label = paste(Cell.Percentage, "%")),
    vjust     = -0.5,
    size      = 5,
    fontface  = "bold",
    position  = position_dodge(0.9),
    color     = "black"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title   = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.y  = element_text(size = 11, colour = "black", face = "bold"),
    legend.text  = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold")
  )

ggsave("Fig5D_CHRNA1_Cell_Cycle_Percentage.png", p5D_bar, width = 7, height = 6, dpi = 600)

cat("FIGURE 5 SCRIPT FINISHED.\n")
