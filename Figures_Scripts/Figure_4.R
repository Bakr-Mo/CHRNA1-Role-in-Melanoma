# =============================================================================
# FIGURE 4C–D: Stage-based Survival and Hallmark GSVA (TCGA-SKCM)
# Panels:
#   C – Kaplan–Meier survival by melanoma stage (TCGA-SKCM)
#   D – Hallmark GSVA heatmap across melanoma stages (highlight myogenesis)
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "survival", "survminer",
  "Seurat", "edgeR", "biomaRt", "org.Hs.eg.db",
  "escape", "GSVA", "pheatmap"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. LOAD TCGA COUNTS AND CLINICAL (ASSUMES SKCMhg38.rda FROM FIGURE 1)
# =============================================================================

load("SKCMhg38.rda")  # 'data' SummarizedExperiment from HTSeq-Counts

all.genes.pc <- read.csv("DEG-All.csv")  # or your combined annotated count table
# all.genes.pc must contain columns: ensembl_gene_id, external_gene_name, and count columns.

clinical <- data@colData
clinical <- as.data.frame(clinical)

# =============================================================================
# 1. SURVIVAL BY AJCC PATHOLOGIC STAGE – PANEL 4C
# =============================================================================

library(survminer)
library(survival)

clinical.w.Stage <- clinical

clinical.w.Stage[clinical.w.Stage == "Stage IA" | clinical.w.Stage == "Stage IB"] <- "Stage I"
clinical.w.Stage[clinical.w.Stage == "Stage IIA" |
                   clinical.w.Stage == "Stage IIB" |
                   clinical.w.Stage == "Stage IIC"] <- "Stage II"
clinical.w.Stage[clinical.w.Stage == "Stage IIIA" |
                   clinical.w.Stage == "Stage IIIB" |
                   clinical.w.Stage == "Stage IIIC"] <- "Stage III"
clinical.w.Stage[clinical.w.Stage == "Not Reported"] <- NA

clinical.w.Stage$deceased <- clinical.w.Stage$vital_status == "Dead"

clinical.w.Stage$overall_survival <- ifelse(
  clinical.w.Stage$deceased,
  clinical.w.Stage$days_to_death,
  clinical.w.Stage$days_to_last_follow_up
)

fit_stage <- survfit(
  Surv(overall_survival, deceased) ~ ajcc_pathologic_stage,
  data = clinical.w.Stage
)

pval_stage <- surv_pvalue(fit_stage, data = clinical.w.Stage)$pval
print(pval_stage)

Survival_plot <- ggsurvplot(
  fit_stage,
  data       = clinical.w.Stage,
  pval       = FALSE,
  risk.table = TRUE,
  size       = 2.5,
  legend.labs = c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV"),
  legend.title = "Melanoma Stage",
  palette      = "lancet",
  font.legend  = c(25, "bold"),
  font.x       = c(25, "bold"),
  font.y       = c(25, "bold"),
  legend       = c(0.8, 0.8)
)

Survival_plot$plot <- Survival_plot$plot +
  annotate(
    "text",
    x      = 2600,
    y      = 0.1,
    label  = "log-rank p = 0.0008",
    cex    = 6,
    vjust  = 0,
    hjust  = 1.1,
    fontface = 4,
    color  = "blue"
  )

ggsave(
  filename = "Fig4C_TCGA_Stage_Survival.png",
  plot     = print(Survival_plot$plot, newpage = FALSE),
  device   = "png",
  width    = 7,
  height   = 5,
  dpi      = 1200
)
dev.off()

# =============================================================================
# 2. HALLMARK GSVA PER STAGE – PANEL 4D
# =============================================================================

# 2.1 Build Seurat object from non-normalized counts --------------------------

TCGA.NOT.NORM <- all.genes.pc
rownames(TCGA.NOT.NORM) <- TCGA.NOT.NORM$ensembl_gene_id
Gene.SYMBOL <- TCGA.NOT.NORM$external_gene_name
TCGA.NOT.NORM <- TCGA.NOT.NORM[, 2:ncol(TCGA.NOT.NORM)]

TCGA.Seurat <- CreateSeuratObject(TCGA.NOT.NORM, project = "TCGA")
TCGA.Seurat <- NormalizeData(TCGA.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)

Cluster.IDs <- clinical.w.Stage$ajcc_pathologic_stage
TCGA.Seurat <- SetIdent(TCGA.Seurat, value = Cluster.IDs)

# 2.2 Differential expression per stage (optional, for marker selection) -----

DEG.TCGA <- FindAllMarkers(
  TCGA.Seurat,
  only.pos       = FALSE,
  test.use       = "DESeq2",
  min.pct        = 0.25,
  logfc.threshold = 2
)

DEG.TCGA.Pos <- DEG.TCGA[DEG.TCGA$avg_log2FC > 0, ]
DEG.TCGA.Neg <- DEG.TCGA[DEG.TCGA$avg_log2FC < 0, ]

top100.TCGA.markers <- DEG.TCGA %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  as.data.frame()

top100.TCGA.markers <- top100.TCGA.markers %>%
  distinct(gene, .keep_all = TRUE)

TCGA.Gene.Symbol <- bitr(
  top100.TCGA.markers$gene,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = "org.Hs.eg.db"
)

top100.TCGA.markers$SYMBOL <- TCGA.Gene.Symbol$SYMBOL

# 2.3 Average expression and GSVA ---------------------------------------------

TCGA.Average <- AverageExpression(TCGA.Seurat, return.seurat = TRUE)
rownames(TCGA.Average[["RNA"]]@data) <- Gene.SYMBOL

TCGA.Average.table <- TCGA.Average[["RNA"]]@data[
  rownames(TCGA.Average[["RNA"]]@data) %in% top100.TCGA.markers$SYMBOL,
]

gene.sets.H <- getGeneSets(library = "H", species = "Homo sapiens")

TCGA.GSVA.Enrich <- gsva(
  TCGA.Average.table,
  gene.sets.H,
  min.sz = 5,
  max.sz = 500,
  kcdf  = "Gaussian"
)

# Reorder columns to Stage0, Stage I, Stage II, Stage III, Stage IV
TCGA.GSVA.Enrich <- TCGA.GSVA.Enrich[, c("Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV")]

write.csv(TCGA.GSVA.Enrich, "Fig4D_TCGA_Stage_Hallmark_GSVA.csv")

ROWnames <- lapply(
  substr(gsub("_", " ", gsub("", "", rownames(TCGA.GSVA.Enrich))), 1, 500),
  function(x) bquote(bold(.(x)))
)
COLnames <- lapply(
  colnames(TCGA.GSVA.Enrich),
  function(x) bquote(bold(.(x)))
)

png("Fig4D_TCGA_Stage_Hallmarks_Heatmap.png", width = 8, height = 10, units = "in", res = 1200)
pheatmap::pheatmap(
  TCGA.GSVA.Enrich,
  fontsize    = 12,
  angle_col   = 90,
  cluster_cols = FALSE,
  labels_col  = as.expression(COLnames),
  labels_row  = as.expression(ROWnames)
)
dev.off()

cat("FIGURE 4C–D SCRIPT FINISHED.\n")
