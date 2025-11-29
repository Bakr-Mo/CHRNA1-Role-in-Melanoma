# =============================================================================
# FIGURE 1: CHRNA1 Expression in Melanoma
# Panels:
#   A – TCGA-SKCM: CHRNA1 expression in Metastatic vs Primary tumors
#   B – GSE65904: CHRNA1 expression in Metastatic vs Primary tumors
#   C – TCGA-SKCM (high CHRNA1): Stage I&II vs Stage III&IV
#   D – TCGA-SKCM (high CHRNA1): N0 vs N3
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "tidyverse", "edgeR", "limma",
  "biomaRt", "ggsignif", "GEOquery", "Biobase"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# PART I – TCGA-SKCM: DOWNLOAD, DE, AND CHRNA1 EXPRESSION
# =============================================================================

# 1. Download HTSeq counts for TCGA-SKCM ---------------------------------------
query.htseq <- GDCquery(
  project      = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

samplesDown <- getResults(query.htseq, cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")  # Primary
dataSmTM <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TM")  # Metastatic

queryDown <- GDCquery(
  project      = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  barcode       = c(dataSmTP, dataSmTM)
)

GDCdownload(queryDown)
prepare <- GDCprepare(queryDown, save = TRUE, save.filename = "SKCMhg38.rda")

load("SKCMhg38.rda")  # loads 'data'

# 2. EdgeR differential expression (Metastatic vs Primary) --------------------
counts <- assay(data)
sample.type <- data$sample_type

DGE <- DGEList(counts, group = sample.type)
group <- factor(sample.type)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
colnames(design)[colnames(design) == "Primary Tumor"] <- "Primary"

keep <- filterByExpr(DGE, design)
DGE <- DGE[keep, , keep.lib.sizes = FALSE]

DGE.normalized <- calcNormFactors(DGE, method = "TMM")

Disp <- estimateDisp(DGE.normalized, design, robust = TRUE)
QL   <- glmQLFit(Disp, design, robust = TRUE)

M.vs.P <- makeContrasts(Metastatic - Primary, levels = design)
TR     <- glmTreat(QL, contrast = M.vs.P, lfc = log2(1.5))

DEG <- topTags(TR, n = Inf, p = 0.05)$table
All.DEG <- QL$table

# 3. Annotate genes via biomaRt -----------------------------------------------
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "m.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

my.attribute <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")

annot.DEG <- getBM(
  attributes = my.attribute,
  filters    = "ensembl_gene_id",
  values     = rownames(DEG),
  mart       = ensembl
)

annot.All.DEG <- getBM(
  attributes = my.attribute,
  filters    = "ensembl_gene_id",
  values     = rownames(All.DEG),
  mart       = ensembl
)

Annot.DEG.combine <- as.data.frame(DEG) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(annot.DEG, "ensembl_gene_id")

Annot.All.DEG.combine <- as.data.frame(All.DEG) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(annot.All.DEG, "ensembl_gene_id")

DEG.PC <- Annot.DEG.combine[Annot.DEG.combine$gene_biotype == "protein_coding", ]
DEG.PC <- DEG.PC %>% distinct(external_gene_name, .keep_all = TRUE)

All.DEG.pc <- Annot.All.DEG.combine[Annot.All.DEG.combine$gene_biotype == "protein_coding", ]
All.DEG.pc <- All.DEG.pc %>% distinct(external_gene_name, .keep_all = TRUE)

rownames(DEG.PC) <- DEG.PC$external_gene_name
DEG.pc <- DEG.PC[, 2:6]

rownames(All.DEG.pc) <- All.DEG.pc$external_gene_name
All.DEG.pc <- All.DEG.pc[, 2:5]

write.csv(DEG.pc,      "DE-genes.csv")
write.csv(All.DEG.pc,  "DEG-All.csv")

# 4. CHRNA1 expression from TR fitted values (Metastatic vs Primary) ---------
DEG.TR.fitted.expr <- as.matrix(TR$fitted.values)

Annot.TR.expr <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = rownames(DEG.TR.fitted.expr),
  mart       = ensembl
)

combined.Annot.TR <- as.data.frame(DEG.TR.fitted.expr) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(Annot.TR.expr, "ensembl_gene_id")

TR.expr.pc <- combined.Annot.TR[combined.Annot.TR$gene_biotype == "protein_coding", ]
TR.expr.pc <- TR.expr.pc %>% distinct(external_gene_name, .keep_all = TRUE)

rownames(TR.expr.pc) <- TR.expr.pc$external_gene_name
TR.expr.pc <- TR.expr.pc[, 2:ncol(TR.expr.pc)]

# Columns are samples; use sample_type from 'data'
colnames(TR.expr.pc) <- data$sample_type
type.TR.expr.pc <- colnames(TR.expr.pc)

transposed.TR.expr.pc <- as.data.frame(t(TR.expr.pc))
transposed.TR.expr.pc$condition <- type.TR.expr.pc
rownames(transposed.TR.expr.pc) <- NULL

# 5. FIGURE 1A – TCGA-SKCM Metastatic vs Primary -----------------------------
p1A <- ggplot(transposed.TR.expr.pc, aes(x = condition, y = CHRNA1, fill = condition)) +
  geom_boxplot(show.legend = FALSE, colour = "black", width = 0.4) +
  scale_fill_brewer(palette = "Set1") +
  ylab(expression(italic("CHRNA1") ~ "Expression Level")) +
  ggtitle("TCGA-SKCM") +
  theme(
    panel.background = element_rect(fill = NA, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgray", size = 0.3),
    axis.title.x  = element_blank(),
    plot.title    = element_text(size = 18, colour = "black", face = "bold", hjust = 0.5),
    axis.title.y  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x   = element_text(size = 13, colour = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, colour = "black", face = "bold")
  ) +
  geom_signif(
    comparisons       = list(c("Metastatic", "Primary Tumor")),
    test              = "t.test",
    map_signif_level  = TRUE,
    col               = "black"
  ) +
  scale_x_discrete(labels = c(
    "Metastatic"    = "Metastatic\n(n = 367)",
    "Primary Tumor" = "Primary\n(n = 103)"
  ))

ggsave("Fig1A_CHRNA1_TCGA_Metastatic_vs_Primary.png", p1A, width = 4, height = 4.5, dpi = 1200)

# =============================================================================
# PART II – TCGA CLINICAL: HIGH CHRNA1, STAGE AND N CLASSIFICATION
# =============================================================================

clinical <- data@colData

# CHRNA1 expression (use same TR expression, averaging if needed)
transposed.TR.expr.pc$barcode <- colnames(TR.expr.pc)
clinical_df <- as.data.frame(clinical)
clinical_df$barcode <- rownames(clinical_df)

merged_expr_clin <- dplyr::left_join(
  clinical_df,
  transposed.TR.expr.pc[, c("barcode", "CHRNA1")],
  by = "barcode"
)

CHRNA1.exp.clin <- data.frame(
  A1.expr = merged_expr_clin$CHRNA1,
  stage   = merged_expr_clin$ajcc_pathologic_stage,
  N       = merged_expr_clin$ajcc_pathologic_n
)

median_value <- median(CHRNA1.exp.clin$A1.expr, na.rm = TRUE)
CHRNA1.exp.clin$level <- ifelse(CHRNA1.exp.clin$A1.expr > median_value, "High", "Low")
CHRNA1.exp.clin <- CHRNA1.exp.clin[CHRNA1.exp.clin$level == "High", ]

# ---- Stage grouping (C) ----
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Stage IA", "Stage IB")]                      <- "Stage I"
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Stage IIA", "Stage IIB", "Stage IIC")]       <- "Stage II"
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC")]    <- "Stage III"
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Not Reported", "Stage 0")]                   <- NA
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Stage I", "Stage II")]                       <- "Stage I & II"
CHRNA1.exp.clin$stage[CHRNA1.exp.clin$stage %in% c("Stage III", "Stage IV")]                     <- "Stage III & IV"

df_stage <- na.omit(CHRNA1.exp.clin[, c("stage", "A1.expr")])

p1C <- ggplot(df_stage, aes(x = stage, y = A1.expr, fill = stage)) +
  geom_boxplot(show.legend = FALSE, colour = "black", width = 0.4) +
  scale_fill_manual(values = c("green", "magenta")) +
  ylab(expression(italic("CHRNA1") ~ "Expression Level")) +
  theme(
    panel.background = element_rect(fill = NA, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgray", size = 0.3),
    axis.title.x  = element_blank(),
    axis.title    = element_text(size = 20, colour = "black"),
    axis.text.x   = element_text(size = 13, colour = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, colour = "black", face = "bold")
  ) +
  geom_signif(
    comparisons      = list(c("Stage I & II", "Stage III & IV")),
    test             = "t.test",
    map_signif_level = TRUE,
    col              = "black"
  )

ggsave("Fig1C_CHRNA1_TCGA_High_Stage.png", p1C, width = 4, height = 4.5, dpi = 1200)

# ---- N classification (D) ----
CHRNA1.exp.clin$N[CHRNA1.exp.clin$N %in% c("N0", "N1", "N2", "N3")] -> CHRNA1.exp.clin$N
df_N <- na.omit(CHRNA1.exp.clin[, c("N", "A1.expr")])

# For the figure they showed N0 vs N3 only:
df_N_filtered <- df_N[df_N$N %in% c("N0", "N3"), ]

p1D <- ggplot(df_N_filtered, aes(x = N, y = A1.expr, fill = N)) +
  geom_boxplot(show.legend = FALSE, colour = "black", width = 0.4) +
  scale_fill_manual(values = c("green", "magenta")) +
  ylab(expression(italic("CHRNA1") ~ "Expression Level")) +
  theme(
    panel.background = element_rect(fill = NA, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgray", size = 0.3),
    axis.title.x  = element_blank(),
    axis.title    = element_text(size = 20, colour = "black"),
    axis.text.x   = element_text(size = 13, colour = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, colour = "black", face = "bold")
  ) +
  geom_signif(
    comparisons      = list(c("N0", "N3")),
    test             = "t.test",
    map_signif_level = TRUE,
    col              = "black"
  )

ggsave("Fig1D_CHRNA1_TCGA_High_N0_vs_N3.png", p1D, width = 4, height = 4.5, dpi = 1200)

# =============================================================================
# PART III – GSE65904: METASTATIC VS PRIMARY (FIG 1B)
# =============================================================================

GSE65904 <- getGEO("GSE65904", GSEMatrix = TRUE)
gse <- GSE65904[[1]]

expr_raw <- Biobase::exprs(gse)
matrix_geo <- edgeR::cpm(expr_raw, prior.count = 2)
sampleInfo <- Biobase::pData(gse)

sampleInfo <- dplyr::select(sampleInfo, "title", "tumor stage:ch1")
sampleInfo <- dplyr::rename(sampleInfo, patient = "title", group = "tumor stage:ch1")

annotation <- Biobase::fData(gse) %>% dplyr::select(Symbol)

rownames(matrix_geo) <- annotation$Symbol
colnames(matrix_geo) <- paste(sampleInfo$group, sep = "-")

matrix_geo <- matrix_geo[, !grepl("^Unknown", colnames(matrix_geo))]
State.geo  <- colnames(matrix_geo)

transposed.geo <- as.data.frame(t(matrix_geo))
transposed.geo$condition <- State.geo
rownames(transposed.geo) <- NULL
transposed.geo <- transposed.geo[, !duplicated(colnames(transposed.geo))]

p1B <- ggplot(transposed.geo, aes(x = condition, y = CHRNA1, fill = condition)) +
  geom_boxplot(show.legend = FALSE, colour = "black", width = 0.4) +
  scale_fill_brewer(palette = "Set1") +
  ylab(expression(italic("CHRNA1") ~ "Expression Level")) +
  ggtitle("GSE65904") +
  theme(
    panel.background = element_rect(fill = NA, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgray", size = 0.3),
    axis.title.x  = element_blank(),
    plot.title    = element_text(size = 18, colour = "black", face = "bold", hjust = 0.5),
    axis.title.y  = element_text(size = 20, colour = "black", face = "bold"),
    axis.text.x   = element_text(size = 13, colour = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, colour = "black", face = "bold")
  ) +
  geom_signif(
    comparisons      = list(c("Metastatic", "Primary")),
    test             = "t.test",
    map_signif_level = TRUE,
    col              = "black",
    y_position       = max(transposed.geo$CHRNA1, na.rm = TRUE) * 1.05
  ) +
  scale_x_discrete(labels = c(
    "Metastatic" = "Metastatic\n(n = 188)",
    "Primary"    = "Primary\n(n = 16)"
  ))

ggsave("Fig1B_CHRNA1_GSE65904_Metastatic_vs_Primary.png", p1B, width = 4, height = 4.5, dpi = 1200)

cat("FIGURE 1 SCRIPT FINISHED.\n")
