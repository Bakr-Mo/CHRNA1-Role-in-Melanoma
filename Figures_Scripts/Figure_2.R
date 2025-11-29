# =============================================================================
# FIGURE 2: CHRN Correlations and EMT / RHO-ROCK Associations
# Panels:
#   A – TCGA-SKCM: nAChR correlation matrix (Primary vs Metastatic; FPKM)
#   B – GSE65904: nAChR correlation matrix
#   C – CHRNA1/CHRNB1/CHRNG vs ZEB1 (TCGA & GSE65904; scatter + r,p)
#   D – CHRNA1/CHRNB1/CHRNG vs RHO/ROCK genes (TCGA; scatter + r,p)
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "tidyverse", "edgeR",
  "biomaRt", "psych", "Hmisc", "ggplot2", "GEOquery", "Biobase",
  "PerformanceAnalytics"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. TCGA-SKCM: DOWNLOAD FPKM AND ANNOTATE (IF NOT ALREADY DONE)
# =============================================================================

# HTSeq-FPKM query
query.htseq.FPKM <- GDCquery(
  project      = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM"
)

samplesDown.FPKM <- getResults(query.htseq.FPKM, cols = c("cases"))

dataSmTP.FPKM <- TCGAquery_SampleTypes(barcode = samplesDown.FPKM, typesample = "TP")
dataSmTM.FPKM <- TCGAquery_SampleTypes(barcode = samplesDown.FPKM, typesample = "TM")

queryDown.FPKM <- GDCquery(
  project      = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM",
  barcode       = c(dataSmTP.FPKM, dataSmTM.FPKM)
)

GDCdownload(queryDown.FPKM)
prepare.FPKM <- GDCprepare(queryDown.FPKM, save = TRUE, save.filename = "SKCM_FPKM.rda")

load("SKCM_FPKM.rda")   # loads 'data'

SKCM.FPKM <- assay(data)
sample_type <- data$sample_type

# Ensembl annotation
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "m.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

SKCM.FPKM.annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = rownames(SKCM.FPKM),
  mart       = ensembl
)

SKCM.FPKM.annot <- as.data.frame(SKCM.FPKM) %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(SKCM.FPKM.annot, "ensembl_gene_id")

SKCM.FPKM.pc <- SKCM.FPKM.annot[SKCM.FPKM.annot$gene_biotype == "protein_coding", ]
SKCM.FPKM.pc <- SKCM.FPKM.pc %>% distinct(external_gene_name, .keep_all = TRUE)

rownames(SKCM.FPKM.pc) <- SKCM.FPKM.pc$external_gene_name
SKCM.FPKM.pc <- SKCM.FPKM.pc[, 2:ncol(SKCM.FPKM.pc)]

# =============================================================================
# 1. nAChR CORRELATION MATRICES (TCGA FPKM) – PANEL 2A
# =============================================================================

SKCM.PC <- SKCM.FPKM.pc
colnames(SKCM.PC) <- sample_type

transp.SKCM <- as.data.frame(t(SKCM.PC))
rownames(transp.SKCM) <- NULL
transp.SKCM$type <- sample_type

chrn_genes <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7",
  "CHRNA9", "CHRNA10", "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4",
  "CHRND", "CHRNE", "CHRNG"
)

# Primary
Primary.Exp <- transp.SKCM[transp.SKCM$type == "Primary Tumor", ]
Prim.nChRs <- Primary.Exp %>%
  dplyr::select(all_of(chrn_genes))

Prim.Corr.nChRs <- psych::corr.test(Prim.nChRs, method = "pearson")

png("Fig2A_TCGA_Primary_CHRN_Corr_FPKM.png", width = 9, height = 8.5, units = "in", res = 1200)
psych::cor.plot(
  Prim.Corr.nChRs$r,
  short   = FALSE,
  numbers = TRUE,
  cex     = 0.7,
  cex.axis = 0.8,
  xlas    = 2,
  xaxis   = 3,
  alpha   = 1,
  stars   = TRUE,
  pval    = Prim.Corr.nChRs$p,
  upper   = FALSE
)
dev.off()

# Metastatic
Metastatic.Exp <- transp.SKCM[transp.SKCM$type == "Metastatic", ]
Metas.nChRs <- Metastatic.Exp %>%
  dplyr::select(all_of(chrn_genes))

Metas.Corr.nChRs <- psych::corr.test(Metas.nChRs, method = "pearson")

png("Fig2A_TCGA_Metastatic_CHRN_Corr_FPKM.png", width = 9, height = 8.5, units = "in", res = 1200)
psych::cor.plot(
  Metas.Corr.nChRs$r,
  short   = FALSE,
  numbers = TRUE,
  cex     = 0.7,
  cex.axis = 0.8,
  xlas    = 2,
  xaxis   = 3,
  alpha   = 1,
  stars   = TRUE,
  pval    = Metas.Corr.nChRs$p
)
dev.off()

# =============================================================================
# 2. GSE65904 nAChR CORRELATION MATRICES – PANEL 2B
# =============================================================================

GSE65904 <- getGEO("GSE65904", GSEMatrix = TRUE)
gse <- GSE65904[[1]]

expr_raw <- Biobase::exprs(gse)
matrix <- edgeR::cpm(expr_raw, prior.count = 2)
sampleInfo <- Biobase::pData(gse)
sampleInfo <- dplyr::select(sampleInfo, "title", "tumor stage:ch1")
sampleInfo <- dplyr::rename(sampleInfo, patient = "title", group = "tumor stage:ch1")

annotation <- Biobase::fData(gse) %>% dplyr::select(Symbol)
rownames(matrix) <- annotation$Symbol
colnames(matrix) <- paste(sampleInfo$group, sep = "-")

matrix <- matrix[, !grepl("^Unknown", colnames(matrix))]
State.geo <- colnames(matrix)

transposed.geo <- as.data.frame(t(matrix))
transposed.geo$condition <- State.geo
rownames(transposed.geo) <- NULL
transposed.geo <- transposed.geo[, !duplicated(colnames(transposed.geo))]

# Primary
Primary.Geo <- transposed.geo[transposed.geo$condition == "Primary", ]
Prim.nChRs.geo <- Primary.Geo %>% dplyr::select(all_of(chrn_genes))
Prim.Corr.nChRs.geo <- psych::corr.test(Prim.nChRs.geo, method = "pearson")

png("Fig2B_GSE65904_Primary_CHRN_Corr.png", width = 9, height = 8.5, units = "in", res = 1200)
psych::cor.plot(
  Prim.Corr.nChRs.geo$r,
  short   = FALSE,
  numbers = TRUE,
  cex     = 0.7,
  cex.axis = 0.8,
  xlas    = 2,
  xaxis   = 3,
  alpha   = 1,
  stars   = TRUE,
  pval    = Prim.Corr.nChRs.geo$p,
  upper   = FALSE
)
dev.off()

# Metastatic
Metastatic.Geo <- transposed.geo[transposed.geo$condition == "Metastatic", ]
Metas.nChRs.geo <- Metastatic.Geo %>% dplyr::select(all_of(chrn_genes))
Metas.Corr.nChRs.geo <- psych::corr.test(Metas.nChRs.geo, method = "pearson")

png("Fig2B_GSE65904_Metastatic_CHRN_Corr.png", width = 9, height = 8.5, units = "in", res = 1200)
psych::cor.plot(
  Metas.Corr.nChRs.geo$r,
  short   = FALSE,
  numbers = TRUE,
  cex     = 0.7,
  cex.axis = 0.8,
  xlas    = 2,
  xaxis   = 3,
  alpha   = 1,
  stars   = TRUE,
  pval    = Metas.Corr.nChRs.geo$p
)
dev.off()

# =============================================================================
# 3. EMT-TFs AND RHO/ROCK vs CHRNs IN METASTATIC TCGA – PANELS 2C, 2D
# =============================================================================

trans.SKCM.FPKM.pc <- as.data.frame(t(SKCM.FPKM.pc))
rownames(trans.SKCM.FPKM.pc) <- NULL
Sample.Type.FPKM <- sample_type
trans.SKCM.FPKM.pc$type <- Sample.Type.FPKM

Metastatic.FPKM <- trans.SKCM.FPKM.pc[trans.SKCM.FPKM.pc$type == "Metastatic", ]
Metastatic.FPKM <- Metastatic.FPKM[, 1:(ncol(Metastatic.FPKM) - 1)]
Metastatic.FPKM.Log <- log2(Metastatic.FPKM + 1)

# Utility plotting function
plot_corr_scatter <- function(df, xgene, ygene, outfile, xlab = xgene, ylab = ygene) {
  fit  <- lm(reformulate(xgene, ygene), data = df)
  pear <- cor(df[[xgene]], df[[ygene]], method = "pearson")
  
  p <- ggplot(df, aes_string(x = xgene, y = ygene)) +
    geom_point(size = 2.5) +
    geom_smooth(formula = y ~ x, method = "lm", se = FALSE, size = 1.5, color = "blue") +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 18, colour = "black", face = "bold.italic"),
      axis.text.x  = element_text(face = "bold", size = 14),
      axis.text.y  = element_text(face = "bold", size = 14)
    ) +
    xlab(xlab) + ylab(ylab) +
    geom_label(
      aes(
        x = min(df[[xgene]], na.rm = TRUE),
        y = max(df[[ygene]], na.rm = TRUE)
      ),
      size      = 5,
      fontface  = "bold.italic",
      label.size = 0,
      colour    = "blue",
      hjust     = 0,
      label     = paste(
        "p =", signif(summary(fit)$coef[2, 4], 2),
        "\nr =", signif(pear, digits = 3)
      )
    )
  
  ggsave(outfile, p, width = 4, height = 3.5, dpi = 1200)
}

# CHRNG vs ZEB1 example (Fig 2C left bottom panel)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNG", "ZEB1",
  "Fig2C_TCGA_CHRNG_ZEB1_FPKM.png",
  xlab = "CHRNG", ylab = "ZEB1"
)

# Additional combinations for CHRNA1/CHRNB1 vs ZEB1
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNA1", "ZEB1",
  "Fig2C_TCGA_CHRNA1_ZEB1_FPKM.png",
  xlab = "CHRNA1", ylab = "ZEB1"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNB1", "ZEB1",
  "Fig2C_TCGA_CHRNB1_ZEB1_FPKM.png",
  xlab = "CHRNB1", ylab = "ZEB1"
)

# RHO/ROCK correlations (Fig 2D)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNA1", "RHOA",
  "Fig2D_TCGA_CHRNA1_RHOA_FPKM.png",
  xlab = "CHRNA1", ylab = "RHOA"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNA1", "RHOC",
  "Fig2D_TCGA_CHRNA1_RHOC_FPKM.png",
  xlab = "CHRNA1", ylab = "RHOC"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNB1", "ROCK1",
  "Fig2D_TCGA_CHRNB1_ROCK1_FPKM.png",
  xlab = "CHRNB1", ylab = "ROCK1"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNB1", "ROCK2",
  "Fig2D_TCGA_CHRNB1_ROCK2_FPKM.png",
  xlab = "CHRNB1", ylab = "ROCK2"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNB1", "RHOB",
  "Fig2D_TCGA_CHRNB1_RHOB_FPKM.png",
  xlab = "CHRNB1", ylab = "RHOB"
)
plot_corr_scatter(
  Metastatic.FPKM.Log, "CHRNG", "ROCK2",
  "Fig2D_TCGA_CHRNG_ROCK2_FPKM.png",
  xlab = "CHRNG", ylab = "ROCK2"
)

# Optional combined Spearman correlation chart (as in your script)
Metas.nChRs.emt <- Metastatic.FPKM.Log %>%
  dplyr::select(
    CHRNA1, CHRNB1, CHRNG, ZEB1,
    RHOA, RHOB, RHOC, ROCK1, ROCK2
  )

png("Fig2D_TCGA_Spearman_Corr_nAChR_EMT_RHO_ROCK.png", width = 9, height = 8, units = "in", res = 1200)
PerformanceAnalytics::chart.Correlation(Metas.nChRs.emt, histogram = TRUE, pch = 19)
dev.off()

cat("FIGURE 2 SCRIPT FINISHED.\n")
