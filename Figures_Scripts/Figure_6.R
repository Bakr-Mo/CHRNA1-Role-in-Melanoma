# =============================================================================
# FIGURE 6: High-CHRNA1 Melanoma vs Rhabdomyosarcomatous Melanoma
# Panels:
#   A – Venn diagram: myogenesis genes enriched in TCGA high-CHRNA1 vs scRNA
#   B – Heatmap: muscle-related genes in
#        - TCGA low-CHRNA1 samples (average)
#        - TCGA high-CHRNA1 samples (average)
#        - Rhabdomyosarcomatous melanoma samples (average)
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "tidyverse", "ComplexHeatmap", "grid", "ggvenn"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. INPUTS
# =============================================================================
# Assumes:
#   - SKCM.FPKM.pc: protein-coding FPKM matrix (genes x TCGA samples),
#                   rownames = gene symbols, columns = TCGA barcodes
#   - "Rhabdo_Melanoma_TPM.csv": TPM expression for rhabdomyosarcomatous melanoma
#       columns: gene, sample1, sample2, ...
#   - "Myo_Genes2.csv": table with columns 'TCGA' and 'scRNA' listing myogenesis genes
#       enriched in TCGA high-CHRNA1 and in scRNA CHRNA1+ clusters.

# If SKCM.FPKM.pc not in memory, load from an RDS you created earlier.

# =============================================================================
# 1. DEFINE HIGH/LOW CHRNA1 GROUPS IN TCGA-SKCM
# =============================================================================

FPKM.Count <- SKCM.FPKM.pc

t.FPKM.Count <- as.data.frame(t(FPKM.Count))
median_value_A1_FPKM <- median(t.FPKM.Count$CHRNA1)
summary(t.FPKM.Count$CHRNA1)

# Use your original cutoffs
t.FPKM.Count$level <- ifelse(
  t.FPKM.Count$CHRNA1 >= 1.4166, "High",
  ifelse(t.FPKM.Count$CHRNA1 <= 0.1474, "Low", "Middle")
)
table(t.FPKM.Count$level)

High.Exp.A1.Samples <- rownames(t.FPKM.Count[t.FPKM.Count$level == "High", ])
Low.Exp.A1.Samples  <- rownames(t.FPKM.Count[t.FPKM.Count$level == "Low", ])

High.Exp.A1.Count <- FPKM.Count[, colnames(FPKM.Count) %in% High.Exp.A1.Samples]
colnames(High.Exp.A1.Count) <- paste(replicate(ncol(High.Exp.A1.Count), "TCGA High Expression CHRNA1"))

Low.Exp.A1.Count <- FPKM.Count[, colnames(FPKM.Count) %in% Low.Exp.A1.Samples]
colnames(Low.Exp.A1.Count) <- paste(replicate(ncol(Low.Exp.A1.Count), "TCGA Low Expression CHRNA1"))

# =============================================================================
# 2. RHABDOMYOSARCOMATOUS MELANOMA TPM AND COMBINED MATRIX – PANEL 6B
# =============================================================================

Rhabdo.Melanoma.TPM <- read.csv("Rhabdo_Melanoma_TPM.csv", header = TRUE, sep = ",")
Rhabdo.Melanoma.TPM <- Rhabdo.Melanoma.TPM %>%
  distinct(gene, .keep_all = TRUE)

rownames(Rhabdo.Melanoma.TPM) <- Rhabdo.Melanoma.TPM$gene
Rhabdo.Melanoma.TPM <- Rhabdo.Melanoma.TPM[, 2:ncol(Rhabdo.Melanoma.TPM)]

colnames(Rhabdo.Melanoma.TPM) <- rep("Rhabdomyosarcomatous Melanoma", ncol(Rhabdo.Melanoma.TPM))

# Match genes between TCGA FPKM and Rhabdo TPM
Matched.Genes <- match(rownames(FPKM.Count), rownames(Rhabdo.Melanoma.TPM))
Rhabdo.Melanoma.TPM <- Rhabdo.Melanoma.TPM[Matched.Genes, ]
rownames(Rhabdo.Melanoma.TPM) <- rownames(FPKM.Count)

Melanoma.AND.Rhabdo <- cbind(High.Exp.A1.Count, Low.Exp.A1.Count, Rhabdo.Melanoma.TPM)

# Averages for each group
Averg.High.A1.FPKM <- rowMeans(Melanoma.AND.Rhabdo[, colnames(Melanoma.AND.Rhabdo) == "TCGA High Expression CHRNA1"])
Averg.Low.A1.FPKM  <- rowMeans(Melanoma.AND.Rhabdo[, colnames(Melanoma.AND.Rhabdo) == "TCGA Low Expression CHRNA1"])
Averg.Rhabdo.TPM   <- rowMeans(Melanoma.AND.Rhabdo[, colnames(Melanoma.AND.Rhabdo) == "Rhabdomyosarcomatous Melanoma"])

Average.Combined <- data.frame(
  "TCGA Low Expression CHRNA1"  = Averg.Low.A1.FPKM,
  "TCGA High Expression CHRNA1" = Averg.High.A1.FPKM,
  "Rhabdomyosarcomatous Melanoma" = Averg.Rhabdo.TPM
)

# Log2 transform for comparable scale
Average.Combined.Log <- log2(Average.Combined + 1e-6)

# =============================================================================
# 3. SUBSET TO MYOGENESIS GENES PRESENT IN TCGA & scRNA – PANEL 6B
# =============================================================================

Myo.Genes.In.TCGA.and.scRNA <- c(
  "CHRNA1", "LAMA2", "CSRP3", "CHRNG", "SLN", "SGCA", "STC2", "TNNI1", "TNNI2", "NCAM1",
  "MYOG", "CHRNB1", "DMPK", "CAV3", "TNNC1", "MYBPH", "MYL4", "MYLPF", "ACTC1", "DES",
  "MYL1", "BIN1", "TNNT2", "TNNT3", "MYF6", "MYH8", "CDH13", "PVALB", "MYL6B",
  "SPARC", "ENO3", "ATP2A1", "TPM2", "TNNT1"
)

Combined.Myo.Genes <- Average.Combined.Log[
  rownames(Average.Combined.Log) %in% Myo.Genes.In.TCGA.and.scRNA,
]
Combined.Myo.Genes <- as.matrix(Combined.Myo.Genes)

colnames(Combined.Myo.Genes) <- c(
  "Low Expressing-CHRNA1\nTCGA-SKCM samples",
  "High Expressing-CHRNA1\nTCGA-SKCM samples",
  "Rhabdomyosarcomatous\nMelanoma"
)

# Heatmap (Panel 6B)
MyoGenes.heatmap <- ComplexHeatmap::Heatmap(
  Combined.Myo.Genes,
  column_title      = NULL,
  name              = "Expression",
  show_row_names    = TRUE,
  show_column_names = TRUE,
  show_column_dend  = FALSE,
  show_row_dend     = FALSE,
  column_names_gp   = grid::gpar(fontsize = 12, fontface = "bold"),
  row_names_gp      = grid::gpar(fontsize = 12, fontface = "bold.italic")
)

png("Fig6B_CHRNA1_Myogenesis_Heatmap_High_Low_Rhabdo.png", width = 4, height = 9, units = "in", res = 1200)
draw(MyoGenes.heatmap)
dev.off()

# =============================================================================
# 4. VENN DIAGRAM: MYOGENESIS GENES IN TCGA VS scRNA – PANEL 6A
# =============================================================================

Myo.Genes <- read.csv("Myo_Genes2.csv", header = TRUE, sep = ",")

# Intersection (for reporting)
inters <- intersect(Myo.Genes$TCGA, Myo.Genes$scRNA)
Intersected.df <- data.frame(`Intersected Genes` = inters)
write.csv(Intersected.df, "Fig6A_intersected_Myo_genes.csv", row.names = FALSE)

# Venn via ggvenn
venn_list <- list(
  "TCGA"  = Myo.Genes$TCGA,
  "scRNA" = Myo.Genes$scRNA
)

p_venn <- ggvenn(
  venn_list,
  set_name_size = 10,
  text_color    = "black",
  text_size     = 6
)

ggsave("Fig6A_Myogenesis_TCGA_scRNA_Venn.png", p_venn, width = 8, height = 6, dpi = 1200)

cat("FIGURE 6 SCRIPT FINISHED.\n")
