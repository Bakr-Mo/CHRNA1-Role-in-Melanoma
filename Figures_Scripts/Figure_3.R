# =============================================================================
# FIGURE 3A–B: Enrichment of CHRNA1-correlated Genes in Metastatic TCGA-SKCM
# Panels:
#   A – GO/BP bar plot (enrichR output)
#   B – GSEA Hallmark myogenesis enrichment curve
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "tidyverse", "edgeR",
  "biomaRt", "Hmisc", "enrichR", "msigdbr", "fgsea", "data.table",
  "ggplot2"
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) BiocManager::install(new_pkgs, ask = FALSE)
lapply(required_pkgs, library, character.only = TRUE)

set.seed(1234)
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. FPKM PREP (reuse from Figure_2.R or load SKCM_FPKM.rda)
# =============================================================================

load("SKCM_FPKM.rda")  # 'data' object with FPKM

SKCM.FPKM <- assay(data)
sample_type <- data$sample_type

# Annotation
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
# 1. CORRELATIONS BETWEEN DEGs (UP-REGULATED) AND CHRNA1 (OPTIONAL SECTION)
# =============================================================================
# This block reproduces your DEGs correlation, but Figure 3 uses A1-only correlations.
# You can skip re-running DEG filtering here if already done.

# Example assumes DEG.pc (from edgeR treat) is already saved:
# load or read "DE-genes.csv" if needed
if (file.exists("DE-genes.csv")) {
  DEG.pc <- read.csv("DE-genes.csv", header = TRUE, row.names = 1)
} else {
  stop("DE-genes.csv not found; run differential expression script first.")
}

transposed.SKCM.FPKM.pc <- as.data.frame(t(SKCM.FPKM.pc))
rownames(transposed.SKCM.FPKM.pc) <- NULL

Sample.Type.FPKM <- sample_type

DEGs.FPKM <- transposed.SKCM.FPKM.pc[, colnames(transposed.SKCM.FPKM.pc) %in% rownames(DEG.pc[DEG.pc$logFC > 0, ])]

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row    = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor    = (cormat)[ut],
    p      = pmat[ut]
  )
}

DEGs.Cor.matrix.FPKM <- Hmisc::rcorr(as.matrix(DEGs.FPKM), type = "pearson")
corANDsig.DEGs.FPKM <- flattenCorrMatrix(DEGs.Cor.matrix.FPKM$r, DEGs.Cor.matrix.FPKM$P)

DEGs.cor.with.A1.FPKM <- corANDsig.DEGs.FPKM[
  corANDsig.DEGs.FPKM$column == "CHRNA1" | corANDsig.DEGs.FPKM$row == "CHRNA1",
]
High.cor.DEGs.A1.FPKM <- DEGs.cor.with.A1.FPKM[
  DEGs.cor.with.A1.FPKM$cor >= 0.3 | DEGs.cor.with.A1.FPKM$cor <= -0.3,
]

write.csv(High.cor.DEGs.A1.FPKM, file = "High.cor.DEGs.A1.FPKM.csv", row.names = FALSE)
Adj.High.cor.DEGs.A1.FPKM <- read.csv("High.cor.DEGs.A1.FPKM.csv", header = TRUE)

# =============================================================================
# 2. CORRELATED GENES TO CHRNA1 IN METASTATIC SAMPLES ONLY (FPKM) – CORE PART
# =============================================================================

SKCM.FPKM.pc.SampleType <- SKCM.FPKM.pc
colnames(SKCM.FPKM.pc.SampleType) <- sample_type

SKCM.FPKM.pc.Metastatic <- t(SKCM.FPKM.pc.SampleType[, colnames(SKCM.FPKM.pc.SampleType) == "Metastatic"])
rownames(SKCM.FPKM.pc.Metastatic) <- NULL

Cor.Mtx.FPKM.Metast <- Hmisc::rcorr(SKCM.FPKM.pc.Metastatic, type = "pearson")
corANDsig.FPKM.Metast <- flattenCorrMatrix(Cor.Mtx.FPKM.Metast$r, Cor.Mtx.FPKM.Metast$P)

Cor.with.A1 <- corANDsig.FPKM.Metast[
  corANDsig.FPKM.Metast$column == "CHRNA1" | corANDsig.FPKM.Metast$row == "CHRNA1",
]

Cor.A1 <- Cor.with.A1[Cor.with.A1$cor >= 0.3 & Cor.A1$p < 0.05, ]
Cor.A1 <- Cor.A1 %>% drop_na()

write.csv(Cor.A1, file = "CorA1_FPKM_metastatic.csv", row.names = FALSE)
Adj.cor.A1 <- read.csv("CorA1_FPKM_metastatic.csv", header = TRUE)

# Prepare gene list ranked by correlation
Adj.cor.A1$gene <- Adj.cor.A1$row
de.genes <- Adj.cor.A1 %>%
  arrange(desc(cor)) %>%
  dplyr::select(gene, cor)

ranks <- deframe(de.genes)

# =============================================================================
# 3. ENRICHMENT FOR CHRNA1-CORRELATED GENES – PANEL 3A (enrichR)
# =============================================================================

library(enrichR)

dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2016", "MSigDB_Hallmark_2020")
gene_symbols <- unique(Adj.cor.A1$gene)

enriched <- enrichr(gene_symbols, dbs)

# GO BP bar plot (panel 3A)
go_bp <- enriched[["GO_Biological_Process_2021"]]

png("Fig3A_GO_BP_CHRNA1_correlated.png", width = 10, height = 6, units = "in", res = 1200)
plotEnrich(
  go_bp,
  orderBy = "P.value",
  numChar = 100,
  showTerms = 10,
  title = "GO: Biological Process",
  xlab  = ""
) +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 15, colour = "black", face = "bold"),
    axis.text.y  = element_text(size = 15, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 12, colour = "black", face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
dev.off()

write.csv(go_bp, "Fig3A_GO_BP_CHRNA1_correlated_table.csv", row.names = FALSE)

# =============================================================================
# 4. HALLMARK GSEA (HALLMARK_MYOGENESIS) – PANEL 3B
# =============================================================================

HALL <- msigdbr(category = "H", species = "Homo sapiens")
fgsea_sets <- HALL %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fwrite(fgseaResTidy, file = "Fig3B_Hallmark_GSEA_CHRNA1_correlated.csv")

# Hallmark MYOGENESIS enrichment curve
png("Fig3B_HALLMARK_MYOGENESIS_CHRNA1.png", width = 6, height = 5, units = "in", res = 1200)
plotEnrichment(
  fgsea_sets[["HALLMARK_MYOGENESIS"]],
  ranks
) +
  labs(title = "HALLMARK_MYOGENESIS") +
  theme_bw() +
  geom_line(color = "red", size = 2) +
  theme(
    plot.title   = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 15, colour = "black", face = "bold"),
    axis.title.y = element_text(size = 15, colour = "black", face = "bold"),
    axis.text.y  = element_text(size = 15, colour = "black", face = "bold"),
    axis.text.x  = element_text(size = 12, colour = "black", face = "bold"),
    plot.margin  = margin(20, 20, 20, 20)
  )
dev.off()

cat("FIGURE 3A–B SCRIPT FINISHED.\n")
