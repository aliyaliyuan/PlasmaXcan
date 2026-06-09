# =============================================================================
# Pseudobulk Aggregation: Hepatocyte_pericentral
# =============================================================================

## Load Packages
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(rhdf5)
library(tibble)

# =============================================================================
# 1. Load & Filter Seurat Object
# =============================================================================

LD <- readRDS("~/Zenodo/processed_seurat_object/MASLD_snRNA_seq_seurat_v4.rds")

# Keep cells with >200 gene features
LD_filtered <- subset(LD, subset = nFeature_RNA > 200)

# Keep genes present in at least 3 cells
counts <- GetAssayData(LD_filtered, assay = "RNA", layer = "counts")
genesin3 <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]
LD_filtered <- LD_filtered[genesin3, ]

cat("Dimensions after cell/gene filtering:\n")
print(dim(LD_filtered))

# =============================================================================
# 2. Filter to Enformer Genes
# =============================================================================

enformer_metadata <- read.csv("/home/aliya/Liver/metadata.csv")
enformer_genes <- enformer_metadata$external_gene_name

seurat_genes <- rownames(LD_filtered)
common_genes <- intersect(seurat_genes, enformer_genes)

cat("Genes in Seurat:         ", length(seurat_genes), "\n")
cat("Genes in Enformer:       ", length(enformer_genes), "\n")
cat("Common genes:            ", length(common_genes), "\n")
cat("Missing from Enformer:   ", length(setdiff(seurat_genes, enformer_genes)), "\n")
cat("Missing from Seurat:     ", length(setdiff(enformer_genes, seurat_genes)), "\n")

# Subset counts to common genes and rebuild Seurat object
counts <- GetAssayData(LD, assay = "RNA", layer = "counts")
common_genes <- intersect(rownames(counts), enformer_genes)
counts_sub <- counts[common_genes, , drop = FALSE]

LD_filteredEnf <- CreateSeuratObject(counts = counts_sub, meta.data = LD@meta.data)

cat("Dimensions after Enformer gene filtering:\n")
print(dim(LD_filteredEnf))
cat("\nCell type counts:\n")
print(table(LD_filteredEnf$Cell_type_detailed))

# =============================================================================
# 3. Subset to Hepatocyte_pericentral
# =============================================================================

hep_obj <- subset(LD_filteredEnf, subset = Cell_type_detailed == "Hepatocyte_pericentral")

cat("\nHepatocyte_pericentral cells: ", ncol(hep_obj), "\n")

# =============================================================================
# 4. Normalize & Aggregate by Donor
# =============================================================================

hep_obj <- NormalizeData(hep_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

IbG <- AggregateExpression(
  hep_obj,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "DonorID",
  verbose = TRUE
)

# =============================================================================
# 5. Write Output
# =============================================================================

IbG_df <- as.data.frame(IbG$RNA)
IbG_df <- tibble::rownames_to_column(IbG_df, var = "gene_name")

dir.create("/home/aliya/Liver/newIbG/", recursive = TRUE, showWarnings = FALSE)

out_path <- "/home/aliya/Liver/newIbG/IbG_Hepatocyte_pericentral.tsv"

write.table(
  IbG_df,
  file = out_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nFile written to:", out_path, "\n")

# =============================================================================
# 6. Verify Output
# =============================================================================

Heppy <- read.delim(out_path)
cat("Dimensions of output (genes x donors):", dim(Heppy), "\n")
head(Heppy)

print("Done")

