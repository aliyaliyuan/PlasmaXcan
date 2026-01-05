## Load Packages

library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(reticulate)

## Convert cRNASeq data to Pseudobulk data
# Read in Seurat object
LD <- readRDS("~/Zenodo/processed_seurat_object/MASLD_snRNA_seq_seurat_v4.rds")

# Filter the count matrix to include only cells expressing more than 200 gene features with each included gene feature present in at least three cells

# Filter cells and genes simultaneously, Keep cells with >200 genes
LD_filtered <- subset(LD, subset = nFeature_RNA > 200)

# Filter genes present in at least 3 cells
# Get counts matrix
counts <- GetAssayData(LD_filtered, assay = "RNA", layer = "counts")

# Identify genes present in at least 3 cells
genesin3 <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]

# Subset the Seurat object to keep only these genes
LD_filtered <- LD_filtered[genesin3, ]

head(LD_filtered)

table(LD_filtered$Cell_type_detailed)

#Remove non-coding sequences by filtering out genes that are absent in the Enformer-derived epigenomic features file 
library(rhdf5)

str(h5ls("/home/aliya/Liver/batch1/HG00096.h5"))

myh5 <- h5read("/home/aliya/Liver/batch1/HG00096.h5",
               name = "chr10_100009946_100009948_predictions")

head(myh5)

dim(myh5) #Should be 5313 4

meta <- read.csv("/home/aliya/Liver/batch1/metadata.csv")

#Debug Print Statements
head(meta)
dim(meta)

dim(LD_filtered)  
head(rownames(LD_filtered)) 
head(colnames(LD_filtered))  

# Load the Enformer metadata
enformer_metadata <- read.csv("/home/aliya/Liver/batch1/metadata.csv")

head(enformer_metadata)
# Get gene names from Seurat object 
seurat_genes <- rownames(LD_filtered)

# Get gene names from Enformer metadata
enformer_genes <- enformer_metadata$external_gene_name

# Find genes present in both
common_genes <- intersect(seurat_genes, enformer_genes)
length(common_genes) #Number of overlapping genes

# Debugging Print Statements
missing_from_enformer <- setdiff(seurat_genes, enformer_genes)
missing_from_seurat <- setdiff(enformer_genes, seurat_genes)

cat("Genes in Seurat:", length(seurat_genes), "\n")
cat("Genes in Enformer:", length(enformer_genes), "\n") 
cat("Common genes:", length(common_genes), "\n")
cat("Missing from Enformer:", length(missing_from_enformer), "\n")
cat("Missing from Seurat:", length(missing_from_seurat), "\n")

# Keep genes found in BOTH Seurat object from liver data AND the Enformer-derived epigenomics data
# Get raw counts
counts <- GetAssayData(LD, assay = "RNA", layer = "counts")

# Keep only the common genes
common_genes <- intersect(rownames(counts), enformer_genes)
counts_sub <- counts[common_genes, , drop = FALSE]

# Build new object with reduced genes but full metadata
LD_filteredEnf <- CreateSeuratObject(counts = counts_sub, meta.data = LD@meta.data)

# Generate individual-by-gene matrix of aggregate mean expression per gene

for (ct in names(celltype_list)) {
  message("Processing cell type: ", ct)
  obj <- celltype_list[[ct]]
  
  # Normalize before aggregation 
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
  
  # Aggregate
  IbG <- AggregateExpression(
    obj,
    assays = "RNA",
    return.seurat = FALSE,
    group.by = "DonorID",
    verbose = TRUE
  )
  
  IbG_mat <- IbG$RNA
  
  # --- Convert to data frame with gene_name column ---
  IbG_df <- as.data.frame(IbG_mat)
  IbG_df <- tibble::rownames_to_column(IbG_df, var = "gene_name")
  
  # --- Write TSV with gene_name explicitly in file ---
  out_path <- paste0("/home/aliya/Liver/batch1/IbG_", ct, ".tsv")
  write.table(
    IbG_df,
    file = out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("File: ", out_path)
}


Heppy <- read.delim('/home/aliya/Liver/batch1/IbG_Hepatocyte_interzonal.tsv', sep = "\t", header = TRUE)


