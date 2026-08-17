# Title: Data Acquisition and Preprocessing for GSE63514
# Description: Download GEO dataset, clean phenotypic data, and extract expression matrix

library(GEOquery)
library(Biobase)

# Fetch dataset from NCBI GEO
gse_id <- "GSE63514"
gse <- getGEO(gse_id, GSEMatrix = TRUE)

if (length(gse) > 1) idx <- 1 else idx <- 1
gse_data <- gse[[idx]]

# Extract expression matrix and phenotype information
exprs_matrix <- exprs(gse_data)
pheno_data <- pData(gse_data)

# Log2 transformation if needed
max_val <- max(exprs_matrix, na.rm = TRUE)
if (max_val > 100) {
  exprs_matrix <- log2(exprs_matrix + 1)
}

# Save processed matrices
saveRDS(exprs_matrix, file = "exprs_matrix.rds")
saveRDS(pheno_data, file = "pheno_data.rds")
