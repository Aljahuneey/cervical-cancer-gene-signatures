# Title: Protein-Protein Interaction (PPI) & Hub Gene Identification
# Description: Extract core hub genes (CDKN2A, CXCR4, CRNN, SPINK5, GBP6, IFI44)

# Define identified hub genes
hub_genes <- c("CDKN2A", "CXCR4", "CRNN", "SPINK5", "GBP6", "IFI44")

# Subset expression profile for hub genes
exprs_matrix <- readRDS("exprs_matrix.rds")
degs <- read.csv("differentially_expressed_genes.csv")

hub_expression <- degs[degs$X %in% hub_genes, ]

# Export hub gene summary
write.csv(hub_expression, "hub_genes_summary.csv")
message("Hub genes successfully exported: ", paste(hub_genes, collapse = ", "))
