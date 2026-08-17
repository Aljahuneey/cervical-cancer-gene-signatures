# Title: Differential Expression Analysis Using limma
# Description: Fit linear models to identify DEGs (|log2FC| > 1 and adj. P < 0.05)

library(limma)

# Load saved preprocessed data
exprs_matrix <- readRDS("exprs_matrix.rds")
pheno_data <- readRDS("pheno_data.rds")

# Define target groups (adjust group names according to your sample labels)
groups <- factor(ifelse(grepl("cancer|tumor", pheno_data$title, ignore.case = TRUE), "Cancer", "Normal"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# Linear modeling via limma
fit <- lmFit(exprs_matrix, design)
contrast_matrix <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract top DEGs
deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)

# Filter 155 DEGs based on threshold (|log2FC| > 1 and adj.P.Val < 0.05)
sig_degs <- subset(deg_results, abs(logFC) > 1 & adj.P.Val < 0.05)

# Export results
write.csv(deg_results, "all_expression_results.csv")
write.csv(sig_degs, "differentially_expressed_genes.csv")
