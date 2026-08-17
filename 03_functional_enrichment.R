# Title: GO and KEGG Pathway Enrichment Profiling
# Description: Conduct functional annotation for identified DEGs

library(clusterProfiler)
library(org.Hs.eg.db)

# Load significant DEGs
degs <- read.csv("differentially_expressed_genes.csv")
gene_list <- degs$X  # Probe IDs / Gene Symbols

# Map Gene Symbols to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = as.character(gene_list),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)

# GO Biological Process Enrichment
go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# KEGG Pathway Enrichment
kegg_results <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pvalueCutoff  = 0.05
)

# Save enrichment outputs
write.csv(as.data.frame(go_results), "GO_enrichment_results.csv")
write.csv(as.data.frame(kegg_results), "KEGG_enrichment_results.csv")
