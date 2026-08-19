suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: export_astro_midas_umap.R <input.rds> <embedding.tsv> <cells.txt>")
}

obj <- readRDS(args[[1]])
emb <- Embeddings(obj, "midas.umap")
write.table(
  data.frame(cell = rownames(emb), emb, check.names = FALSE),
  file = args[[2]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(rownames(emb), con = args[[3]])
