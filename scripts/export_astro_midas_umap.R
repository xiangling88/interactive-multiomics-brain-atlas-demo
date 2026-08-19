suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: export_astro_midas_umap.R <input.rds> <input_meta.txt> <embedding.tsv> <aligned_meta.txt>")
}

rds_path <- args[[1]]
meta_path <- args[[2]]
embedding_path <- args[[3]]
aligned_meta_path <- args[[4]]

obj <- readRDS(rds_path)
print(obj)
print(Reductions(obj))

astro_umap <- Embeddings(obj, reduction = "midas.umap")
print(dim(astro_umap))
print(head(rownames(astro_umap)))
print(colnames(astro_umap))

rds_cells <- rownames(astro_umap)
cat("RDS cells:", length(rds_cells), "\n")
cat("RDS unique cells:", length(unique(rds_cells)), "\n")
if (length(rds_cells) != length(unique(rds_cells))) {
  stop("Duplicated RDS midas.umap cell names")
}

astro_meta <- read.delim(
  meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)
cat("Astro meta dimensions:", paste(dim(astro_meta), collapse = " x "), "\n")
print(colnames(astro_meta))
print(head(astro_meta))

raw_meta_cells <- rownames(astro_meta)
candidate_cells <- list(
  rownames_raw = raw_meta_cells,
  rownames_cell_prefix = paste0("cell_", raw_meta_cells)
)
for (col in c("cell", "cell_id", "barcode", "Cell", "Unnamed: 0", "Unnamed..0", "X")) {
  if (col %in% colnames(astro_meta)) {
    candidate_cells[[col]] <- as.character(astro_meta[[col]])
    candidate_cells[[paste0(col, "_cell_prefix")]] <- paste0("cell_", as.character(astro_meta[[col]]))
  }
}

matches <- vapply(candidate_cells, function(x) sum(rds_cells %in% x), numeric(1))
print(matches)
best_name <- names(which.max(matches))
best_cells <- candidate_cells[[best_name]]
cat("Selected Astro meta cell column/rule:", best_name, "\n")
cat("Matched:", max(matches), "\n")
cat("Unmatched:", length(rds_cells) - max(matches), "\n")
if (max(matches) != length(rds_cells)) {
  missing <- rds_cells[!(rds_cells %in% best_cells)]
  print(head(missing, 20))
  stop("Astrocyte meta did not match RDS midas.umap cells at 100%")
}

idx <- match(rds_cells, best_cells)
if (any(is.na(idx))) {
  stop("NA match index after 100% membership check")
}

astro_meta_new <- astro_meta[idx, , drop = FALSE]
astro_meta_new <- data.frame(
  cell_id = rds_cells,
  astro_meta_new,
  check.names = FALSE
)
if (!all(as.character(astro_meta_new[["cell_id"]]) == rds_cells)) {
  stop("Astrocyte cell order mismatch after reorder")
}
if (!("celltype_astro_my" %in% colnames(astro_meta_new))) {
  stop("celltype_astro_my not found in astrocyte metadata")
}
print(table(astro_meta_new[["celltype_astro_my"]], useNA = "ifany"))

astro_umap_df <- data.frame(
  cell_id = rds_cells,
  UMAP1 = astro_umap[, 1],
  UMAP2 = astro_umap[, 2],
  check.names = FALSE
)
write.table(
  astro_umap_df,
  file = embedding_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  astro_meta_new,
  file = aligned_meta_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[PASS] Astrocyte RDS midas.umap and metadata alignment\n")
