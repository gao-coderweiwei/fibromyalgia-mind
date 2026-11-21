# Purpose: 01 prepare specificity composites.
# Inputs: MIND edge tables and masks.
# Outputs: Reference means and direct contrasts.
# Manuscript: S2.

args <- commandArgs(trailingOnly = TRUE)
stage_dir <- normalizePath(if (length(args) > 0) args[[1]] else getwd())
root_dir <- Sys.getenv("FM_WORK")
network_dir <- Sys.getenv("MIND_STAGE")

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
wide_path <- file.path(network_dir, "outputs", "tables", "mind_edges_wide.tsv")
if (!file.exists(wide_path)) wide_path <- paste0(wide_path, ".gz")
wide <- read_tsv(wide_path)
edge_dictionary <- read_tsv(file.path(network_dir, "outputs", "tables", "mind_edge_dictionary.tsv"))
masks <- read_tsv(file.path(network_dir, "outputs", "masks", "analysis_network_masks.tsv"))
metric_masks <- data.frame(
  metric = c("within_SalVentAttn", "within_Default", "within_DorsAttn", "Default_DorsAttn", "SalVentAttn_SomMot", "Default_SomMot", "DorsAttn_SomMot", "SalVentAttn_Default", "SalVentAttn_DorsAttn"),
  mask_name = c("within SalVentAttn", "within Default", "within DorsAttn", "Default-DorsAttn", "SalVentAttn-SomMot", "Default-SomMot", "DorsAttn-SomMot", "SalVentAttn-Default", "SalVentAttn-DorsAttn"),
  stringsAsFactors = FALSE
)

mask_columns <- function(mask_name) {
  selected <- masks[masks$mask_name == mask_name, , drop = FALSE]
  if (nrow(selected) != 1) stop("Expected one mask row for ", mask_name)
  edge_ids <- as.integer(strsplit(selected$edge_ids[[1]], ",", fixed = TRUE)[[1]])
  if (length(edge_ids) != selected$edge_count[[1]]) stop("Mask edge count mismatch for ", mask_name)
  if (!all(edge_ids %in% as.integer(edge_dictionary$edge_id))) stop("Mask edge identifier not found in edge dictionary")
  columns <- sprintf("edge_%05d", edge_ids)
  if (!all(columns %in% names(wide))) stop("Mask edge column not found in wide table")
  columns
}

output <- wide[c("cohort", "subject_id", "group", "sex", "age")]
output$group_code <- ifelse(output$group == "FM", 1L, ifelse(output$group == "HC", 0L, NA_integer_))
if (anyNA(output$group_code)) stop("Unexpected group encoding")
for (index in seq_len(nrow(metric_masks))) {
  output[[metric_masks$metric[[index]]]] <- rowMeans(wide[mask_columns(metric_masks$mask_name[[index]])])
}
output$SalVentAttn_Default_minus_Default_DorsAttn <- output$SalVentAttn_Default - output$Default_DorsAttn
output$SalVentAttn_DorsAttn_minus_Default_DorsAttn <- output$SalVentAttn_DorsAttn - output$Default_DorsAttn
output$SalVentAttn_DorsAttn_minus_SalVentAttn_Default <- output$SalVentAttn_DorsAttn - output$SalVentAttn_Default
output <- output[c("cohort", "subject_id", "group", "group_code", "sex", "age", metric_masks$metric, "SalVentAttn_Default_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_SalVentAttn_Default")]
write.table(output, file.path(stage_dir, "outputs", "network_specificity_composite_values.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

edge_counts <- data.frame(metric = metric_masks$metric, mask_name = metric_masks$mask_name, edge_count = sapply(metric_masks$mask_name, function(name) masks$edge_count[masks$mask_name == name]))
edge_counts <- rbind(edge_counts, data.frame(metric = c("SalVentAttn_Default_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_SalVentAttn_Default"), mask_name = c("SalVentAttn-Default minus Default-DorsAttn", "SalVentAttn-DorsAttn minus Default-DorsAttn", "SalVentAttn-DorsAttn minus SalVentAttn-Default"), edge_count = NA_integer_))
write.table(edge_counts, file.path(stage_dir, "outputs", "network_specificity_edge_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(parameter = c("model", "permutation_method", "permutations", "seed", "family_A", "family_B", "family_C"), value = c("metric ~ group_FM + age", "freedman_lane", Sys.getenv("R_PERMUTATIONS"), Sys.getenv("PERM_SEED"), "3 within-network composites", "4 cross-network reference composites", "3 composite differences")), file.path(stage_dir, "config", "analysis_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
