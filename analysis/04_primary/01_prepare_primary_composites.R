# Purpose: Prepare the primary cross-network MIND measures.
# Inputs: MIND edge tables and masks.
# Outputs: Primary means and demographics.
# Manuscript: Primary results.

args <- commandArgs(trailingOnly = TRUE)
stage_dir <- normalizePath(if (length(args) > 0) args[[1]] else file.path(dirname(sys.frame(1)$ofile), ".."))
root_dir <- Sys.getenv("FM_WORK")
network_dir <- Sys.getenv("MIND_STAGE")

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
wide_path <- file.path(network_dir, "outputs", "tables", "mind_edges_wide.tsv")
if (!file.exists(wide_path)) wide_path <- paste0(wide_path, ".gz")
wide <- read_tsv(wide_path)
edge_dictionary <- read_tsv(file.path(network_dir, "outputs", "tables", "mind_edge_dictionary.tsv"))
masks <- read_tsv(file.path(network_dir, "outputs", "masks", "analysis_network_masks.tsv"))

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

default_columns <- mask_columns("SalVentAttn-Default")
dorsal_columns <- mask_columns("SalVentAttn-DorsAttn")
output <- wide[c("cohort", "subject_id", "group", "sex", "age")]
output$group_code <- ifelse(output$group == "FM", 1L, ifelse(output$group == "HC", 0L, NA_integer_))
if (anyNA(output$group_code)) stop("Unexpected group encoding")
output$SalVentAttn_Default <- rowMeans(wide[default_columns])
output$SalVentAttn_DorsAttn <- rowMeans(wide[dorsal_columns])
output <- output[c("cohort", "subject_id", "group", "group_code", "sex", "age", "SalVentAttn_Default", "SalVentAttn_DorsAttn")]
write.table(output, file.path(stage_dir, "outputs", "primary_composite_values.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

demographics <- do.call(rbind, lapply(c("FM", "HC"), function(group_name) {
  ages <- as.numeric(output$age[output$group == group_name])
  data.frame(group = group_name, n = length(ages), age_mean = mean(ages), age_sd = sd(ages), age_min = min(ages), age_max = max(ages))
}))
write.table(demographics, file.path(stage_dir, "outputs", "group_demographics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(parameter = c("model", "permutation_method", "permutations", "seed", "multiple_testing"), value = c("metric ~ group_FM + age", "freedman_lane", Sys.getenv("R_PERMUTATIONS"), Sys.getenv("PERM_SEED"), "Holm across two primary endpoints")), file.path(stage_dir, "config", "analysis_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
