# Purpose: 03 write specificity summary.
# Inputs: Reference metrics and models.
# Outputs: Family-wise results table.
# Manuscript: S2.

args <- commandArgs(trailingOnly = TRUE)
stage_dir <- normalizePath(if (length(args) > 0) args[[1]] else getwd())
root_dir <- Sys.getenv("FM_WORK")
network_dir <- Sys.getenv("MIND_STAGE")

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
values <- read_tsv(file.path(stage_dir, "outputs", "network_specificity_composite_values.tsv"))
edge_counts <- read_tsv(file.path(stage_dir, "outputs", "network_specificity_edge_counts.tsv"))
masks <- read_tsv(file.path(network_dir, "outputs", "masks", "analysis_network_masks.tsv"))
edge_dictionary <- read_tsv(file.path(network_dir, "outputs", "tables", "mind_edge_dictionary.tsv"))
wide_path <- file.path(network_dir, "outputs", "tables", "mind_edges_wide.tsv")
if (!file.exists(wide_path)) wide_path <- paste0(wide_path, ".gz")
wide <- read_tsv(wide_path)
edge_columns <- grep("^edge_[0-9]{5}$", names(wide), value = TRUE)
specification <- data.frame(
  analysis_family = c(rep("within_network", 3), rep("cross_network_reference", 4), rep("direct_composite_difference", 3)),
  metric = c("within_SalVentAttn", "within_Default", "within_DorsAttn", "Default_DorsAttn", "SalVentAttn_SomMot", "Default_SomMot", "DorsAttn_SomMot", "SalVentAttn_Default_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_SalVentAttn_Default"),
  stringsAsFactors = FALSE
)
hedges_g <- function(x_fm, x_hc) {
  n_fm <- length(x_fm)
  n_hc <- length(x_hc)
  pooled_sd <- sqrt(((n_fm - 1) * var(x_fm) + (n_hc - 1) * var(x_hc)) / (n_fm + n_hc - 2))
  (1 - 3 / (4 * (n_fm + n_hc) - 9)) * (mean(x_fm) - mean(x_hc)) / pooled_sd
}
rows <- lapply(seq_len(nrow(specification)), function(index) {
  item <- specification[index, ]
  fit_data <- readRDS(file.path(stage_dir, "outputs", paste0("model_", item$metric, ".rds")))
  group_row <- fit_data$summary_table["group_FM", , drop = FALSE]
  metric_values <- as.numeric(values[[item$metric]])
  fm <- metric_values[values$group == "FM"]
  hc <- metric_values[values$group == "HC"]
  edge_count <- edge_counts$edge_count[edge_counts$metric == item$metric]
  data.frame(analysis_family = item$analysis_family, metric = item$metric, edge_count = if (length(edge_count) == 1) edge_count else NA_integer_, n_fm = length(fm), n_hc = length(hc), mean_fm = mean(fm), sd_fm = sd(fm), mean_hc = mean(hc), sd_hc = sd(hc), beta_group = group_row[["Estimate"]], se_group = group_row[["Std. Error"]], t_group = group_row[["t value"]], p_two_sided_permutation = group_row[["resampled Pr(>|t|)"]], p_holm_family = NA_real_, beta_standardized = group_row[["Estimate"]] * sd(values$group_code) / sd(metric_values), hedges_g = hedges_g(fm, hc), model_formula = "metric ~ group_FM + age", permutation_method = "freedman_lane", permutations = as.integer(Sys.getenv("R_PERMUTATIONS")), seed = as.integer(Sys.getenv("PERM_SEED")), software_package = "permuco", software_version = fit_data$package_version, stringsAsFactors = FALSE)
})
results <- do.call(rbind, rows)
for (family_name in unique(results$analysis_family)) {
  index <- results$analysis_family == family_name
  results$p_holm_family[index] <- p.adjust(results$p_two_sided_permutation[index], method = "holm")
}
write.table(results, file.path(stage_dir, "outputs", "network_specificity_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

required_masks <- edge_counts$mask_name[!is.na(edge_counts$edge_count)]
mask_ids <- unlist(lapply(required_masks, function(name) as.integer(strsplit(masks$edge_ids[masks$mask_name == name][[1]], ",", fixed = TRUE)[[1]])))
checks <- data.frame(
  check = c("subject_count", "group_counts", "edge_count_per_subject", "all_required_masks_present", "mask_edge_ids_present_in_edge_table", "within_network_edge_counts", "cross_network_reference_edge_counts", "difference_metrics_created", "composite_values_finite", "age_finite", "permutations_50000", "seed_recorded", "holm_correction_by_family_applied"),
  status = c(nrow(values) == 66, sum(values$group == "FM") == 33 && sum(values$group == "HC") == 33, length(edge_columns) == 19900 && all(rowSums(is.finite(as.matrix(wide[edge_columns]))) == 19900), all(required_masks %in% masks$mask_name), all(mask_ids %in% as.integer(edge_dictionary$edge_id)), all(!is.na(edge_counts$edge_count[edge_counts$metric %in% c("within_SalVentAttn", "within_Default", "within_DorsAttn")])), all(!is.na(edge_counts$edge_count[edge_counts$metric %in% c("Default_DorsAttn", "SalVentAttn_SomMot", "Default_SomMot", "DorsAttn_SomMot")])), all(c("SalVentAttn_Default_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_Default_DorsAttn", "SalVentAttn_DorsAttn_minus_SalVentAttn_Default") %in% names(values)), all(is.finite(as.matrix(values[specification$metric]))), all(is.finite(as.numeric(values$age))), all(results$permutations == as.integer(Sys.getenv("R_PERMUTATIONS"))), all(results$seed == as.integer(Sys.getenv("PERM_SEED"))), all(is.finite(results$p_holm_family))),
  detail = c(as.character(nrow(values)), "FM=33; HC=33", "19900 finite upper-triangle edges per subject", "true", "true", "three masks", "four masks", "three metrics", "all metrics", "all participants", "permuco lmperm np=50000", Sys.getenv("PERM_SEED"), "Holm within A/B/C", "not run", "read-only inputs"),
  stringsAsFactors = FALSE
)
checks$status <- ifelse(checks$status, "passed", "failed")
write.table(checks, file.path(stage_dir, "outputs", "analysis_checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
packages <- data.frame(package = c("permuco", "permute", "Rcpp"), version = c(as.character(packageVersion("permuco")), as.character(packageVersion("permute")), as.character(packageVersion("Rcpp"))))
write.table(packages, file.path(stage_dir, "software", "R_package_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(software = c("Rscript", "R", "permuco_library"), version_or_path = c(Sys.which("Rscript"), R.version.string, paste(.libPaths(), collapse = .Platform$path.sep))), file.path(stage_dir, "software", "software_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
