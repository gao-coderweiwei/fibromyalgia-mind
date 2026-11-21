# Purpose: 03 write summary tables.
# Inputs: Primary model RDS and means.
# Outputs: Primary results and checks.
# Manuscript: Primary results.

args <- commandArgs(trailingOnly = TRUE)
stage_dir <- normalizePath(if (length(args) > 0) args[[1]] else file.path(dirname(sys.frame(1)$ofile), ".."))
root_dir <- Sys.getenv("FM_WORK")
network_dir <- Sys.getenv("MIND_STAGE")

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
values <- read_tsv(file.path(stage_dir, "outputs", "primary_composite_values.tsv"))
wide_path <- file.path(network_dir, "outputs", "tables", "mind_edges_wide.tsv")
if (!file.exists(wide_path)) wide_path <- paste0(wide_path, ".gz")
wide <- read_tsv(wide_path)
edge_columns <- grep("^edge_[0-9]{5}$", names(wide), value = TRUE)
masks <- read_tsv(file.path(network_dir, "outputs", "masks", "analysis_network_masks.tsv"))
edge_dictionary <- read_tsv(file.path(network_dir, "outputs", "tables", "mind_edge_dictionary.tsv"))
specification <- read_tsv(file.path(Sys.getenv("ATLAS_STAGE"), "config", "analysis_specification.tsv"))

endpoint_map <- data.frame(network_pair = c("SalVentAttn-Default", "SalVentAttn-DorsAttn"), metric_column = c("SalVentAttn_Default", "SalVentAttn_DorsAttn"), stringsAsFactors = FALSE)
hedges_g <- function(x_fm, x_hc) {
  n_fm <- length(x_fm)
  n_hc <- length(x_hc)
  pooled_sd <- sqrt(((n_fm - 1) * var(x_fm) + (n_hc - 1) * var(x_hc)) / (n_fm + n_hc - 2))
  correction <- 1 - 3 / (4 * (n_fm + n_hc) - 9)
  correction * (mean(x_fm) - mean(x_hc)) / pooled_sd
}

rows <- lapply(seq_len(nrow(endpoint_map)), function(index) {
  endpoint <- endpoint_map[index, ]
  fit_data <- readRDS(file.path(stage_dir, "outputs", paste0("model_", endpoint$metric_column, ".rds")))
  table <- fit_data$summary_table
  group_row <- table["group_FM", , drop = FALSE]
  metric <- as.numeric(values[[endpoint$metric_column]])
  fm <- metric[values$group == "FM"]
  hc <- metric[values$group == "HC"]
  mask_row <- masks[masks$mask_name == endpoint$network_pair, , drop = FALSE]
  data.frame(
    network_pair = endpoint$network_pair,
    edge_count = as.integer(mask_row$edge_count),
    n_fm = length(fm), n_hc = length(hc),
    mean_fm = mean(fm), sd_fm = sd(fm), mean_hc = mean(hc), sd_hc = sd(hc),
    beta_group = group_row[["Estimate"]], se_group = group_row[["Std. Error"]], t_group = group_row[["t value"]],
    p_two_sided_permutation = group_row[["resampled Pr(>|t|)"]], p_holm = NA_real_,
    beta_standardized = group_row[["Estimate"]] * sd(values$group_code) / sd(metric),
    hedges_g = hedges_g(fm, hc), model_formula = "metric ~ group_FM + age", permutation_method = "freedman_lane",
    permutations = as.integer(Sys.getenv("R_PERMUTATIONS")), seed = as.integer(Sys.getenv("PERM_SEED")), software_package = "permuco", software_version = fit_data$package_version,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)
results$p_holm <- p.adjust(results$p_two_sided_permutation, method = "holm")
write.table(results, file.path(stage_dir, "outputs", "primary_composite_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

mask_present <- all(endpoint_map$network_pair %in% masks$mask_name)
mask_ids_present <- all(unlist(lapply(endpoint_map$network_pair, function(name) as.integer(strsplit(masks$edge_ids[masks$mask_name == name][[1]], ",", fixed = TRUE)[[1]]))) %in% as.integer(edge_dictionary$edge_id))
checks <- data.frame(
  check = c("subject_count", "group_counts", "edge_count_per_subject", "primary_masks_present", "mask_edge_ids_present_in_edge_table", "SalVentAttn_Default_edge_count", "SalVentAttn_DorsAttn_edge_count", "composite_values_finite", "age_finite", "permutations_50000", "seed_recorded", "holm_correction_applied"),
  status = c(nrow(values) == 66, sum(values$group == "FM") == 33 && sum(values$group == "HC") == 33, length(edge_columns) == 19900 && all(rowSums(is.finite(as.matrix(wide[edge_columns]))) == 19900), mask_present, mask_ids_present, results$edge_count[results$network_pair == "SalVentAttn-Default"] == 1012, results$edge_count[results$network_pair == "SalVentAttn-DorsAttn"] == 572, all(is.finite(as.matrix(values[, endpoint_map$metric_column]))), all(is.finite(as.numeric(values$age))), all(results$permutations == as.integer(Sys.getenv("R_PERMUTATIONS"))), all(results$seed == as.integer(Sys.getenv("PERM_SEED"))), all(is.finite(results$p_holm))),
  detail = c(as.character(nrow(values)), "FM=33; HC=33", "19900 finite upper-triangle edges per subject", as.character(mask_present), as.character(mask_ids_present), as.character(results$edge_count[results$network_pair == "SalVentAttn-Default"]), as.character(results$edge_count[results$network_pair == "SalVentAttn-DorsAttn"]), "both composites", "all participants", "permuco lmperm np=50000", Sys.getenv("PERM_SEED"), "two endpoints", "not run", "read-only inputs"),
  stringsAsFactors = FALSE
)
checks$status <- ifelse(checks$status, "passed", "failed")
write.table(checks, file.path(stage_dir, "outputs", "analysis_checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

package_rows <- data.frame(package = c("permuco", "permute", "Rcpp"), version = c(as.character(packageVersion("permuco")), as.character(packageVersion("permute")), as.character(packageVersion("Rcpp"))))
write.table(package_rows, file.path(stage_dir, "software", "R_package_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
software_rows <- data.frame(software = c("Rscript", "R", "permuco_library"), version_or_path = c(Sys.which("Rscript"), R.version.string, file.path(stage_dir, "software", "R")))
write.table(software_rows, file.path(stage_dir, "software", "software_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
