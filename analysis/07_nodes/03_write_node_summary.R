# Purpose: 03 write node summary.
# Inputs: Parcel metrics and models.
# Outputs: BH-adjusted parcel results.
# Manuscript: S4.

args <- commandArgs(trailingOnly = TRUE)
stage_dir <- normalizePath(if (length(args) > 0) args[[1]] else getwd())
root_dir <- Sys.getenv("FM_WORK")
network_dir <- Sys.getenv("MIND_STAGE")
nbs_nodes_path <- file.path(Sys.getenv("NBS_STAGE"), "outputs", "nbs_component_nodes.tsv")

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
values <- read_tsv(file.path(stage_dir, "outputs", "node_contribution_values.tsv"))
edge_counts <- read_tsv(file.path(stage_dir, "outputs", "node_contribution_edge_counts.tsv"))
edge_dictionary <- read_tsv(file.path(network_dir, "outputs", "tables", "mind_edge_dictionary.tsv"))
masks <- read_tsv(file.path(network_dir, "outputs", "masks", "analysis_network_masks.tsv"))
pair_counts <- read_tsv(file.path(network_dir, "outputs", "masks", "network_pair_edge_counts.tsv"))
nbs_nodes <- read_tsv(nbs_nodes_path)
pairs <- c("SalVentAttn-Default", "SalVentAttn-DorsAttn")

hedges_g <- function(x_fm, x_hc) {
  n_fm <- length(x_fm)
  n_hc <- length(x_hc)
  pooled_sd <- sqrt(((n_fm - 1) * var(x_fm) + (n_hc - 1) * var(x_hc)) / (n_fm + n_hc - 2))
  (1 - 3 / (4 * (n_fm + n_hc) - 9)) * (mean(x_fm) - mean(x_hc)) / pooled_sd
}

metric_keys <- unique(edge_counts[c("network_pair", "node", "parcel", "network", "edge_count")])
rows <- lapply(seq_len(nrow(metric_keys)), function(index) {
  item <- metric_keys[index, ]
  metric_values <- values[values$network_pair == item$network_pair & values$node == item$node, , drop = FALSE]
  safe_pair <- gsub("-", "_", item$network_pair, fixed = TRUE)
  fit_data <- readRDS(file.path(stage_dir, "outputs", paste0("model_", safe_pair, "_node_", item$node, ".rds")))
  group_row <- fit_data$summary_table["group_FM", , drop = FALSE]
  fm <- metric_values$node_metric[metric_values$group == "FM"]
  hc <- metric_values$node_metric[metric_values$group == "HC"]
  overlap <- any(nbs_nodes$network_pair == item$network_pair & nbs_nodes$direction == "FM_gt_HC" & nbs_nodes$node == item$node)
  data.frame(network_pair = item$network_pair, node = item$node, parcel = item$parcel, network = item$network, edge_count = item$edge_count, n_fm = length(fm), n_hc = length(hc), mean_fm = mean(fm), sd_fm = sd(fm), mean_hc = mean(hc), sd_hc = sd(hc), beta_group = group_row[["Estimate"]], se_group = group_row[["Std. Error"]], t_group = group_row[["t value"]], p_two_sided_permutation = group_row[["resampled Pr(>|t|)"]], p_fdr_module = NA_real_, beta_standardized = group_row[["Estimate"]] * sd(metric_values$group_code) / sd(metric_values$node_metric), hedges_g = hedges_g(fm, hc), model_formula = "node_metric ~ group_FM + age", permutation_method = "freedman_lane", permutations = as.integer(Sys.getenv("R_PERMUTATIONS")), seed = as.integer(Sys.getenv("PERM_SEED")), software_package = "permuco", software_version = fit_data$package_version, in_fm_gt_hc_nbs_component = overlap, stringsAsFactors = FALSE)
})
results <- do.call(rbind, rows)
for (network_pair in pairs) {
  index <- results$network_pair == network_pair
  results$p_fdr_module[index] <- p.adjust(results$p_two_sided_permutation[index], method = "BH")
}
write.table(results, file.path(stage_dir, "outputs", "node_contribution_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

overlap_summary <- do.call(rbind, lapply(pairs, function(network_pair) {
  pair_results <- results[results$network_pair == network_pair, , drop = FALSE]
  data.frame(network_pair = network_pair, n_node_metrics = nrow(pair_results), n_fdr_significant = sum(pair_results$p_fdr_module < 0.05), n_nodes_in_fm_gt_hc_nbs_component = sum(pair_results$in_fm_gt_hc_nbs_component), n_fdr_significant_in_fm_gt_hc_nbs_component = sum(pair_results$p_fdr_module < 0.05 & pair_results$in_fm_gt_hc_nbs_component), stringsAsFactors = FALSE)
}))
write.table(overlap_summary, file.path(stage_dir, "outputs", "node_nbs_overlap_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

pair_edge_match <- vapply(pairs, function(network_pair) {
  nrow(edge_dictionary[edge_dictionary$network_pair == network_pair, , drop = FALSE]) == pair_counts$edge_count[pair_counts$mask_name == network_pair] &&
    nrow(edge_dictionary[edge_dictionary$network_pair == network_pair, , drop = FALSE]) == masks$edge_count[masks$mask_name == network_pair]
}, logical(1))
check_values <- c(
  nrow(unique(values[c("cohort", "subject_id")])) == 66L,
  sum(unique(values[c("subject_id", "group")])$group == "FM") == 33L && sum(unique(values[c("subject_id", "group")])$group == "HC") == 33L,
  all(pairs %in% unique(values$network_pair)) && all(pairs %in% unique(results$network_pair)),
  nrow(values) == 66L * nrow(edge_counts) && all(vapply(split(values, paste(values$network_pair, values$node, sep = "_")), nrow, integer(1)) == 66L),
  all(edge_counts$edge_count == vapply(seq_len(nrow(edge_counts)), function(i) sum(edge_dictionary$network_pair == edge_counts$network_pair[[i]] & (edge_dictionary$node_i == edge_counts$node[[i]] | edge_dictionary$node_j == edge_counts$node[[i]])), integer(1))),
  all(pair_edge_match),
  all(is.finite(values$node_metric)),
  all(is.finite(values$age)),
  all(results$permutations == as.integer(Sys.getenv("R_PERMUTATIONS"))),
  all(results$seed == as.integer(Sys.getenv("PERM_SEED"))),
  all(is.finite(results$p_fdr_module)),
  all(c("network_pair", "direction", "node") %in% names(nbs_nodes)) && any(nbs_nodes$direction == "FM_gt_HC")
)
checks <- data.frame(check = c("subject_count", "group_counts", "required_network_pairs_present", "node_metrics_created", "node_counts_dynamic", "edge_ids_present_in_edge_table", "node_metric_values_finite", "age_finite", "permutations_50000", "seed_recorded", "fdr_by_module_applied", "fm_gt_hc_component_nodes_available"), status = ifelse(check_values, "passed", "failed"), detail = c("66", "FM=33; HC=33", "two pairs", paste(nrow(edge_counts), "node metrics"), "derived from edge dictionary", "dictionary and masks agree", "all values", "all participants", "permuco lmperm np=50000", Sys.getenv("PERM_SEED"), "BH within each network pair", "FM_gt_HC component-node membership"), stringsAsFactors = FALSE)
write.table(checks, file.path(stage_dir, "outputs", "analysis_checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
packages <- data.frame(package = c("permuco", "permute", "Rcpp"), version = c(as.character(packageVersion("permuco")), as.character(packageVersion("permute")), as.character(packageVersion("Rcpp"))))
write.table(packages, file.path(stage_dir, "software", "R_package_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(software = c("Rscript", "R", "permuco_library"), version_or_path = c(Sys.which("Rscript"), R.version.string, Sys.getenv("R_LIBS"))), file.path(stage_dir, "software", "software_versions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
