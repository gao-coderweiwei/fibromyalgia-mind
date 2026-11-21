# Purpose: 01 prepare node metrics.
# Inputs: MIND edges, masks, and NBS nodes.
# Outputs: Parcel means and edge counts.
# Manuscript: S4.

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
pair_counts <- read_tsv(file.path(network_dir, "outputs", "masks", "network_pair_edge_counts.tsv"))

pairs <- c("SalVentAttn-Default", "SalVentAttn-DorsAttn")
if (!all(pairs %in% masks$mask_name) || !all(pairs %in% pair_counts$mask_name)) stop("Required network-pair masks are missing")
if (!all(pairs %in% edge_dictionary$network_pair)) stop("Required network pairs are missing from the edge dictionary")

metadata <- wide[c("cohort", "subject_id", "group", "sex", "age")]
metadata$group_code <- ifelse(metadata$group == "FM", 1L, ifelse(metadata$group == "HC", 0L, NA_integer_))
if (anyNA(metadata$group_code)) stop("Unexpected group coding")

value_rows <- list()
count_rows <- list()
row_index <- 1L
for (network_pair in pairs) {
  pair_edges <- edge_dictionary[edge_dictionary$network_pair == network_pair, , drop = FALSE]
  expected_count <- pair_counts$edge_count[pair_counts$mask_name == network_pair]
  mask_count <- masks$edge_count[masks$mask_name == network_pair]
  if (nrow(pair_edges) != expected_count || nrow(pair_edges) != mask_count) stop("Network-pair edge count mismatch")
  nodes <- sort(unique(c(pair_edges$node_i, pair_edges$node_j)))
  for (node in nodes) {
    incident <- pair_edges[pair_edges$node_i == node | pair_edges$node_j == node, , drop = FALSE]
    if (nrow(incident) == 0L) stop("Node has no incident network-pair edges")
    node_network <- unique(c(incident$network_i[incident$node_i == node], incident$network_j[incident$node_j == node]))
    node_parcel <- unique(c(incident$parcel_i[incident$node_i == node], incident$parcel_j[incident$node_j == node]))
    if (length(node_network) != 1L || length(node_parcel) != 1L) stop("Node annotation is inconsistent")
    edge_columns <- sprintf("edge_%05d", incident$edge_id)
    if (!all(edge_columns %in% names(wide))) stop("Edge identifier is absent from the wide table")
    metrics <- rowMeans(wide[edge_columns])
    value_rows[[row_index]] <- data.frame(metadata, network_pair = network_pair, node = node, parcel = node_parcel, network = node_network, edge_count = nrow(incident), node_metric = metrics, stringsAsFactors = FALSE)
    count_rows[[row_index]] <- data.frame(network_pair = network_pair, node = node, parcel = node_parcel, network = node_network, edge_count = nrow(incident), stringsAsFactors = FALSE)
    row_index <- row_index + 1L
  }
}
values <- do.call(rbind, value_rows)
values <- values[c("cohort", "subject_id", "group", "group_code", "sex", "age", "network_pair", "node", "parcel", "network", "edge_count", "node_metric")]
edge_counts <- do.call(rbind, count_rows)
write.table(values, file.path(stage_dir, "outputs", "node_contribution_values.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_counts, file.path(stage_dir, "outputs", "node_contribution_edge_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
parameters <- data.frame(parameter = c("model", "permutation_method", "permutations", "seed", "multiple_testing", "network_pairs"), value = c("node_metric ~ group_FM + age", "freedman_lane", Sys.getenv("R_PERMUTATIONS"), Sys.getenv("PERM_SEED"), "BH-FDR within network pair", paste(pairs, collapse = "; ")))
write.table(parameters, file.path(stage_dir, "config", "analysis_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
