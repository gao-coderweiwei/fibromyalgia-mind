# Purpose: 02 run node models.
# Inputs: Parcel metrics and age.
# Outputs: Parcel permutation model RDS.
# Manuscript: S4.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: 02_run_node_models.R <stage_dir> <network_pair> <node>")
stage_dir <- normalizePath(args[[1]])
network_pair <- args[[2]]
node <- as.integer(args[[3]])
library(permuco)

values <- read.delim(file.path(stage_dir, "outputs", "node_contribution_values.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
model_rows <- values[values$network_pair == network_pair & values$node == node, , drop = FALSE]
if (nrow(model_rows) == 0L) stop("Requested node metric is absent")
model_data <- data.frame(node_metric = as.numeric(model_rows$node_metric), group_FM = as.integer(model_rows$group_code), age = as.numeric(model_rows$age))
if (any(!is.finite(as.matrix(model_data)))) stop("Model data contain non-finite values")
set.seed(as.integer(Sys.getenv("PERM_SEED")))
fit <- permuco::lmperm(node_metric ~ group_FM + age, data = model_data, np = as.integer(Sys.getenv("R_PERMUTATIONS")), method = "freedman_lane")
safe_pair <- gsub("-", "_", network_pair, fixed = TRUE)
saveRDS(list(network_pair = network_pair, node = node, summary_table = summary(fit), package_version = as.character(packageVersion("permuco"))), file.path(stage_dir, "outputs", paste0("model_", safe_pair, "_node_", node, ".rds")))
