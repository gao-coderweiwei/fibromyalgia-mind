# Purpose: 02 run primary models.
# Inputs: Primary means and age.
# Outputs: Permutation model RDS.
# Manuscript: Primary results.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: 02_run_primary_models.R <stage_dir> <metric_column>")
stage_dir <- normalizePath(args[[1]])
metric_column <- args[[2]]
library(permuco)

values <- read.delim(file.path(stage_dir, "outputs", "primary_composite_values.tsv"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
if (!metric_column %in% names(values)) stop("Unknown metric column")
model_data <- data.frame(metric = as.numeric(values[[metric_column]]), group_FM = as.integer(values$group_code), age = as.numeric(values$age))
if (any(!is.finite(as.matrix(model_data)))) stop("Model data contain non-finite values")
set.seed(as.integer(Sys.getenv("PERM_SEED")))
fit <- permuco::lmperm(metric ~ group_FM + age, data = model_data, np = as.integer(Sys.getenv("R_PERMUTATIONS")), method = "freedman_lane")
summary_table <- summary(fit)
saveRDS(list(metric_column = metric_column, fit = fit, summary_table = summary_table, package_version = as.character(packageVersion("permuco"))), file.path(stage_dir, "outputs", paste0("model_", metric_column, ".rds")))
