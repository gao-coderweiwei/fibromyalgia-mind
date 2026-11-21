# Purpose: Primary robustness.
# Inputs: Primary values and MIND edges; FreeSurfer eTIV.
# Outputs: Model, influence, and eTIV tables.
# Manuscript: Primary results and S5A.

stage <- Sys.getenv("ROBUST_STAGE")
project <- Sys.getenv("FM_WORK")
primary <- Sys.getenv("PRIMARY_STAGE")
networks <- Sys.getenv("MIND_STAGE")
freesurfer <- Sys.getenv("FS_DIR")

read_tsv <- function(path) read.delim(path, check.names = FALSE)
write_tsv <- function(data, name) {
  write.table(data, file.path(stage, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

values <- read_tsv(file.path(primary, "outputs/primary_composite_values.tsv"))
original <- read_tsv(file.path(primary, "outputs/primary_composite_results.tsv"))
edges <- read_tsv(file.path(networks, "outputs/tables/mind_edges_wide.tsv"))
edge_columns <- grep("^edge_[0-9]{5}$", names(edges), value = TRUE)
stopifnot(
  length(edge_columns) == 19900L, !anyDuplicated(values$subject_id),
  !anyDuplicated(edges$subject_id), setequal(values$subject_id, edges$subject_id)
)
values$global_mind <- rowMeans(edges[match(values$subject_id, edges$subject_id), edge_columns])
values$group_FM <- values$group_code
values$eTIV_l <- vapply(values$subject_id, function(subject) {
  measure <- grep("^# Measure EstimatedTotalIntraCranialVol,",
    readLines(file.path(freesurfer, subject, "stats/aseg.stats")),
    value = TRUE
  )
  stopifnot(length(measure) == 1L)
  as.numeric(strsplit(measure, ",", fixed = TRUE)[[1]][4]) / 1e6
}, numeric(1))
stopifnot(all(is.finite(values$global_mind)), all(is.finite(values$eTIV_l)))

hedges_g <- function(fm, hc) {
  df <- length(fm) + length(hc) - 2L
  pooled <- sqrt(((length(fm) - 1) * var(fm) + (length(hc) - 1) * var(hc)) / df)
  (1 - 3 / (4 * df - 1)) * (mean(fm) - mean(hc)) / pooled
}

effect_size <- function(data) {
  fm <- data$metric[data$group_FM == 1]
  hc <- data$metric[data$group_FM == 0]
  set.seed(as.integer(Sys.getenv("BOOT_SEED")))
  draws <- replicate(as.integer(Sys.getenv("BOOT_DRAWS")), hedges_g(sample(fm, replace = TRUE), sample(hc, replace = TRUE)))
  c(hedges_g(fm, hc), unname(quantile(draws, c(.025, .975))))
}

fit_group <- function(analysis, endpoint, data, outcome, covariates, original_p = NULL,
                      report_g = FALSE) {
  data$metric <- data[[outcome]]
  formula <- reformulate(c("group_FM", covariates), response = "metric")
  model <- lm(formula, data = data, na.action = na.fail)
  interval <- confint(model, "group_FM")
  if (is.null(original_p)) {
    set.seed(as.integer(Sys.getenv("PERM_SEED")))
    # Permute reduced-model residuals before refitting the full model.
    test <- permuco::lmperm(formula, data = data, np = as.integer(Sys.getenv("R_PERMUTATIONS")), method = "freedman_lane")
    p_value <- summary(test)["group_FM", "resampled Pr(>|t|)"]
  } else {
    p_value <- original_p
  }
  effect <- if (report_g) effect_size(data) else rep(NA_real_, 3)
  data.frame(analysis, endpoint,
    n_fm = sum(data$group_FM == 1),
    n_hc = sum(data$group_FM == 0), beta = unname(coef(model)["group_FM"]),
    ci95_low = interval[1], ci95_high = interval[2], p_permutation = p_value,
    g = effect[1], g_ci95_low = effect[2], g_ci95_high = effect[3]
  )
}

endpoints <- original$network_pair
metrics <- sub("-", "_", endpoints, fixed = TRUE)
rows <- lapply(seq_along(metrics), function(index) {
  metric <- metrics[index]
  endpoint <- endpoints[index]
  data <- values
  edge_count <- original$edge_count[index]
  data$global_excluding_pair <- (length(edge_columns) * data$global_mind -
    edge_count * data[[metric]]) / (length(edge_columns) - edge_count)
  rbind(
    fit_group("primary", endpoint, data, metric, "age",
      original$p_two_sided_permutation[index],
      report_g = TRUE
    ),
    fit_group("adjust_eTIV", endpoint, data, metric, c("age", "eTIV_l")),
    fit_group(
      "adjust_global_excluding_pair", endpoint, data, metric,
      c("age", "global_excluding_pair")
    )
  )
})

default_metric <- "SalVentAttn_Default"
controls <- which(values$group_FM == 0)
low_control <- controls[which.min(values[[default_metric]][controls])]
influence_model <- lm(SalVentAttn_Default ~ group_FM + age, data = values)
excluded <- fit_group("omit_low_control", endpoints[1], values[-low_control, ],
  default_metric, "age",
  report_g = TRUE
)
full_beta <- unname(coef(influence_model)["group_FM"])
influence <- data.frame(
  subject_id = values$subject_id[low_control],
  mind = values[[default_metric]][low_control],
  cooks_distance = unname(cooks.distance(influence_model)[low_control]),
  dfbetas_group = unname(dfbetas(influence_model)[low_control, "group_FM"]),
  beta_change_percent = 100 * (excluded$beta / full_beta - 1)
)

results <- do.call(rbind, c(rows, list(
  fit_group("global", "Whole cortex", values, "global_mind", "age", report_g = TRUE),
  excluded
)))
results$p_holm <- ave(results$p_permutation, results$analysis, FUN = function(p_values) {
  if (length(p_values) > 1L) p.adjust(p_values, "holm") else NA_real_
})
replicated <- results[results$analysis == "primary", ]
stopifnot(
  isTRUE(all.equal(replicated$beta, original$beta_group)),
  isTRUE(all.equal(replicated$g, original$hedges_g)),
  isTRUE(all.equal(replicated$p_holm, original$p_holm))
)

etiv_fm <- values$eTIV_l[values$group_FM == 1]
etiv_hc <- values$eTIV_l[values$group_FM == 0]
etiv_test <- t.test(etiv_fm, etiv_hc)
etiv <- data.frame(
  unit = "litres", n_fm = length(etiv_fm), n_hc = length(etiv_hc),
  mean_fm = mean(etiv_fm), sd_fm = sd(etiv_fm),
  mean_hc = mean(etiv_hc), sd_hc = sd(etiv_hc),
  difference = mean(etiv_fm) - mean(etiv_hc),
  ci95_low = etiv_test$conf.int[1], ci95_high = etiv_test$conf.int[2],
  p_welch = etiv_test$p.value
)

write_tsv(results, "model_results.tsv")
write_tsv(influence, "low_control_influence.tsv")
write_tsv(etiv, "etiv_summary.tsv")
