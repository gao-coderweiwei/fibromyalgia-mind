# Purpose: Analyze surface morphometry.
# Inputs: Parcel morphometry, MIND means, and eTIV.
# Outputs: Morphometry sensitivity and auxiliary regional tables.
# Manuscript: S5C.

stage <- Sys.getenv("SURFACE_STAGE")
project <- Sys.getenv("FM_WORK")
primary <- Sys.getenv("PRIMARY_STAGE")
freesurfer <- Sys.getenv("FS_DIR")

read_tsv <- function(path) read.delim(path, check.names = FALSE)
write_tsv <- function(data, name) {
  write.table(data, file.path(stage, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

mind <- read_tsv(file.path(primary, "outputs/primary_composite_values.tsv"))
mind$group_FM <- mind$group_code
mind$eTIV_l <- vapply(mind$subject_id, function(subject) {
  line <- grep("^# Measure EstimatedTotalIntraCranialVol,",
    readLines(file.path(freesurfer, subject, "stats/aseg.stats")),
    value = TRUE
  )
  as.numeric(strsplit(line, ",", fixed = TRUE)[[1]][4]) / 1e6
}, numeric(1))
measures <- read_tsv(file.path(stage, "parcel_measures.tsv"))
stopifnot(
  setequal(mind$subject_id, measures$subject_id), length(unique(measures$parcel)) == 200L,
  all(table(measures$subject_id, measures$parcel) == 1L)
)
regional_data <- merge(measures, mind[c("subject_id", "group_FM", "age", "eTIV_l")],
  by = "subject_id", sort = FALSE
)
features <- c("thickness_mm", "area_mm2", "volume_mm3")

group_effect <- function(model) {
  group <- coef(summary(model))["group_FM", ]
  interval <- confint(model, "group_FM")
  data.frame(
    n = nobs(model), beta_fm_minus_hc = unname(group["Estimate"]),
    ci95_low = interval[1], ci95_high = interval[2], p_ols = unname(group["Pr(>|t|)"])
  )
}

region_results <- do.call(rbind, lapply(split(regional_data, regional_data$parcel), function(data) {
  do.call(rbind, lapply(features, function(feature) {
    covariates <- c("group_FM", "age", if (feature == "volume_mm3") "eTIV_l")
    fit <- lm(reformulate(covariates, response = feature), data = data, na.action = na.fail)
    outcome <- data[[feature]]
    data.frame(
      parcel = data$parcel[1], hemisphere = data$hemisphere[1], network = data$network[1],
      feature, mean_fm = mean(outcome[data$group_FM == 1]),
      mean_hc = mean(outcome[data$group_FM == 0]), group_effect(fit)
    )
  }))
}))
region_results$q_bh_within_feature <- ave(region_results$p_ols, region_results$feature,
  FUN = function(p) p.adjust(p, "BH")
)

network_means <- aggregate(thickness_mm ~ subject_id + network, measures, mean)
network_sums <- aggregate(cbind(area_mm2, volume_mm3) ~ subject_id + network, measures, sum)
network_data <- merge(network_means, network_sums, by = c("subject_id", "network"))

targets <- c("Default", "DorsAttn")
sensitivity <- do.call(rbind, lapply(targets, function(target) {
  endpoint <- paste("SalVentAttn", target, sep = "_")
  data <- mind
  data$metric <- data[[endpoint]]
  baseline <- lm(metric ~ group_FM + age, data = data, na.action = na.fail)
  baseline_beta <- unname(coef(baseline)["group_FM"])
  salience <- network_data[network_data$network == "SalVentAttn", ]
  comparison <- network_data[network_data$network == target, ]

  do.call(rbind, lapply(features, function(feature) {
    data$salience_value <- salience[[feature]][match(data$subject_id, salience$subject_id)]
    data$target_value <- comparison[[feature]][match(data$subject_id, comparison$subject_id)]
    covariates <- c(
      "group_FM", "age", "salience_value", "target_value",
      if (feature == "volume_mm3") "eTIV_l"
    )
    formula <- reformulate(covariates, response = "metric")
    fit <- lm(formula, data = data, na.action = na.fail)
    set.seed(as.integer(Sys.getenv("PERM_SEED")))
    # Permute reduced-model residuals before refitting the full model.
    test <- permuco::lmperm(formula, data = data, np = as.integer(Sys.getenv("R_PERMUTATIONS")), method = "freedman_lane")
    effect <- group_effect(fit)
    data.frame(endpoint, feature, baseline_beta, effect,
      beta_change_percent = 100 * (effect$beta_fm_minus_hc / baseline_beta - 1),
      p_permutation = summary(test)["group_FM", "resampled Pr(>|t|)"]
    )
  }))
}))
sensitivity$p_holm <- p.adjust(sensitivity$p_permutation, "holm")
sensitivity$p_ols <- NULL

write_tsv(region_results, "regional_group_effects.tsv")
write_tsv(sensitivity, "mind_morphometry_sensitivity.tsv")
