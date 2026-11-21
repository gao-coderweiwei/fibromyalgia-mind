# Purpose: Qc clinical.
# Inputs: Primary values, BIDS phenotype, and MRIQC.
# Outputs: Participant, IQM, sensitivity, and clinical tables.
# Manuscript: Table 1; S1, S5A, S7.

stage <- Sys.getenv("CLINICAL_STAGE")
project <- Sys.getenv("FM_WORK")
dataset <- Sys.getenv("BIDS_DIR")
primary <- Sys.getenv("PRIMARY_STAGE")

read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = c("", "NA", "n/a"))
write_tsv <- function(data, name) {
  write.table(data, file.path(stage, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
align <- function(data, keys, subject_ids) {
  stopifnot(!anyDuplicated(keys), setequal(keys, subject_ids))
  data[match(subject_ids, keys), , drop = FALSE]
}

mind <- read_tsv(file.path(primary, "outputs/primary_composite_values.tsv"))
participants <- read_tsv(file.path(dataset, "phenotype/participants.tsv"))
clinical <- read_tsv(file.path(dataset, "phenotype/clinical.tsv"))
mriqc <- read_tsv(file.path(Sys.getenv("MRIQC_OUT"), "group_T1w.tsv"))
participants <- align(participants, participants$participant_id, mind$subject_id)
clinical <- align(clinical, clinical$participant_id, mind$subject_id)
mriqc <- align(mriqc, sub("_T1w$", "", mriqc$bids_name), mind$subject_id)
stopifnot(
  all(mind$group == participants$group), all(mind$age == participants$age),
  all(clinical$gp == ifelse(mind$group == "FM", 0, 1))
)

fm_only <- function(measure) ifelse(mind$group == "FM", measure, NA_real_)
clinical_data <- data.frame(
  age_years = mind$age,
  education_years = clinical[["1_yr_study"]],
  bmi = participants$bmi,
  pain_vas = participants$vas_pain_initial,
  ham_d = clinical[["6_hamd_total"]],
  ham_a = clinical[["7_hama_total"]],
  current_major_depression = clinical[["5_mdda"]],
  past_major_depression = clinical[["5_mddp"]],
  widespread_pain_index = clinical[["3_wpi_sc"]],
  symptom_severity = clinical[["3_sss_ts"]],
  years_since_fm_diagnosis = fm_only(participants$fm_diagnosis_years),
  fiq = fm_only(clinical[["12_ts"]]),
  daily_medication_count = fm_only(clinical[["1_no_meds_d"]]),
  any_daily_medication = fm_only(as.integer(clinical[["1_no_meds_d"]] > 0)),
  any_crisis_medication = fm_only(as.integer(clinical$no_meds_c > 0))
)
qc_data <- mriqc[c("cjv", "cnr", "efc", "fber", "fwhm_avg", "snr_total")]
qc_data$fber[which(qc_data$fber == -1)] <- NA_real_

describe <- function(data, binary = character(), fm_only_fields = character(),
                     mean_sd = character(), digits = 2L) {
  rows <- lapply(names(data), function(variable) {
    do.call(rbind, lapply(c("FM", "HC"), function(group) {
      eligible <- group == "FM" || !variable %in% fm_only_fields
      measures <- data[[variable]][mind$group == group]
      observed <- measures[!is.na(measures)]
      count <- length(observed)
      statistic <- if (variable %in% binary) {
        "count (%)"
      } else if (variable %in% mean_sd) {
        "mean (SD)"
      } else {
        "median [Q1, Q3]"
      }
      summary <- if (!eligible) {
        "Not applicable"
      } else if (!count) {
        "Not observed"
      } else if (variable %in% binary) {
        sprintf("%d (%.1f%%)", sum(observed == 1), 100 * mean(observed == 1))
      } else if (variable %in% mean_sd) {
        sprintf("%.2f (%.2f)", mean(observed), sd(observed))
      } else {
        quartiles <- quantile(observed, c(.25, .5, .75))
        sprintf("%.*f [%.*f, %.*f]", digits, quartiles[2], digits, quartiles[1], digits, quartiles[3])
      }
      data.frame(variable, group, eligible,
        n_observed = if (eligible) count else NA_integer_,
        n_missing = if (eligible) length(measures) - count else NA_integer_,
        statistic, summary, minimum = if (count) min(observed) else NA_real_,
        maximum = if (count) max(observed) else NA_real_
      )
    }))
  })
  do.call(rbind, rows)
}

binary_fields <- c(
  "current_major_depression", "past_major_depression",
  "any_daily_medication", "any_crisis_medication"
)
fm_fields <- c(
  "years_since_fm_diagnosis", "fiq", "daily_medication_count",
  "any_daily_medication", "any_crisis_medication"
)

fit_effect <- function(data, term, covariates) {
  model_formula <- reformulate(c(term, covariates), response = "metric")
  fit <- lm(model_formula, data = data, na.action = na.fail)
  interval <- confint(fit, term)
  set.seed(as.integer(Sys.getenv("PERM_SEED")))
  permutation <- permuco::lmperm(model_formula,
    data = data, np = as.integer(Sys.getenv("R_PERMUTATIONS")),
    method = "freedman_lane"
  )
  data.frame(
    n = nrow(data), beta = unname(coef(fit)[term]), ci95_low = interval[1],
    ci95_high = interval[2],
    p_permutation = summary(permutation)[term, "resampled Pr(>|t|)"]
  )
}

metrics <- c("SalVentAttn_Default", "SalVentAttn_DorsAttn")
qc_models <- do.call(rbind, lapply(metrics, function(metric) {
  data <- data.frame(
    metric = mind[[metric]], age = mind$age,
    group_FM = mind$group_code, cjv = qc_data$cjv
  )
  data <- data[complete.cases(data), ]
  baseline <- lm(metric ~ group_FM + age, data = data)
  data.frame(
    endpoint = metric, covariate = "cjv", n_fm = sum(data$group_FM == 1),
    n_hc = sum(data$group_FM == 0), baseline_beta = unname(coef(baseline)["group_FM"]),
    fit_effect(data, "group_FM", c("age", "cjv"))
  )
}))
qc_models$p_holm <- p.adjust(qc_models$p_permutation, "holm")

increments <- c(
  pain_vas = 10, fiq = 10, years_since_fm_diagnosis = 5,
  ham_d = 10, ham_a = 10, any_daily_medication = 1
)
clinical_models <- do.call(rbind, lapply(metrics, function(metric) {
  do.call(rbind, lapply(names(increments), function(variable) {
    data <- data.frame(
      metric = mind[[metric]], age = mind$age,
      clinical_predictor = clinical_data[[variable]] / increments[[variable]]
    )
    data <- data[mind$group == "FM" & complete.cases(data), ]
    data.frame(
      endpoint = metric, predictor = variable, increment = increments[[variable]],
      fit_effect(data, "clinical_predictor", "age")
    )
  }))
}))
clinical_models$q_bh <- p.adjust(clinical_models$p_permutation, "BH")

write_tsv(
  describe(clinical_data, binary_fields, fm_fields, c("age_years", "education_years")),
  "participant_characteristics.tsv"
)
write_tsv(describe(qc_data, digits = 3L), "mriqc_descriptive.tsv")
write_tsv(qc_models, "mriqc_sensitivity.tsv")
write_tsv(clinical_models, "clinical_associations.tsv")
