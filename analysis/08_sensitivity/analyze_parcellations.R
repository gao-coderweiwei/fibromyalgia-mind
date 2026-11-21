# Purpose: Analyze parcellations.
# Inputs: Alternative means and Schaefer 200 results.
# Outputs: Parcellation effect table.
# Manuscript: S5B.

stage <- Sys.getenv("PARCELLATION_STAGE")
project <- Sys.getenv("FM_WORK")
primary <- Sys.getenv("PRIMARY_STAGE")

read_tsv <- function(path) read.delim(path, check.names = FALSE)
baseline <- read_tsv(file.path(primary, "outputs/primary_composite_values.tsv"))
original <- read_tsv(file.path(primary, "outputs/primary_composite_results.tsv"))
intervals <- read_tsv(file.path(Sys.getenv("ROBUST_STAGE"), "model_results.tsv"))
alternative <- read_tsv(file.path(stage, "alternative_composites.tsv"))
endpoints <- c("SalVentAttn_Default", "SalVentAttn_DorsAttn")
stopifnot(
  nrow(alternative) == 264L,
  !anyDuplicated(alternative[c("subject_id", "parcels", "endpoint")]),
  setequal(baseline$subject_id, alternative$subject_id)
)

hedges_g <- function(fm, hc) {
  degrees <- length(fm) + length(hc) - 2
  pooled <- sqrt(((length(fm) - 1) * var(fm) + (length(hc) - 1) * var(hc)) / degrees)
  (1 - 3 / (4 * degrees - 1)) * (mean(fm) - mean(hc)) / pooled
}

estimate <- function(parcels, endpoint) {
  if (parcels == 200) {
    prior <- original[original$network_pair == sub("_", "-", endpoint), ]
    primary_interval <- intervals[intervals$analysis == "primary" &
      intervals$endpoint == prior$network_pair, ]
    return(data.frame(parcels, endpoint,
      n_edges = prior$edge_count, n = prior$n_fm + prior$n_hc,
      beta = prior$beta_group, beta_low = primary_interval$ci95_low,
      beta_high = primary_interval$ci95_high, g = prior$hedges_g,
      g_low = primary_interval$g_ci95_low, g_high = primary_interval$g_ci95_high,
      p_permutation = prior$p_two_sided_permutation, p_holm = prior$p_holm,
      family = "original_2_endpoints"
    ))
  }
  subset <- alternative[alternative$parcels == parcels & alternative$endpoint == endpoint, ]
  stopifnot(setequal(subset$subject_id, baseline$subject_id), length(unique(subset$n_edges)) == 1L)
  data <- data.frame(
    metric = subset$mind[match(baseline$subject_id, subset$subject_id)],
    group_FM = baseline$group_code, age = baseline$age
  )
  model <- lm(metric ~ group_FM + age, data = data, na.action = na.fail)
  limits <- confint(model, "group_FM")
  set.seed(as.integer(Sys.getenv("PERM_SEED")))
  test <- permuco::lmperm(metric ~ group_FM + age,
    data = data, np = as.integer(Sys.getenv("R_PERMUTATIONS")),
    method = "freedman_lane"
  )
  fm <- data$metric[data$group_FM == 1]
  hc <- data$metric[data$group_FM == 0]
  set.seed(as.integer(Sys.getenv("BOOT_SEED")))
  bootstrap <- replicate(as.integer(Sys.getenv("BOOT_DRAWS")), hedges_g(sample(fm, replace = TRUE), sample(hc, replace = TRUE)))
  g_limits <- quantile(bootstrap, c(.025, .975), names = FALSE)
  data.frame(parcels, endpoint,
    n_edges = subset$n_edges[1], n = nrow(data),
    beta = unname(coef(model)["group_FM"]), beta_low = limits[1], beta_high = limits[2],
    g = hedges_g(fm, hc), g_low = g_limits[1], g_high = g_limits[2],
    p_permutation = summary(test)["group_FM", "resampled Pr(>|t|)"],
    p_holm = NA_real_, family = "alternative_4_tests"
  )
}

results <- do.call(rbind, lapply(endpoints, function(endpoint) {
  do.call(rbind, lapply(c(100, 200, 300), estimate, endpoint = endpoint))
}))
alternative_rows <- results$parcels != 200
results$p_holm[alternative_rows] <- p.adjust(results$p_permutation[alternative_rows], "holm")
write.table(results, file.path(stage, "parcellation_effects.tsv"),
  sep = "\t",
  quote = FALSE, row.names = FALSE
)

forest_panel <- function(effect, lower, upper, scale, label, show_labels = FALSE) {
  positions <- c(6, 5, 4, 2, 1, 0)
  original_rows <- results$parcels == 200
  colors <- ifelse(results$endpoint == endpoints[1], "#007A78", "#AD4B42")
  limits <- range(c(0, results[[lower]], results[[upper]]) * scale)
  plot(NA, xlim = limits, ylim = c(-.55, 7.3), xlab = label, ylab = "", axes = FALSE)
  axis(1, cex.axis = .86)
  if (show_labels) axis(2, at = positions, labels = results$parcels, las = 1, tick = FALSE)
  abline(v = 0, lty = 2, col = "#909090")
  segments(results[[lower]] * scale, positions, results[[upper]] * scale, positions,
    col = colors, lwd = 1.45
  )
  points(results[[effect]] * scale, positions,
    pch = ifelse(original_rows, 18, 16),
    cex = ifelse(original_rows, 1.35, .88), col = colors
  )
  text(limits[1], c(6.8, 2.8), c("SalVentAttn-Default", "SalVentAttn-DorsAttn"),
    adj = 0, font = 2, cex = .82
  )
}

pdf(file.path(stage, "parcellation_forest.pdf"), width = 7.3, height = 4.1, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4.7, 3.2, .6, .7), oma = c(.7, 0, 0, 0))
forest_panel(
  "beta", "beta_low", "beta_high", 1000,
  expression("Age-adjusted difference (" * beta %*% 10^3 * ")"), TRUE
)
forest_panel("g", "g_low", "g_high", 1, "Unadjusted Hedges' g")
mtext("Bars: pointwise 95% CIs. Diamonds: original 200-parcel analysis.",
  side = 1,
  outer = TRUE, cex = .75, line = -.2
)
invisible(dev.off())
