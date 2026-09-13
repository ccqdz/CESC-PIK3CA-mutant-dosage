# ================================================================
# TCGA-CESC PIK3CA mutant dosage study
# FINAL LOCKED REPRODUCIBLE ANALYSIS SCRIPT
# Target journal: Gynecologic Oncology
# Repository-ready version: relative paths only
# ================================================================
#
# This script starts from the five validated patient-level RDS files:
#   master_tcga_cesc.rds
#   pik_final_tcga_cesc.rds
#   hpv_patient_tcga_cesc.rds
#   msk_sample_master.rds
#   msk_pik_mut.rds
#
# It reproduces the locked analyses, manuscript tables, main figures,
# supplementary figures, QC checkpoints, and sessionInfo().
#
# IMPORTANT LOCKED DEFINITIONS
# 1. TCGA tumor_vaf is stored as percentage (0-100): VAF = tumor_vaf / 100.
# 2. Mutant CN = VAF * [purity * total_CN + 2*(1-purity)] / purity.
# 3. FAM = mutant_CN / total_CN.
# 4. WGD multiple-testing family: 5 prespecified endpoints.
# 5. 3q multiple-testing family: 4 adjusted mutant-dosage endpoints.
# 6. Total CN and mutant CN 3q models adjust for ploidy + WGD.
# 7. Mutant CN/ploidy 3q model adjusts for WGD only.
# 8. FAM 3q model adjusts for ploidy + WGD.
# 9. Bulk sequencing does NOT directly phase PIK3CA mutation to major/minor homolog.
# 10. MSK CNA=2 is reported as a discrete high-level CNA call, not as a
#     continuous ABSOLUTE-equivalent amplification measurement.
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260912)

# ----------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------
required_pkgs <- c("data.table", "ggplot2", "patchwork", "scales", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first with install.packages()."
  )
}

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)

# ----------------------------------------------------------------
# 1. Paths
# ----------------------------------------------------------------
# Repository mode:
#   Rscript analysis/01_locked_analysis.R
#
# Optional custom repository root:
#   Rscript analysis/01_locked_analysis.R /path/to/repository
#
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  PROJECT_DIR <- normalizePath(args[1], mustWork = TRUE)
} else {
  PROJECT_DIR <- normalizePath(getwd(), mustWork = TRUE)
}

INPUT_DIR  <- file.path(PROJECT_DIR, "data", "derived")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
TAB_DIR    <- file.path(OUTPUT_DIR, "tables")

for (d in c(OUTPUT_DIR, FIG_DIR, TAB_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

MASTER_RDS     <- file.path(INPUT_DIR, "master_tcga_cesc.rds")
PIK_RDS        <- file.path(INPUT_DIR, "pik_final_tcga_cesc.rds")
HPV_RDS        <- file.path(INPUT_DIR, "hpv_patient_tcga_cesc.rds")
MSK_MASTER_RDS <- file.path(INPUT_DIR, "msk_sample_master.rds")
MSK_MUT_RDS    <- file.path(INPUT_DIR, "msk_pik_mut.rds")

required_files <- c(MASTER_RDS, PIK_RDS, HPV_RDS, MSK_MASTER_RDS, MSK_MUT_RDS)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required RDS files:\n", paste(missing_files, collapse = "\n"))
}

# ----------------------------------------------------------------
# 2. Load inputs
# ----------------------------------------------------------------
master            <- as.data.table(readRDS(MASTER_RDS))
pik_final         <- as.data.table(readRDS(PIK_RDS))
hpv_patient       <- as.data.table(readRDS(HPV_RDS))
msk_sample_master <- as.data.table(readRDS(MSK_MASTER_RDS))
msk_pik_mut       <- as.data.table(readRDS(MSK_MUT_RDS))

# ----------------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------------
pick_col <- function(dat, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    if (is.null(label)) label <- paste(candidates, collapse = " / ")
    stop(
      "Required column not found for ", label, ".\nTried: ",
      paste(candidates, collapse = ", "), "\nAvailable columns:\n",
      paste(names(dat), collapse = ", ")
    )
  }
  NA_character_
}

safe_wilcox <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 3 || length(unique(g[ok])) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x[ok] ~ g[ok], exact = FALSE)$p.value)
}

med_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(median = NA_real_, q1 = NA_real_, q3 = NA_real_))
  c(
    median = median(x),
    q1 = unname(quantile(x, 0.25)),
    q3 = unname(quantile(x, 0.75))
  )
}

extract_term <- function(model, term, endpoint, model_label) {
  sm <- coef(summary(model))
  if (!(term %in% rownames(sm))) {
    return(data.table(
      endpoint = endpoint, model = model_label,
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_
    ))
  }
  ci <- suppressMessages(confint(model, parm = term))
  data.table(
    endpoint = endpoint,
    model = model_label,
    beta = unname(sm[term, "Estimate"]),
    se = unname(sm[term, "Std. Error"]),
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    p = unname(sm[term, "Pr(>|t|)"])
  )
}

fmt_p <- function(p, prefix = "P") {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(paste0(prefix, "<0.001"))
  paste0(prefix, "=", sprintf("%.3f", p))
}

near_equal <- function(x, target, tol) {
  isTRUE(all.equal(as.numeric(x), as.numeric(target), tolerance = tol, check.attributes = FALSE))
}

# ----------------------------------------------------------------
# 4. Locked cohort QC
# ----------------------------------------------------------------
stopifnot(
  nrow(master) == 178,
  nrow(pik_final) == 47,
  nrow(msk_sample_master) == 177,
  sum(msk_sample_master$pik3ca_mut, na.rm = TRUE) == 44,
  sum(msk_sample_master$pik3ca_cna == 2, na.rm = TRUE) == 4
)

cat("===== LOCKED COHORT QC =====\n")
cat("TCGA-CESC =", nrow(master), "\n")
cat("TCGA PIK3CA-mutant =", nrow(pik_final), "\n")
cat("MSK =", nrow(msk_sample_master), "\n")
cat("MSK PIK3CA-mutant =", sum(msk_sample_master$pik3ca_mut, na.rm = TRUE), "\n")
cat("MSK CNA=2 =", sum(msk_sample_master$pik3ca_cna == 2, na.rm = TRUE), "\n\n")

# ----------------------------------------------------------------
# 5. Standardize TCGA variables
# ----------------------------------------------------------------
pid_col    <- pick_col(pik_final, c("patient_id", "Patient_ID", "case_id"), TRUE, "patient ID")
vaf_col    <- pick_col(pik_final, c("tumor_vaf", "Tumor_VAF", "t_vaf"), TRUE, "tumor VAF")
purity_col <- pick_col(pik_final, c("purity", "Purity", "absolute_purity"), TRUE, "purity")
ploidy_col <- pick_col(pik_final, c("ploidy", "Ploidy", "absolute_ploidy"), TRUE, "ploidy")
tcn_col    <- pick_col(pik_final, c("total_cn", "Modal_Total_CN", "modal_total_cn"), TRUE, "local total CN")
wgd_col    <- pick_col(pik_final, c("wgd_positive", "WGD", "wgd"), TRUE, "WGD")
q3_col     <- pick_col(pik_final, c("chr3q_arm_call", "chr3q_gain", "gain_3q"), TRUE, "3q gain")

major_col <- pick_col(
  pik_final,
  c("major_cn", "Modal_HSCN_1", "modal_hscn_1", "hscn1"),
  required = FALSE
)
minor_col <- pick_col(
  pik_final,
  c("minor_cn", "Modal_HSCN_2", "modal_hscn_2", "hscn2"),
  required = FALSE
)
gd_col <- pick_col(
  pik_final,
  c("genome_doublings", "genome_doubling_count", "GD", "gd", "wgd_rounds"),
  required = FALSE
)

pik_final[, patient_id_std := as.character(get(pid_col))]
pik_final[, tumor_vaf_std := as.numeric(get(vaf_col))]
pik_final[, purity_std := as.numeric(get(purity_col))]
pik_final[, ploidy_std := as.numeric(get(ploidy_col))]
pik_final[, total_cn_std := as.numeric(get(tcn_col))]

# WGD -> 0/1
wraw <- pik_final[[wgd_col]]
if (is.logical(wraw)) {
  pik_final[, wgd01 := as.integer(get(wgd_col))]
} else if (is.numeric(wraw) || is.integer(wraw)) {
  pik_final[, wgd01 := fifelse(is.na(get(wgd_col)), NA_integer_, as.integer(as.numeric(get(wgd_col)) > 0))]
} else {
  tmp <- tolower(trimws(as.character(wraw)))
  pik_final[, wgd01 := fifelse(
    tmp %in% c("1", "true", "yes", "positive", "wgd+", "gd1", "gd2"), 1L,
    fifelse(tmp %in% c("0", "false", "no", "negative", "wgd-", "gd0"), 0L, NA_integer_)
  )]
}
pik_final[, wgd_group := factor(wgd01, levels = c(0, 1), labels = c("WGD−", "WGD+"))]

# 3q -> 0/1
qraw <- pik_final[[q3_col]]
if (is.logical(qraw)) {
  pik_final[, gain3q01 := as.integer(get(q3_col))]
} else if (is.numeric(qraw) || is.integer(qraw)) {
  pik_final[, gain3q01 := fifelse(is.na(get(q3_col)), NA_integer_, as.integer(as.numeric(get(q3_col)) > 0))]
} else {
  tmp <- tolower(trimws(as.character(qraw)))
  pik_final[, gain3q01 := fifelse(
    tmp %in% c("1", "true", "yes", "gain", "3q+", "positive"), 1L,
    fifelse(tmp %in% c("0", "false", "no", "neutral", "3q-", "negative"), 0L, NA_integer_)
  )]
}
pik_final[, gain3q_group := factor(gain3q01, levels = c(0, 1), labels = c("3q−", "3q+"))]

# ----------------------------------------------------------------
# 6. Recalculate dosage from source fields
# ----------------------------------------------------------------
pik_final[, vaf_fraction := tumor_vaf_std / 100]
pik_final[, mutant_cn :=
            vaf_fraction * (purity_std * total_cn_std + 2 * (1 - purity_std)) / purity_std]
pik_final[, fam := mutant_cn / total_cn_std]
pik_final[, cn_over_ploidy := total_cn_std / ploidy_std]
pik_final[, mutant_cn_over_ploidy := mutant_cn / ploidy_std]

cat("===== DOSAGE QC =====\n")
cat("Median mutant CN =", median(pik_final$mutant_cn, na.rm = TRUE), "\n")
cat("Median FAM =", median(pik_final$fam, na.rm = TRUE), "\n")
cat("FAM range =", paste(range(pik_final$fam, na.rm = TRUE), collapse = " to "), "\n\n")

# ----------------------------------------------------------------
# 7. WGD analysis
# ----------------------------------------------------------------
wgd_endpoints <- data.table(
  endpoint = c("Total CN", "CN / ploidy", "Mutant CN", "Mutant CN / ploidy", "FAM"),
  variable = c("total_cn_std", "cn_over_ploidy", "mutant_cn", "mutant_cn_over_ploidy", "fam")
)

wgd_summary <- rbindlist(lapply(seq_len(nrow(wgd_endpoints)), function(i) {
  ep <- wgd_endpoints$endpoint[i]
  vv <- wgd_endpoints$variable[i]
  x0 <- pik_final[wgd01 == 0, get(vv)]
  x1 <- pik_final[wgd01 == 1, get(vv)]
  s0 <- med_iqr(x0)
  s1 <- med_iqr(x1)
  data.table(
    endpoint = ep,
    variable = vv,
    WGD_minus_n = sum(is.finite(x0)),
    WGD_minus_median = s0["median"],
    WGD_minus_q1 = s0["q1"],
    WGD_minus_q3 = s0["q3"],
    WGD_plus_n = sum(is.finite(x1)),
    WGD_plus_median = s1["median"],
    WGD_plus_q1 = s1["q1"],
    WGD_plus_q3 = s1["q3"],
    p = safe_wilcox(pik_final[[vv]], pik_final$wgd_group)
  )
}))
wgd_summary[, fdr := p.adjust(p, method = "BH")]

print(wgd_summary)
fwrite(wgd_summary, file.path(TAB_DIR, "Table_S1_WGD_endpoints.csv"))

# ----------------------------------------------------------------
# 8. Allele-specific descriptive analysis
# ----------------------------------------------------------------
allele_summary <- NULL
allele_proximity <- NULL
major_gain_result <- NULL

if (!is.na(major_col) && !is.na(minor_col)) {
  pik_final[, major_cn_std := as.numeric(get(major_col))]
  pik_final[, minor_cn_std := as.numeric(get(minor_col))]

  allele_summary <- pik_final[, .(
    n = .N,
    major_median = median(major_cn_std, na.rm = TRUE),
    major_q1 = quantile(major_cn_std, 0.25, na.rm = TRUE),
    major_q3 = quantile(major_cn_std, 0.75, na.rm = TRUE),
    minor_median = median(minor_cn_std, na.rm = TRUE),
    minor_q1 = quantile(minor_cn_std, 0.25, na.rm = TRUE),
    minor_q3 = quantile(minor_cn_std, 0.75, na.rm = TRUE)
  ), by = wgd_group]

  pik_final[, close_major := abs(mutant_cn - major_cn_std) <= 0.5]
  pik_final[, close_minor := abs(mutant_cn - minor_cn_std) <= 0.5]
  pik_final[, allele_proximity := fifelse(
    close_major & close_minor, "Both",
    fifelse(close_major, "Major only", fifelse(close_minor, "Minor only", "Neither"))
  )]

  allele_proximity <- pik_final[, .N, by = allele_proximity]
  allele_proximity[, allele_proximity := factor(
    allele_proximity,
    levels = c("Both", "Major only", "Minor only", "Neither")
  )]
  setorder(allele_proximity, allele_proximity)

  fwrite(allele_summary, file.path(TAB_DIR, "Table_S4_allele_specific_summary.csv"))
  fwrite(allele_proximity, file.path(TAB_DIR, "Table_S4_mutantCN_allele_proximity_descriptive.csv"))

  # Formal GD0/GD1 retained-major comparison only when true GD count exists.
  if (!is.na(gd_col)) {
    pik_final[, gd_count := as.integer(round(as.numeric(get(gd_col))))]
    gd01 <- pik_final[gd_count %in% c(0L, 1L)]
    gd01[, expected_retained_major := fifelse(gd_count == 0L, 1, 2)]
    gd01[, major_gain_above_expected := major_cn_std > expected_retained_major]
    gt <- table(GD = gd01$gd_count, MajorGain = gd01$major_gain_above_expected)
    if (length(dim(gt)) == 2 && all(dim(gt) == c(2, 2))) {
      ft <- fisher.test(gt)
      major_gain_result <- data.table(
        p = ft$p.value,
        odds_ratio = unname(ft$estimate),
        ci_low = ft$conf.int[1],
        ci_high = ft$conf.int[2]
      )
      fwrite(major_gain_result, file.path(TAB_DIR, "Table_S4_WGD_aware_major_gain.csv"))
    }
  }
}

# ----------------------------------------------------------------
# 9. 3q analysis
# ----------------------------------------------------------------
pik_3q <- pik_final[!is.na(gain3q01)]
stopifnot(nrow(pik_3q) == 40, sum(pik_3q$gain3q01 == 0) == 14, sum(pik_3q$gain3q01 == 1) == 26)

q3_endpoints <- data.table(
  endpoint = c("Total CN", "CN / ploidy", "Mutant CN", "Mutant CN / ploidy", "FAM"),
  variable = c("total_cn_std", "cn_over_ploidy", "mutant_cn", "mutant_cn_over_ploidy", "fam")
)

gain3q_summary <- rbindlist(lapply(seq_len(nrow(q3_endpoints)), function(i) {
  ep <- q3_endpoints$endpoint[i]
  vv <- q3_endpoints$variable[i]
  x0 <- pik_3q[gain3q01 == 0, get(vv)]
  x1 <- pik_3q[gain3q01 == 1, get(vv)]
  s0 <- med_iqr(x0)
  s1 <- med_iqr(x1)
  data.table(
    endpoint = ep,
    variable = vv,
    q3_minus_n = sum(is.finite(x0)),
    q3_minus_median = s0["median"],
    q3_minus_q1 = s0["q1"],
    q3_minus_q3 = s0["q3"],
    q3_plus_n = sum(is.finite(x1)),
    q3_plus_median = s1["median"],
    q3_plus_q1 = s1["q1"],
    q3_plus_q3 = s1["q3"],
    wilcox_p = safe_wilcox(pik_3q[[vv]], pik_3q$gain3q_group)
  )
}))
fwrite(gain3q_summary, file.path(TAB_DIR, "Table_3q_raw_endpoints.csv"))

wgd_ploidy_spearman <- suppressWarnings(
  cor.test(pik_3q$wgd01, pik_3q$ploidy_std, method = "spearman", exact = FALSE)
)

# Locked models
fit_total_base <- lm(total_cn_std ~ ploidy_std + wgd01, data = pik_3q)
fit_total_full <- lm(total_cn_std ~ ploidy_std + wgd01 + gain3q01, data = pik_3q)
fit_total_ploidy3q <- lm(total_cn_std ~ ploidy_std + gain3q01, data = pik_3q)
fit_total_wgd3q <- lm(total_cn_std ~ wgd01 + gain3q01, data = pik_3q)

fit_mut_base <- lm(mutant_cn ~ ploidy_std + wgd01, data = pik_3q)
fit_mut_full <- lm(mutant_cn ~ ploidy_std + wgd01 + gain3q01, data = pik_3q)
fit_mutnorm <- lm(mutant_cn_over_ploidy ~ wgd01 + gain3q01, data = pik_3q)
fit_fam_base <- lm(fam ~ ploidy_std + wgd01, data = pik_3q)
fit_fam_full <- lm(fam ~ ploidy_std + wgd01 + gain3q01, data = pik_3q)

results_3q <- rbindlist(list(
  extract_term(fit_total_full, "gain3q01", "Total CN", "ploidy + WGD + 3q"),
  extract_term(fit_mut_full, "gain3q01", "Mutant CN", "ploidy + WGD + 3q"),
  extract_term(fit_mutnorm, "gain3q01", "Mutant CN / ploidy", "WGD + 3q"),
  extract_term(fit_fam_full, "gain3q01", "FAM", "ploidy + WGD + 3q")
))
results_3q[, fdr := p.adjust(p, method = "BH")]
fwrite(results_3q, file.path(TAB_DIR, "Table_3q_adjusted_endpoints.csv"))

model_fit_summary <- data.table(
  model = c(
    "Total: ploidy + WGD",
    "Total: ploidy + WGD + 3q",
    "Total: ploidy + 3q",
    "Total: WGD + 3q",
    "Mutant: ploidy + WGD",
    "Mutant: ploidy + WGD + 3q",
    "Mutant/ploidy: WGD + 3q",
    "FAM: ploidy + WGD",
    "FAM: ploidy + WGD + 3q"
  ),
  R2 = c(
    summary(fit_total_base)$r.squared,
    summary(fit_total_full)$r.squared,
    summary(fit_total_ploidy3q)$r.squared,
    summary(fit_total_wgd3q)$r.squared,
    summary(fit_mut_base)$r.squared,
    summary(fit_mut_full)$r.squared,
    summary(fit_mutnorm)$r.squared,
    summary(fit_fam_base)$r.squared,
    summary(fit_fam_full)$r.squared
  ),
  adj_R2 = c(
    summary(fit_total_base)$adj.r.squared,
    summary(fit_total_full)$adj.r.squared,
    summary(fit_total_ploidy3q)$adj.r.squared,
    summary(fit_total_wgd3q)$adj.r.squared,
    summary(fit_mut_base)$adj.r.squared,
    summary(fit_mut_full)$adj.r.squared,
    summary(fit_mutnorm)$adj.r.squared,
    summary(fit_fam_base)$adj.r.squared,
    summary(fit_fam_full)$adj.r.squared
  ),
  AIC = c(
    AIC(fit_total_base), AIC(fit_total_full), AIC(fit_total_ploidy3q), AIC(fit_total_wgd3q),
    AIC(fit_mut_base), AIC(fit_mut_full), AIC(fit_mutnorm), AIC(fit_fam_base), AIC(fit_fam_full)
  )
)
fwrite(model_fit_summary, file.path(TAB_DIR, "Table_S2_model_fit_summary.csv"))

nested_tests <- data.table(
  comparison = c("Total CN: +3q", "Mutant CN: +3q", "FAM: +3q"),
  p = c(
    anova(fit_total_base, fit_total_full)$`Pr(>F)`[2],
    anova(fit_mut_base, fit_mut_full)$`Pr(>F)`[2],
    anova(fit_fam_base, fit_fam_full)$`Pr(>F)`[2]
  )
)
fwrite(nested_tests, file.path(TAB_DIR, "Table_S2_nested_model_tests.csv"))

# ----------------------------------------------------------------
# 10. Influence / sensitivity analysis
# ----------------------------------------------------------------
cook <- cooks.distance(fit_total_full)
cook_cutoff <- 4 / nobs(fit_total_full)
model_rows <- as.integer(rownames(model.frame(fit_total_full)))

influence_table <- data.table(
  model_row = seq_along(cook),
  original_row = model_rows,
  cooks_distance = as.numeric(cook),
  influential = as.numeric(cook) > cook_cutoff
)
influence_table[, patient_id := pik_3q$patient_id_std[original_row]]
setorder(influence_table, -cooks_distance)
fwrite(influence_table, file.path(TAB_DIR, "Table_S3_influence_diagnostics.csv"))

influential_ids <- influence_table[influential == TRUE, patient_id]
sensitivity_result <- NULL
if (length(influential_ids) > 0) {
  sens_dat <- pik_3q[!patient_id_std %in% influential_ids]
  fit_sens <- lm(total_cn_std ~ ploidy_std + wgd01 + gain3q01, data = sens_dat)
  sensitivity_result <- extract_term(
    fit_sens, "gain3q01", "Total CN", "Sensitivity excluding Cook > 4/n"
  )
  fwrite(sensitivity_result, file.path(TAB_DIR, "Table_S3_sensitivity_excluding_influential.csv"))
}

# ----------------------------------------------------------------
# 11. MSK external validation
# ----------------------------------------------------------------
tcga_n <- nrow(master)
tcga_mut_n <- nrow(pik_final)
msk_n <- nrow(msk_sample_master)
msk_mut_n <- sum(msk_sample_master$pik3ca_mut, na.rm = TRUE)

prev_tab <- matrix(
  c(tcga_mut_n, tcga_n - tcga_mut_n, msk_mut_n, msk_n - msk_mut_n),
  nrow = 2, byrow = TRUE,
  dimnames = list(Cohort = c("TCGA-CESC", "MSK"), Status = c("Mutant", "Wild-type"))
)
prev_fisher <- fisher.test(prev_tab)

msk_prevalence <- data.table(
  cohort = c("TCGA-CESC", "MSK"),
  mutant_n = c(tcga_mut_n, msk_mut_n),
  total_n = c(tcga_n, msk_n)
)
msk_prevalence[, prevalence := mutant_n / total_n]
fwrite(msk_prevalence, file.path(TAB_DIR, "Table_MSK_TCGA_prevalence.csv"))

# MSK mutation x CNA=2
msk_mut_cna_tab <- table(
  Mutation = msk_sample_master$pik3ca_mut,
  CNA2 = msk_sample_master$pik3ca_cna == 2
)
msk_cna_fisher <- fisher.test(msk_mut_cna_tab)

# Four states
if ("pik3ca_state" %in% names(msk_sample_master)) {
  msk_state <- msk_sample_master[, .N, by = pik3ca_state]
  setnames(msk_state, "pik3ca_state", "state")
} else {
  tmp <- copy(msk_sample_master)
  tmp[, state := fifelse(
    pik3ca_mut & pik3ca_cna == 2, "Mutation + CNA=2",
    fifelse(pik3ca_mut, "Mutation only", fifelse(pik3ca_cna == 2, "CNA=2 only", "Neither"))
  )]
  msk_state <- tmp[, .N, by = state]
}
state_levels <- c("Neither", "Mutation only", "Mutation + CNA=2", "CNA=2 only")
msk_state[, state := factor(state, levels = state_levels)]
msk_state[, pct := 100 * N / sum(N)]
setorder(msk_state, state)
fwrite(msk_state, file.path(TAB_DIR, "Table_MSK_PIK3CA_states.csv"))

# Mutation spectrum
protein_col <- pick_col(
  msk_pik_mut,
  c("protein_change", "HGVSp_Short", "hgvsp_short", "HGVSp", "hgvsp", "Protein_Change", "proteinChange"),
  required = FALSE
)
mutation_spectrum <- NULL
if (!is.na(protein_col)) {
  tmp <- copy(msk_pik_mut)
  tmp[, protein_change_std := sub("^p\\.", "", as.character(get(protein_col)))]
  events <- tmp[, .N, by = protein_change_std]
  events[, mutation_group := fifelse(
    protein_change_std %in% c("E545K", "E542K", "E453K", "E726K"),
    protein_change_std, "Other"
  )]
  mutation_spectrum <- events[, .(N = sum(N)), by = mutation_group]
  mutation_spectrum[, pct := 100 * N / sum(N)]
  mutation_spectrum[, mutation_group := factor(
    mutation_group,
    levels = c("E545K", "E542K", "E453K", "E726K", "Other")
  )]
  setorder(mutation_spectrum, mutation_group)
  fwrite(mutation_spectrum, file.path(TAB_DIR, "Table_S6_MSK_mutation_spectrum.csv"))
}

# Grouped histology
if ("histology_group" %in% names(msk_sample_master)) {
  msk_hist_group <- msk_sample_master[, .(
    N = .N,
    mutant_n = sum(pik3ca_mut, na.rm = TRUE)
  ), by = histology_group]
  msk_hist_group[, prevalence := mutant_n / N]
  msk_hist_group[, histology_group := factor(
    histology_group,
    levels = c("SCC", "Endocervical adenocarcinoma", "Adenosquamous", "Gastric-type", "Other")
  )]
  setorder(msk_hist_group, histology_group)
  fwrite(msk_hist_group, file.path(TAB_DIR, "Table_S6_MSK_histology_grouped.csv"))
} else {
  msk_hist_group <- NULL
}

# ----------------------------------------------------------------
# 12. HPV exploratory analyses
# ----------------------------------------------------------------
master_pid_col <- pick_col(master, c("patient_id", "Patient_ID", "case_id"), TRUE, "master patient ID")
hpv_pid_col <- pick_col(hpv_patient, c("patient_id", "Patient_ID", "case_id", "hpvsrc_sample_id"), TRUE, "HPV patient ID")

master_link <- copy(master)
master_link[, patient_id_std := as.character(get(master_pid_col))]
hpv_sub <- copy(hpv_patient)
hpv_sub[, patient_id_std := as.character(get(hpv_pid_col))]
hpv_sub[, pik3ca_mut := patient_id_std %in% pik_final$patient_id_std]

mwgd <- pick_col(master_link, c("wgd_positive", "WGD", "wgd"), required = FALSE)
mq3 <- pick_col(master_link, c("chr3q_arm_call", "chr3q_gain", "gain_3q"), required = FALSE)
merge_dt <- unique(master_link[, .(patient_id_std)])
if (!is.na(mwgd)) {
  x <- unique(master_link[, .(patient_id_std, wgd_for_hpv = get(mwgd))])
  merge_dt <- merge(merge_dt, x, by = "patient_id_std", all.x = TRUE)
}
if (!is.na(mq3)) {
  x <- unique(master_link[, .(patient_id_std, gain3q_for_hpv = get(mq3))])
  merge_dt <- merge(merge_dt, x, by = "patient_id_std", all.x = TRUE)
}
hpv_sub <- merge(hpv_sub, merge_dt, by = "patient_id_std", all.x = TRUE)

hpv_type_col <- pick_col(
  hpv_sub,
  c("hpv_type_standard", "hpvsrc_final_type_hierarchical", "hpvsrc_hpv_types"),
  required = FALSE
)
hpv_genotype_summary <- NULL
hpv_genotype_test <- NULL
if (!is.na(hpv_type_col)) {
  hpv_sub[, hpv_type_raw := as.character(get(hpv_type_col))]
  hpv_sub[, hpv_group := fifelse(
    grepl("HPV16", hpv_type_raw, ignore.case = TRUE), "HPV16",
    fifelse(grepl("HPV18", hpv_type_raw, ignore.case = TRUE), "HPV18", "Other")
  )]
  hpv_genotype_summary <- hpv_sub[, .(
    N = .N,
    PIK3CA_mut_n = sum(pik3ca_mut, na.rm = TRUE),
    PIK3CA_mut_pct = 100 * mean(pik3ca_mut, na.rm = TRUE)
  ), by = hpv_group]
  hpv_genotype_test <- fisher.test(table(hpv_sub$hpv_group, hpv_sub$pik3ca_mut))
  fwrite(hpv_genotype_summary, file.path(TAB_DIR, "Table_S5_HPV_genotype.csv"))
}

# Secondary integration: use explicit aliases first, then regex fallback.
integration_col <- pick_col(
  hpv_sub,
  c(
    "hpvsrc_secondary_integration", "hpvsrc_secondary_integration_status",
    "secondary_integration", "secondary_integration_status", "hpv_integration_secondary"
  ),
  required = FALSE
)
if (is.na(integration_col)) {
  z <- grep("secondary.*integration|integration.*secondary", names(hpv_sub), value = TRUE, ignore.case = TRUE)
  if (length(z)) integration_col <- z[1]
}

hpv_integration_summary <- NULL
hpv_integration_tests <- NULL
if (!is.na(integration_col)) {
  hpv_sub[, integration_raw := trimws(as.character(get(integration_col)))]
  hpv_int <- hpv_sub[tolower(integration_raw) %in% c("yes", "no")]
  hpv_int[, integration_yes := tolower(integration_raw) == "yes"]

  hpv_integration_summary <- hpv_int[, .(
    N = .N,
    PIK3CA_mut_n = sum(pik3ca_mut),
    PIK3CA_mut_pct = 100 * mean(pik3ca_mut)
  ), by = integration_yes]

  ft1 <- fisher.test(table(hpv_int$integration_yes, hpv_int$pik3ca_mut))
  hpv_integration_tests <- data.table(
    endpoint = "PIK3CA mutation",
    p = ft1$p.value,
    odds_ratio = unname(ft1$estimate),
    ci_low = ft1$conf.int[1],
    ci_high = ft1$conf.int[2]
  )

  if ("wgd_for_hpv" %in% names(hpv_int)) {
    z <- hpv_int[!is.na(wgd_for_hpv)]
    z[, wgd_yes := as.numeric(wgd_for_hpv) > 0]
    ft <- fisher.test(table(z$integration_yes, z$wgd_yes))
    hpv_integration_tests <- rbind(hpv_integration_tests, data.table(
      endpoint = "WGD", p = ft$p.value, odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2]
    ), fill = TRUE)
  }

  if ("gain3q_for_hpv" %in% names(hpv_int)) {
    z <- hpv_int[!is.na(gain3q_for_hpv)]
    z[, gain3q_yes := as.numeric(gain3q_for_hpv) > 0]
    ft <- fisher.test(table(z$integration_yes, z$gain3q_yes))
    hpv_integration_tests <- rbind(hpv_integration_tests, data.table(
      endpoint = "3q gain", p = ft$p.value, odds_ratio = unname(ft$estimate),
      ci_low = ft$conf.int[1], ci_high = ft$conf.int[2]
    ), fill = TRUE)
  }

  fwrite(hpv_integration_summary, file.path(TAB_DIR, "Table_S5_HPV_integration_summary.csv"))
  fwrite(hpv_integration_tests, file.path(TAB_DIR, "Table_S5_HPV_integration_tests.csv"))
}

# ----------------------------------------------------------------
# 13. Main Table 1 export
# ----------------------------------------------------------------
main_table1 <- rbindlist(list(
  data.table(
    domain = "WGD",
    endpoint = wgd_summary$endpoint,
    estimate_or_median = paste0(
      sprintf("%.3f", wgd_summary$WGD_minus_median), " vs ",
      sprintf("%.3f", wgd_summary$WGD_plus_median)
    ),
    ci = NA_character_,
    p = wgd_summary$p,
    fdr = wgd_summary$fdr
  ),
  data.table(
    domain = "3q adjusted",
    endpoint = results_3q$endpoint,
    estimate_or_median = sprintf("β=%.3f", results_3q$beta),
    ci = sprintf("%.3f to %.3f", results_3q$ci_low, results_3q$ci_high),
    p = results_3q$p,
    fdr = results_3q$fdr
  )
), fill = TRUE)
fwrite(main_table1, file.path(TAB_DIR, "Table1_core_results.csv"))

# ----------------------------------------------------------------
# 14. Locked QC assertions
# ----------------------------------------------------------------
# These tolerances are intentionally tight enough to detect accidental changes
# while allowing harmless floating-point differences.
stopifnot(
  near_equal(median(pik_final$mutant_cn, na.rm = TRUE), 1.3368, 5e-4),
  near_equal(median(pik_final$fam, na.rm = TRUE), 0.4229, 5e-4),
  sum(pik_final$wgd01 == 0, na.rm = TRUE) == 35,
  sum(pik_final$wgd01 == 1, na.rm = TRUE) == 12
)

expected_wgd_p <- c(0.0006750165, 0.8643675965, 0.0009468488, 0.3602408128, 0.1680616635)
stopifnot(max(abs(wgd_summary$p - expected_wgd_p), na.rm = TRUE) < 1e-7)

expected_3q_beta <- c(1.36204621, 0.38579593, 0.20192228, -0.07465301)
expected_3q_p <- c(1.002918e-05, 7.495209e-02, 5.008082e-02, 1.121375e-01)
stopifnot(
  max(abs(results_3q$beta - expected_3q_beta), na.rm = TRUE) < 1e-6,
  max(abs(results_3q$p - expected_3q_p), na.rm = TRUE) < 1e-6,
  abs(prev_fisher$p.value - 0.8080991) < 1e-6,
  abs(msk_cna_fisher$p.value - 0.04800117) < 1e-6
)

# ----------------------------------------------------------------
# 15. Figure style — final higher-saturation palette
# ----------------------------------------------------------------
COL_BLUE  <- "#2F78AD"
COL_BLUE2 <- "#4E9BC3"
COL_TEAL  <- "#3C9B91"
COL_ROSE  <- "#DF4F73"
COL_GREY  <- "#B9C3CC"
COL_GREEN <- "#55A98A"
COL_DARK  <- "#2B2B2B"
COL_LIGHT <- "#F3F5F7"
COL_BLUE_FILL  <- "#E7F0F7"
COL_GREEN_FILL <- "#E4F2ED"
COL_ROSE_FILL  <- "#FBE5EB"
COL_GREY_FILL  <- "#EEF1F4"

base_theme <- theme_classic(base_size = 11) +
  theme(
    text = element_text(color = COL_DARK),
    axis.text = element_text(color = COL_DARK),
    axis.title = element_text(color = COL_DARK),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, margin = margin(b = 6)),
    legend.title = element_blank(),
    legend.position = "top",
    plot.margin = margin(8, 10, 8, 8)
  )

# ----------------------------------------------------------------
# 16. Figure 1 — study design
# Solid arrowheads are drawn manually to avoid hollow-arrow rendering.
# Arrow lengths are shortened slightly so heads stop cleanly before box interiors.
# Panel B label is deliberately outside the first box.
# ----------------------------------------------------------------
fig1_path_pdf <- file.path(FIG_DIR, "Figure1_study_design.pdf")
fig1_path_png <- file.path(FIG_DIR, "Figure1_study_design.png")

make_figure1 <- function(filename, png = FALSE) {
  if (png) {
    png(filename, width = 3300, height = 1755, res = 300)
  } else {
    pdf(filename, width = 11, height = 5.85, useDingbats = FALSE)
  }
  grid.newpage()

  draw_box <- function(cx, cy, w, h, label, fill, border = COL_DARK, fs = 10) {
    grid.roundrect(
      x = unit(cx, "npc"), y = unit(cy, "npc"),
      width = unit(w, "npc"), height = unit(h, "npc"),
      r = unit(0.02, "npc"),
      gp = gpar(fill = fill, col = border, lwd = 1.5)
    )
    grid.text(
      label, x = unit(cx, "npc"), y = unit(cy, "npc"),
      gp = gpar(col = COL_DARK, fontsize = fs, fontface = "plain")
    )
  }

  draw_arrow <- function(x1, y1, x2, y2, col = "#565656") {
    ah <- 0.008
    aw <- 0.007
    if (abs(y2 - y1) < 1e-9) {
      sgn <- sign(x2 - x1)
      grid.lines(
        x = unit(c(x1, x2 - sgn * ah), "npc"),
        y = unit(c(y1, y2), "npc"),
        gp = gpar(col = col, lwd = 1.6)
      )
      grid.polygon(
        x = unit(c(x2, x2 - sgn * ah, x2 - sgn * ah), "npc"),
        y = unit(c(y2, y2 + aw, y2 - aw), "npc"),
        gp = gpar(fill = col, col = col)
      )
    } else {
      sgn <- sign(y2 - y1)
      grid.lines(
        x = unit(c(x1, x2), "npc"),
        y = unit(c(y1, y2 - sgn * ah), "npc"),
        gp = gpar(col = col, lwd = 1.6)
      )
      grid.polygon(
        x = unit(c(x2, x2 + aw, x2 - aw), "npc"),
        y = unit(c(y2, y2 - sgn * ah, y2 - sgn * ah), "npc"),
        gp = gpar(fill = col, col = col)
      )
    }
  }

  # Fixed four-column geometry. All boxes in rows 2-4 have identical width.
  cx <- c(0.15, 0.40, 0.65, 0.90)
  bw <- 0.20
  top_y <- 0.84
  mid_y <- 0.57
  low_y <- 0.31
  b_y <- 0.095
  gap <- 0.008

  grid.text("A", x = unit(0.028, "npc"), y = unit(0.955, "npc"),
            gp = gpar(fontsize = 15, fontface = "bold"))

  draw_box(cx[1], top_y, bw, 0.14, "TCGA-CESC\n178 tumors", COL_LIGHT)
  draw_box(cx[2], top_y, bw, 0.14, "PIK3CA-mutant\n47 tumors", COL_BLUE_FILL, COL_BLUE)
  draw_box(cx[3], top_y, bw, 0.14, "WGD status\n35 WGD− / 12 WGD+", COL_ROSE_FILL, COL_ROSE)
  draw_box(cx[4], top_y, bw, 0.14, "3q evaluable\n40 tumors", COL_GREEN_FILL, COL_GREEN)
  for (i in 1:3) {
    draw_arrow(cx[i] + bw/2 + gap, top_y, cx[i+1] - bw/2 - gap, top_y)
  }

  draw_box(cx[1], mid_y, bw, 0.14, "Mutation VAF\n(% → fraction)", COL_BLUE_FILL, COL_BLUE)
  draw_box(cx[2], mid_y, bw, 0.14, "ABSOLUTE\npurity + ploidy", COL_BLUE_FILL, COL_BLUE)
  draw_box(cx[3], mid_y, bw, 0.14, "Local PIK3CA\ntotal CN", COL_ROSE_FILL, COL_ROSE)
  draw_box(cx[4], mid_y, bw, 0.14, "Allele-specific\nmajor/minor CN", COL_GREY_FILL, COL_GREY)

  draw_box(cx[1], low_y, bw, 0.14, "Mutant CN\nm = VAF[ρC+2(1−ρ)]/ρ", COL_BLUE_FILL, COL_BLUE)
  draw_box(cx[2], low_y, bw, 0.14, "Ploidy-normalized\nCN and mutant CN", COL_GREEN_FILL, COL_GREEN)
  draw_box(cx[3], low_y, bw, 0.14, "FAM\nmutant CN / total CN", COL_ROSE_FILL, COL_ROSE)
  draw_box(cx[4], low_y, bw, 0.14, "WGD-aware\ncontext", COL_GREY_FILL, COL_GREY)

  # Equal-length vertical connectors.
  for (xx in cx) {
    draw_arrow(xx, mid_y - 0.07 - gap, xx, low_y + 0.07 + gap)
  }
  # Equal-length horizontal connectors.
  for (i in 1:3) {
    draw_arrow(cx[i] + bw/2 + gap, low_y, cx[i+1] - bw/2 - gap, low_y)
  }

  grid.text("B", x = unit(0.020, "npc"), y = unit(b_y, "npc"),
            gp = gpar(fontsize = 15, fontface = "bold"))
  draw_box(cx[1], b_y, bw, 0.11, "Genome-wide\nWGD / ploidy", COL_BLUE_FILL, COL_BLUE)
  draw_box(cx[2], b_y, bw, 0.11, "Regional\nchromosome 3q gain", COL_GREEN_FILL, COL_GREEN)
  draw_box(cx[3], b_y, bw, 0.11, "Local locus\nPIK3CA total CN", COL_ROSE_FILL, COL_ROSE)
  draw_box(cx[4], b_y, bw, 0.11, "Mutant-specific\nmutant CN / FAM", COL_GREY_FILL, COL_GREY)
  for (i in 1:3) {
    draw_arrow(cx[i] + bw/2 + gap, b_y, cx[i+1] - bw/2 - gap, b_y)
  }

  dev.off()
}
make_figure1(fig1_path_pdf, png = FALSE)
make_figure1(fig1_path_png, png = TRUE)

# ----------------------------------------------------------------
# 17. Figure 2 — WGD and PIK3CA dosage
# Statistics are placed in subtitles to prevent title/text overlap.
# P values <0.001 are displayed conventionally as P<0.001.
# ----------------------------------------------------------------
fig2_dat <- copy(wgd_summary)
fig2_dat[, endpoint := factor(
  endpoint,
  levels = c("Total CN", "Mutant CN", "CN / ploidy", "Mutant CN / ploidy", "FAM")
)]

make_wgd_panel <- function(ep) {
  z <- fig2_dat[endpoint == ep]
  dd <- data.table(
    group = factor(c("WGD−", "WGD+"), levels = c("WGD−", "WGD+")),
    median = c(z$WGD_minus_median, z$WGD_plus_median)
  )
  ggplot(dd, aes(x = group, y = median, fill = group)) +
    geom_col(width = 0.62, color = COL_DARK, linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", median)), vjust = -0.40, size = 3.0) +
    scale_fill_manual(values = c("WGD−" = COL_BLUE2, "WGD+" = COL_ROSE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      x = NULL,
      y = if (ep == "Total CN") "Median" else NULL,
      title = ep,
      subtitle = paste0(fmt_p(z$p, "P"), "\n", fmt_p(z$fdr, "FDR"))
    ) +
    base_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8.8, margin = margin(b = 8))
    )
}

p2 <- wrap_plots(lapply(levels(fig2_dat$endpoint), make_wgd_panel), nrow = 1) +
  plot_annotation(
    title = "Whole-genome doubling increases absolute, but not relative, PIK3CA dosage",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
  )

ggsave(file.path(FIG_DIR, "Figure2_WGD_PIK3CA_dosage.pdf"), p2, width = 12.3, height = 3.8)
ggsave(file.path(FIG_DIR, "Figure2_WGD_PIK3CA_dosage.tiff"), p2, width = 12.3, height = 3.8, dpi = 600, compression = "lzw")

# ----------------------------------------------------------------
# 18. Figure 3 — chromosome 3q gain and PIK3CA dosage
# Dot/interval layout avoids overlapping labels and uses standard P/FDR notation.
# ----------------------------------------------------------------
raw_plot_dat <- gain3q_summary[
  endpoint %in% c("Total CN", "Mutant CN", "Mutant CN / ploidy", "FAM")
]
raw_plot_dat[, endpoint := factor(
  endpoint,
  levels = rev(c("Total CN", "Mutant CN", "Mutant CN / ploidy", "FAM"))
)]
raw_plot_dat[, p_label := vapply(wilcox_p, fmt_p, character(1), prefix = "P")]
raw_plot_dat[, label_x := pmax(q3_minus_median, q3_plus_median) + fifelse(
  pmax(q3_minus_median, q3_plus_median) > 1, 0.12, 0.07
)]

p3a <- ggplot(raw_plot_dat, aes(y = endpoint)) +
  geom_segment(
    aes(x = q3_minus_median, xend = q3_plus_median, yend = endpoint),
    color = "#A9B0B5", linewidth = 0.8
  ) +
  geom_point(aes(x = q3_minus_median), color = COL_BLUE2, size = 3.1, shape = 21, fill = COL_BLUE2, stroke = 0.45) +
  geom_point(aes(x = q3_plus_median), color = COL_ROSE, size = 3.1, shape = 21, fill = COL_ROSE, stroke = 0.45) +
  geom_text(aes(x = label_x, label = p_label), hjust = 0, size = 3.0) +
  coord_cartesian(xlim = c(0.20, 3.70), clip = "off") +
  labs(x = "Median endpoint value", y = NULL, title = "A  Raw 3q− versus 3q+ summaries") +
  base_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 11.5),
    plot.margin = margin(8, 20, 8, 8)
  )

coef_dt <- copy(results_3q)
coef_dt[, endpoint := factor(endpoint, levels = rev(c("Total CN", "Mutant CN", "Mutant CN / ploidy", "FAM")))]
coef_dt[, fdr_label := vapply(fdr, fmt_p, character(1), prefix = "FDR")]
coef_dt[, label_x := ci_high + 0.08]

p3b <- ggplot(coef_dt, aes(y = endpoint, x = beta)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#8F969B", linewidth = 0.5) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = endpoint), color = COL_DARK, linewidth = 0.9) +
  geom_point(color = COL_ROSE, fill = COL_ROSE, shape = 21, stroke = 0.45, size = 3.1) +
  geom_text(aes(x = label_x, label = fdr_label), hjust = 0, size = 3.0) +
  coord_cartesian(xlim = c(-0.28, 2.30), clip = "off") +
  labs(x = "Adjusted 3q coefficient (95% CI)", y = NULL, title = "B  Endpoint-specific adjusted models") +
  base_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 11.5),
    plot.margin = margin(8, 24, 8, 8)
  )

fig3 <- p3a | p3b + plot_layout(widths = c(1.02, 1))
ggsave(file.path(FIG_DIR, "Figure3_3q_PIK3CA_dosage.pdf"), fig3, width = 12, height = 5.32)
ggsave(file.path(FIG_DIR, "Figure3_3q_PIK3CA_dosage.tiff"), fig3, width = 12, height = 5.32, dpi = 600, compression = "lzw")

# ----------------------------------------------------------------
# 19. Figure 4 — external genomic validation
# ----------------------------------------------------------------
p4a <- ggplot(msk_prevalence, aes(x = cohort, y = prevalence, fill = cohort)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d/%d\n(%.1f%%)", mutant_n, total_n, 100 * prevalence)), vjust = -0.35, size = 3.6) +
  scale_fill_manual(values = c("TCGA-CESC" = COL_BLUE, "MSK" = COL_TEAL)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.36), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "PIK3CA mutation prevalence", subtitle = paste0("Fisher P=", fmt_p(prev_fisher$p.value))) +
  base_theme + theme(legend.position = "none")

if (!is.null(mutation_spectrum)) {
  mut_cols <- c("E545K" = COL_ROSE, "E542K" = COL_BLUE, "E453K" = COL_BLUE2, "E726K" = COL_TEAL, "Other" = COL_GREY)
  p4b <- ggplot(mutation_spectrum, aes(x = mutation_group, y = N, fill = mutation_group)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = paste0(N, "\n(", sprintf("%.0f%%", pct), ")")), vjust = -0.3, size = 3.3) +
    scale_fill_manual(values = mut_cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Mutation events") +
    base_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  fig4 <- p4a | p4b + plot_annotation(title = "External genomic validation in MSK", tag_levels = "A")
} else {
  fig4 <- p4a
}

ggsave(file.path(FIG_DIR, "Figure4_external_validation.pdf"), fig4, width = 10, height = 4.8)
ggsave(file.path(FIG_DIR, "Figure4_external_validation.tiff"), fig4, width = 10, height = 4.8, dpi = 600, compression = "lzw")

# ----------------------------------------------------------------
# 20. Figure 5 — MSK biological context
# ----------------------------------------------------------------
if (!is.null(msk_hist_group)) {
  hist_cols <- c(
    "SCC" = COL_ROSE,
    "Endocervical adenocarcinoma" = COL_BLUE,
    "Adenosquamous" = COL_BLUE2,
    "Gastric-type" = COL_GREEN,
    "Other" = COL_GREY
  )
  p5a <- ggplot(msk_hist_group, aes(x = histology_group, y = prevalence, fill = histology_group)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = sprintf("%d/%d\n(%.1f%%)", mutant_n, N, 100 * prevalence)), vjust = -0.30, size = 3.1) +
    scale_fill_manual(values = hist_cols, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "PIK3CA mutation prevalence") +
    base_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

  state_cols <- c("Neither" = COL_GREY, "Mutation only" = COL_BLUE, "Mutation + CNA=2" = COL_ROSE, "CNA=2 only" = COL_GREEN)
  p5b <- ggplot(msk_state, aes(x = "MSK", y = pct, fill = state)) +
    geom_col(width = 0.55) + coord_flip() +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    geom_text(
      data = msk_state[state %in% c("Neither", "Mutation only")],
      aes(label = sprintf("%d (%.1f%%)", N, pct)),
      position = position_stack(vjust = 0.5), size = 3.4
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "MSK cohort composition") +
    base_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")

  fig5 <- p5a / p5b + plot_layout(heights = c(1.25, 0.75)) +
    plot_annotation(title = "MSK biological context of PIK3CA alterations", tag_levels = "A")
  ggsave(file.path(FIG_DIR, "Figure5_MSK_context.pdf"), fig5, width = 10, height = 7.8)
  ggsave(file.path(FIG_DIR, "Figure5_MSK_context.tiff"), fig5, width = 10, height = 7.8, dpi = 600, compression = "lzw")
}

# ----------------------------------------------------------------
# 21. Supplementary Figure S1 — genomic context and model fit
# ----------------------------------------------------------------
pS1a <- ggplot(pik_3q, aes(x = ploidy_std, y = total_cn_std, fill = gain3q_group, shape = gain3q_group)) +
  geom_point(size = 3.0, alpha = 0.8, color = "white", stroke = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = COL_DARK) +
  facet_wrap(~wgd_group) +
  scale_fill_manual(values = c("3q−" = COL_GREY, "3q+" = COL_ROSE)) +
  labs(x = "Tumor ploidy", y = "PIK3CA total CN", title = "A  Genomic context") +
  base_theme

fit_plot <- model_fit_summary[1:4]
fit_plot[, model_short := factor(
  c("Ploidy + WGD", "Ploidy + WGD + 3q", "Ploidy + 3q", "WGD + 3q"),
  levels = c("Ploidy + WGD", "Ploidy + WGD + 3q", "Ploidy + 3q", "WGD + 3q")
)]
pS1b <- ggplot(fit_plot, aes(x = model_short, y = R2, fill = model_short)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = paste0("R²=", sprintf("%.3f", R2), "\nAIC=", sprintf("%.1f", AIC))), vjust = -0.25, size = 3.0) +
  scale_fill_manual(values = c(COL_GREY, COL_ROSE, COL_BLUE, COL_GREEN)) +
  scale_y_continuous(limits = c(0, 0.72), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "R²", title = "B  Total-CN model fit") +
  base_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

figS1 <- pS1a | pS1b

ggsave(file.path(FIG_DIR, "SupplementaryFigureS1_genomic_context_model_fit.pdf"), figS1, width = 11, height = 4.8)
ggsave(file.path(FIG_DIR, "SupplementaryFigureS1_genomic_context_model_fit.tiff"), figS1, width = 11, height = 4.8, dpi = 600, compression = "lzw")

# ----------------------------------------------------------------
# 22. Supplementary Figure S2 — allele-specific context
# ----------------------------------------------------------------
if (!is.na(major_col) && !is.na(minor_col)) {
  allele_long <- rbindlist(list(
    pik_final[, .(patient_id_std, wgd_group, allele = "Major CN", value = major_cn_std)],
    pik_final[, .(patient_id_std, wgd_group, allele = "Minor CN", value = minor_cn_std)]
  ))
  pS2a <- ggplot(allele_long, aes(x = wgd_group, y = value, fill = wgd_group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA) +
    geom_jitter(width = 0.10, size = 1.3, alpha = 0.65) +
    facet_wrap(~allele) +
    scale_fill_manual(values = c("WGD−" = COL_BLUE2, "WGD+" = COL_ROSE)) +
    labs(x = NULL, y = "Allele-specific copy number", title = "A  Allele-specific PIK3CA copy-number context") +
    base_theme + theme(legend.position = "none")

  prox <- copy(allele_proximity)
  pS2b <- ggplot(prox, aes(x = allele_proximity, y = N, fill = allele_proximity)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = N), vjust = -0.3, size = 3.3) +
    scale_fill_manual(values = c("Both" = COL_GREY, "Major only" = COL_BLUE, "Minor only" = COL_GREEN, "Neither" = COL_ROSE)) +
    labs(x = NULL, y = "Tumors", title = "B  Mutant-CN proximity (descriptive only)") +
    base_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

  figS2 <- pS2a | pS2b
  ggsave(file.path(FIG_DIR, "SupplementaryFigureS2_allele_specific_context.pdf"), figS2, width = 10, height = 4.5)
  ggsave(file.path(FIG_DIR, "SupplementaryFigureS2_allele_specific_context.tiff"), figS2, width = 10, height = 4.5, dpi = 600, compression = "lzw")
}

# ----------------------------------------------------------------
# 23. Save locked analysis objects and software versions
# ----------------------------------------------------------------
saveRDS(
  list(
    master = master,
    pik_final = pik_final,
    hpv_patient = hpv_patient,
    msk_sample_master = msk_sample_master,
    msk_pik_mut = msk_pik_mut,
    wgd_summary = wgd_summary,
    allele_summary = allele_summary,
    allele_proximity = allele_proximity,
    major_gain_result = major_gain_result,
    gain3q_summary = gain3q_summary,
    results_3q = results_3q,
    model_fit_summary = model_fit_summary,
    nested_tests = nested_tests,
    influence_table = influence_table,
    sensitivity_result = sensitivity_result,
    msk_prevalence = msk_prevalence,
    msk_state = msk_state,
    mutation_spectrum = mutation_spectrum,
    msk_hist_group = msk_hist_group,
    hpv_genotype_summary = hpv_genotype_summary,
    hpv_genotype_test = hpv_genotype_test,
    hpv_integration_summary = hpv_integration_summary,
    hpv_integration_tests = hpv_integration_tests
  ),
  file.path(OUTPUT_DIR, "FINAL_locked_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))

# ----------------------------------------------------------------
# 24. Human-readable QC report
# ----------------------------------------------------------------
qc <- c(
  "TCGA-CESC PIK3CA FINAL LOCKED QC REPORT",
  "========================================",
  paste0("TCGA cohort: ", nrow(master)),
  paste0("PIK3CA-mutant: ", nrow(pik_final)),
  paste0("Median mutant CN: ", sprintf("%.6f", median(pik_final$mutant_cn, na.rm = TRUE))),
  paste0("Median FAM: ", sprintf("%.6f", median(pik_final$fam, na.rm = TRUE))),
  paste0("WGD− / WGD+: ", sum(pik_final$wgd01 == 0), " / ", sum(pik_final$wgd01 == 1)),
  paste0("3q evaluable / 3q+: ", nrow(pik_3q), " / ", sum(pik_3q$gain3q01 == 1)),
  "",
  "WGD endpoints:",
  capture.output(print(wgd_summary[, .(endpoint, WGD_minus_median, WGD_plus_median, p, fdr)])),
  "",
  "Adjusted 3q endpoints:",
  capture.output(print(results_3q[, .(endpoint, beta, ci_low, ci_high, p, fdr)])),
  "",
  "Influential observations:",
  capture.output(print(influence_table[influential == TRUE])),
  "",
  paste0("TCGA vs MSK mutation Fisher P: ", format(prev_fisher$p.value, digits = 8)),
  paste0("MSK mutation vs CNA=2 Fisher P: ", format(msk_cna_fisher$p.value, digits = 8)),
  "",
  "ALL LOCKED NUMERICAL ASSERTIONS PASSED."
)
writeLines(qc, file.path(OUTPUT_DIR, "FINAL_QC_REPORT.txt"))

cat("\n============================================================\n")
cat("FINAL ANALYSIS COMPLETE — ALL LOCKED QC ASSERTIONS PASSED\n")
cat("Output directory:\n", OUTPUT_DIR, "\n", sep = "")
cat("============================================================\n")
