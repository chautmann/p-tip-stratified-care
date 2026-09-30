#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 13: Concordance and Sensitivity Analyses
#
# This script reproduces analyses of reliable clinical improvement
# according to decision-rule treatment classification and three
# sensitivity analyses.
#
# Specifically, it:
# (1) describes reliable clinical improvement among cases receiving
#     the recommended intensity, receiving a non-recommended intensity,
#     or having no predicted clinically meaningful advantage of HIT
#     over LIT;
# (2) compares reliable improvement between cases receiving the
#     recommended versus non-recommended intensity;
# (3) compares cases with no predicted clinically meaningful advantage
#     versus cases receiving a non-recommended intensity;
# (4) conducts Sensitivity Analysis 1 by repeating the primary
#     recommended versus non-recommended comparison while adjusting
#     for prespecified baseline covariates; and
# (5) conducts Sensitivity Analysis 2 by repeating the adjusted
#     concordance analysis after restricting HIT recipients receiving
#     the recommended intensity to maximum treatment durations of
#     6, 9, and 12 months; and
# (6) conducts Sensitivity Analysis 3 using the FBB-SSV Oppositional
#     Defiant Disorder scale as a continuous outcome.
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library(tidyverse)
library(broom)

load(file.path("01_data", "roc_noimp.RData"))
load(file.path("01_data", "ptip_hit_fil0.RData"))
load(file.path("01_data", "ptip_lit_fil0.RData"))


#######
# 2. Quality checks
#######

# Recommendation-group counts
roc_noimp %>%
  count(cut_of)

# Check uniqueness of synthetic analysis IDs
if (anyDuplicated(roc_noimp$id) > 0) {
  stop("ERROR: id is not unique in roc_noimp.")
}

# Check that the original participant ID is available
if (!"original_id" %in% names(roc_noimp)) {
  stop("ERROR: original_id is missing from roc_noimp.")
}

# Check uniqueness of original participant IDs among HIT cases
if (
  anyDuplicated(
    roc_noimp %>%
    filter(intervention == "hit") %>%
    pull(original_id)
  ) > 0
) {
  stop("ERROR: original_id is not unique among HIT cases in roc_noimp.")
}

# Check uniqueness of participant IDs in the original HIT dataset
if (anyDuplicated(ptip_hit_fil0$id) > 0) {
  stop("ERROR: id is not unique in ptip_hit_fil0.")
}

if (anyDuplicated(ptip_lit_fil0$id) > 0) {
  stop("ERROR: id is not unique in ptip_lit_fil0.")
}


#######
# 3. Descriptive results:
# Reliable clinical improvement by decision-rule treatment classification
#######

group_results <- roc_noimp %>%
  group_by(cut_of) %>%
  summarise(
    n = n(),
    improved = sum(rci == 1),
    pct_improved = 100 * mean(rci == 1),
    .groups = "drop"
  )

group_results


#######
# 4. Primary concordance comparison:
# Recommended vs. non-recommended treatment intensity
#######

dat_recommended <- roc_noimp %>%
  filter(
    cut_of %in% c(
      "Received non-optimal",
      "Received optimal"
    )
  ) %>%
  mutate(
    rci = factor(
      rci,
      levels = c(0, 1)
    ),
    cut_of = factor(
      cut_of,
      levels = c(
        "Received non-optimal",
        "Received optimal"
      )
    )
  )


#######
# 4a. Descriptive results for primary comparison
#######

primary_results <- dat_recommended %>%
  group_by(cut_of) %>%
  summarise(
    n = n(),
    improved = sum(rci == 1),
    pct_improved = 100 * mean(rci == 1),
    .groups = "drop"
  )

primary_results


#######
# 4b. Unadjusted logistic regression
#######

model_recommended <- glm(
  rci ~ cut_of,
  data = dat_recommended,
  family = binomial
)

or_recommended <- broom::tidy(
  model_recommended,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    term == "cut_ofReceived optimal"
  ) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

or_recommended
or_recommended %>%
  mutate(
    across(c(estimate, conf.low, conf.high), ~ sprintf("%.2f", .x))
  )

#######
# 5. Additional comparison:
# No predicted clinically meaningful advantage vs.
# non-recommended treatment intensity
#######

dat_no_advantage <- roc_noimp %>%
  filter(
    cut_of %in% c(
      "Received non-optimal",
      "No difference"
    )
  ) %>%
  mutate(
    rci = factor(
      rci,
      levels = c(0, 1)
    ),
    cut_of = factor(
      cut_of,
      levels = c(
        "Received non-optimal",
        "No difference"
      )
    )
  )


#######
# 5a. Unadjusted logistic regression
#######

model_no_advantage <- glm(
  rci ~ cut_of,
  data = dat_no_advantage,
  family = binomial
)

or_no_advantage <- broom::tidy(
  model_no_advantage,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    term == "cut_ofNo difference"
  ) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

or_no_advantage
or_no_advantage %>%
  mutate(
    across(c(estimate, conf.low, conf.high), ~ sprintf("%.2f", .x))
  )

#######
# 6. Sensitivity Analysis 1:
# Covariate-adjusted comparison of recommended vs.
# non-recommended treatment intensity
#
# Adjustment variables:
# - Baseline CBCL Aggressive Behavior
# - Child age
# - Single-parent status
# - Maternal education
# - Paternal education
# - Number of children in the household
#######

model_recommended_adjusted <- glm(
  rci ~ cut_of +
    cbcl_av_baseline +
    age +
    single_parent +
    casmin_km +
    casmin_kv +
    number_children,
  data = dat_recommended,
  family = binomial
)

or_recommended_adjusted <- broom::tidy(
  model_recommended_adjusted,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    term == "cut_ofReceived optimal"
  ) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

or_recommended_adjusted
or_recommended_adjusted %>%
  mutate(
    across(c(estimate, conf.low, conf.high), ~ sprintf("%.2f", .x))
  )

#######
# 7. Sensitivity Analysis 2:
# Treatment-duration-restricted concordance analyses
#
# This analysis examines whether the association between
# receiving the recommended intensity and reliable clinical
# improvement could be explained by longer treatment duration
# among HIT recipients.
#
# Treatment duration is added for HIT cases using the preserved
# original participant ID.
#######

hit_duration <- ptip_hit_fil0 %>%
  transmute(
    original_id = as.character(id),
    intervention = "hit",
    ther_duration
  )

roc_duration <- roc_noimp %>%
  mutate(
    original_id = as.character(original_id),
    intervention = as.character(intervention)
  ) %>%
  left_join(
    hit_duration,
    by = c(
      "original_id",
      "intervention"
    )
  )


#######
# 7a. Quality checks after merge
#######

# Recommendation-group counts should remain unchanged
roc_duration %>%
  count(cut_of)

# Check availability of treatment-duration information
# among cases receiving the recommended intensity
roc_duration %>%
  filter(
    cut_of == "Received optimal"
  ) %>%
  summarise(
    n = n(),
    n_duration_available = sum(!is.na(ther_duration)),
    n_duration_missing = sum(is.na(ther_duration))
  )


#######
# 7b. Define treatment-duration cutoffs
#
# ther_duration is measured in weeks.
# Months are converted using 4.345 weeks per month.
#######

weeks_per_month <- 4.345

cut_6m_weeks  <-  6 * weeks_per_month
cut_9m_weeks  <-  9 * weeks_per_month
cut_12m_weeks <- 12 * weeks_per_month


#######
# 7c. Function to restrict analyses by HIT treatment duration
#
# Cases receiving a non-recommended intensity received LIT and
# therefore remain in each analysis.
#
# Cases receiving the recommended intensity received HIT and
# are retained only when treatment duration is at or below the
# specified cutoff.
#######

restrict_by_duration <- function(data, cutoff_weeks) {
  
  data %>%
    filter(
      cut_of %in% c(
        "Received non-optimal",
        "Received optimal"
      )
    ) %>%
    filter(
      cut_of == "Received non-optimal" |
        (
          cut_of == "Received optimal" &
            !is.na(ther_duration) &
            ther_duration <= cutoff_weeks
        )
    ) %>%
    mutate(
      rci = factor(
        rci,
        levels = c(0, 1)
      ),
      cut_of = factor(
        cut_of,
        levels = c(
          "Received non-optimal",
          "Received optimal"
        )
      )
    )
}


dat_6m <- restrict_by_duration(
  roc_duration,
  cut_6m_weeks
)

dat_9m <- restrict_by_duration(
  roc_duration,
  cut_9m_weeks
)

dat_12m <- restrict_by_duration(
  roc_duration,
  cut_12m_weeks
)


#######
# 7d. Sample sizes by duration cutoff
#######

dat_6m %>%
  count(cut_of)

dat_9m %>%
  count(cut_of)

dat_12m %>%
  count(cut_of)


#######
# 7e. Descriptive reliable-improvement results
#######

summarise_improvement <- function(data, duration_label) {
  
  data %>%
    group_by(cut_of) %>%
    summarise(
      improved = sum(rci == 1, na.rm = TRUE),
      total = n(),
      percent = 100 * improved / total,
      .groups = "drop"
    ) %>%
    mutate(
      duration = duration_label,
      n_N = paste0(
        improved,
        "/",
        total
      )
    ) %>%
    select(
      duration,
      cut_of,
      n_N,
      percent
    )
}


tab_6m <- summarise_improvement(
  dat_6m,
  "≤ 6"
)

tab_9m <- summarise_improvement(
  dat_9m,
  "≤ 9"
)

tab_12m <- summarise_improvement(
  dat_12m,
  "≤ 12"
)

tab_all <- bind_rows(
  tab_6m,
  tab_9m,
  tab_12m
)

tab_all


#######
# 7f. Function for covariate-adjusted logistic regression
#
# The same covariate adjustment as in Sensitivity Analysis 1
# is used.
#######

run_adjusted_or <- function(data, duration_label) {
  
  fit <- glm(
    rci ~ cut_of +
      cbcl_av_baseline +
      age +
      single_parent +
      casmin_km +
      casmin_kv +
      number_children,
    data = data,
    family = binomial
  )
  
  broom::tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == "cut_ofReceived optimal"
    ) %>%
    mutate(
      duration = duration_label
    ) %>%
    select(
      duration,
      term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}


#######
# 7g. ≤ 6 months
#
# The manuscript reports that the adjusted effect was not
# estimable because only five cases receiving the recommended
# intensity remained at this cutoff.
#
# Therefore, no adjusted model is fitted for this cutoff.
#######

dat_6m %>%
  count(cut_of)


#######
# 7h. ≤ 9 months
#######

or_9m <- run_adjusted_or(
  dat_9m,
  "≤ 9"
)

or_9m


#######
# 7i. ≤ 12 months
#######

or_12m <- run_adjusted_or(
  dat_12m,
  "≤ 12"
)

or_12m


#######
# 7j. Combine adjusted odds ratios
#######

or_duration <- bind_rows(
  or_9m,
  or_12m
)

or_duration


#######
# 8. Supplementary Table S7
#######

or_table <- or_duration %>%
  transmute(
    duration,
    OR = sprintf("%.2f", estimate),
    CI = paste0(
      "[",
      sprintf("%.2f", conf.low),
      ", ",
      sprintf("%.2f", conf.high),
      "]"
    )
  )

table_s7 <- tab_all %>%
  mutate(
    percent = sprintf("%.1f", percent)
  ) %>%
  pivot_wider(
    names_from = cut_of,
    values_from = c(
      n_N,
      percent
    )
  ) %>%
  left_join(
    or_table,
    by = "duration"
  ) %>%
  mutate(
    duration = factor(
      duration,
      levels = c(
        "≤ 6",
        "≤ 9",
        "≤ 12"
      )
    )
  ) %>%
  arrange(
    duration
  )


#######
# Final Supplementary Table S7
#######

table_s7

#######
# 9. Sensitivity Analysis 3:
# FBB-SSV Oppositional Defiant Disorder outcome
#
# This analysis examines whether the concordance findings generalize
# to a related continuous measure of conduct problems.
#
# FBB-SSV data are linked to the held-out test sample using the
# preserved original participant ID and treatment group.
#######

ssv_all <- bind_rows(
  ptip_hit_fil0 %>%
    transmute(
      original_id = as.character(id),
      intervention = "hit",
      ssv_odd_baseline,
      ssv_odd_post
    ),

  ptip_lit_fil0 %>%
    transmute(
      original_id = as.character(id),
      intervention = "lit",
      ssv_odd_baseline,
      ssv_odd_post
    )
)


#######
# 9a. Quality checks and merge
#######

if (
  anyDuplicated(
    ssv_all %>%
      select(
        original_id,
        intervention
      )
  ) > 0
) {
  stop("ERROR: Participant ID is not unique within treatment group in SSV data.")
}

n_before_ssv <- nrow(roc_noimp)

roc_ssv <- roc_noimp %>%
  mutate(
    original_id = as.character(original_id),
    intervention = as.character(intervention)
  ) %>%
  left_join(
    ssv_all,
    by = c(
      "original_id",
      "intervention"
    )
  )

if (nrow(roc_ssv) != n_before_ssv) {
  stop("ERROR: Row count changed after FBB-SSV merge.")
}


#######
# 9b. Prepare analysis dataset
#
# Restrict to cases with a recommendation for HIT and complete
# FBB-SSV baseline and posttreatment data.
#######

analysis_ssv <- roc_ssv %>%
  filter(
    cut_of %in% c(
      "Received non-optimal",
      "Received optimal"
    )
  ) %>%
  mutate(
    recommended = if_else(
      cut_of == "Received optimal",
      1,
      0
    )
  ) %>%
  filter(
    !is.na(ssv_odd_baseline),
    !is.na(ssv_odd_post)
  )


#######
# 9c. Sample-size checks
#######

# Primary concordance test sample
roc_ssv %>%
  filter(
    cut_of %in% c(
      "Received non-optimal",
      "Received optimal"
    )
  ) %>%
  summarise(
    n = n()
  )

# FBB-SSV sensitivity-analysis sample
analysis_ssv %>%
  summarise(
    n = n()
  )

# Sample size by concordance group
analysis_ssv %>%
  count(
    cut_of
  )


#######
# 9d. ANCOVA model
#
# Posttreatment FBB-SSV ODD symptoms are predicted from concordance
# with the treatment recommendation while adjusting for baseline
# FBB-SSV ODD symptoms and the same baseline covariates used in
# Sensitivity Analysis 1.
#######

model_ssv <- lm(
  ssv_odd_post ~ recommended +
    ssv_odd_baseline +
    age +
    single_parent +
    casmin_km +
    casmin_kv +
    number_children,
  data = analysis_ssv
)

ssv_results <- broom::tidy(
  model_ssv,
  conf.int = TRUE
)

ssv_results


#######
# 9e. Pooled baseline FBB-SSV standard deviation
#######

sd_pooled_ssv <- analysis_ssv %>%
  group_by(
    recommended
  ) %>%
  summarise(
    sd = sd(
      ssv_odd_baseline,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  summarise(
    sd_pooled = sqrt(
      (
        (n[1] - 1) * sd[1]^2 +
          (n[2] - 1) * sd[2]^2
      ) /
        (sum(n) - 2)
    )
  ) %>%
  pull(
    sd_pooled
  )

sd_pooled_ssv


#######
# 9f. Standardized effect size
#
# Cohen's d is calculated by dividing the adjusted regression
# coefficient and its confidence limits by the pooled baseline
# FBB-SSV standard deviation.
#######

d_ssv <- ssv_results %>%
  filter(
    term == "recommended"
  ) %>%
  mutate(
    cohens_d = estimate / sd_pooled_ssv,
    d_low = conf.low / sd_pooled_ssv,
    d_high = conf.high / sd_pooled_ssv
  ) %>%
  select(
    estimate,
    conf.low,
    conf.high,
    p.value,
    cohens_d,
    d_low,
    d_high
  )

d_ssv
