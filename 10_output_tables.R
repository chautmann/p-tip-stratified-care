#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 10: Output Tables
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)
library(tidyverse)
library(writexl)

load(file.path("01_data", "split_ptip.RData"))


#######
# 2. Table 1
# Baseline demographic and clinical characteristics
# of the total sample and the LIT and HIT samples
#######


#######
# 2.1 Prepare data
#######

lit_hit <- bind_rows(
  lit = split_ptip$lit_noimp$split$data,
  hit = split_ptip$hit_noimp$split$data,
  .id = "grp"
) %>%
  select(-id)

lit <- split_ptip$lit_noimp$split$data


#######
# 2.2 Define categorical levels to be reported
#######

categorical_level <- c(
  sex = "female",
  medication = "yes",
  mother = "biological mother",
  single_parent = "yes"
)


#######
# 2.3 Helper function:
# continuous variables -> M (SD)
# categorical variables -> n (%)
#######

format_desc <- function(data, var) {
  
  x <- data[[var]]
  
  if (is.numeric(x)) {
    
    sprintf(
      "%.2f (%.2f)",
      mean(x, na.rm = TRUE),
      sd(x, na.rm = TRUE)
    )
    
  } else {
    
    value <- categorical_level[[var]]
    
    if (is.null(value)) {
      stop(
        paste0(
          "No reporting category defined for variable: ",
          var
        )
      )
    }
    
    n <- sum(
      x == value,
      na.rm = TRUE
    )
    
    N <- sum(
      !is.na(x)
    )
    
    sprintf(
      "%d (%.1f)",
      n,
      100 * n / N
    )
  }
}


#######
# 2.4 Format p-values
#######

format_p <- function(p) {
  
  case_when(
    p < .001 ~ "< .001",
    p < .01  ~ "< .01",
    p < .05  ~ "< .05",
    TRUE     ~ sprintf("%.3f", p)
  )
}


#######
# 2.5 Compare LIT and HIT
#######

get_p <- function(data, var) {
  
  x <- data[[var]]
  
  if (is.numeric(x)) {
    
    p <- data %>%
      transmute(
        grp = factor(grp),
        out = .data[[var]]
      ) %>%
      t_test(
        out ~ grp,
        order = c("hit", "lit")
      ) %>%
      pull(p_value)
    
  } else {
    
    tab <- table(
      data$grp,
      x
    )
    
    p <- chisq.test(tab)$p.value
  }
  
  format_p(p)
}


#######
# 2.6 Variables and table order
#######

variables <- c(
  "sex",
  "age",
  "cbcl_ad_baseline",
  "cbcl_ap_baseline",
  "cbcl_av_baseline",
  "medication",
  "mother",
  "single_parent",
  "casmin_km",
  "casmin_kv",
  "number_children"
)


#######
# 2.7 Create Table 1
#######

table_1 <- tibble(
  var = variables
) %>%
  mutate(
    
    type = if_else(
      map_lgl(
        var,
        ~ is.numeric(lit_hit[[.x]])
      ),
      1,
      0
    ),
    
    tot = map_chr(
      var,
      ~ format_desc(
        lit_hit,
        .x
      )
    ),
    
    lit = map_chr(
      var,
      ~ format_desc(
        filter(
          lit_hit,
          grp == "lit"
        ),
        .x
      )
    ),
    
    hit = map_chr(
      var,
      ~ format_desc(
        filter(
          lit_hit,
          grp == "hit"
        ),
        .x
      )
    ),
    
    p = map_chr(
      var,
      ~ get_p(
        lit_hit,
        .x
      )
    )
  )


#######
# 2.8 Save Table 1
#######

write_xlsx(
  table_1,
  file.path(
    "03_output",
    "table_1.xlsx"
  )
)


#######
# 3. Supplementary Table S1
# Baseline demographic and clinical characteristics
# across the three LIT trials
#######


#######
# 3.1 Prepare data
#######

lit_tab <- lit %>%
  mutate(
    study = case_when(
      str_starts(
        str_to_lower(id),
        "adopt"
      ) ~ "ADOPT",
      
      str_starts(
        str_to_lower(id),
        "esca"
      ) ~ "ESCAadol",
      
      str_starts(
        str_to_lower(id),
        "wash"
      ) ~ "WASH",
      
      TRUE ~ NA_character_
    )
  )


#######
# 3.2 Helper functions
#######

mean_sd <- function(x) {
  
  sprintf(
    "%.2f (%.2f)",
    mean(
      x,
      na.rm = TRUE
    ),
    sd(
      x,
      na.rm = TRUE
    )
  )
}


n_pct <- function(x, value) {
  
  n <- sum(
    x == value,
    na.rm = TRUE
  )
  
  N <- sum(
    !is.na(x)
  )
  
  sprintf(
    "%d (%.1f)",
    n,
    100 * n / N
  )
}


#######
# 3.3 Function to create one table column
#######

make_column <- function(data) {
  
  c(
    "Child",
    
    n_pct(
      data$sex,
      "female"
    ),
    
    mean_sd(
      data$age
    ),
    
    mean_sd(
      data$cbcl_ad_baseline
    ),
    
    mean_sd(
      data$cbcl_ap_baseline
    ),
    
    mean_sd(
      data$cbcl_av_baseline
    ),
    
    n_pct(
      data$medication,
      "yes"
    ),
    
    "Family",
    
    n_pct(
      data$mother,
      "biological mother"
    ),
    
    n_pct(
      data$single_parent,
      "yes"
    ),
    
    mean_sd(
      data$casmin_km
    ),
    
    mean_sd(
      data$casmin_kv
    ),
    
    mean_sd(
      data$number_children
    )
  )
}


#######
# 3.4 Create Supplementary Table S1
#######

table_s1 <- tibble(
  
  measure = c(
    "Child",
    "Sex, female, n (%)",
    "Age, M (SD)",
    "CBCL Anxious/Depressed, M (SD)",
    "CBCL Attention Problems, M (SD)",
    "CBCL Aggressive Behavior, M (SD)",
    "Medication, yes, n (%)",
    
    "Family",
    "Biological mother, yes, n (%)",
    "Single-parent status, yes, n (%)",
    "Maternal education, M (SD)",
    "Paternal education, M (SD)",
    "No. of children in household, M (SD)"
  ),
  
  `Total LIT\n(N = 375)` =
    make_column(
      lit_tab
    ),
  
  `ADOPT\n(n = 183)` =
    make_column(
      filter(
        lit_tab,
        study == "ADOPT"
      )
    ),
  
  `ESCAadol\n(n = 39)` =
    make_column(
      filter(
        lit_tab,
        study == "ESCAadol"
      )
    ),
  
  `WASH\n(n = 153)` =
    make_column(
      filter(
        lit_tab,
        study == "WASH"
      )
    )
)


#######
# 3.5 Save Supplementary Table S1
#######

write_xlsx(
  table_s1,
  file.path(
    "03_output",
    "table_s1.xlsx"
  )
)


#######
# 4. Supplementary Tables S2 and S3
# Bivariate correlations among baseline predictors
# and posttreatment outcome
#
# S2 = LIT sample
# S3 = HIT sample
#######


#######
# 4.1 Function to create correlation table
#######

make_cor_table <- function(data, group) {
  
  cor_dat <- data %>%
    filter(
      grp == group
    ) %>%
    transmute(
      outcome         = outcome,
      sex             = as.numeric(sex == "female"),
      age             = age,
      cbcl_ad         = cbcl_ad_baseline,
      cbcl_ap         = cbcl_ap_baseline,
      cbcl_av         = cbcl_av_baseline,
      medication      = as.numeric(medication == "yes"),
      mother          = as.numeric(mother == "biological mother"),
      single_parent   = as.numeric(single_parent == "yes"),
      maternal_edu    = casmin_km,
      paternal_edu    = casmin_kv,
      number_children = number_children
    )
  
  
  labels <- c(
    "Outcome: CBCL Aggressive Behavior",
    "Sex, female = 1",
    "Age",
    "CBCL Anxious/Depressed",
    "CBCL Attention Problems",
    "CBCL Aggressive Behavior",
    "Medication, yes = 1",
    "Biological mother, yes = 1",
    "Single-parent status, yes = 1",
    "Maternal education",
    "Paternal education",
    "No. of children in household"
  )
  
  
  #######
  # Correlation and significance stars
  #######
  
  cor_cell <- function(x, y) {
    
    dat <- tibble(
      x,
      y
    ) %>%
      drop_na()
    
    test <- cor.test(
      dat$x,
      dat$y,
      method = "pearson"
    )
    
    r <- unname(
      test$estimate
    )
    
    p <- test$p.value
    
    stars <- case_when(
      p < .001 ~ "***",
      p < .01  ~ "**",
      p < .05  ~ "*",
      TRUE     ~ ""
    )
    
    r <- sprintf(
      "%.2f",
      r
    ) %>%
      str_remove(
        "^0"
      ) %>%
      str_replace(
        "^-0",
        "–"
      )
    
    paste0(
      r,
      stars
    )
  }
  
  
  #######
  # Create empty correlation table
  #######
  
  n_var <- ncol(
    cor_dat
  )
  
  tab <- matrix(
    "",
    nrow = n_var,
    ncol = n_var - 1
  )
  
  
  #######
  # Fill lower triangle
  #######
  
  for (i in 2:n_var) {
    
    for (j in 1:(i - 1)) {
      
      tab[i, j] <- cor_cell(
        cor_dat[[i]],
        cor_dat[[j]]
      )
    }
  }
  
  
  #######
  # Convert to tibble
  #######
  
  out <- as_tibble(
    tab,
    .name_repair = "minimal"
  )
  
  names(out) <- 1:(n_var - 1)
  
  out %>%
    mutate(
      nr = paste0(
        row_number(),
        "."
      ),
      measure = labels,
      .before = 1
    )
}


#######
# 4.2 Create Supplementary Table S2
# LIT sample
#######

table_s2_lit <- make_cor_table(
  lit_hit,
  "lit"
)


#######
# 4.3 Create Supplementary Table S3
# HIT sample
#######

table_s3_hit <- make_cor_table(
  lit_hit,
  "hit"
)


#######
# 4.4 Save Supplementary Tables S2 and S3
#######

write_xlsx(
  table_s2_lit,
  file.path(
    "03_output",
    "table_s2_lit.xlsx"
  )
)

write_xlsx(
  table_s3_hit,
  file.path(
    "03_output",
    "table_s3_hit.xlsx"
  )
)