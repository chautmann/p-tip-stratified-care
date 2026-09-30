#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 7: Model Performance Table
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(writexl)


#######
# 2. Load fitted model results
#######

model_files <- c(
  "lm_model.RData",
  "lasso_model.RData",
  "decision_tree_model.RData",
  "random_forest_model.RData",
  "gbm_model.RData"
)

walk(
  model_files,
  ~ load(
    file.path(
      "02_models",
      .x
    ),
    envir = globalenv()
  )
)

models <- str_remove(
  model_files,
  "\\.RData$"
)


#######
# 3. Extract and combine model-performance results
#######

make_tab <- function(model, suffix) {
  
  temp <- get(model)
  
  lit <- temp[[paste0("lit_", suffix)]]$out_tab %>%
    rename_with(
      ~ paste0("lit_", .x),
      -c(imp, model, grp)
    ) %>%
    select(
      -grp
    )
  
  hit <- temp[[paste0("hit_", suffix)]]$out_tab %>%
    rename_with(
      ~ paste0("hit_", .x),
      -c(imp, model, grp)
    ) %>%
    select(
      -grp
    )
  
  left_join(
    lit,
    hit,
    by = join_by(
      imp,
      model
    )
  ) %>%
    mutate(
      model = case_when(
        str_detect(model, "lm")            ~ "Linear Regression",
        str_detect(model, "lasso")         ~ "Lasso Regression",
        str_detect(model, "decision_tree") ~ "Decision Tree",
        str_detect(model, "random_forest") ~ "Random Forest",
        str_detect(model, "gbm")           ~ "GBM"
      )
    )
}


#######
# 4. Create complete-case and single-imputation tables
#######

tab_noimp <- map_dfr(
  models,
  make_tab,
  suffix = "noimp"
) %>%
  mutate(
    analysis = "Complete-case",
    .before = 1
  )


tab_imp <- map_dfr(
  models,
  make_tab,
  suffix = "imp"
) %>%
  mutate(
    analysis = "Single imputation",
    .before = 1
  )


#######
# 5. Combine and arrange results
#######

model_order <- c(
  "Linear Regression",
  "Lasso Regression",
  "Decision Tree",
  "Random Forest",
  "GBM"
)


model_performance <- bind_rows(
  tab_noimp,
  tab_imp
) %>%
  mutate(
    model = factor(
      model,
      levels = model_order
    )
  ) %>%
  arrange(
    analysis,
    model
  ) %>%
  select(
    analysis,
    model,
    lit_n,
    lit_r2_train,
    lit_r2_test,
    lit_rmse,
    hit_n,
    hit_r2_train,
    hit_r2_test,
    hit_rmse
  )


#######
# 6. Save model-performance table
#######

write_xlsx(
  model_performance,
  file.path(
    "03_output",
    "model_performance.xlsx"
  )
)