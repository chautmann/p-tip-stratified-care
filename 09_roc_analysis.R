#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 9: ROC Analysis
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(writexl)
library(pROC)

load(file.path("01_data", "test_pred_noimp.RData"))
load(file.path("01_data", "test_pred_imp.RData"))
load(file.path("01_data", "train_pred_noimp.RData"))
load(file.path("01_data", "train_pred_imp.RData"))


#######
# 2. Define function for ROC analysis
#######

fun_roc <- function(imp) {
  
  # imp <- 0
  
  
  #######
  # 2.1 Set imputation condition and extract data
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  test_pred <- get(
    paste0("test_pred_", imp_name)
  )
  
  train_pred <- get(
    paste0("train_pred_", imp_name)
  )
  
  
  #######
  # 2.2 ROC analysis in training sample
  #
  # Higher PAI values indicate greater predicted advantage
  # for HIT and are expected to be associated with reliable
  # clinical improvement (RCI = 1).
  #######
  
  roc_train <- pROC::roc(
    response = train_pred$rci,
    predictor = train_pred$pai,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  
  #######
  # 2.3 Select optimal cut-off using Youden index
  #######
  
  cut_fit <- pROC::coords(
    roc_train,
    "best",
    ret = c(
      "threshold",
      "sensitivity",
      "specificity"
    ),
    best.method = "youden"
  )
  
  cut_of <- cut_fit %>%
    pull(threshold)
  
  
  #######
  # 2.4 Apply cut-off to training and test samples
  #######
  
  test_pred <- test_pred %>%
    mutate(
      prog_intervention = case_when(
        pai <= cut_of ~ "Below Cut-Of",
        pai > cut_of  ~ "Above Cut-Of"
      ),
      
      cut_of = case_when(
        intervention == "hit" &
          prog_intervention == "Above Cut-Of" ~ "Received optimal",
        
        intervention == "lit" &
          prog_intervention == "Above Cut-Of" ~ "Received non-optimal",
        
        prog_intervention == "Below Cut-Of" ~ "No difference"
      )
    )
  
  train_pred <- train_pred %>%
    mutate(
      prog_intervention = case_when(
        pai <= cut_of ~ "Below Cut-Of",
        pai > cut_of  ~ "Above Cut-Of"
      )
    )
  
  
  #######
  # 2.5 Save classified test data
  #######
  
  temp_name <- paste0(
    "roc_",
    imp_name
  )
  
  assign(
    temp_name,
    test_pred
  )
  
  save(
    list = temp_name,
    file = file.path(
      "01_data",
      paste0(
        temp_name,
        ".RData"
      )
    )
  )
  
  
  #######
  # 2.6 Sensitivity and specificity
  # Test sample
  #######
  
  cnt <- test_pred %>%
    count(
      rci,
      prog_intervention
    )
  
  true_positive <- cnt %>%
    filter(
      rci == 1 &
        prog_intervention == "Above Cut-Of"
    ) %>%
    pull(n)
  
  false_negative <- cnt %>%
    filter(
      rci == 1 &
        prog_intervention == "Below Cut-Of"
    ) %>%
    pull(n)
  
  sensitivity_test <-
    true_positive /
    (true_positive + false_negative)
  
  true_negative <- cnt %>%
    filter(
      rci == 0 &
        prog_intervention == "Below Cut-Of"
    ) %>%
    pull(n)
  
  false_positive <- cnt %>%
    filter(
      rci == 0 &
        prog_intervention == "Above Cut-Of"
    ) %>%
    pull(n)
  
  specificity_test <-
    true_negative /
    (true_negative + false_positive)
  
  
  #######
  # 2.7 Sensitivity and specificity
  # Training sample
  #######
  
  cnt <- train_pred %>%
    count(
      rci,
      prog_intervention
    )
  
  true_positive <- cnt %>%
    filter(
      rci == 1 &
        prog_intervention == "Above Cut-Of"
    ) %>%
    pull(n)
  
  false_negative <- cnt %>%
    filter(
      rci == 1 &
        prog_intervention == "Below Cut-Of"
    ) %>%
    pull(n)
  
  sensitivity_train <-
    true_positive /
    (true_positive + false_negative)
  
  true_negative <- cnt %>%
    filter(
      rci == 0 &
        prog_intervention == "Below Cut-Of"
    ) %>%
    pull(n)
  
  false_positive <- cnt %>%
    filter(
      rci == 0 &
        prog_intervention == "Above Cut-Of"
    ) %>%
    pull(n)
  
  specificity_train <-
    true_negative /
    (true_negative + false_positive)
  
  
  #######
  # 2.8 AUC and 95% CI
  # Training and test samples
  #######
  
  roc_test <- pROC::roc(
    response = test_pred$rci,
    predictor = test_pred$pai,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_train <- as.numeric(
    pROC::auc(roc_train)
  )
  
  auc_test <- as.numeric(
    pROC::auc(roc_test)
  )
  
  auc_train_ci <- pROC::ci.auc(
    roc_train,
    method = "delong"
  )
  
  auc_test_ci <- pROC::ci.auc(
    roc_test,
    method = "delong"
  )
  
  
  #######
  # 2.9 Create output table
  #######
  
  tab_out <- tibble(
    data = c(
      "train",
      "test"
    ),
    
    auc = c(
      sprintf(
        "%.2f [%.2f, %.2f]",
        auc_train,
        auc_train_ci[1],
        auc_train_ci[3]
      ),
      
      sprintf(
        "%.2f [%.2f, %.2f]",
        auc_test,
        auc_test_ci[1],
        auc_test_ci[3]
      )
    ),
    
    sensitivity = c(
      sprintf(
        "%.2f",
        sensitivity_train
      ),
      
      sprintf(
        "%.2f",
        sensitivity_test
      )
    ),
    
    specificity = c(
      sprintf(
        "%.2f",
        specificity_train
      ),
      
      sprintf(
        "%.2f",
        specificity_test
      )
    )
  )
  
  
  #######
  # 2.10 Save output table
  #######
  
  table_name <- paste0(
    "auc_sens_spec_",
    imp_name
  )
  
  write_xlsx(
    tab_out,
    file.path(
      "03_output",
      paste0(
        table_name,
        ".xlsx"
      )
    )
  )
  
  
  #######
  # 2.11 ROC figure
  #######
  
  roc_data <- tibble(
    Sensitivity = rev(
      roc_test$sensitivities
    ),
    
    Specificity = rev(
      roc_test$specificities
    )
  )
  
  auc_lab <- paste0(
    "AUC = ",
    sprintf(
      "%.2f",
      auc_test
    ),
    ", 95% CI [",
    sprintf(
      "%.2f",
      auc_test_ci[1]
    ),
    ", ",
    sprintf(
      "%.2f",
      auc_test_ci[3]
    ),
    "]"
  )
  
  gg <- ggplot(
    roc_data,
    aes(
      x = 1 - Specificity,
      y = Sensitivity
    )
  ) +
    geom_line(
      color = "#e41a1c",
      linewidth = 1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "#377eb8",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = auc_lab,
      size = 8,
      hjust = 0,
      vjust = 1
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme(
      plot.title = element_text(
        size = rel(2.5)
      ),
      axis.text = element_text(
        size = rel(2)
      ),
      axis.title = element_text(
        size = rel(2.25)
      )
    )
  
  
  #######
  # 2.12 Save ROC figure
  #######
  
  gg_name <- paste0(
    "figure2_",
    imp_name
  )
  
  ggsave(
    file.path(
      "03_output",
      "figures",
      paste0(
        gg_name,
        ".png"
      )
    ),
    plot = gg,
    width = 16,
    height = 13.4
  )
  
  
  #######
  # 2.13 Display key results
  #######
  
  cat(
    "\n",
    "Analysis: ", imp_name, "\n",
    "PAI threshold: ", cut_of, "\n",
    "Train AUC: ", auc_train, "\n",
    "Test AUC: ", auc_test, "\n",
    "Test AUC 95% CI: ",
    auc_test_ci[1], " - ", auc_test_ci[3], "\n",
    "Test sensitivity: ", sensitivity_test, "\n",
    "Test specificity: ", specificity_test, "\n",
    sep = ""
  )
}


#######
# 3. Run analyses
#######

fun_roc(0)

fun_roc(1)