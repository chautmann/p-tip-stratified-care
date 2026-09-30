#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 8: Personalized Predictions Using Final LASSO Models
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)
library(tidyverse)

load(file.path("01_data", "prediction_data.RData"))
load(file.path("01_data", "split_ptip.RData"))


#######
# 2. Reliable-change parameter
#
# RCI denominator used in the original analyses.
# The value is defined explicitly here to avoid dependency
# on a separate binary RData file.
#######

sdiff <- 3.216558


#######
# 3. Define prediction function
#######

fun_pred_test <- function(imp) {
  
  # imp <- 0
  
  
  #######
  # 3.1 Set imputation condition
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  data <- prediction_data[[imp_name]] %>%
    mutate(
      baseline = cbcl_av_baseline
    )
  
  
  #######
  # 3.2 Load final LIT and HIT LASSO models
  #######
  
  algo_list <- c("lit", "hit") %>%
    set_names() %>%
    map(\(grp) {
      
      model_name <- paste0(
        "lasso_",
        grp,
        "_",
        imp_name
      )
      
      load(
        file.path(
          "02_models",
          paste0(model_name, ".RData")
        ),
        envir = .GlobalEnv
      )
      
      get(
        model_name,
        envir = .GlobalEnv
      )
    })
  
  
  #######
  # 3.3 Define bootstrap prediction function
  #     for one treatment group
  #######
  
  fun_pred <- function(algo_grp) {
    
    # algo_grp <- "lit"
    
    algo <- algo_list[[algo_grp]]
    
    
    #######
    # Extract corresponding training data
    #######
    
    train_data_name <- paste0(
      algo_grp,
      "_",
      imp_name
    )
    
    train <- split_ptip[[train_data_name]]$train
    
    
    #######
    # Generate bootstrap samples
    #######
    
    set.seed(9112001)
    
    bootstrap_resamples <- bootstraps(
      train,
      times = 1000
    )
    
    
    #######
    # Refit the selected LASSO model specification
    # in each bootstrap sample and generate predictions
    # for all individuals in the held-out test sample
    #######
    
    bootstrap_resamples <- bootstrap_resamples %>%
      mutate(
        predictions = map(
          splits,
          function(split) {
            
            fit <- algo %>%
              fit(
                data = analysis(split)
              )
            
            predict(
              fit,
              new_data = data
            ) %>%
              mutate(
                .row = row_number()
              )
          }
        )
      )
    
    
    #######
    # Define prediction variable names
    #######
    
    pred_var_name_mean <- paste0(
      algo_grp,
      "_pred"
    )
    
    pred_var_name_loci <- paste0(
      algo_grp,
      "_loci"
    )
    
    pred_var_name_upci <- paste0(
      algo_grp,
      "_upci"
    )
    
    
    #######
    # Summarize bootstrap predictions
    #
    # The mean across the 1,000 bootstrap refits is used
    # as the individual point prediction.
    #
    # The 2.5th and 97.5th percentiles describe the
    # bootstrap uncertainty distribution of the predicted
    # outcome for each individual.
    #######
    
    bootstrap_fin <- bootstrap_resamples %>%
      unnest(predictions) %>%
      group_by(.row) %>%
      summarize(
        !!pred_var_name_mean := mean(.pred),
        
        !!pred_var_name_loci := quantile(
          .pred,
          probs = 0.025
        ),
        
        !!pred_var_name_upci := quantile(
          .pred,
          probs = 0.975
        ),
        
        .groups = "drop"
      )
    
    
    #######
    # Assemble prediction data
    #######
    
    pred_data <- bind_cols(
      data,
      bootstrap_fin
    ) %>%
      select(
        id,
        intervention,
        baseline,
        outcome,
        starts_with(algo_grp)
      )
    
    return(pred_data)
  }
  
  
  #######
  # 3.4 Generate bootstrap-based LIT and HIT predictions
  #######
  
  pred_lit <- fun_pred("lit")
  
  pred_hit <- fun_pred("hit")
  
  pred <- pred_lit %>%
    left_join(
      pred_hit,
      by = join_by(
        id,
        intervention,
        baseline,
        outcome
      )
    )
  
  pred <- data %>%
    select(
      id,
      everything()
    ) %>%
    left_join(
      pred,
      by = join_by(
        id,
        outcome,
        intervention,
        baseline
      )
    )
  
  
  #######
  # 3.5 Calculate MPI, PAI, and reliable clinical improvement
  #
  # MPI = mean predicted outcome across LIT and HIT
  #
  # PAI = predicted outcome under LIT minus predicted
  #       outcome under HIT
  #
  # Positive PAI values indicate a predicted advantage
  # for HIT because lower outcome scores are more favorable.
  #######
  
  pred <- pred %>%
    mutate(
      mpi = (
        lit_pred + hit_pred
      ) / 2,
      
      pai =
        lit_pred - hit_pred,
      
      rci_nr = (
        outcome - cbcl_av_baseline
      ) / sdiff,
      
      rci = if_else(
        rci_nr <= -1.96,
        1,
        0
      )
    )
  
  
  #######
  # 3.6 Generate training-sample predictions
  #
  # Training predictions are obtained directly from the
  # final fitted LIT and HIT LASSO models. These predictions
  # are subsequently used to derive the PAI threshold in
  # the ROC analysis.
  #######
  
  train_df <- bind_rows(
    split_ptip[[paste0(
      "lit_",
      imp_name
    )]]$train,
    
    split_ptip[[paste0(
      "hit_",
      imp_name
    )]]$train
  )
  
  train_pred <- bind_cols(
    
    predict(
      algo_list$lit,
      new_data = train_df
    ) %>%
      rename(
        lit_pred = .pred
      ),
    
    predict(
      algo_list$hit,
      new_data = train_df
    ) %>%
      rename(
        hit_pred = .pred
      ),
    
    train_df %>%
      select(
        cbcl_av_baseline,
        outcome
      )
  )
  
  train_pred <- train_pred %>%
    mutate(
      mpi = (
        lit_pred + hit_pred
      ) / 2,
      
      pai =
        lit_pred - hit_pred,
      
      rci_nr = (
        outcome - cbcl_av_baseline
      ) / sdiff,
      
      rci = if_else(
        rci_nr <= -1.96,
        1,
        0
      )
    )
  
  
  #######
  # 3.7 Save test-sample predictions
  #######
  
  save_name <- paste0(
    "test_pred_",
    imp_name
  )
  
  assign(
    save_name,
    pred
  )
  
  save(
    list = save_name,
    file = file.path(
      "01_data",
      paste0(
        save_name,
        ".RData"
      )
    )
  )
  
  
  #######
  # 3.8 Save training-sample predictions
  #######
  
  save_name <- paste0(
    "train_pred_",
    imp_name
  )
  
  assign(
    save_name,
    train_pred
  )
  
  save(
    list = save_name,
    file = file.path(
      "01_data",
      paste0(
        save_name,
        ".RData"
      )
    )
  )
}


#######
# 4. Run prediction function
# Complete-case and single-imputation analyses
#######

fun_pred_test(0)

fun_pred_test(1)