#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 4: Decision Tree Model
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)

load(file.path(getwd(), "01_data", "split_ptip.RData"))


#######
# 2. Define modeling function
#######

fun_dec_tree <- function(grp, imp) {
  
  set.seed(24120000)
  
  # grp <- "lit"
  # imp <- 0
  
  # Define imputation condition and dataset name
  imp_name <- if_else(imp == 1, "imp", "noimp")
  
  df_name <- paste0(grp, "_", imp_name)
  
  # Extract data split and training/test data
  temp <- split_ptip[[df_name]]
  
  split <- temp$split
  
  split$data <- split$data %>%
    select(-id)
  
  train <- temp$train %>%
    select(-id)
  
  # Initialize model output object
  decision_tree_model <- list()
  
  
  #######
  # 2.1 Specify and tune decision tree model
  #######
  
  decision_tree_model$tune_spec <- decision_tree(
    cost_complexity = tune(),
    tree_depth = tune(),
    min_n = tune()
  ) %>%
    set_engine("rpart") %>%
    set_mode("regression")
  
  decision_tree_model$tree_grid <- grid_max_entropy(
    cost_complexity(),
    tree_depth(),
    min_n(),
    size = 200
  )
  
  set.seed(201157)
  
  decision_tree_model$folds <- vfold_cv(train)
  
  set.seed(11901)
  
  decision_tree_model$workflow <- workflow() %>%
    add_model(decision_tree_model$tune_spec) %>%
    add_formula(outcome ~ .)
  
  decision_tree_model$tree_res <-
    decision_tree_model$workflow %>%
    tune_grid(
      resamples = decision_tree_model$folds,
      grid = decision_tree_model$tree_grid
    )
  
  decision_tree_model$best_tree <- decision_tree_model$tree_res %>%
    select_best(metric = "rmse")
  
  decision_tree_model$final_decision_tree <-
    decision_tree_model$workflow %>%
    finalize_workflow(decision_tree_model$best_tree)
  
  decision_tree_model$fit <-
    decision_tree_model$final_decision_tree %>%
    last_fit(split)
  
  decision_tree_model$fin_fit <- decision_tree_model$fit %>%
    collect_metrics()
  
  
  #######
  # 2.2 Training-set performance
  #######
  
  decision_tree_model$train_fit <-
    decision_tree_model$final_decision_tree %>%
    fit(train)
  
  # Generate predictions for the training set
  decision_tree_model$train_pred <- predict(
    decision_tree_model$train_fit,
    new_data = train
  ) %>%
    bind_cols(train) %>%
    select(outcome, .pred)
  
  # Calculate training-set metrics
  decision_tree_model$train_fit_metrics <-
    decision_tree_model$train_pred %>%
    metrics(
      truth = outcome,
      estimate = .pred
    )
  
  
  #######
  # 2.3 Bootstrap 95% confidence intervals for model performance
  #######
  
  set.seed(79110)
  
  # Training-set bootstrap
  decision_tree_model$train_r2_boot <- bootstraps(
    decision_tree_model$train_pred,
    times = 2000
  ) %>%
    mutate(
      r2 = map_dbl(
        splits,
        ~ rsq_vec(
          analysis(.x)$outcome,
          analysis(.x)$.pred
        )
      )
    )
  
  decision_tree_model$train_r2_ci <- quantile(
    decision_tree_model$train_r2_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  # Test-set predictions
  decision_tree_model$test_pred <- decision_tree_model$fit %>%
    collect_predictions() %>%
    select(outcome, .pred)
  
  # Calculate R2 and RMSE for each bootstrap test sample
  decision_tree_model$test_metrics_boot <- bootstraps(
    decision_tree_model$test_pred,
    times = 2000
  ) %>%
    mutate(
      r2 = map_dbl(
        splits,
        ~ rsq_vec(
          truth = analysis(.x)$outcome,
          estimate = analysis(.x)$.pred
        )
      ),
      rmse = map_dbl(
        splits,
        ~ rmse_vec(
          truth = analysis(.x)$outcome,
          estimate = analysis(.x)$.pred
        )
      )
    )
  
  # 95% CI for test-set R2
  decision_tree_model$test_r2_ci <- quantile(
    decision_tree_model$test_metrics_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  # 95% CI for test-set RMSE
  decision_tree_model$test_rmse_ci <- quantile(
    decision_tree_model$test_metrics_boot$rmse,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  #######
  # 2.4 Format model performance for output table
  #######
  
  decision_tree_model$test_fit <- tibble(
    .metric = c("rmse", "rsq"),
    .estimate = c(
      rmse_vec(
        decision_tree_model$test_pred$outcome,
        decision_tree_model$test_pred$.pred
      ),
      rsq_vec(
        decision_tree_model$test_pred$outcome,
        decision_tree_model$test_pred$.pred
      )
    ),
    lower = c(
      decision_tree_model$test_rmse_ci[["2.5%"]],
      decision_tree_model$test_r2_ci[["2.5%"]]
    ),
    upper = c(
      decision_tree_model$test_rmse_ci[["97.5%"]],
      decision_tree_model$test_r2_ci[["97.5%"]]
    )
  )
  
  decision_tree_model$train_fit <- tibble(
    .metric = "rsq",
    .estimate = decision_tree_model$train_fit_metrics %>%
      filter(.metric == "rsq") %>%
      pull(.estimate),
    lower = decision_tree_model$train_r2_ci[["2.5%"]],
    upper = decision_tree_model$train_r2_ci[["97.5%"]]
  )
  
  format_metric <- function(data, metric, estimate = ".estimate") {
    
    estimate_value <- data %>%
      filter(.metric == metric) %>%
      pull(all_of(estimate))
    
    lower_value <- data %>%
      filter(.metric == metric) %>%
      pull(lower)
    
    upper_value <- data %>%
      filter(.metric == metric) %>%
      pull(upper)
    
    sprintf(
      "%.2f [%.2f, %.2f]",
      estimate_value,
      lower_value,
      upper_value
    )
  }
  
  r2_train <- format_metric(
    decision_tree_model$train_fit,
    "rsq"
  )
  
  r2_test <- format_metric(
    decision_tree_model$test_fit,
    "rsq"
  )
  
  rmse <- format_metric(
    decision_tree_model$test_fit,
    "rmse"
  )
  
  
  #######
  # 2.5 Create summary output
  #######
  
  decision_tree_model$out_tab <- tibble(
    grp = grp,
    imp = imp_name,
    n = temp$n,
    model = "3_decision_tree",
    r2_train = r2_train,
    r2_test = r2_test,
    rmse = rmse
  )
  
  return(decision_tree_model)
}


#######
# 3. Prepare analysis map for treatment group and imputation condition
#
# grp = treatment group: lit, hit
# imp = imputation: 1 = yes, 0 = no
#######

map_tab <- crossing(
  grp = factor(
    c("lit", "hit"),
    levels = c("lit", "hit")
  ),
  imp = c(1, 0)
) %>%
  mutate(
    name = paste(
      grp,
      if_else(imp == 1, "imp", "noimp"),
      sep = "_"
    )
  )


#######
# 4. Run decision tree models across analysis conditions
#######

decision_tree_model <- pmap(
  list(
    map_tab$grp,
    map_tab$imp
  ),
  fun_dec_tree
) %>%
  set_names(map_tab$name)


#######
# 5. Save final model results
#######

save(
  decision_tree_model,
  file = file.path(
    getwd(),
    "02_models",
    "decision_tree_model.RData"
  )
)