#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 6: GBM Model
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)
library(finetune)

load(file.path(getwd(), "01_data", "split_ptip.RData"))


#######
# 2. Define modeling function
#######

fun_gbm <- function(grp, imp) {
  
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
  gbm_model <- list()
  
  
  #######
  # 2.1 Data preprocessing and model specification
  #######
  
  gbm_model$recipe <- recipe(outcome ~ ., data = train) %>%
    step_dummy(all_nominal_predictors())
  
  gbm_model$model <- boost_tree(
    trees = tune(),
    min_n = tune(),
    mtry = tune(),
    learn_rate = 0.1
  ) %>%
    set_engine("xgboost") %>%
    set_mode("regression")
  
  gbm_model$workflow <- workflow(
    gbm_model$recipe,
    gbm_model$model
  )
  
  
  #######
  # 2.2 Tune GBM parameters
  #######
  
  # Create cross-validation folds
  set.seed(24120000)
  
  gbm_model$folds <- vfold_cv(train)
  
  # Tune model parameters
  set.seed(2390)
  
  gbm_model$tuning <- tune_race_anova(
    gbm_model$workflow,
    resamples = gbm_model$folds,
    grid = 35,
    control = control_race(verbose_elim = TRUE)
  )
  
  
  #######
  # 2.3 Finalize model
  #######
  
  gbm_model$final_gbm <- gbm_model$workflow %>%
    finalize_workflow(
      select_best(
        gbm_model$tuning,
        metric = "rmse"
      )
    )
  
  
  #######
  # 2.4 Fit final model and evaluate test-set performance
  #######
  
  set.seed(345)
  
  gbm_model$fit <- gbm_model$final_gbm %>%
    last_fit(split)
  
  gbm_model$fin_fit <- gbm_model$fit %>%
    collect_metrics()
  
  
  #######
  # 2.5 Training-set performance
  #######
  
  # Use the same seed as for the final test-set fit so that the
  # separate training-set GBM fit is reproducible
  set.seed(345)
  
  gbm_model$train_fit <- gbm_model$final_gbm %>%
    fit(train)
  
  # Generate predictions for the training set
  gbm_model$train_pred <- predict(
    gbm_model$train_fit,
    new_data = train
  ) %>%
    bind_cols(train) %>%
    select(outcome, .pred)
  
  # Calculate training-set metrics
  gbm_model$train_fit_metrics <- gbm_model$train_pred %>%
    metrics(
      truth = outcome,
      estimate = .pred
    )
  
  
  #######
  # 2.6 Bootstrap 95% confidence intervals for model performance
  #######
  
  set.seed(79110)
  
  # Training-set bootstrap
  gbm_model$train_r2_boot <- bootstraps(
    gbm_model$train_pred,
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
  
  gbm_model$train_r2_ci <- quantile(
    gbm_model$train_r2_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  # Test-set predictions
  gbm_model$test_pred <- gbm_model$fit %>%
    collect_predictions() %>%
    select(outcome, .pred)
  
  # Calculate R2 and RMSE for each bootstrap test sample
  gbm_model$test_metrics_boot <- bootstraps(
    gbm_model$test_pred,
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
  gbm_model$test_r2_ci <- quantile(
    gbm_model$test_metrics_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  # 95% CI for test-set RMSE
  gbm_model$test_rmse_ci <- quantile(
    gbm_model$test_metrics_boot$rmse,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  #######
  # 2.7 Format model performance for output table
  #######
  
  gbm_model$fin_fit <- tibble(
    .metric = c("rmse", "rsq"),
    .estimate = c(
      rmse_vec(
        gbm_model$test_pred$outcome,
        gbm_model$test_pred$.pred
      ),
      rsq_vec(
        gbm_model$test_pred$outcome,
        gbm_model$test_pred$.pred
      )
    ),
    lower = c(
      gbm_model$test_rmse_ci[["2.5%"]],
      gbm_model$test_r2_ci[["2.5%"]]
    ),
    upper = c(
      gbm_model$test_rmse_ci[["97.5%"]],
      gbm_model$test_r2_ci[["97.5%"]]
    )
  )
  
  gbm_model$train_fit <- tibble(
    .metric = "rsq",
    .estimate = gbm_model$train_fit_metrics %>%
      filter(.metric == "rsq") %>%
      pull(.estimate),
    lower = gbm_model$train_r2_ci[["2.5%"]],
    upper = gbm_model$train_r2_ci[["97.5%"]]
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
    gbm_model$train_fit,
    "rsq"
  )
  
  r2_test <- format_metric(
    gbm_model$fin_fit,
    "rsq"
  )
  
  rmse <- format_metric(
    gbm_model$fin_fit,
    "rmse"
  )
  
  
  #######
  # 2.8 Create summary output
  #######
  
  gbm_model$out_tab <- tibble(
    grp = grp,
    imp = imp_name,
    n = temp$n,
    model = "5_gbm",
    r2_train = r2_train,
    r2_test = r2_test,
    rmse = rmse
  )
  
  return(gbm_model)
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
# 4. Run GBM models across analysis conditions
#######

gbm_model <- pmap(
  list(
    map_tab$grp,
    map_tab$imp
  ),
  fun_gbm
) %>%
  set_names(map_tab$name)


#######
# 5. Save final model results
#######

save(
  gbm_model,
  file = file.path(
    getwd(),
    "02_models",
    "gbm_model.RData"
  )
)