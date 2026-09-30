#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 5: Random Forest Model
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

fun_ran_for <- function(grp, imp) {
  
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
  random_forest_model <- list()
  
  
  #######
  # 2.1 Specify random forest model
  #######
  
  cores <- parallel::detectCores()
  
  random_forest_model$model <- rand_forest(
    mtry = tune(),
    min_n = tune(),
    trees = tune()
  ) %>%
    set_engine(
      "ranger",
      num.threads = cores
    ) %>%
    set_mode("regression")
  
  random_forest_model$recipe <- recipe(
    outcome ~ .,
    data = train
  )
  
  random_forest_model$workflow <- workflow() %>%
    add_model(random_forest_model$model) %>%
    add_recipe(random_forest_model$recipe)
  
  
  #######
  # 2.2 Tune random forest parameters
  #######
  
  random_forest_model$mtry_final <- finalize(
    mtry(),
    select(train, -outcome)
  )
  
  random_forest_model$tree_grid <- grid_max_entropy(
    mtry = random_forest_model$mtry_final,
    min_n(),
    trees(),
    size = 200
  )
  
  # Create cross-validation folds
  set.seed(201157)
  
  random_forest_model$folds <- vfold_cv(
    train,
    strata = outcome
  )
  
  # Tune model parameters
  set.seed(11901)
  
  random_forest_model$tuning <-
    random_forest_model$workflow %>%
    tune_grid(
      random_forest_model$folds,
      grid = random_forest_model$tree_grid,
      control = control_grid(
        save_pred = TRUE
      ),
      metrics = metric_set(
        rmse,
        rsq
      )
    )
  
  
  #######
  # 2.3 Finalize model
  #######
  
  random_forest_model$lowest_rmse <-
    random_forest_model$tuning %>%
    select_best(metric = "rmse")
  
  random_forest_model$update_model <- rand_forest(
    mtry = random_forest_model$lowest_rmse$mtry,
    min_n = random_forest_model$lowest_rmse$min_n,
    trees = random_forest_model$lowest_rmse$trees
  ) %>%
    set_engine(
      "ranger",
      num.threads = cores,
      importance = "impurity"
    ) %>%
    set_mode("regression")
  
  random_forest_model$final_random_forest <-
    random_forest_model$workflow %>%
    update_model(random_forest_model$update_model)
  
  
  #######
  # 2.4 Fit final model and evaluate test-set performance
  #######
  
  set.seed(345)
  
  random_forest_model$fit <-
    random_forest_model$final_random_forest %>%
    last_fit(split)
  
  random_forest_model$fin_fit <-
    random_forest_model$fit %>%
    collect_metrics()
  
  
  #######
  # 2.5 Training-set performance
  #######
  
  # Use the same seed as for the final test-set fit so that the
  # separate training-set random forest is reproducible
  set.seed(345)
  
  random_forest_model$train_fit <-
    random_forest_model$final_random_forest %>%
    fit(train)
  
  # Generate predictions for the training set
  random_forest_model$train_pred <- predict(
    random_forest_model$train_fit,
    new_data = train
  ) %>%
    bind_cols(train) %>%
    select(outcome, .pred)
  
  # Calculate training-set metrics
  random_forest_model$train_fit_metrics <-
    random_forest_model$train_pred %>%
    metrics(
      truth = outcome,
      estimate = .pred
    )
  
  
  #######
  # 2.6 Bootstrap 95% confidence intervals for model performance
  #######
  
  set.seed(79110)
  
  # Training-set bootstrap
  random_forest_model$train_r2_boot <- bootstraps(
    random_forest_model$train_pred,
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
  
  random_forest_model$train_r2_ci <- quantile(
    random_forest_model$train_r2_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  # Test-set predictions
  random_forest_model$test_pred <-
    random_forest_model$fit %>%
    collect_predictions() %>%
    select(outcome, .pred)
  
  # Calculate R2 and RMSE for each bootstrap test sample
  random_forest_model$test_metrics_boot <- bootstraps(
    random_forest_model$test_pred,
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
  random_forest_model$test_r2_ci <- quantile(
    random_forest_model$test_metrics_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  # 95% CI for test-set RMSE
  random_forest_model$test_rmse_ci <- quantile(
    random_forest_model$test_metrics_boot$rmse,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  #######
  # 2.7 Format model performance for output table
  #######
  
  random_forest_model$fin_fit <- tibble(
    .metric = c(
      "rmse",
      "rsq"
    ),
    .estimate = c(
      rmse_vec(
        random_forest_model$test_pred$outcome,
        random_forest_model$test_pred$.pred
      ),
      rsq_vec(
        random_forest_model$test_pred$outcome,
        random_forest_model$test_pred$.pred
      )
    ),
    lower = c(
      random_forest_model$test_rmse_ci[["2.5%"]],
      random_forest_model$test_r2_ci[["2.5%"]]
    ),
    upper = c(
      random_forest_model$test_rmse_ci[["97.5%"]],
      random_forest_model$test_r2_ci[["97.5%"]]
    )
  )
  
  random_forest_model$train_fit <- tibble(
    .metric = "rsq",
    .estimate = random_forest_model$train_fit_metrics %>%
      filter(.metric == "rsq") %>%
      pull(.estimate),
    lower = random_forest_model$train_r2_ci[["2.5%"]],
    upper = random_forest_model$train_r2_ci[["97.5%"]]
  )
  
  format_metric <- function(
    data,
    metric,
    estimate = ".estimate"
  ) {
    
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
    random_forest_model$train_fit,
    "rsq"
  )
  
  r2_test <- format_metric(
    random_forest_model$fin_fit,
    "rsq"
  )
  
  rmse <- format_metric(
    random_forest_model$fin_fit,
    "rmse"
  )
  
  
  #######
  # 2.8 Create summary output
  #######
  
  random_forest_model$out_tab <- tibble(
    grp = grp,
    imp = imp_name,
    n = temp$n,
    model = "4_random_forest",
    r2_train = r2_train,
    r2_test = r2_test,
    rmse = rmse
  )
  
  return(random_forest_model)
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
# 4. Run random forest models across analysis conditions
#######

random_forest_model <- pmap(
  list(
    map_tab$grp,
    map_tab$imp
  ),
  fun_ran_for
) %>%
  set_names(map_tab$name)


#######
# 5. Save final model results
#######

save(
  random_forest_model,
  file = file.path(
    getwd(),
    "02_models",
    "random_forest_model.RData"
  )
)