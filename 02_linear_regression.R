#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 2: Linear Regression Model
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("tidymodels")

load(file.path(getwd(), "01_data", "split_ptip.RData"))


#######
# 2. Define modeling function
#######

fun_lm <- function(grp, imp){
  
  set.seed(24120000)
  
  # grp <- "lit"
  # imp <- 0
  
  # Define imputation condition and dataset name
  imp_name <- if_else(imp == 1, "imp", "noimp")
  
  df_name <- paste0(grp, "_", imp_name)
  
  # Extract data split and training/test data
  temp <- split_ptip[[df_name]]
  
  split <- temp$split
  train <- temp$train
  
  # Initialize model output object
  lm_model <- list()
  
  
  #######
  # 2.1 Data preprocessing
  #######
  
  lm_model$recipe <- recipe(outcome ~ ., data = train) %>%
    step_rm(id) %>%
    step_zv(all_predictors()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  
  #######
  # 2.2 Specify linear regression model
  #######
  
  lm_model$model <- linear_reg() %>%
    set_engine("lm")
  
  
  #######
  # 2.3 Create modeling workflow
  #######
  
  lm_model$workflow <- workflow() %>%
    add_model(lm_model$model) %>%
    add_recipe(lm_model$recipe)
  
  
  #######
  # 2.4 Fit model and evaluate test-set performance
  #######
  
  lm_model$fit <- last_fit(
    lm_model$workflow,
    split
  )
  
  # Extract predictions for the test set
  lm_model$pred <- lm_model$fit %>%
    collect_predictions()
  
  
  #######
  # 2.5 Bootstrap confidence intervals for test-set performance
  #######
  
  set.seed(25653)
  
  boot_metrics <- bootstraps(
    lm_model$pred,
    times = 1000
  ) %>%
    mutate(
      metrics = map(
        splits,
        ~ tibble(
          rmse = rmse_vec(
            analysis(.x)$outcome,
            analysis(.x)$.pred
          ),
          r2 = rsq_vec(
            analysis(.x)$outcome,
            analysis(.x)$.pred
          )
        )
      )
    ) %>%
    select(metrics) %>%
    unnest(metrics)
  
  
  # Calculate point estimates and 95% bootstrap confidence intervals
  lm_model$fin_fit <- tibble(
    .metric = c("rmse", "rsq"),
    .estimate = c(
      rmse_vec(lm_model$pred$outcome, lm_model$pred$.pred),
      rsq_vec(lm_model$pred$outcome, lm_model$pred$.pred)
    ),
    lower = c(
      quantile(boot_metrics$rmse, 0.025),
      quantile(boot_metrics$r2, 0.025)
    ),
    upper = c(
      quantile(boot_metrics$rmse, 0.975),
      quantile(boot_metrics$r2, 0.975)
    )
  )
  
  
  #######
  # 2.6 Training-set performance
  #######
  
  # Fit the model to the full training set
  lm_model$train_model <- lm_model$workflow %>%
    fit(train)
  
  # Generate predictions for the training set
  lm_model$train_pred <- predict(
    lm_model$train_model,
    new_data = train
  ) %>%
    bind_cols(train) %>%
    select(outcome, .pred)
  
  # Calculate training-set R2
  lm_model$train_fit_metrics <- lm_model$train_pred %>%
    metrics(
      truth = outcome,
      estimate = .pred
    )
  
  
  #######
  # 2.7 Bootstrap confidence interval for training-set R2
  #######
  
  set.seed(25653)
  
  lm_model$train_r2_boot <- bootstraps(
    lm_model$train_pred,
    times = 1000
  ) %>%
    mutate(
      r2 = map_dbl(
        splits,
        ~ rsq_vec(
          truth = analysis(.x)$outcome,
          estimate = analysis(.x)$.pred
        )
      )
    )
  
  lm_model$train_r2_ci <- quantile(
    lm_model$train_r2_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  lm_model$train_fit <- tibble(
    .metric = "rsq",
    .estimate = lm_model$train_fit_metrics %>%
      filter(.metric == "rsq") %>%
      pull(.estimate),
    lower = lm_model$train_r2_ci[["2.5%"]],
    upper = lm_model$train_r2_ci[["97.5%"]]
  )
  
  
  #######
  # 2.8 Format model performance for output table
  #######
  
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
  
  r2_train <- format_metric(lm_model$train_fit, "rsq")
  r2_test  <- format_metric(lm_model$fin_fit, "rsq")
  rmse     <- format_metric(lm_model$fin_fit, "rmse")
  
  
  #######
  # 2.9 Create summary output
  #######
  
  lm_model$out_tab <- tibble(
    grp = grp,
    imp = imp_name,
    n = temp$n,
    model = "1_lm",
    r2_train = r2_train,
    r2_test = r2_test,
    rmse = rmse
  )
  
  return(lm_model)
}


#######
# 3. Prepare analysis map for treatment group and imputation condition
#
# grp = treatment group: lit, hit
# imp = imputation: 1 = yes, 0 = no
#######

map_tab <- crossing(
  grp = factor(c("lit", "hit"), levels = c("lit", "hit")),
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
# 4. Run linear regression models across analysis conditions
#######

lm_model <- pmap(
  list(map_tab$grp, map_tab$imp),
  fun_lm
) %>%
  set_names(map_tab$name)


#######
# 5. Save model results
#######

save(
  lm_model,
  file = file.path(
    getwd(),
    "02_models",
    "lm_model.RData"
  )
)