#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 3: Lasso Regression Model
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)
library(vip)

load(file.path(getwd(), "01_data", "split_ptip.RData"))


#######
# 1.1 Variable labels for variable-importance figures
#######

meta_tab <- tibble::tribble(
  ~old,                ~new,                              ~nr,
  "age",               "Age",                               5,
  "casmin_km",         "Education (mother)",                10,
  "casmin_kv",         "Education (father)",                11,
  "medication",        "Medication",                         6,
  "mother",            "Biological mother",                  7,
  "number_children",   "Number of children",                 9,
  "sex",               "Sex",                                4,
  "single_parent",     "Single parent",                      8,
  "cbcl_ad_baseline",  "CBCL Anxious/depressed",             1,
  "cbcl_ap_baseline",  "CBCL Attention problems",            2,
  "cbcl_av_baseline",  "CBCL Aggressive behaviour",          3,
  "outcome",           "Outcome: CBCL Aggressive Behaviour", 0
)


#######
# 2. Define modeling function
#######

fun_lasso <- function(grp, imp) {
  
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
  lasso_model <- list()
  
  
  #######
  # 2.1 Data preprocessing
  #######
  
  lasso_model$recipe <- recipe(outcome ~ ., data = train) %>%
    update_role(id, new_role = "id") %>%
    step_zv(all_predictors()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  lasso_model$workflow <- workflow() %>%
    add_recipe(lasso_model$recipe)
  
  
  #######
  # 2.2 Tune lasso parameters
  #######
  
  set.seed(201157)
  
  lasso_model$cv_folds <- vfold_cv(train)
  
  lasso_model$tune_spec <- linear_reg(
    penalty = tune(),
    mixture = 1
  ) %>%
    set_engine("glmnet")
  
  lasso_model$lambda_grid <- grid_max_entropy(
    penalty(),
    size = 200
  )
  
  set.seed(2390)
  
  lasso_model$lasso_grid <- tune_grid(
    lasso_model$workflow %>%
      add_model(lasso_model$tune_spec),
    resamples = lasso_model$cv_folds,
    grid = lasso_model$lambda_grid
  )
  
  
  #######
  # 2.3 Finalize model
  #######
  
  lasso_model$lowest_rmse <- lasso_model$lasso_grid %>%
    select_best(metric = "rmse")
  
  lasso_model$final_lasso <- finalize_workflow(
    lasso_model$workflow %>%
      add_model(lasso_model$tune_spec),
    lasso_model$lowest_rmse
  )
  
  
  #######
  # 2.4 Fit final model and evaluate test-set performance
  #######
  
  lasso_model$fit <- last_fit(
    lasso_model$final_lasso,
    split
  )
  
  lasso_model$fin_fit <- lasso_model$fit %>%
    collect_metrics() %>%
    select(.metric, .estimate)
  
  
  #######
  # 2.5 Training-set performance
  #######
  
  lasso_model$train_fit <- lasso_model$final_lasso %>%
    fit(train)
  
  # Generate predictions for the training set
  lasso_model$train_pred <- predict(
    lasso_model$train_fit,
    new_data = train
  ) %>%
    bind_cols(train) %>%
    select(outcome, .pred)
  
  # Calculate training-set metrics
  lasso_model$train_fit_metrics <- lasso_model$train_pred %>%
    metrics(
      truth = outcome,
      estimate = .pred
    )
  
  
  #######
  # 2.6 Bootstrap 95% confidence intervals for model performance
  #######
  
  set.seed(79110)
  
  # Training-set bootstrap
  lasso_model$train_r2_boot <- bootstraps(
    lasso_model$train_pred,
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
  
  lasso_model$train_r2_ci <- quantile(
    lasso_model$train_r2_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  # Test-set predictions
  lasso_model$test_pred <- lasso_model$fit %>%
    collect_predictions() %>%
    select(outcome, .pred)
  
  # Calculate R2 and RMSE for each bootstrap test sample
  lasso_model$test_metrics_boot <- bootstraps(
    lasso_model$test_pred,
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
  lasso_model$test_r2_ci <- quantile(
    lasso_model$test_metrics_boot$r2,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  # 95% CI for test-set RMSE
  lasso_model$test_rmse_ci <- quantile(
    lasso_model$test_metrics_boot$rmse,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  
  #######
  # 2.7 Save final fitted model
  #######
  
  save_model <- lasso_model$fit %>%
    extract_workflow(.)
  
  model_name <- paste0(
    "lasso_",
    grp,
    "_",
    imp_name
  )
  
  assign(model_name, save_model)
  
  save(
    list = model_name,
    file = file.path(
      getwd(),
      "02_models",
      paste0(model_name, ".RData")
    )
  )
  
  
  #######
  # 2.8 Figure S3:
  # Observed versus predicted values for the final lasso models
  # in the held-out test sample
  #######
  
  # Use held-out test-set predictions
  pred_graph_df <- lasso_model$test_pred
  
  
  # Correlation -------------------------------------------------------------
  
  cor_fit <- cor.test(
    pred_graph_df$outcome,
    pred_graph_df$.pred
  )
  
  p_label <- case_when(
    cor_fit$p.value < .001 ~ "< .001",
    cor_fit$p.value < .01  ~ "< .01",
    cor_fit$p.value < .05  ~ "< .05",
    TRUE                   ~ paste0(
      "= ",
      sprintf("%.3f", cor_fit$p.value)
    )
  )
  
  annotate_label <- sprintf(
    "r = %.2f, p %s",
    unname(cor_fit$estimate),
    p_label
  )
  
  
  # Annotation position -----------------------------------------------------
  
  x_pos <- min(
    pred_graph_df$outcome,
    na.rm = TRUE
  ) +
    0.10 * diff(
      range(
        pred_graph_df$outcome,
        na.rm = TRUE
      )
    )
  
  y_pos <- max(
    pred_graph_df$.pred,
    na.rm = TRUE
  ) -
    0.05 * diff(
      range(
        pred_graph_df$.pred,
        na.rm = TRUE
      )
    )
  
  
  # Plot --------------------------------------------------------------------
  
  gg_plot <- ggplot(
    pred_graph_df,
    aes(
      x = outcome,
      y = .pred
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.4
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "#377eb8",
      linetype = "dashed",
      linewidth = 1.75
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = annotate_label,
      hjust = 0,
      size = 7.5
    ) +
    labs(
      x = "Observed",
      y = "Predicted"
    ) +
    coord_cartesian(
      xlim = c(0, 38),
      ylim = c(0, 38)
    ) +
    theme(
      plot.title = element_text(
        size = rel(2.5),
        hjust = 0.5
      ),
      axis.title = element_text(
        size = rel(2.25)
      ),
      axis.text = element_text(
        size = rel(2)
      )
    )
  
  
  # Save --------------------------------------------------------------------
  
  ggsave(
    file.path(
      "03_output",
      "figures",
      paste0(
        "figureS3_",
        grp,
        "_",
        imp_name,
        ".png"
      )
    ),
    plot = gg_plot,
    width = 16,
    height = 13.4
  )
  
  
  #######
  # 2.9 Figure S4:
  # Variable importance patterns for the final lasso regression models
  #######
  
  gg_plot <- lasso_model$final_lasso %>%
    fit(train) %>%
    extract_fit_parsnip() %>%
    vi(
      lambda = lasso_model$lowest_rmse$penalty
    ) %>%
    left_join(
      meta_tab %>%
        mutate(
          old = case_when(
            old == "medication"    ~ "medication_yes",
            old == "mother"        ~ "mother_biological.mother",
            old == "sex"           ~ "sex_female",
            old == "single_parent" ~ "single_parent_yes",
            TRUE                   ~ old
          )
        ),
      join_by(Variable == old)
    ) %>%
    mutate(
      Variable = new
    ) %>%
    select(
      -new,
      -nr
    ) %>%
    mutate(
      Importance = abs(Importance),
      Variable = forcats::fct_reorder(
        Variable,
        Importance
      )
    ) %>%
    ggplot(
      aes(
        x = Importance,
        y = Variable,
        fill = Sign
      )
    ) +
    geom_col() +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    theme(
      plot.title = element_text(
        size = rel(2.5),
        hjust = 0.5
      ),
      axis.title = element_text(
        size = rel(2.25)
      ),
      axis.text = element_text(
        size = rel(2)
      ),
      legend.title = element_text(
        size = rel(2.25)
      ),
      legend.text = element_text(
        size = rel(2)
      ),
      legend.position = "bottom"
    ) +
    scale_fill_manual(
      values = c(
        "POS" = "#377eb8",
        "NEG" = "#e41a1c"
      )
    )
  
  figure_name <- paste0(
    "figureS4_",
    grp,
    "_",
    imp_name
  )
  
  ggsave(
    file.path(
      getwd(),
      "03_output",
      "figures",
      paste0(
        figure_name,
        ".png"
      )
    ),
    plot = gg_plot,
    width = 16,
    height = 13.4
  )
  
  
  #######
  # 2.10 Format model performance for output table
  #######
  
  lasso_model$fin_fit <- tibble(
    .metric = c(
      "rmse",
      "rsq"
    ),
    .estimate = c(
      rmse_vec(
        lasso_model$test_pred$outcome,
        lasso_model$test_pred$.pred
      ),
      rsq_vec(
        lasso_model$test_pred$outcome,
        lasso_model$test_pred$.pred
      )
    ),
    lower = c(
      lasso_model$test_rmse_ci[["2.5%"]],
      lasso_model$test_r2_ci[["2.5%"]]
    ),
    upper = c(
      lasso_model$test_rmse_ci[["97.5%"]],
      lasso_model$test_r2_ci[["97.5%"]]
    )
  )
  
  lasso_model$train_fit <- tibble(
    .metric = "rsq",
    .estimate = lasso_model$train_fit_metrics %>%
      filter(.metric == "rsq") %>%
      pull(.estimate),
    lower = lasso_model$train_r2_ci[["2.5%"]],
    upper = lasso_model$train_r2_ci[["97.5%"]]
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
    lasso_model$train_fit,
    "rsq"
  )
  
  r2_test <- format_metric(
    lasso_model$fin_fit,
    "rsq"
  )
  
  rmse <- format_metric(
    lasso_model$fin_fit,
    "rmse"
  )
  
  
  #######
  # 2.11 Create summary output
  #######
  
  lasso_model$out_tab <- tibble(
    grp = grp,
    imp = imp_name,
    n = temp$n,
    model = "2_lasso",
    r2_train = r2_train,
    r2_test = r2_test,
    rmse = rmse
  )
  
  return(lasso_model)
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
      if_else(
        imp == 1,
        "imp",
        "noimp"
      ),
      sep = "_"
    )
  )


#######
# 4. Run lasso models across analysis conditions
#######

lasso_model <- pmap(
  list(
    map_tab$grp,
    map_tab$imp
  ),
  fun_lasso
) %>%
  set_names(map_tab$name)


#######
# 5. Save final model results
#######

save(
  lasso_model,
  file = file.path(
    getwd(),
    "02_models",
    "lasso_model.RData"
  )
)