#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 12: Predictive Performance of Lasso Models in HIT Across
# Subsamples With Varying Therapy Lengths
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidymodels)
library(lubridate)
library(writexl)

load(
  file.path(
    "01_data",
    "ptip_hit_fil0.RData"
  )
)


#######
# 2. Define lasso function
#######

fun_lasso <- function(imp, month) {
  
  # imp <- 0
  # month <- 12
  
  
  #######
  # 2.1 Preprocess data
  #######
  
  df <- ptip_hit_fil0 %>%
    mutate(
      outcome = cbcl_av_post
    ) %>%
    select(
      id,
      ther_duration,
      age,
      casmin_km,
      casmin_kv,
      medication,
      mother,
      number_children,
      sex,
      single_parent,
      cbcl_ad_baseline,
      cbcl_ap_baseline,
      cbcl_av_baseline,
      outcome
    ) %>%
    mutate_if(
      is.factor,
      ~ factor(.)
    )
  
  
  #######
  # 2.2 Restrict sample by maximum therapy duration
  #######
  
  dur <- duration(
    month,
    units = "month"
  )
  
  dur_weeks <- as.numeric(
    dur,
    units = "weeks"
  )
  
  df <- df %>%
    filter(
      ther_duration <= dur_weeks
    ) %>%
    select(
      -ther_duration
    )
  
  
  #######
  # 2.3 Split data and handle missing values
  #######
  
  set.seed(16592)
  
  split <- initial_split(
    df,
    strata = outcome
  )
  
  train <- training(split)
  test <- testing(split)
  
  
  if (imp == 1) {
    
    split <- initial_split(
      df,
      strata = outcome
    )
    
    train <- training(split)
    train_id <- train$id
    
    test <- testing(split)
    test_id <- test$id
    
    
    data_imp <- missForest::missForest(
      train %>%
        select(-id) %>%
        as.data.frame()
    )
    
    train <- data_imp$ximp %>%
      tibble() %>%
      mutate(
        id = train_id
      ) %>%
      select(
        id,
        everything()
      )
    
    
    data_imp <- missForest::missForest(
      test %>%
        select(-id) %>%
        as.data.frame()
    )
    
    test <- data_imp$ximp %>%
      tibble() %>%
      mutate(
        id = test_id
      ) %>%
      select(
        id,
        everything()
      )
    
    
    split$data <- bind_rows(
      train,
      test
    ) %>%
      arrange(id)
    
    train <- train %>%
      arrange(id)
    
    test <- test %>%
      arrange(id)
    
    n <- nrow(split$data)
    
  } else if (imp == 0) {
    
    df <- df %>%
      drop_na()
    
    split <- initial_split(
      df,
      strata = outcome
    )
    
    train <- training(split)
    test <- testing(split)
    
    n <- nrow(split$data)
  }
  
  
  #######
  # 2.4 Recipe
  #######
  
  rec <- recipe(
    outcome ~ .,
    data = train
  ) %>%
    update_role(
      id,
      new_role = "id"
    ) %>%
    step_zv(
      all_predictors()
    ) %>%
    step_dummy(
      all_nominal_predictors()
    ) %>%
    step_normalize(
      all_numeric_predictors()
    )
  
  
  #######
  # 2.5 Specify workflow
  #######
  
  wf <- workflow() %>%
    add_recipe(rec)
  
  
  #######
  # 2.6 Tune lasso parameters using 10-fold cross-validation
  #######
  
  set.seed(201157)
  
  cv_folds <- vfold_cv(
    train,
    v = 10
  )
  
  
  tune_spec <- linear_reg(
    penalty = tune(),
    mixture = 1
  ) %>%
    set_engine("glmnet")
  
  
  lambda_grid <- grid_max_entropy(
    penalty(),
    size = 200
  )
  
  
  doParallel::registerDoParallel()
  
  
  set.seed(2390)
  
  lasso_grid <- tune_grid(
    wf %>%
      add_model(tune_spec),
    resamples = cv_folds,
    grid = lambda_grid
  )
  
  
  #######
  # 2.7 Finalize and evaluate model
  #######
  
  lowest_rmse <- lasso_grid %>%
    select_best(
      metric = "rmse"
    )
  
  
  final_lasso <- finalize_workflow(
    wf %>%
      add_model(tune_spec),
    lowest_rmse
  )
  
  
  fit <- last_fit(
    final_lasso,
    split
  )
  
  
  test_pred <- fit %>%
    collect_predictions() %>%
    select(
      outcome,
      .pred
    )
  
  
  # Point estimates
  test_r2 <- rsq_vec(
    truth = test_pred$outcome,
    estimate = test_pred$.pred
  )
  
  test_rmse <- rmse_vec(
    truth = test_pred$outcome,
    estimate = test_pred$.pred
  )
  
  
  #######
  # 2.8 Bootstrap 95% confidence intervals for test-set performance
  #######
  
  set.seed(79110)
  
  test_metrics_boot <- bootstraps(
    test_pred,
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
  
  
  test_r2_ci <- quantile(
    test_metrics_boot$r2,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )
  
  
  test_rmse_ci <- quantile(
    test_metrics_boot$rmse,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )
  
  
  #######
  # 2.9 Create output
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  
  out_tab <- tibble(
    imp = imp_name,
    month = month,
    n = n,
    
    r2 = test_r2,
    r2_lower = test_r2_ci[["2.5%"]],
    r2_upper = test_r2_ci[["97.5%"]],
    
    rmse = test_rmse,
    rmse_lower = test_rmse_ci[["2.5%"]],
    rmse_upper = test_rmse_ci[["97.5%"]]
  )
  
  
  return(out_tab)
}


#######
# 3. Run analyses across maximum therapy durations
#
# imp = imputation:
#   0 = complete-case analysis
#   1 = single imputation
#######

map_tab <- crossing(
  imp = 0,
  month = c(
    6,
    9,
    12,
    15,
    18,
    21,
    24,
    27,
    30,
    33,
    36
  )
)


tab <- pmap(
  list(
    map_tab$imp,
    map_tab$month
  ),
  fun_lasso
) %>%
  reduce(bind_rows)


tab_unimp <- tab %>%
  filter(
    imp == "noimp"
  ) %>%
  select(
    -imp
  )


#######
# 4. Format output table
#######

tab_output <- tab_unimp %>%
  mutate(
    r2_ci = sprintf(
      "%.2f [%.2f, %.2f]",
      r2,
      r2_lower,
      r2_upper
    ),
    
    rmse_ci = sprintf(
      "%.2f [%.2f, %.2f]",
      rmse,
      rmse_lower,
      rmse_upper
    )
  ) %>%
  select(
    month,
    n,
    r2,
    r2_lower,
    r2_upper,
    rmse,
    rmse_lower,
    rmse_upper,
    r2_ci,
    rmse_ci
  )


write_xlsx(
  tab_output,
  file.path(
    "03_output",
    "hit_duration_lasso.xlsx"
  )
)


#######
# 5. Figure S2
# Distribution of therapy duration in the HIT primary analytic sample
#######

dur_36 <- duration(
  36,
  units = "month"
)

dur_36_weeks <- as.numeric(
  dur_36,
  units = "weeks"
)


df_gg <- ptip_hit_fil0 %>%
  mutate(
    outcome = cbcl_av_post
  ) %>%
  select(
    id,
    ther_duration,
    age,
    casmin_km,
    casmin_kv,
    medication,
    mother,
    number_children,
    sex,
    single_parent,
    cbcl_ad_baseline,
    cbcl_ap_baseline,
    cbcl_av_baseline,
    outcome
  ) %>%
  filter(
    ther_duration <= dur_36_weeks
  ) %>%
  drop_na() %>%
  mutate(
    ther_duration = ther_duration / 4.345
  )


stats <- df_gg %>%
  summarise(
    mean_duration = mean(
      ther_duration
    ),
    sd_duration = sd(
      ther_duration
    )
  )


gg_plot <- ggplot(
  df_gg,
  aes(
    x = ther_duration
  )
) +
  geom_histogram(
    binwidth = 1,
    color = "grey30",
    fill = "grey85"
  ) +
  geom_vline(
    xintercept = stats$mean_duration,
    linewidth = 0.8
  ) +
  labs(
    x = "Therapy Duration (Months)",
    y = "Frequency"
  )


ggsave(
  file.path(
    "03_output",
    "figures",
    "figureS2.png"
  ),
  plot = gg_plot,
  width = 16,
  height = 9
)