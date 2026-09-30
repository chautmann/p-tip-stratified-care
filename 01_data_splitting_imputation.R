#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 1: Data Splitting and Imputation
#######


#######
# 1. Setup
# Clear workspace, set working directory, and load packages
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("tidymodels")
library("missForest")
library("tidyverse")


#######
# 2. Define function for data splitting and imputation
#######

fun_split_imp <- function(grp, imp){
  # grp <- "lit"
  # imp <- 1
  
  # 2.1 Load data
  df_name <- paste0("ptip_", grp, "_fil0")
  
  load(paste0(getwd(), "/01_data/", df_name,".RData"))
  
  df <- get(df_name)
  
  # 2.2 Remove variables not used in the analysis
  if("ther_duration" %in% names(df)){
    df <- df %>%
      select(-ther_duration)
  }
  
  # 2.3 Prepare analysis variables
  df <- df %>%
    mutate(outcome = cbcl_av_post) %>%
    select(id, age, casmin_km, casmin_kv, medication, mother, number_children,
           sex, single_parent, cbcl_ad_baseline, cbcl_ap_baseline,
           cbcl_av_baseline, outcome) %>%
    mutate_if(is.factor, ~ factor(.))
  
  set.seed(16592)
  
  
  if(imp==1){
    # 2.4 Split and impute data
    split <- initial_split(df,
                           strata = outcome)
    
    train <- training(split)
    
    train_id <- train$id
    
    test <- testing(split)
    
    test_id <- test$id
    
    
    data_imp <- missForest(train %>%
                             select(-id) %>%
                             as.data.frame())
    
    train <- data_imp$ximp %>%
      tibble() %>%
      mutate(id = train_id) %>%
      select(id, everything())
    
    data_imp <- missForest(test %>%
                             select(-id) %>%
                             as.data.frame())
    
    test <- data_imp$ximp %>%
      tibble() %>%
      mutate(id = test_id) %>%
      select(id, everything())
    
    split$data <- bind_rows(train, test) %>%
      arrange(id)
    
    train <- train %>%
      arrange(id) 
    
    test <- test %>%
      arrange(id) 
    
    n <- nrow(split$data)
    
  }else if(imp==0){
    # 2.5 Split complete-case data
    df <- df %>%
      drop_na() 
    
    split <- initial_split(df,
                           strata = outcome)
    
    train <- training(split) 
    
    test <- testing(split)
    
    n <- nrow(split$data)
  }
  
  out <- list()
  
  out$split <- split
  
  out$train <- train
  
  out$test <- test
  
  out$n <- n
  
  return(out)}


#######
# 3. Prepare analysis map for treatment group and imputation condition
#######

map_tab <- crossing(
  tibble(grp = factor(c("lit", "hit"), levels = c("lit", "hit"))),
  tibble(imp = c(1,0))
) %>%
  mutate(nr = c(1:n()),
         imp_name = case_when(imp==1 ~ "imp",
                              imp==0 ~ "noimp")) %>%
  unite("name", c(grp, imp_name), remove = F)


split_ptip <- map(map_tab$nr, function(x){
  # x<-2
  
  grp <- map_tab %>% 
    filter(nr==x) %>%
    pull(grp)
  
  imp <- map_tab %>% 
    filter(nr==x) %>%
    pull(imp)
  
  ret <- fun_split_imp(grp, imp)
  
  return(ret)}) %>%
  set_names(map_tab$name)


#######
# 4. Save split datasets
#######

save(split_ptip, file = paste0(getwd(), "/01_data/split_ptip.RData"))


#######
# 5. Create data for prediction analyses
#######

fun_pred_data <- function(imp) {
  
  names(split_ptip) %>%
    str_subset(paste0("_", imp)) %>%
    map_dfr(\(x) {
      
      split_ptip[[x]]$test %>%
        mutate(
          intervention = str_remove(x, paste0("_", imp, "$"))
        )
      
    }) %>%
    mutate(
      intervention = factor(intervention),
      original_id = id,
      id = paste0(intervention, "_", row_number())
    )
}


prediction_data <- list(
  noimp = fun_pred_data("noimp"),
  imp   = fun_pred_data("imp")
)


#######
# 6. Save prediction data
#######

save(
  prediction_data,
  file = file.path(getwd(), "01_data", "prediction_data.RData")
)