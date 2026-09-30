#######
# This R script was developed by Felix Oswald and Christopher Hautmann
# as part of the reproducible analysis workflow for the manuscript
# “Stratified Care for Child and Adolescent Conduct Problems:
# A Data-Driven Decision Rule for Initial Treatment Intensity.”
#
# Script 11: Figures
#######


#######
# 1. Setup
# Clear workspace, set working directory, load packages, and load data
#######

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

load(file.path("01_data", "test_pred_noimp.RData"))
load(file.path("01_data", "test_pred_imp.RData"))
load(file.path("01_data", "roc_noimp.RData"))
load(file.path("01_data", "roc_imp.RData"))

# Ensure output directory exists
dir.create(
  file.path("03_output", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)


#######
# 2. Figure 1
# Distribution of Personalized Advantage Index (PAI) Scores
#######

plot_pai_distribution <- function(imp) {
  
  # imp <- 0
  
  
  #######
  # 2.1 Set imputation condition and extract data
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  df_name <- paste0(
    "test_pred_",
    imp_name
  )
  
  pred <- get(df_name)
  
  
  #######
  # 2.2 Prepare data for plot
  #######
  
  df_gg <- pred %>%
    filter(
      intervention != "tau"
    ) %>%
    select(
      pai
    ) %>%
    mutate(
      pai = round(pai * 2) / 2,
      grid_color = "Personalized Advantage Index"
    )
  
  df_gg <- df_gg %>%
    count(
      pai,
      grid_color
    ) %>%
    mutate(
      pct = prop.table(n) * 100
    )
  
  
  #######
  # 2.3 Create figure
  #######
  
  gg_plot <- df_gg %>%
    ggplot(
      aes(
        x = pai,
        y = pct
      )
    ) +
    geom_bar(
      stat = "identity",
      position = "identity",
      width = 0.4
    ) +
    labs(
      y = "% of Sample",
      x = "Personalized Advantage Index"
    ) +
    scale_x_continuous(
      limits = c(
        round(min(df_gg$pai)) - 1,
        round(max(df_gg$pai)) + 1
      ),
      breaks = seq(
        round(min(df_gg$pai)) - 1,
        round(max(df_gg$pai)) + 1,
        by = 1
      )
    ) +
    scale_y_continuous(
      limits = c(
        0,
        max(df_gg$pct) + 1
      ),
      labels = scales::percent_format(
        scale = 1
      ),
      breaks = seq(
        round(min(df_gg$pct)),
        round(max(df_gg$pct)) + 1,
        by = 2
      )
    ) +
    theme(
      plot.title = element_text(
        size = rel(2.5),
        hjust = 0.5
      ),
      axis.title = element_text(
        size = rel(2)
      ),
      axis.text = element_text(
        size = rel(1.75)
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = rel(1.75)
      )
    )
  
  
  #######
  # 2.4 Save figure
  #######
  
  gg_name <- paste0(
    "figure1_",
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
    plot = gg_plot,
    width = 16,
    height = 13.4
  )
}


#######
# Run Figure 1
#######

plot_pai_distribution(0)

plot_pai_distribution(1)


#######
# 3. Figure 3
# Relationship Between the Personalized Advantage Index (PAI)
# and the Mean Prognostic Index (MPI)
#######

plot_pai_mpi <- function(imp) {
  
  # imp <- 0
  
  
  #######
  # 3.1 Set imputation condition and extract data
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  df_name <- paste0(
    "test_pred_",
    imp_name
  )
  
  pred <- get(df_name)
  
  
  #######
  # 3.2 Calculate correlation
  #######
  
  cor_fit <- cor.test(
    pred$mpi,
    pred$pai
  )
  
  p_value <- cor_fit$p.value
  
  p_value <- case_when(
    p_value >= .05 ~ paste0(
      "= ",
      round(
        p_value,
        digits = 3
      )
    ),
    p_value < .05 & p_value >= .01 ~ "< .05",
    p_value < .01 & p_value >= .001 ~ "< .01",
    p_value < .001 ~ "< .001"
  )
  
  
  #######
  # 3.3 Create axis and correlation labels
  #######
  
  mean_mpi <- sprintf(
    "%.2f",
    mean(pred$mpi)
  )
  
  sd_mpi <- sprintf(
    "%.2f",
    sd(pred$mpi)
  )
  
  mean_pai <- sprintf(
    "%.2f",
    mean(pred$pai)
  )
  
  sd_pai <- sprintf(
    "%.2f",
    sd(pred$pai)
  )
  
  x_label <- bquote(
    Mean~Prognostic~Index~
      (italic(M) == .(mean_mpi) *
         "," ~ italic(SD) == .(sd_mpi))
  )
  
  y_label <- bquote(
    Personalized~Advantage~Index~
      (italic(M) == .(mean_pai) *
         "," ~ italic(SD) == .(sd_pai))
  )
  
  annotate_label <- paste0(
    "r = ",
    round(
      cor_fit$estimate,
      digits = 2
    ),
    ", p ",
    p_value
  )
  
  
  #######
  # 3.4 Create figure
  #######
  
  gg_plot <- pred %>%
    ggplot(
      aes(
        x = mpi,
        y = pai
      )
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 1.75
    ) +
    geom_point(
      size = 3,
      shape = 16
    ) +
    labs(
      x = x_label,
      y = y_label
    ) +
    annotate(
      "text",
      x = mean(pred$mpi),
      y = mean(pred$pai) + 2 * sd(pred$pai),
      label = annotate_label,
      size = rel(7.5)
    ) +
    theme(
      plot.title = element_text(
        size = rel(2.5),
        hjust = 0.5
      ),
      axis.title = element_text(
        size = rel(2)
      ),
      axis.text = element_text(
        size = rel(1.75)
      ),
      legend.title = element_text(
        size = rel(2)
      ),
      legend.text = element_text(
        size = rel(1.75)
      )
    )
  
  
  #######
  # 3.5 Save figure
  #######
  
  gg_name <- paste0(
    "figure3_",
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
    plot = gg_plot,
    width = 16,
    height = 13.4
  )
}


#######
# Run Figure 3
#######

plot_pai_mpi(0)

plot_pai_mpi(1)


#######
# 4. Figure 4
# Reliable Clinical Improvement by
# Decision-Rule Treatment Classification
#######

plot_rci_classification <- function(imp) {
  
  # imp <- 0
  
  
  #######
  # 4.1 Set imputation condition and extract data
  #######
  
  imp_name <- case_when(
    imp == 1 ~ "imp",
    imp == 0 ~ "noimp"
  )
  
  df_name <- paste0(
    "roc_",
    imp_name
  )
  
  pred <- get(df_name)
  
  
  #######
  # 4.2 Prepare data for plot
  #######
  
  plot_dat <- pred %>%
    mutate(
      group = case_when(
        cut_of == "No difference" ~
          "No Predicted\nClinical Advantage",
        
        cut_of == "Received non-optimal" ~
          "Received\nNon-Recommended Intensity",
        
        cut_of == "Received optimal" ~
          "Received\nRecommended Intensity"
      ),
      
      group = factor(
        group,
        levels = c(
          "No Predicted\nClinical Advantage",
          "Received\nNon-Recommended Intensity",
          "Received\nRecommended Intensity"
        )
      )
    ) %>%
    count(
      group,
      rci
    ) %>%
    group_by(
      group
    ) %>%
    mutate(
      n_total = sum(n),
      pct = 100 * n / n_total
    ) %>%
    ungroup() %>%
    filter(
      rci == 1
    ) %>%
    mutate(
      label = sprintf(
        "%d/%d (%.1f%%)",
        n,
        n_total,
        pct
      )
    )
  
  
  #######
  # 4.3 Unadjusted odds ratio
  #######
  
  or_dat <- pred %>%
    filter(
      cut_of %in% c(
        "Received non-optimal",
        "Received optimal"
      )
    ) %>%
    mutate(
      rci = factor(rci),
      
      cut_of = factor(
        cut_of,
        levels = c(
          "Received non-optimal",
          "Received optimal"
        )
      )
    )
  
  fit_or <- glm(
    rci ~ cut_of,
    data = or_dat,
    family = binomial
  )
  
  or_res <- broom::tidy(
    fit_or,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      str_detect(
        term,
        "^cut_of"
      )
    )
  
  or_label <- sprintf(
    "OR = %.2f, 95%% CI [%.2f, %.2f]",
    or_res$estimate,
    or_res$conf.low,
    or_res$conf.high
  )
  
  
  #######
  # 4.4 Create figure
  #######
  
  y_bracket <- max(plot_dat$pct) + 6
  
  y_top <- y_bracket + 3
  
  gg_plot <- ggplot(
    plot_dat,
    aes(
      x = group,
      y = pct,
      fill = group
    )
  ) +
    geom_col(
      width = .9
    ) +
    
    # n/N (%) above bars
    geom_text(
      aes(
        label = label
      ),
      vjust = -0.5,
      size = 4.5
    ) +
    
    # Bracket: vertical left
    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = y_bracket,
      yend = y_top
    ) +
    
    # Bracket: horizontal
    annotate(
      "segment",
      x = 2,
      xend = 3,
      y = y_top,
      yend = y_top
    ) +
    
    # Bracket: vertical right
    annotate(
      "segment",
      x = 3,
      xend = 3,
      y = y_bracket,
      yend = y_top
    ) +
    
    # Odds-ratio label
    annotate(
      "text",
      x = 2.5,
      y = y_top + 1,
      label = or_label,
      fontface = "italic",
      size = 4.2,
      vjust = 0
    ) +
    
    scale_fill_manual(
      values = c(
        "grey85",
        "grey65",
        "grey30"
      )
    ) +
    
    scale_y_continuous(
      limits = c(
        0,
        y_top + 7
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(
      base_size = 15
    ) +
    
    theme(
      legend.position = "none",
      
      axis.text.x = element_text(
        size = 13,
        lineheight = .9
      ),
      
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      panel.grid.minor = element_blank()
    )
  
  
  #######
  # 4.5 Save figure
  #######
  
  gg_name <- paste0(
    "figure4_",
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
    plot = gg_plot,
    width = 16,
    height = 13.4
  )
}


#######
# Run Figure 4
#######

plot_rci_classification(0)

plot_rci_classification(1)