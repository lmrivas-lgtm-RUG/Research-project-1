# ============================================================
# 03_soil_indicator_models.R
# Individual soil-indicator analyses
# RP1 - Luciano Rivas
# ============================================================

# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(broom)
library(car)


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

input_file <- "data/processed/final_soil_variables_table_updated.csv"

results_dir <-
  "results/soil_indicator_models"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load prepared dataset
# ------------------------------------------------------------

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)


stopifnot(
  nrow(dat) == 34,
  n_distinct(dat$Plot) == 34
)


# ------------------------------------------------------------
# 5. Standardise model variables
# ------------------------------------------------------------

management_levels <- c(
  "Diversified",
  "Intermediate",
  "Intensive"
)


dat <- dat |>
  mutate(
    management_type = factor(
      `Management type`,
      levels = management_levels
    ),
    
    clay_silt =
      as.numeric(
        ClaySilt_Average_percent
      ),
    
    log_godwit_density =
      log1p(
        TerritorialsPerHa
      )
  )


stopifnot(
  !any(
    is.na(dat$management_type)
  )
)


# ------------------------------------------------------------
# 6. Define soil indicators
# ------------------------------------------------------------

indicator_specs <- tribble(
  
  ~variable,                  ~indicator,
  
  "SOM_Average_g_kg_soil",   "SOM",
  "POM_Average_g_kg_soil",   "POM",
  "MAOM_mean_g_kg",          "MAOM",
  "Moisture_percent_DS",     "Soil moisture",
  "BD_0_10_mean",            "Bulk density",
  "pH_mean",                 "pH",
  "aryl_sul_mean",           "Arylsulfatase",
  "ac_fos_mean",             "Acid phosphatase",
  "beta_glu_mean",           "Beta-glucosidase",
  "urease_mean",             "Urease",
  "Shannon",                 "Shannon",
  "Shannon_SRS",             "Shannon-SRS"
)


# ------------------------------------------------------------
# 7. Significance function
# ------------------------------------------------------------

significance_label <- function(p) {
  
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ ".",
    TRUE ~ "ns"
  )
}


# ============================================================
# PART A
# Soil indicator ~ clay+silt + management
# ============================================================


# ------------------------------------------------------------
# 8. Fit models explaining individual soil indicators
# ------------------------------------------------------------

fit_soil_response_model <- function(variable) {
  
  model_data <- dat |>
    transmute(
      response =
        as.numeric(
          .data[[variable]]
        ),
      
      clay_silt,
      
      management_type
    ) |>
    drop_na()
  
  
  lm(
    response ~
      clay_silt +
      management_type,
    
    data = model_data
  )
}


soil_variance_models <- map(
  indicator_specs$variable,
  fit_soil_response_model
)


names(
  soil_variance_models
) <- indicator_specs$indicator


# ------------------------------------------------------------
# 9. Partition explained variation
# ------------------------------------------------------------

extract_variance_partition <- function(
    model,
    indicator_name
) {
  
  type2 <- car::Anova(
    model,
    type = 2
  ) |>
    broom::tidy() |>
    filter(
      term != "Residuals"
    )
  
  
  model_r2 <-
    summary(model)$r.squared
  
  
  predictor_ss <-
    sum(
      type2$sumsq
    )
  
  
  explained <- type2 |>
    mutate(
      
      indicator =
        indicator_name,
      
      component =
        case_when(
          term == "clay_silt" ~
            "Clay-Silt",
          
          term == "management_type" ~
            "Management type"
        ),
      
      # Rescale the Type-II sums of squares so that
      # predictor contributions sum to model R².
      r2_share =
        model_r2 *
        (
          sumsq /
            predictor_ss
        ),
      
      percent_explained =
        r2_share * 100,
      
      p_value =
        p.value,
      
      significance =
        significance_label(
          p_value
        )
    ) |>
    
    select(
      indicator,
      component,
      percent_explained,
      p_value,
      significance
    )
  
  
  unexplained <- tibble(
    indicator =
      indicator_name,
    
    component =
      "Unexplained",
    
    percent_explained =
      (1 - model_r2) * 100,
    
    p_value =
      NA_real_,
    
    significance =
      NA_character_
  )
  
  
  bind_rows(
    explained,
    unexplained
  )
}


soil_variance_partition <-
  imap_dfr(
    soil_variance_models,
    extract_variance_partition
  )


# ------------------------------------------------------------
# 10. Soil-indicator model summary
# ------------------------------------------------------------

extract_soil_model_summary <- function(
    model,
    indicator_name
) {
  
  type2 <- car::Anova(
    model,
    type = 2
  ) |>
    broom::tidy()
  
  
  f_stat <-
    summary(model)$fstatistic
  
  
  tibble(
    
    indicator =
      indicator_name,
    
    n =
      nobs(model),
    
    R2 =
      summary(model)$r.squared,
    
    adjusted_R2 =
      summary(model)$adj.r.squared,
    
    AIC =
      AIC(model),
    
    overall_F =
      unname(
        f_stat[1]
      ),
    
    overall_p =
      pf(
        f_stat[1],
        f_stat[2],
        f_stat[3],
        lower.tail = FALSE
      ),
    
    clay_silt_p =
      type2 |>
      filter(
        term == "clay_silt"
      ) |>
      pull(
        p.value
      ),
    
    management_p =
      type2 |>
      filter(
        term == "management_type"
      ) |>
      pull(
        p.value
      )
  )
}


soil_variance_model_summary <-
  imap_dfr(
    soil_variance_models,
    extract_soil_model_summary
  )


# ============================================================
# PART B
# Godwit density ~ individual soil indicator
# ============================================================


# ------------------------------------------------------------
# 11. Function to fit individual Godwit models
# ------------------------------------------------------------

fit_godwit_indicator_model <- function(
    variable
) {
  
  model_data <- dat |>
    transmute(
      
      Plot,
      
      godwit_density =
        TerritorialsPerHa,
      
      log_godwit_density,
      
      soil_indicator =
        as.numeric(
          .data[[variable]]
        ),
      
      clay_silt,
      
      management_type
    ) |>
    
    drop_na(
      log_godwit_density,
      soil_indicator,
      clay_silt,
      management_type
    )
  
  
  lm(
    log_godwit_density ~
      soil_indicator +
      clay_silt +
      management_type,
    
    data = model_data
  )
}


soil_godwit_models <- map(
  indicator_specs$variable,
  fit_godwit_indicator_model
)


names(
  soil_godwit_models
) <- indicator_specs$indicator


# ------------------------------------------------------------
# 12. Extract Godwit-model results
# ------------------------------------------------------------

extract_godwit_results <- function(
    model,
    indicator_name
) {
  
  coefficients <-
    broom::tidy(
      model,
      conf.int = TRUE
    )
  
  
  type2 <-
    car::Anova(
      model,
      type = 2
    ) |>
    broom::tidy()
  
  
  indicator_row <-
    coefficients |>
    filter(
      term == "soil_indicator"
    )
  
  
  clay_row <-
    type2 |>
    filter(
      term == "clay_silt"
    )
  
  
  management_row <-
    type2 |>
    filter(
      term == "management_type"
    )
  
  
  indicator_anova <-
    type2 |>
    filter(
      term == "soil_indicator"
    )
  
  
  # Total sum of squares of the response
  response_values <-
    model$model$log_godwit_density
  
  
  total_ss <-
    sum(
      (
        response_values -
          mean(response_values)
      )^2
    )
  
  
  # Semi-partial R²:
  # unique proportion of total response variance
  # attributed to the soil indicator.
  semi_partial_R2 <-
    indicator_anova$sumsq /
    total_ss
  
  
  tibble(
    
    indicator =
      indicator_name,
    
    estimate =
      indicator_row$estimate,
    
    std_error =
      indicator_row$std.error,
    
    conf_low =
      indicator_row$conf.low,
    
    conf_high =
      indicator_row$conf.high,
    
    p_value =
      indicator_row$p.value,
    
    significance =
      significance_label(
        indicator_row$p.value
      ),
    
    semi_partial_R2 =
      semi_partial_R2,
    
    clay_silt_p =
      clay_row$p.value,
    
    management_p =
      management_row$p.value,
    
    R2 =
      summary(model)$r.squared,
    
    adjusted_R2 =
      summary(model)$adj.r.squared,
    
    AIC =
      AIC(model),
    
    n =
      nobs(model)
  )
}


godwit_individual_indicator_results <-
  imap_dfr(
    soil_godwit_models,
    extract_godwit_results
  ) |>
  
  arrange(
    p_value
  )


# ------------------------------------------------------------
# 13. Identify significant soil indicators
# ------------------------------------------------------------

godwit_significant_indicators <-
  godwit_individual_indicator_results |>
  
  filter(
    p_value < 0.05
  )


# ------------------------------------------------------------
# 14. Save results
# ------------------------------------------------------------

write_csv(
  soil_variance_partition,
  file.path(
    results_dir,
    "soil_indicator_variance_partition.csv"
  )
)


write_csv(
  soil_variance_model_summary,
  file.path(
    results_dir,
    "soil_indicator_variance_model_summary.csv"
  )
)


write_csv(
  godwit_individual_indicator_results,
  file.path(
    results_dir,
    "godwit_individual_soil_indicator_models.csv"
  )
)


write_csv(
  godwit_significant_indicators,
  file.path(
    results_dir,
    "godwit_significant_soil_indicators.csv"
  )
)


# ------------------------------------------------------------
# 15. Summary
# ------------------------------------------------------------

cat(
  "\nIndividual soil-indicator analyses complete.\n"
)


cat(
  "\nSignificant soil indicators associated with territorial Godwit density:\n"
)


print(
  godwit_significant_indicators |>
    select(
      indicator,
      estimate,
      p_value,
      semi_partial_R2,
      n
    )
)


cat(
  "\nResults saved to:\n",
  results_dir,
  "\n"
)
