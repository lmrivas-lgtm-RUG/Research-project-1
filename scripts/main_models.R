# ============================================================
# 04_main_models.R
# Main statistical models M1-M6
# RP1 - Luciano Rivas
# ============================================================


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(tibble)
library(broom)
library(car)
library(emmeans)


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

input_file <- "data/processed/final_soil_variables_table_updated.csv"

aes_file <- "data/raw/AES_final_table_2023.csv"

results_dir <-
  "results/main_models"


dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load soil + SHI dataset
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
# 5. Load AES income dataset
# ------------------------------------------------------------

aes_lookup <- read_csv(
  aes_file,
  show_col_types = FALSE
) |>
  transmute(
    
    field_id =
      as.numeric(Field_ID),
    
    AES_subsidies_per_ha =
      as.numeric(
        AES_income_EUR_per_ha
      )
  ) |>
  distinct(
    field_id,
    .keep_all = TRUE
  )


if (
  anyDuplicated(
    aes_lookup$field_id
  ) > 0
) {
  stop(
    "Duplicated field IDs found in AES dataset."
  )
}


# ------------------------------------------------------------
# 6. Prepare common analysis dataset
# ------------------------------------------------------------

management_levels <- c(
  "Diversified",
  "Intermediate",
  "Intensive"
)


analysis_data <- dat |>
  mutate(
    
    field_id =
      as.numeric(Plot),
    
    clay_silt =
      as.numeric(
        ClaySilt_Average_percent
      ),
    
    management_type =
      factor(
        `Management type`,
        levels = management_levels
      )
  ) |>
  
  left_join(
    aes_lookup,
    by = "field_id"
  ) |>
  
  mutate(
    
    # Fields without AES payments are assigned zero income.
    AES_subsidies_per_ha =
      coalesce(
        AES_subsidies_per_ha,
        0
      ),
    
    log1p_TerritorialsPerHa =
      log1p(
        TerritorialsPerHa
      )
  )


stopifnot(
  nrow(analysis_data) == 34,
  n_distinct(analysis_data$field_id) == 34
)


if (
  any(
    is.na(
      analysis_data$management_type
    )
  )
) {
  stop(
    "Unexpected management category detected."
  )
}


# ------------------------------------------------------------
# 7. Helper functions
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


model_n <- function(fit) {
  
  nrow(
    stats::model.frame(fit)
  )
}


# ============================================================
# M1
# Response: SHI
# Predictors: clay_silt + management_type
# ============================================================

m1_data <- analysis_data |>
  filter(
    complete.cases(
      SHI,
      clay_silt,
      management_type
    )
  )


m1 <- lm(
  SHI ~
    clay_silt +
    management_type,
  
  data = m1_data
)


# ============================================================
# COMMON GODWIT DATASET
# ============================================================

# M2, M3, M5 and M6 use exactly the same observations.
# This allows direct comparison of their AIC values.

godwit_model_data <- analysis_data |>
  filter(
    complete.cases(
      TerritorialsPerHa,
      SHI,
      AES_subsidies_per_ha,
      clay_silt,
      management_type
    )
  )


# ============================================================
# M2
# Response: log1p(TerritorialsPerHa)
# Predictors: clay_silt + management_type
# ============================================================

m2 <- lm(
  log1p(TerritorialsPerHa) ~
    clay_silt +
    management_type,
  
  data = godwit_model_data
)


# ============================================================
# M3
# Response: log1p(TerritorialsPerHa)
# Predictors: SHI + clay_silt + management_type
# ============================================================

m3 <- lm(
  log1p(TerritorialsPerHa) ~
    SHI +
    clay_silt +
    management_type,
  
  data = godwit_model_data
)


# ============================================================
# M4
# Response: SHI
# Predictors: AES subsidies + clay_silt + management_type
# ============================================================

m4_data <- analysis_data |>
  filter(
    complete.cases(
      SHI,
      AES_subsidies_per_ha,
      clay_silt,
      management_type
    )
  )


m4 <- lm(
  SHI ~
    AES_subsidies_per_ha +
    clay_silt +
    management_type,
  
  data = m4_data
)


# ============================================================
# M5
# Response: log1p(TerritorialsPerHa)
# Predictors: AES subsidies + clay_silt + management_type
# ============================================================

m5 <- lm(
  log1p(TerritorialsPerHa) ~
    AES_subsidies_per_ha +
    clay_silt +
    management_type,
  
  data = godwit_model_data
)


# ============================================================
# M6
# Response: log1p(TerritorialsPerHa)
# Predictors: SHI + AES subsidies + clay_silt + management_type
# ============================================================

m6 <- lm(
  log1p(TerritorialsPerHa) ~
    SHI +
    AES_subsidies_per_ha +
    clay_silt +
    management_type,
  
  data = godwit_model_data
)


# ------------------------------------------------------------
# 8. Verify that M1-M6 are lm objects
# ------------------------------------------------------------

stopifnot(
  inherits(m1, "lm"),
  inherits(m2, "lm"),
  inherits(m3, "lm"),
  inherits(m4, "lm"),
  inherits(m5, "lm"),
  inherits(m6, "lm")
)


# ------------------------------------------------------------
# 9. Verify sample sizes
# ------------------------------------------------------------

godwit_model_n <- c(
  M2 = model_n(m2),
  M3 = model_n(m3),
  M5 = model_n(m5),
  M6 = model_n(m6)
)


# All candidate Godwit models must use the same fields
stopifnot(
  length(
    unique(
      godwit_model_n
    )
  ) == 1
)


# ------------------------------------------------------------
# 10. Model formulas
# ------------------------------------------------------------

model_formulas <- c(
  
  M1 =
    "SHI ~ clay_silt + management_type",
  
  M2 =
    "log1p(TerritorialsPerHa) ~ clay_silt + management_type",
  
  M3 =
    "log1p(TerritorialsPerHa) ~ SHI + clay_silt + management_type",
  
  M4 =
    "SHI ~ AES_subsidies_per_ha + clay_silt + management_type",
  
  M5 =
    "log1p(TerritorialsPerHa) ~ AES_subsidies_per_ha + clay_silt + management_type",
  
  M6 =
    "log1p(TerritorialsPerHa) ~ SHI + AES_subsidies_per_ha + clay_silt + management_type"
)


# ============================================================
# MODEL COEFFICIENTS
# ============================================================


# ------------------------------------------------------------
# 11. Function to extract coefficients
# ------------------------------------------------------------

extract_coefficients <- function(
    fit,
    model_name
) {
  
  broom::tidy(
    fit,
    conf.int = TRUE
  ) |>
    mutate(
      
      model =
        model_name,
      
      significance =
        significance_label(
          p.value
        ),
      
      .before = 1
    )
}


main_model_coefficients <- bind_rows(
  
  extract_coefficients(
    m1,
    "M1"
  ),
  
  extract_coefficients(
    m2,
    "M2"
  ),
  
  extract_coefficients(
    m3,
    "M3"
  ),
  
  extract_coefficients(
    m4,
    "M4"
  ),
  
  extract_coefficients(
    m5,
    "M5"
  ),
  
  extract_coefficients(
    m6,
    "M6"
  )
)


# ============================================================
# TYPE-II ANOVA
# ============================================================


# ------------------------------------------------------------
# 12. Function to extract Type-II tests
# ------------------------------------------------------------

extract_type2 <- function(
    fit,
    model_name
) {
  
  car::Anova(
    fit,
    type = 2
  ) |>
    broom::tidy() |>
    mutate(
      
      model =
        model_name,
      
      significance =
        significance_label(
          p.value
        ),
      
      .before = 1
    )
}


main_model_type2_anova <- bind_rows(
  
  extract_type2(
    m1,
    "M1"
  ),
  
  extract_type2(
    m2,
    "M2"
  ),
  
  extract_type2(
    m3,
    "M3"
  ),
  
  extract_type2(
    m4,
    "M4"
  ),
  
  extract_type2(
    m5,
    "M5"
  ),
  
  extract_type2(
    m6,
    "M6"
  )
)


# ============================================================
# MODEL METRICS
# ============================================================


# ------------------------------------------------------------
# 13. Function to extract model metrics
# ------------------------------------------------------------

extract_metrics <- function(
    fit,
    model_name
) {
  
  s <- summary(fit)
  
  tibble(
    
    model =
      model_name,
    
    formula =
      unname(
        model_formulas[
          model_name
        ]
      ),
    
    n =
      model_n(
        fit
      ),
    
    R2 =
      s$r.squared,
    
    adjusted_R2 =
      s$adj.r.squared,
    
    AIC =
      stats::AIC(
        fit
      )
  )
}


main_model_metrics <- bind_rows(
  
  extract_metrics(
    m1,
    "M1"
  ),
  
  extract_metrics(
    m2,
    "M2"
  ),
  
  extract_metrics(
    m3,
    "M3"
  ),
  
  extract_metrics(
    m4,
    "M4"
  ),
  
  extract_metrics(
    m5,
    "M5"
  ),
  
  extract_metrics(
    m6,
    "M6"
  )
)


# ============================================================
# GODWIT MODEL COMPARISON
# ============================================================


# ------------------------------------------------------------
# 14. Compare M2, M3, M5 and M6
# ------------------------------------------------------------

godwit_model_comparison <-
  main_model_metrics |>
  
  filter(
    model %in%
      c(
        "M2",
        "M3",
        "M5",
        "M6"
      )
  ) |>
  
  arrange(
    AIC
  ) |>
  
  mutate(
    delta_AIC =
      AIC -
      min(AIC)
  ) |>
  
  select(
    model,
    formula,
    AIC,
    delta_AIC,
    R2,
    adjusted_R2,
    n
  )


# ------------------------------------------------------------
# 15. Nested comparison: M2 vs M3
# ------------------------------------------------------------

m2_vs_m3_nested <-
  anova(
    m2,
    m3
  ) |>
  broom::tidy()


# ============================================================
# ESTIMATED MARGINAL MEANS
# ============================================================


# ------------------------------------------------------------
# 16. M1 management contrasts
# ------------------------------------------------------------

m1_emmeans <-
  emmeans(
    m1,
    ~ management_type
  )


m1_emmeans_table <-
  broom::tidy(
    m1_emmeans
  )


m1_pairwise <-
  pairs(
    m1_emmeans,
    adjust = "tukey"
  ) |>
  broom::tidy()


# ------------------------------------------------------------
# 17. M3 management contrasts
# ------------------------------------------------------------

m3_emmeans <-
  emmeans(
    m3,
    ~ management_type
  )


m3_emmeans_table <-
  broom::tidy(
    m3_emmeans
  )


m3_pairwise <-
  pairs(
    m3_emmeans,
    adjust = "tukey"
  ) |>
  broom::tidy()


# ------------------------------------------------------------
# 18. M4 management contrasts
# ------------------------------------------------------------

m4_emmeans <-
  emmeans(
    m4,
    ~ management_type
  )


m4_emmeans_table <-
  broom::tidy(
    m4_emmeans
  )


m4_pairwise <-
  pairs(
    m4_emmeans,
    adjust = "tukey"
  ) |>
  broom::tidy()


# ------------------------------------------------------------
# 19. M5 management contrasts
# ------------------------------------------------------------

m5_emmeans <-
  emmeans(
    m5,
    ~ management_type
  )


m5_emmeans_table <-
  broom::tidy(
    m5_emmeans
  )


m5_pairwise <-
  pairs(
    m5_emmeans,
    adjust = "tukey"
  ) |>
  broom::tidy()


# ------------------------------------------------------------
# 20. M6 management contrasts
# ------------------------------------------------------------

m6_emmeans <-
  emmeans(
    m6,
    ~ management_type
  )


m6_emmeans_table <-
  broom::tidy(
    m6_emmeans
  )


m6_pairwise <-
  pairs(
    m6_emmeans,
    adjust = "tukey"
  ) |>
  broom::tidy()


# ============================================================
# STANDARDIZED GODWIT MODELS
# ============================================================


# ------------------------------------------------------------
# 21. Standardize continuous variables
# ------------------------------------------------------------

godwit_standardized_data <-
  godwit_model_data |>
  mutate(
    
    z_godwit =
      as.numeric(
        scale(
          log1p(
            TerritorialsPerHa
          )
        )
      ),
    
    z_SHI =
      as.numeric(
        scale(
          SHI
        )
      ),
    
    z_AES_subsidies_per_ha =
      as.numeric(
        scale(
          AES_subsidies_per_ha
        )
      ),
    
    z_clay_silt =
      as.numeric(
        scale(
          clay_silt
        )
      )
  )


# ------------------------------------------------------------
# 22. Standardized M3
# ------------------------------------------------------------

m3_std <- lm(
  z_godwit ~
    z_SHI +
    z_clay_silt +
    management_type,
  
  data =
    godwit_standardized_data
)


# ------------------------------------------------------------
# 23. Standardized M5
# ------------------------------------------------------------

m5_std <- lm(
  z_godwit ~
    z_AES_subsidies_per_ha +
    z_clay_silt +
    management_type,
  
  data =
    godwit_standardized_data
)


# ------------------------------------------------------------
# 24. Standardized M6
# ------------------------------------------------------------

m6_std <- lm(
  z_godwit ~
    z_SHI +
    z_AES_subsidies_per_ha +
    z_clay_silt +
    management_type,
  
  data =
    godwit_standardized_data
)


# ------------------------------------------------------------
# 25. Standardized coefficient table
# ------------------------------------------------------------

godwit_standardized_coefficients <- bind_rows(
  
  broom::tidy(
    m3_std,
    conf.int = TRUE
  ) |>
    mutate(
      model = "M3",
      .before = 1
    ),
  
  broom::tidy(
    m5_std,
    conf.int = TRUE
  ) |>
    mutate(
      model = "M5",
      .before = 1
    ),
  
  broom::tidy(
    m6_std,
    conf.int = TRUE
  ) |>
    mutate(
      model = "M6",
      .before = 1
    )
  
) |>
  filter(
    term != "(Intercept)"
  )


# ============================================================
# M6 DIAGNOSTICS
# ============================================================


# ------------------------------------------------------------
# 26. Residual diagnostics
# ------------------------------------------------------------

m6_residual_diagnostics <-
  godwit_model_data |>
  transmute(
    
    field_id,
    
    management_type,
    
    TerritorialsPerHa,
    
    SHI,
    
    AES_subsidies_per_ha,
    
    clay_silt,
    
    fitted =
      fitted(
        m6
      ),
    
    residual =
      residuals(
        m6
      ),
    
    standardized_residual =
      rstandard(
        m6
      ),
    
    studentized_residual =
      rstudent(
        m6
      ),
    
    leverage =
      hatvalues(
        m6
      ),
    
    cooks_distance =
      cooks.distance(
        m6
      )
  )


# ------------------------------------------------------------
# 27. Influential observations
# ------------------------------------------------------------

n_m6 <-
  model_n(
    m6
  )


p_m6 <-
  length(
    coef(
      m6
    )
  )


m6_influential_observations <-
  m6_residual_diagnostics |>
  mutate(
    
    cooks_cutoff =
      4 / n_m6,
    
    leverage_cutoff =
      2 * p_m6 / n_m6,
    
    studentized_residual_cutoff =
      2,
    
    flag_cooks =
      cooks_distance >
      cooks_cutoff,
    
    flag_leverage =
      leverage >
      leverage_cutoff,
    
    flag_studentized_residual =
      abs(
        studentized_residual
      ) >
      studentized_residual_cutoff,
    
    influential =
      flag_cooks |
      flag_leverage |
      flag_studentized_residual
  ) |>
  filter(
    influential
  )


# ============================================================
# SAVE OUTPUTS
# ============================================================


# ------------------------------------------------------------
# 28. Model coefficients
# ------------------------------------------------------------

write_csv(
  main_model_coefficients,
  file.path(
    results_dir,
    "main_model_coefficients.csv"
  )
)


# ------------------------------------------------------------
# 29. Type-II ANOVA
# ------------------------------------------------------------

write_csv(
  main_model_type2_anova,
  file.path(
    results_dir,
    "main_model_type2_anova.csv"
  )
)


# ------------------------------------------------------------
# 30. Model metrics
# ------------------------------------------------------------

write_csv(
  main_model_metrics,
  file.path(
    results_dir,
    "main_model_metrics.csv"
  )
)


# ------------------------------------------------------------
# 31. Godwit candidate-model comparison
# ------------------------------------------------------------

write_csv(
  godwit_model_comparison,
  file.path(
    results_dir,
    "godwit_model_comparison.csv"
  )
)


write_csv(
  m2_vs_m3_nested,
  file.path(
    results_dir,
    "m2_vs_m3_nested_comparison.csv"
  )
)


# ------------------------------------------------------------
# 32. Estimated marginal means
# ------------------------------------------------------------

write_csv(
  m1_emmeans_table,
  file.path(
    results_dir,
    "m1_management_emmeans.csv"
  )
)


write_csv(
  m3_emmeans_table,
  file.path(
    results_dir,
    "m3_management_emmeans.csv"
  )
)


write_csv(
  m4_emmeans_table,
  file.path(
    results_dir,
    "m4_management_emmeans.csv"
  )
)


write_csv(
  m5_emmeans_table,
  file.path(
    results_dir,
    "m5_management_emmeans.csv"
  )
)


write_csv(
  m6_emmeans_table,
  file.path(
    results_dir,
    "m6_management_emmeans.csv"
  )
)


# ------------------------------------------------------------
# 33. Pairwise management comparisons
# ------------------------------------------------------------

write_csv(
  m1_pairwise,
  file.path(
    results_dir,
    "m1_management_pairwise.csv"
  )
)


write_csv(
  m3_pairwise,
  file.path(
    results_dir,
    "m3_management_pairwise.csv"
  )
)


write_csv(
  m4_pairwise,
  file.path(
    results_dir,
    "m4_management_pairwise.csv"
  )
)


write_csv(
  m5_pairwise,
  file.path(
    results_dir,
    "m5_management_pairwise.csv"
  )
)


write_csv(
  m6_pairwise,
  file.path(
    results_dir,
    "m6_management_pairwise.csv"
  )
)


# ------------------------------------------------------------
# 34. Standardized coefficients
# ------------------------------------------------------------

write_csv(
  godwit_standardized_coefficients,
  file.path(
    results_dir,
    "godwit_standardized_coefficients.csv"
  )
)


# ------------------------------------------------------------
# 35. M6 diagnostics
# ------------------------------------------------------------

write_csv(
  m6_residual_diagnostics,
  file.path(
    results_dir,
    "m6_residual_diagnostics.csv"
  )
)


write_csv(
  m6_influential_observations,
  file.path(
    results_dir,
    "m6_influential_observations.csv"
  )
)


# ------------------------------------------------------------
# 36. Dataset used for Godwit candidate models
# ------------------------------------------------------------

write_csv(
  
  godwit_model_data |>
    select(
      field_id,
      management_type,
      TerritorialsPerHa,
      SHI,
      AES_subsidies_per_ha,
      clay_silt
    ),
  
  file.path(
    results_dir,
    "godwit_model_dataset.csv"
  )
)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat(
  "\nMain models M1-M6 completed successfully.\n"
)


cat(
  "\nSample sizes:\n"
)

cat(
  "M1:",
  model_n(m1),
  "\n"
)

cat(
  "M2:",
  model_n(m2),
  "\n"
)

cat(
  "M3:",
  model_n(m3),
  "\n"
)

cat(
  "M4:",
  model_n(m4),
  "\n"
)

cat(
  "M5:",
  model_n(m5),
  "\n"
)

cat(
  "M6:",
  model_n(m6),
  "\n"
)


cat(
  "\nModel metrics:\n"
)

print(
  main_model_metrics
)


cat(
  "\nGodwit candidate-model comparison:\n"
)

print(
  godwit_model_comparison
)


cat(
  "\nM6 coefficients:\n"
)

print(
  main_model_coefficients |>
    filter(
      model == "M6"
    )
)


cat(
  "\nResults saved to:\n",
  results_dir,
  "\n"
)
