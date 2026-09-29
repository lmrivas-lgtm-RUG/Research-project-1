# ============================================================
# 02_build_SHI.R
# Construct the Soil Health Index (SHI)
# ============================================================


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(tibble)


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

input_file <- "data/processed/final_soil_variables_table.csv"

output_file <- "data/processed/final_soil_variables_table_updated.csv"

results_dir <- "results/SHI"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load prepared dataset
# ------------------------------------------------------------

final_soil_variables_table <- read_csv(
  input_file,
  show_col_types = FALSE
)

stopifnot(
  nrow(final_soil_variables_table) == 34
)


# ------------------------------------------------------------
# 5. Define candidate indicators
# ------------------------------------------------------------

full_indicator_names <- c(
  "SOM",
  "POM",
  "MAOM",
  "Soil moisture",
  "Bulk density",
  "pH",
  "Arylsulfatase",
  "Acid phosphatase",
  "Beta-glucosidase",
  "Urease",
  "Shannon",
  "Shannon_SRS"
)


# Indicators retained after correlation screening
selected_pca_indicators <- c(
  "SOM",
  "POM",
  "Soil moisture",
  "Arylsulfatase",
  "pH",
  "Beta-glucosidase",
  "Urease",
  "Shannon_SRS",
  "Shannon"
)


# Final minimum data set
final_mds_indicators <- c(
  "SOM",
  "Beta-glucosidase",
  "Shannon",
  "Urease"
)


# ------------------------------------------------------------
# 6. Scoring function
# ------------------------------------------------------------

more_is_better_score <- function(x) {
  
  rng <- range(
    x,
    na.rm = TRUE
  )
  
  if (
    !all(is.finite(rng)) ||
    isTRUE(all.equal(rng[1], rng[2]))
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }
  
  (x - rng[1]) /
    (rng[2] - rng[1])
}


# ------------------------------------------------------------
# 7. Prepare candidate indicator dataset
# ------------------------------------------------------------

corr_source <- final_soil_variables_table |>
  transmute(
    SOM = as.numeric(SOM_Average_g_kg_soil),
    POM = as.numeric(POM_Average_g_kg_soil),
    MAOM = as.numeric(MAOM_mean_g_kg),
    `Soil moisture` = as.numeric(Moisture_percent_DS),
    `Bulk density` = as.numeric(BD_0_10_mean),
    pH = as.numeric(pH_mean),
    Arylsulfatase = as.numeric(aryl_sul_mean),
    `Acid phosphatase` = as.numeric(ac_fos_mean),
    `Beta-glucosidase` = as.numeric(beta_glu_mean),
    Urease = as.numeric(urease_mean),
    Shannon = as.numeric(Shannon),
    Shannon_SRS = as.numeric(Shannon_SRS)
  )


stopifnot(
  all(
    vapply(
      corr_source,
      is.numeric,
      logical(1)
    )
  )
)


# ------------------------------------------------------------
# 8. Pearson correlation screening
# ------------------------------------------------------------

all_indicator_pearson_matrix <- cor(
  corr_source,
  use = "pairwise.complete.obs",
  method = "pearson"
)


high_correlation_pairs <- combn(
  full_indicator_names,
  2,
  simplify = FALSE
) |>
  lapply(
    function(pair) {
      
      tibble(
        var1 = pair[1],
        var2 = pair[2],
        r = all_indicator_pearson_matrix[
          pair[1],
          pair[2]
        ]
      )
      
    }
  ) |>
  bind_rows() |>
  mutate(
    abs_r = abs(r)
  ) |>
  filter(
    !is.na(r),
    abs_r >= 0.70
  ) |>
  arrange(
    desc(abs_r)
  )


# ------------------------------------------------------------
# 9. Prepare PCA dataset
# ------------------------------------------------------------

shi_pca_selected_data <- final_soil_variables_table |>
  transmute(
    Plot = Plot,
    
    management_type = `Management type`,
    
    SOM = as.numeric(
      SOM_Average_g_kg_soil
    ),
    
    POM = as.numeric(
      POM_Average_g_kg_soil
    ),
    
    `Soil moisture` = as.numeric(
      Moisture_percent_DS
    ),
    
    Arylsulfatase = as.numeric(
      aryl_sul_mean
    ),
    
    pH = as.numeric(
      pH_mean
    ),
    
    `Beta-glucosidase` = as.numeric(
      beta_glu_mean
    ),
    
    Urease = as.numeric(
      urease_mean
    ),
    
    Shannon_SRS = as.numeric(
      Shannon_SRS
    ),
    
    Shannon = as.numeric(
      Shannon
    )
  ) |>
  
  filter(
    if_all(
      all_of(
        selected_pca_indicators
      ),
      ~ !is.na(.)
    )
  )


shi_pca_input <- shi_pca_selected_data |>
  select(
    all_of(
      selected_pca_indicators
    )
  )


stopifnot(
  all(
    vapply(
      shi_pca_input,
      is.numeric,
      logical(1)
    )
  )
)


# ------------------------------------------------------------
# 10. Principal Component Analysis
# ------------------------------------------------------------

shi_pca_new <- prcomp(
  shi_pca_input,
  center = TRUE,
  scale. = TRUE
)


eigenvalues <- shi_pca_new$sdev^2


percent_variance <- (
  eigenvalues /
    sum(eigenvalues)
) * 100


cumulative_variance <- cumsum(
  percent_variance
)


# ------------------------------------------------------------
# 11. PCA summary tables
# ------------------------------------------------------------

shi_pca_eigen_table <- tibble(
  pc = paste0(
    "PC",
    seq_along(eigenvalues)
  ),
  
  eigenvalue = eigenvalues,
  
  percent_variance = percent_variance,
  
  cumulative_variance = cumulative_variance,
  
  retained_pc = (
    eigenvalue > 1 &
      percent_variance > 5
  )
)


shi_pca_loading_table <- shi_pca_new$rotation |>
  as.data.frame() |>
  rownames_to_column(
    var = "indicator"
  )


shi_pca_field_scores <- bind_cols(
  
  shi_pca_selected_data |>
    select(
      Plot,
      management_type
    ),
  
  as.data.frame(
    shi_pca_new$x
  )
)


# ------------------------------------------------------------
# 12. Define final Minimum Data Set
# ------------------------------------------------------------

shi_mds_indicator_map <- tibble(
  
  Indicator = c(
    "SOM",
    "Beta-glucosidase",
    "Shannon",
    "Urease"
  ),
  
  principal_component = c(
    "PC1",
    "PC1",
    "PC2",
    "PC2"
  ),
  
  scoring_function = "More is better"
)


# ------------------------------------------------------------
# 13. Calculate PC weights
# ------------------------------------------------------------

retained_pc_table <- shi_pca_eigen_table |>
  
  filter(
    retained_pc
  ) |>
  
  transmute(
    
    principal_component = pc,
    
    pc_variance_explained =
      percent_variance,
    
    pc_weight =
      percent_variance /
      sum(percent_variance)
  )


# ------------------------------------------------------------
# 14. Calculate indicator weights
# ------------------------------------------------------------

shi_mds_scoring_weighting_table <-
  shi_mds_indicator_map |>
  
  left_join(
    
    shi_pca_loading_table |>
      
      pivot_longer(
        cols = starts_with("PC"),
        names_to = "principal_component",
        values_to = "loading"
      ) |>
      
      mutate(
        abs_loading = abs(loading)
      ) |>
      
      select(
        indicator,
        principal_component,
        abs_loading
      ),
    
    by = c(
      "Indicator" = "indicator",
      "principal_component"
    )
  ) |>
  
  left_join(
    retained_pc_table,
    by = "principal_component"
  ) |>
  
  group_by(
    principal_component
  ) |>
  
  mutate(
    indicator_weight_raw =
      pc_weight *
      abs_loading /
      sum(abs_loading)
  ) |>
  
  ungroup() |>
  
  mutate(
    `Indicator weight` =
      indicator_weight_raw /
      sum(indicator_weight_raw)
  ) |>
  
  transmute(
    Indicator,
    
    `Principal component` =
      principal_component,
    
    `Absolute loading` =
      abs_loading,
    
    `PC variance explained` =
      pc_variance_explained,
    
    `Indicator weight`,
    
    `Scoring function` =
      scoring_function
  )


# Check that weights sum to 1
stopifnot(
  abs(
    sum(
      shi_mds_scoring_weighting_table$
        `Indicator weight`
    ) - 1
  ) < 1e-10
)


# ------------------------------------------------------------
# 15. Final SHI weights
# ------------------------------------------------------------

# These weights reproduce the final SHI used in the analysis.

shi_weights <- c(
  SOM = 0.323,
  Beta_glucosidase = 0.280,
  Shannon = 0.205,
  Urease = 0.192
)


# Check that final weights sum to 1
stopifnot(
  abs(
    sum(shi_weights) - 1
  ) < 1e-10
)


# ------------------------------------------------------------
# 16. Score indicators and calculate SHI
# ------------------------------------------------------------

final_soil_variables_table_updated <-
  final_soil_variables_table |>
  
  mutate(
    
    SOM_score =
      more_is_better_score(
        as.numeric(
          SOM_Average_g_kg_soil
        )
      ),
    
    Beta_glucosidase_score =
      more_is_better_score(
        as.numeric(
          beta_glu_mean
        )
      ),
    
    Shannon_score =
      more_is_better_score(
        as.numeric(
          Shannon
        )
      ),
    
    Urease_score =
      more_is_better_score(
        as.numeric(
          urease_mean
        )
      ),
    
    SHI = if_else(
      
      is.na(SOM_Average_g_kg_soil) |
        is.na(beta_glu_mean) |
        is.na(Shannon) |
        is.na(urease_mean),
      
      NA_real_,
      
      shi_weights[["SOM"]] *
        SOM_score +
        
        shi_weights[["Beta_glucosidase"]] *
        Beta_glucosidase_score +
        
        shi_weights[["Shannon"]] *
        Shannon_score +
        
        shi_weights[["Urease"]] *
        Urease_score
    )
  )


# ------------------------------------------------------------
# 17. Create field-level SHI table
# ------------------------------------------------------------

shi_by_field <-
  final_soil_variables_table_updated |>
  
  transmute(
    
    field_identifier = Plot,
    
    management_type =
      `Management type`,
    
    SOM_score,
    
    Beta_glucosidase_score,
    
    Shannon_score,
    
    Urease_score,
    
    SHI
  )


# ------------------------------------------------------------
# 18. Integrity checks
# ------------------------------------------------------------

# Dataset should still contain 34 fields
stopifnot(
  nrow(
    final_soil_variables_table_updated
  ) == 34
)


# Field IDs should remain unique
stopifnot(
  n_distinct(
    final_soil_variables_table_updated$Plot
  ) == 34
)


# SHI values must remain between 0 and 1
stopifnot(
  
  all(
    final_soil_variables_table_updated$
      SHI[
        !is.na(
          final_soil_variables_table_updated$SHI
        )
      ] >= 0
  ),
  
  all(
    final_soil_variables_table_updated$
      SHI[
        !is.na(
          final_soil_variables_table_updated$SHI
        )
      ] <= 1
  )
)


# ------------------------------------------------------------
# 19. Save processed dataset
# ------------------------------------------------------------

write_csv(
  final_soil_variables_table_updated,
  output_file
)


# ------------------------------------------------------------
# 20. Save key SHI results
# ------------------------------------------------------------

write_csv(
  shi_by_field,
  file.path(
    results_dir,
    "shi_by_field.csv"
  )
)


write_csv(
  shi_mds_scoring_weighting_table,
  file.path(
    results_dir,
    "shi_mds_scoring_weighting_table.csv"
  )
)


write_csv(
  shi_pca_eigen_table,
  file.path(
    results_dir,
    "shi_pca_eigen_table.csv"
  )
)


write_csv(
  shi_pca_loading_table,
  file.path(
    results_dir,
    "shi_pca_loading_table.csv"
  )
)


write_csv(
  shi_pca_field_scores,
  file.path(
    results_dir,
    "shi_pca_field_scores.csv"
  )
)


write_csv(
  high_correlation_pairs,
  file.path(
    results_dir,
    "shi_high_correlation_pairs.csv"
  )
)


# ------------------------------------------------------------
# 21. Summary
# ------------------------------------------------------------

cat(
  "\nSHI construction complete.\n"
)

cat(
  "\nPCA variance explained:\n"
)

print(
  shi_pca_eigen_table
)


cat(
  "\nFinal MDS and weights:\n"
)

print(
  shi_mds_scoring_weighting_table
)


cat(
  "\nNumber of missing SHI values:",
  sum(
    is.na(
      final_soil_variables_table_updated$SHI
    )
  ),
  "\n"
)


cat(
  "\nUpdated dataset saved to:\n",
  output_file,
  "\n"
)
