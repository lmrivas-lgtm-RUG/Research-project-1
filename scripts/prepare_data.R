# ============================================================
# 01_prepare_data.R
# Prepare field-level dataset for subsequent analyses
# RP1 - Luciano Rivas
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(readr)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------

soil_file <- "data/raw/final_soil_variables_table.csv"

fungi_file <- "data/raw/Complete data set - Fungi diversity.csv"

output_file <- "data/processed/final_soil_variables_table.csv"

# ------------------------------------------------------------
# 3. Check input files
# ------------------------------------------------------------

if (!file.exists(soil_file)) {
  stop("Soil dataset not found: ", soil_file)
}

if (!file.exists(fungi_file)) {
  stop("Fungal diversity dataset not found: ", fungi_file)
}


# ------------------------------------------------------------
# 4. Create output directory
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Load field-level soil dataset
# ------------------------------------------------------------

final_soil_variables_table <- read_csv(
  soil_file,
  show_col_types = FALSE
)


# Check expected structure
required_soil_columns <- c(
  "Plot",
  "TerritorialsPerHa"
)

missing_soil_columns <- setdiff(
  required_soil_columns,
  names(final_soil_variables_table)
)

if (length(missing_soil_columns) > 0) {
  stop(
    "Missing columns in soil dataset: ",
    paste(missing_soil_columns, collapse = ", ")
  )
}


# Study contains 34 fields
stopifnot(nrow(final_soil_variables_table) == 34)

# Preserve original field order
original_plot_order <- final_soil_variables_table$Plot


# ------------------------------------------------------------
# 6. Load fungal diversity dataset
# ------------------------------------------------------------

fungi_diversity <- read_csv(
  fungi_file,
  show_col_types = FALSE
)


required_fungi_columns <- c(
  "Plot",
  "Shannon",
  "Shannon_SRS"
)

missing_fungi_columns <- setdiff(
  required_fungi_columns,
  names(fungi_diversity)
)

if (length(missing_fungi_columns) > 0) {
  stop(
    "Missing columns in fungal diversity dataset: ",
    paste(missing_fungi_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 7. Calculate field-level fungal diversity
# ------------------------------------------------------------

fungi_field_means <- fungi_diversity |>
  mutate(
    # Convert replicate IDs such as 1.1, 1.2, 1.3
    # into their corresponding field number
    field = as.numeric(
      sub("\\..*$", "", as.character(Plot))
    ),
    
    # "ns" indicates unavailable sequencing data
    Shannon = na_if(Shannon, "ns"),
    Shannon_SRS = na_if(Shannon_SRS, "ns"),
    
    # Convert diversity metrics to numeric
    Shannon = as.numeric(Shannon),
    Shannon_SRS = as.numeric(Shannon_SRS)
  ) |>
  group_by(field) |>
  summarise(
    Shannon = mean(Shannon, na.rm = TRUE),
    Shannon_SRS = mean(Shannon_SRS, na.rm = TRUE),
    .groups = "drop"
  )


# Check that fungal data produced one row per field
if (anyDuplicated(fungi_field_means$field) > 0) {
  stop("Duplicated field IDs found after fungal aggregation.")
}


# ------------------------------------------------------------
# 8. Merge fungal diversity with soil dataset
# ------------------------------------------------------------

final_soil_variables_table <- final_soil_variables_table |>
  mutate(.row_id = row_number()) |>
  
  # Remove previous Shannon columns if they already exist
  select(
    -any_of(c("Shannon", "Shannon_SRS"))
  ) |>
  
  # Add field-level fungal diversity
  left_join(
    fungi_field_means,
    by = c("Plot" = "field")
  ) |>
  
  # Godwit-density data were unavailable for fields 14-19
  mutate(
    TerritorialsPerHa = if_else(
      Plot %in% 14:19,
      NA_real_,
      as.numeric(TerritorialsPerHa)
    )
  ) |>
  
  # Restore original field order
  arrange(.row_id) |>
  select(-.row_id)


# ------------------------------------------------------------
# 9. Integrity checks
# ------------------------------------------------------------

# Dataset must still contain 34 fields
stopifnot(
  nrow(final_soil_variables_table) == 34
)

# Each field must occur only once
stopifnot(
  n_distinct(final_soil_variables_table$Plot) == 34
)

# Original field order must be preserved
stopifnot(
  identical(
    final_soil_variables_table$Plot,
    original_plot_order
  )
)

# Fields 14-19 must have missing godwit-density data
stopifnot(
  all(
    is.na(
      final_soil_variables_table$TerritorialsPerHa[
        final_soil_variables_table$Plot %in% 14:19
      ]
    )
  )
)


# ------------------------------------------------------------
# 10. Summary check
# ------------------------------------------------------------

missing_summary <- final_soil_variables_table |>
  summarise(
    missing_Shannon = sum(is.na(Shannon)),
    missing_Shannon_SRS = sum(is.na(Shannon_SRS)),
    missing_TerritorialsPerHa = sum(is.na(TerritorialsPerHa))
  )

print(missing_summary)


# ------------------------------------------------------------
# 11. Save processed dataset
# ------------------------------------------------------------

write_csv(
  final_soil_variables_table,
  output_file
)

message(
  "Data preparation complete."
)

message(
  "Processed dataset saved to: ",
  output_file
)
