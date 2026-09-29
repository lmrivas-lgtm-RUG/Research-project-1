# ============================================================
# 05_figures_tables.R
# Publication figures and summary tables
# RP1 - Luciano Rivas
# ============================================================


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggcorrplot)
library(factoextra)
library(scales)
library(patchwork)


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

data_file <- "data/processed/final_soil_variables_table_updated.csv"

soil_variance_file <-
  "results/soil_indicator_models/soil_indicator_variance_partition.csv"

model_coefficients_file <-
  "results/main_models/main_model_coefficients.csv"

model_anova_file <-
  "results/main_models/main_model_type2_anova.csv"

model_metrics_file <-
  "results/main_models/main_model_metrics.csv"

model_comparison_file <-
  "results/main_models/godwit_model_comparison.csv"

standardized_coefficients_file <-
  "results/main_models/godwit_standardized_coefficients.csv"

godwit_dataset_file <-
  "results/main_models/godwit_model_dataset.csv"

pca_eigen_file <-
  "results/SHI/shi_pca_eigen_table.csv"

mds_file <-
  "results/SHI/shi_mds_scoring_weighting_table.csv"

aes_file <- "data/raw/AES_final_table_2023.csv"


# ------------------------------------------------------------
# 4. Output folders
# ------------------------------------------------------------

main_figures_dir <-
  "figures/main"

supp_figures_dir <-
  "figures/supplementary"

tables_dir <-
  "results/publication_tables"


dir.create(
  main_figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supp_figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Load datasets
# ------------------------------------------------------------

dat <- read_csv(
  data_file,
  show_col_types = FALSE
)


soil_variance <- read_csv(
  soil_variance_file,
  show_col_types = FALSE
)


model_coefficients <- read_csv(
  model_coefficients_file,
  show_col_types = FALSE
)


model_anova <- read_csv(
  model_anova_file,
  show_col_types = FALSE
)


model_metrics <- read_csv(
  model_metrics_file,
  show_col_types = FALSE
)


model_comparison <- read_csv(
  model_comparison_file,
  show_col_types = FALSE
)


standardized_coefficients <- read_csv(
  standardized_coefficients_file,
  show_col_types = FALSE
)


godwit_data <- read_csv(
  godwit_dataset_file,
  show_col_types = FALSE
)


pca_eigen <- read_csv(
  pca_eigen_file,
  show_col_types = FALSE
)


mds_table <- read_csv(
  mds_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 6. Common settings
# ------------------------------------------------------------

management_levels <- c(
  "Diversified",
  "Intermediate",
  "Intensive"
)


management_colors <- c(
  Diversified = "#2E7D32",
  Intermediate = "#F9A825",
  Intensive = "#C62828"
)


dat <- dat |>
  mutate(
    management_type =
      factor(
        `Management type`,
        levels = management_levels
      ),
    
    clay_silt =
      as.numeric(
        ClaySilt_Average_percent
      )
  )


format_p_value <- function(p) {
  
  case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0(
      "p = ",
      formatC(
        p,
        format = "f",
        digits = 3
      )
    )
  )
}


save_figure <- function(
    plot_object,
    filename,
    width,
    height
) {
  
  ggsave(
    filename =
      file.path(
        main_figures_dir,
        paste0(
          filename,
          ".png"
        )
      ),
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    dpi =
      600,
    
    bg =
      "white"
  )
  
  
  ggsave(
    filename =
      file.path(
        main_figures_dir,
        paste0(
          filename,
          ".pdf"
        )
      ),
    
    plot =
      plot_object,
    
    width =
      width,
    
    height =
      height,
    
    bg =
      "white"
  )
}


# ============================================================
# FIGURE 1
# Pearson correlations among SHI candidate indicators
# ============================================================


corr_source <- dat |>
  transmute(
    
    SOM =
      as.numeric(
        SOM_Average_g_kg_soil
      ),
    
    POM =
      as.numeric(
        POM_Average_g_kg_soil
      ),
    
    MAOM =
      as.numeric(
        MAOM_mean_g_kg
      ),
    
    `Soil moisture` =
      as.numeric(
        Moisture_percent_DS
      ),
    
    `Bulk density` =
      as.numeric(
        BD_0_10_mean
      ),
    
    pH =
      as.numeric(
        pH_mean
      ),
    
    Arylsulfatase =
      as.numeric(
        aryl_sul_mean
      ),
    
    `Acid phosphatase` =
      as.numeric(
        ac_fos_mean
      ),
    
    `Beta-glucosidase` =
      as.numeric(
        beta_glu_mean
      ),
    
    Urease =
      as.numeric(
        urease_mean
      ),
    
    Shannon =
      as.numeric(
        Shannon
      ),
    
    Shannon_SRS =
      as.numeric(
        Shannon_SRS
      )
  )


pearson_matrix <- cor(
  corr_source,
  use = "pairwise.complete.obs",
  method = "pearson"
)


figure_pearson <- ggcorrplot(
  pearson_matrix,
  
  hc.order = FALSE,
  
  method = "square",
  
  type = "lower",
  
  lab = TRUE,
  
  lab_size = 3,
  
  show.legend = TRUE,
  
  outline.color = "white",
  
  colors = c(
    "#B2182B",
    "white",
    "#2166AC"
  )
) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 74,
        hjust = 1,
        vjust = 1
      ),
    
    axis.text.y =
      element_text(
        size = 11
      ),
    
    legend.title =
      element_blank()
  )


ggsave(
  filename =
    file.path(
      supp_figures_dir,
      "shi_candidate_pearson_correlations.png"
    ),
  
  plot =
    figure_pearson,
  
  width =
    7.5,
  
  height =
    6.5,
  
  dpi =
    600,
  
  bg =
    "white"
)


# ============================================================
# FIGURE 2
# PCA biplot
# ============================================================


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


pca_data <- dat |>
  transmute(
    
    management_type,
    
    SOM =
      as.numeric(
        SOM_Average_g_kg_soil
      ),
    
    POM =
      as.numeric(
        POM_Average_g_kg_soil
      ),
    
    `Soil moisture` =
      as.numeric(
        Moisture_percent_DS
      ),
    
    Arylsulfatase =
      as.numeric(
        aryl_sul_mean
      ),
    
    pH =
      as.numeric(
        pH_mean
      ),
    
    `Beta-glucosidase` =
      as.numeric(
        beta_glu_mean
      ),
    
    Urease =
      as.numeric(
        urease_mean
      ),
    
    Shannon_SRS =
      as.numeric(
        Shannon_SRS
      ),
    
    Shannon =
      as.numeric(
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


pca_input <- pca_data |>
  select(
    all_of(
      selected_pca_indicators
    )
  )


pca_fit <- prcomp(
  pca_input,
  center = TRUE,
  scale. = TRUE
)


pc1_percent <-
  pca_eigen$percent_variance[1]

pc2_percent <-
  pca_eigen$percent_variance[2]


figure_pca <- factoextra::fviz_pca_biplot(
  
  pca_fit,
  
  geom.ind =
    "point",
  
  habillage =
    pca_data$management_type,
  
  pointsize =
    2,
  
  alpha.ind =
    0.75,
  
  repel =
    TRUE,
  
  label =
    "var"
) +
  
  scale_color_manual(
    values =
      management_colors,
    
    name =
      "Management type"
  ) +
  
  labs(
    
    x =
      paste0(
        "PC1 (",
        round(
          pc1_percent,
          1
        ),
        "%)"
      ),
    
    y =
      paste0(
        "PC2 (",
        round(
          pc2_percent,
          1
        ),
        "%)"
      )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  theme_classic() +
  
  theme(
    legend.position =
      "right",
    
    axis.text =
      element_text(
        color = "black"
      )
  )


ggsave(
  filename =
    file.path(
      supp_figures_dir,
      "shi_pca_biplot.png"
    ),
  
  plot =
    figure_pca,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    600,
  
  bg =
    "white"
)


# ============================================================
# FIGURE 3
# PCA scree plot
# ============================================================


scree_data <- pca_eigen |>
  
  mutate(
    
    pc =
      factor(
        pc,
        levels = pc
      ),
    
    variance_label =
      paste0(
        round(
          percent_variance,
          1
        ),
        "%"
      ),
    
    retained =
      eigenvalue > 1
  )


figure_scree <- ggplot(
  scree_data,
  aes(
    x = pc,
    y = eigenvalue,
    group = 1
  )
) +
  
  geom_line(
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      fill = retained
    ),
    shape = 21,
    size = 3.2
  ) +
  
  geom_text(
    aes(
      label = variance_label
    ),
    vjust = -0.7,
    size = 3.2
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  scale_fill_manual(
    values = c(
      `TRUE` = "#2E7D32",
      `FALSE` = "grey75"
    ),
    
    name =
      "Eigenvalue > 1"
  ) +
  
  labs(
    x =
      "Principal component",
    
    y =
      "Eigenvalue"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text =
      element_text(
        color = "black"
      )
  )


ggsave(
  filename =
    file.path(
      supp_figures_dir,
      "shi_pca_scree_plot.png"
    ),
  
  plot =
    figure_scree,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    600,
  
  bg =
    "white"
)


# ============================================================
# FIGURE 4
# Field-level SHI
# ============================================================


shi_field_grid_data <- dat |>
  
  transmute(
    
    Plot =
      as.integer(
        Plot
      ),
    
    SHI,
    
    field_label =
      as.character(
        as.integer(
          Plot
        )
      ),
    
    shi_label =
      if_else(
        is.na(SHI),
        "NA",
        scales::number(
          SHI,
          accuracy = 0.01
        )
      )
  ) |>
  
  arrange(
    Plot
  ) |>
  
  mutate(
    
    panel_index =
      row_number(),
    
    col =
      (
        (
          panel_index - 1L
        ) %% 6L
      ) + 1L,
    
    row =
      (
        (
          panel_index - 1L
        ) %/% 6L
      ) + 1L,
    
    row =
      max(row) -
      row +
      1L
  )


figure_shi_fields <- ggplot(
  
  shi_field_grid_data,
  
  aes(
    x = col,
    y = row,
    fill = SHI
  )
  
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8,
    width = 0.96,
    height = 0.96
  ) +
  
  geom_text(
    aes(
      label = field_label
    ),
    vjust = -0.15,
    size = 4.4,
    fontface = "bold"
  ) +
  
  geom_text(
    aes(
      label = shi_label
    ),
    vjust = 1.35,
    size = 3.5
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    limits = c(0, 1),
    na.value = "grey85",
    name = "Soil Health Index"
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  
  coord_equal() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_void(
    base_size = 12
  ) +
  
  theme(
    legend.position =
      "bottom"
  )


save_figure(
  figure_shi_fields,
  "field_level_SHI",
  9,
  6.5
)


# ============================================================
# FIGURE 5
# SHI vs clay+silt and management
# ============================================================


m1_clay <- model_coefficients |>
  filter(
    model == "M1",
    term == "clay_silt"
  )


m1_management <- model_anova |>
  filter(
    model == "M1",
    term == "management_type"
  )


m1_metrics <- model_metrics |>
  filter(
    model == "M1"
  )


shi_model_data <- dat |>
  filter(
    complete.cases(
      SHI,
      clay_silt,
      management_type
    )
  )


shi_annotation <- paste(
  
  paste0(
    "Clay+silt: b = ",
    formatC(
      m1_clay$estimate,
      format = "f",
      digits = 4
    ),
    ", ",
    format_p_value(
      m1_clay$p.value
    )
  ),
  
  paste0(
    "Management: ",
    format_p_value(
      m1_management$p.value
    )
  ),
  
  paste0(
    "Adj. R² = ",
    formatC(
      m1_metrics$adjusted_R2,
      format = "f",
      digits = 3
    )
  ),
  
  sep = "\n"
)


figure_shi_texture <- ggplot(
  
  shi_model_data,
  
  aes(
    x = clay_silt,
    y = SHI,
    color = management_type
  )
  
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  geom_smooth(
    
    data =
      shi_model_data,
    
    aes(
      x = clay_silt,
      y = SHI,
      group = 1
    ),
    
    inherit.aes =
      FALSE,
    
    method =
      "lm",
    
    se =
      TRUE,
    
    color =
      "black",
    
    fill =
      "grey80"
  ) +
  
  scale_color_manual(
    values =
      management_colors,
    
    name =
      "Management type"
  ) +
  
  labs(
    x =
      "Clay+silt (%)",
    
    y =
      "Soil Health Index (SHI)"
  ) +
  
  annotate(
    
    "text",
    
    x =
      min(
        shi_model_data$clay_silt,
        na.rm = TRUE
      ),
    
    y =
      max(
        shi_model_data$SHI,
        na.rm = TRUE
      ),
    
    label =
      shi_annotation,
    
    hjust =
      0,
    
    vjust =
      1,
    
    size =
      3.3
  ) +
  
  theme_classic() +
  
  theme(
    legend.position =
      "right",
    
    axis.text =
      element_text(
        color = "black"
      )
  )


save_figure(
  figure_shi_texture,
  "SHI_vs_clay_silt_management",
  7.5,
  5.5
)


# ============================================================
# FIGURE 6
# Explained variation in individual soil indicators
# ============================================================


indicator_order <- c(
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
  "Shannon-SRS"
)


variance_plot_data <- soil_variance |>
  
  mutate(
    
    indicator =
      factor(
        indicator,
        levels =
          rev(
            indicator_order
          )
      ),
    
    component =
      factor(
        component,
        levels = c(
          "Clay-Silt",
          "Management type",
          "Unexplained"
        )
      )
  ) |>
  
  arrange(
    indicator,
    component
  ) |>
  
  group_by(
    indicator
  ) |>
  
  mutate(
    
    xmin =
      lag(
        cumsum(
          percent_explained
        ),
        default = 0
      ),
    
    xmax =
      cumsum(
        percent_explained
      ),
    
    x_mid =
      (
        xmin +
          xmax
      ) / 2
  ) |>
  
  ungroup()


variance_sig_data <- variance_plot_data |>
  
  filter(
    component != "Unexplained",
    significance %in%
      c(
        "***",
        "**",
        "*",
        "."
      )
  )


figure_variance <- ggplot() +
  
  geom_rect(
    
    data =
      variance_plot_data,
    
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin =
        as.numeric(
          indicator
        ) - 0.18,
      ymax =
        as.numeric(
          indicator
        ) + 0.18,
      fill = component
    ),
    
    color =
      "white",
    
    linewidth =
      0.7
  ) +
  
  geom_text(
    
    data =
      variance_sig_data,
    
    aes(
      x = x_mid,
      y =
        as.numeric(
          indicator
        ) + 0.27,
      label =
        significance
    ),
    
    vjust =
      0,
    
    size =
      4
  ) +
  
  scale_fill_manual(
    
    values = c(
      "Clay-Silt" = "#8D6E63",
      "Management type" = "#2E7D32",
      "Unexplained" = "#CCCCCC"
    ),
    
    name =
      NULL
  ) +
  
  scale_x_continuous(
    
    limits =
      c(
        0,
        100
      ),
    
    breaks =
      c(
        0,
        25,
        50,
        75,
        100
      ),
    
    expand =
      c(
        0,
        0
      )
  ) +
  
  scale_y_continuous(
    
    breaks =
      seq_along(
        levels(
          variance_plot_data$indicator
        )
      ),
    
    labels =
      levels(
        variance_plot_data$indicator
      ),
    
    expand =
      expansion(
        add = c(
          0.6,
          0.8
        )
      )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    x =
      "Percentage of total variation (%)",
    
    y =
      NULL
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    legend.position =
      "right",
    
    axis.line.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.text.y =
      element_text(
        size = 11
      ),
    
    plot.margin =
      margin(
        10,
        18,
        10,
        5.5
      )
  )


save_figure(
  figure_variance,
  "soil_indicator_explained_variation",
  9,
  6.5
)


# ============================================================
# FIGURE 7
# Godwit density vs SHI - M3
# ============================================================


godwit_data <- godwit_data |>
  
  mutate(
    management_type =
      factor(
        management_type,
        levels =
          management_levels
      ),
    
    log1p_TerritorialsPerHa =
      log1p(
        TerritorialsPerHa
      )
  )


m3_shi <- model_coefficients |>
  filter(
    model == "M3",
    term == "SHI"
  )


m3_metrics <- model_metrics |>
  filter(
    model == "M3"
  )


godwit_annotation <- paste(
  
  paste0(
    "SHI: b = ",
    formatC(
      m3_shi$estimate,
      format = "f",
      digits = 3
    ),
    ", ",
    format_p_value(
      m3_shi$p.value
    )
  ),
  
  paste0(
    "Adj. R² = ",
    formatC(
      m3_metrics$adjusted_R2,
      format = "f",
      digits = 3
    )
  ),
  
  sep = "\n"
)


figure_godwit_shi <- ggplot(
  
  godwit_data,
  
  aes(
    x = SHI,
    y = log1p_TerritorialsPerHa,
    color = management_type
  )
  
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  geom_smooth(
    
    data =
      godwit_data,
    
    aes(
      x = SHI,
      y = log1p_TerritorialsPerHa,
      group = 1
    ),
    
    inherit.aes =
      FALSE,
    
    method =
      "lm",
    
    se =
      TRUE,
    
    color =
      "black",
    
    fill =
      "grey80"
  ) +
  
  scale_color_manual(
    values =
      management_colors,
    
    name =
      "Management type"
  ) +
  
  labs(
    x =
      "Soil Health Index (SHI)",
    
    y =
      "log1p territorial Godwit density"
  ) +
  
  annotate(
    
    "text",
    
    x =
      min(
        godwit_data$SHI,
        na.rm = TRUE
      ),
    
    y =
      max(
        godwit_data$log1p_TerritorialsPerHa,
        na.rm = TRUE
      ),
    
    label =
      godwit_annotation,
    
    hjust =
      0,
    
    vjust =
      1,
    
    size =
      3.3
  ) +
  
  theme_classic() +
  
  theme(
    legend.position =
      "right",
    
    axis.text =
      element_text(
        color = "black"
      )
  )


save_figure(
  figure_godwit_shi,
  "godwit_density_vs_SHI",
  7.5,
  5.5
)


# ============================================================
# FIGURE 8
# Standardized coefficient forest plot
# M3, M5 and M6
# ============================================================

forest_data <- standardized_coefficients |>
  
  filter(
    model %in% c(
      "M3",
      "M5",
      "M6"
    )
  ) |>
  
  mutate(
    
    model = case_when(
      model == "M3" ~ "SHI model",
      model == "M5" ~ "AES model",
      model == "M6" ~ "Integrated model",
      TRUE ~ model
    ),
    
    predictor = case_when(
      term == "z_SHI" ~
        "Soil Health Index",
      
      term == "z_AES_subsidies_per_ha" ~
        "AES subsidies per ha",
      
      term == "z_clay_silt" ~
        "Clay+silt",
      
      term == "management_typeIntermediate" ~
        "Intermediate vs Diversified",
      
      term == "management_typeIntensive" ~
        "Intensive vs Diversified",
      
      TRUE ~ term
    ),
    
    predictor = factor(
      predictor,
      levels = rev(
        c(
          "Soil Health Index",
          "AES subsidies per ha",
          "Clay+silt",
          "Intermediate vs Diversified",
          "Intensive vs Diversified"
        )
      )
    ),
    
    model = factor(
      model,
      levels = c(
        "SHI model",
        "AES model",
        "Integrated model"
      )
    )
  )


model_colours <- c(
  "SHI model" = "#F8766D",
  "AES model" = "#00BA38",
  "Integrated model" = "#619CFF"
)


figure_forest <- ggplot(
  forest_data,
  aes(
    x = estimate,
    y = predictor,
    colour = model
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    position = position_dodge(width = 0.58),
    width = 0,
    linewidth = 0.7,
    orientation = "y"
  ) +
  
  geom_point(
    position = position_dodge(width = 0.58),
    size = 2.6
  ) +
  
  scale_colour_manual(
    values = model_colours,
    name = "Model"
  ) +
  
  labs(
    x = "Standardized coefficient estimate",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black")
  )


save_figure(
  figure_forest,
  "godwit_standardized_coefficients",
  8,
  5.5
)

# ============================================================
# FIGURE 9
# AES subsidies by management and AES practice
# ============================================================


aes_raw <- read_csv(
  aes_file,
  show_col_types = FALSE
) |>
  
  transmute(
    
    field_id =
      as.numeric(
        Field_ID
      ),
    
    AES =
      AES_practice,
    
    AES_subsidies_per_ha =
      as.numeric(
        AES_income_EUR_per_ha
      )
  ) |>
  
  distinct(
    field_id,
    .keep_all = TRUE
  ) |>
  
  left_join(
    
    dat |>
      transmute(
        field_id =
          as.numeric(
            Plot
          ),
        
        management_type
      ),
    
    by =
      "field_id"
  )


# Management panel
aes_management_data <- aes_raw |>
  
  filter(
    !is.na(
      AES_subsidies_per_ha
    ),
    
    management_type %in%
      c(
        "Diversified",
        "Intermediate"
      )
  ) |>
  
  mutate(
    management_type =
      factor(
        management_type,
        levels =
          c(
            "Diversified",
            "Intermediate"
          )
      )
  )


panel_a <- ggplot(
  
  aes_management_data,
  
  aes(
    x = management_type,
    y = AES_subsidies_per_ha,
    fill = management_type
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  
  geom_jitter(
    width = 0.18,
    size = 2,
    alpha = 0.85
  ) +
  
  scale_fill_manual(
    values = c(
      Diversified = "#2E7D32",
      Intermediate = "#F9A825"
    )
  ) +
  
  labs(
    x = NULL,
    y = "AES subsidies per ha (€)"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position =
      "none"
  )


# AES-practice panel
aes_practice_levels <- c(
  "NestProtection",
  "Mosaics",
  "DelayedMowing",
  "HerbRich"
)


aes_practice_data <- aes_raw |>
  
  filter(
    !is.na(
      AES_subsidies_per_ha
    ),
    
    AES %in%
      aes_practice_levels
  ) |>
  
  mutate(
    AES =
      factor(
        AES,
        levels =
          aes_practice_levels
      )
  )


panel_b <- ggplot(
  
  aes_practice_data,
  
  aes(
    x = AES,
    y = AES_subsidies_per_ha,
    fill = AES
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.85
  ) +
  
  scale_fill_manual(
    
    values = c(
      NestProtection = "#F8766D",
      Mosaics = "#7CAE00",
      DelayedMowing = "#00BFC4",
      HerbRich = "#C77CFF"
    ),
    
    drop = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic() +
  
  theme(
    
    legend.position =
      "none",
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


figure_aes <-
  panel_a +
  panel_b +
  patchwork::plot_annotation(
    tag_levels = "A"
  )


save_figure(
  figure_aes,
  "AES_subsidies_distribution",
  11,
  5
)


# ============================================================
# PUBLICATION TABLES
# ============================================================


# ------------------------------------------------------------
# Table: MDS indicators and weights
# ------------------------------------------------------------

table_mds <- mds_table |>
  
  transmute(
    
    Indicator,
    
    Weight =
      round(
        `Indicator weight`,
        3
      ),
    
    Scoring =
      `Scoring function`
  )


write_csv(
  table_mds,
  file.path(
    tables_dir,
    "table_MDS_indicators_weights.csv"
  )
)


# ------------------------------------------------------------
# Table: Godwit candidate-model comparison
# ------------------------------------------------------------

table_model_comparison <- model_comparison |>
  
  select(
    model,
    AIC,
    delta_AIC,
    adjusted_R2,
    n
  ) |>
  
  mutate(
    
    AIC =
      round(
        AIC,
        3
      ),
    
    delta_AIC =
      round(
        delta_AIC,
        3
      ),
    
    adjusted_R2 =
      round(
        adjusted_R2,
        3
      )
  )


write_csv(
  table_model_comparison,
  file.path(
    tables_dir,
    "table_godwit_model_comparison.csv"
  )
)


# ------------------------------------------------------------
# Table: Integrated model M6
# ------------------------------------------------------------

table_m6 <- model_coefficients |>
  
  filter(
    model == "M6"
  ) |>
  
  select(
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    p.value,
    significance
  )


write_csv(
  table_m6,
  file.path(
    tables_dir,
    "table_M6_coefficients.csv"
  )
)


# ============================================================
# COMPLETION
# ============================================================

cat(
  "\nFigures and publication tables created successfully.\n"
)

cat(
  "\nMain figures saved to:\n",
  main_figures_dir,
  "\n"
)

cat(
  "\nSupplementary figures saved to:\n",
  supp_figures_dir,
  "\n"
)

cat(
  "\nPublication tables saved to:\n",
  tables_dir,
  "\n"
)
