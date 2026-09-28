# Research-project-1
Scritps and data source for project analysis

# Soil Health and Black-tailed Godwit in Dutch Dairy Grasslands

This repository contains the data and statistical analysis  developed the first Research project within the MSc in Ecology and Evolution, at the University of Groningen (RUG). 

## Project overview

The study investigates whether field-level soil health, agricultural management, and Agri-Environmental Schemes (AES) are associated with territorial densities of Black-tailed Godwits (*Limosa limosa*) in dairy grasslands of South-West Friesland, the Netherlands.

A total of 34 grassland fields were studied across three management categories depending on the origin of their subsidies, vegetation structure, and management practices:

- **Diversified**
- **Intermediate**
- **Intensive**

The analysis integrates soil biological, chemical, and physical indicators with management types, AES subsidies, and territorial godwit density.

## Main objectives

1. Construct a context-specific **Soil Health Index (SHI)**.
2. Assess relationships between soil conditions, management type, and territorial godwit density.
3. Evaluate whether AES payments align with soil health and godwit occurrence.
4. Compare the explanatory contribution of soil health, management, soil texture, and AES.

## Soil Health Index

The SHI was developed using a Minimum Data Set approach based on Pearson's correlation analysis and Principal Component Analysis (PCA).

The final index included:

- Soil Organic Matter (SOM)
- β-glucosidase activity
- Fungal Shannon diversity
- Urease activity

Soil texture, expressed as clay+silt content, was included separately as an environmental covariate.

## Statistical analysis

Analyses were conducted in **R** and included:

- Exploratory soil-indicator analyses
- Pearson correlations
- Principal Component Analysis
- Soil Health Index construction
- Linear models
- Model comparison using AIC
- Standardized coefficient comparisons

Territorial godwit density was analysed using a `log1p` transformation to account for zero observations and right-skewed distributions.

## Repository structure

```text
├── data/              # Input and processed datasets
├── scripts/           # R scripts used for analysis
├── figures/           # Figures generated from the analyses
├── results/           # Model outputs and summary tables
└── README.md
