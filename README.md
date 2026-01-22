# 🦜 Parrot Project 

Welcome to the **Parrot Project** repository! This project contains scripts and data for the analysis of parrot microbiome samples, focusing on data processing, diversity analysis, and predictive modeling.

We are trying to see if we can identify parrots who have been smuggled from the wild into captivity.

---

## Repository Structure

All analysis scripts are located in the `MelissaAnalysis/` folder. The repository is organized as follows:

### Active Scripts (Current Pipeline)

| Script | Description |
|--------|-------------|
| **01.2_EY_name_qc_process.R** | Processes raw sequence data, performs quality control filtering, creates filtered phyloseq objects |
| **02.2_EY_split_18_classes.R** | Splits 18S data into different eukaryotic classes |
| **03.2_EY_alpha_diversity.R** | Analyzes alpha diversity between groups |
| **04.2_EY_beta_diversity.R** | Analyzes beta diversity between groups |
| **05_EY_ANCOMBC.R** | ANCOMBC2 differential abundance analysis |
| **05_EY_LDA.R** | LDA differential analysis |
| **05_EY_Combine_LDA_ANCOMBC.R** | Combines LDA and ANCOMBC results |
| **05.2_pcoa_levels.R** | PCoA analysis at different taxonomic levels |
| **06_EY_grp_prediction_models.R** | Prediction models for group classification |
| **06_EY_grp_prediction_models_18s.R** | 18S-specific prediction models |
| **06.2_qc_group_pred_models.R** | Quality control for prediction models |
| **paper_plots.Rmd** | Generates all figures for the paper |
| **utils.R** | Utility functions used across scripts |

### Binary Analysis Scripts (Captive vs Wild Comparison)

These scripts perform binary comparisons between Captive and Wild (free ranging) groups, excluding "Wild, seized from traffickers" samples:

| Script | Description |
|--------|-------------|
| **`binary_scripts/02.2_EY_split_classes.R`** | Flexible version that accepts a list of phyloseq objects and processes both 16S and 18S data, splitting into classes (plants, microeukaryotes, nohost) |
| **`binary_scripts/05_LDA_BINARY.R`** | LDA differential analysis for binary comparison (Captive vs Wild, excluding seized samples). Accepts a list of phyloseq objects. |
| **`binary_scripts/06_ANCOMBC_BINARY.R`** | ANCOMBC2 differential abundance analysis for binary comparison (Captive vs Wild, excluding seized samples). Accepts a list of phyloseq objects. |
| **`binary_scripts/07_Combine_LDA_ANCOMBC_BINARY.R`** | Combines binary LDA and ANCOMBC results for visualization. Processes all available dataset/taxonomic level combinations and saves input data as CSVs. |

**Note:** Each binary script has a corresponding `*_EXAMPLE.R` file showing usage.

### Archived Scripts

- **`melissa_scripts/`** - Original scripts by Mel (archived for reference)
- **`old_scripts/`** - Older versions and tutorial scripts (archived)

---

## Analysis Flow

The following diagram shows how scripts and data files are connected in the main analysis:

```
Raw Data
  |
  |-- 00_metadata/          (Sample metadata)
  |-- 00_sequence_data/     (Raw sequence files: 16S, 18S)
  |
  |
[01.2_EY_name_qc_process.R]
  |
  |-- Processes raw sequence data
  |-- Removes problematic samples (S082895, S082806)
  |-- Filters host DNA, chloroplasts, mitochondria
  |-- Performs quality control filtering
  |-- Creates coverage analysis plots
  |-- Outputs to: 01.1_qc_checks/
  |   * phyloseq_otu16m_mat_filt_20_thr10.rds
  |   * phyloseq_otu16n_mat_filt_20_thr10.rds
  |   * phyloseq_otu18m_nohost_mat_filt_20_thr10.rds
  |   * phyloseq_otu18n_nohost_mat_filt_20_thr10.rds
  |
  |
  +---> [02.2_EY_split_18_classes.R] --> 02_split_18_classes/
  |     (Splits 18S into nohost, plants, microeuk)
  |
  +---> [03.2_EY_alpha_diversity.R] --> 03_alpha_diversity/
  |     (Alpha diversity analysis)
  |
  +---> [04.2_EY_beta_diversity.R] --> 04_beta_diversity/
  |     (Beta diversity analysis)
  |
  +---> [05_EY_ANCOMBC.R] --> 05_taxa_driving_groups/
  |     * ancom_16m_16n_18m_18m_trimmed.RData
  |
  +---> [05_EY_LDA.R] --> 05_taxa_driving_groups/LDA/
  |     * lda4_16m_16n_18m_18n.RData
  |
  +---> [05_EY_Combine_LDA_ANCOMBC.R]
  |     (Combines results from ANCOMBC and LDA)
  |
  +---> [05.2_pcoa_levels.R]
  |     (PCoA analysis at different taxonomic levels)
  |
  +---> [06_EY_grp_prediction_models.R] --> 06_grp_prediction_models/
  |     (16S prediction models)
  |
  +---> [06_EY_grp_prediction_models_18s.R]
  |     (18S prediction models)
  |
  +---> [06.2_qc_group_pred_models.R] --> 06.2_qc_group_pred_models/
  |     (QC for prediction models)
  |
  |
  +---> [binary_scripts/] (Binary Analysis Pipeline)
  |     |
  |     +---> [02.2_EY_split_classes.R] --> 02_split_18_classes/
  |     |     (Flexible class splitting for lists of phyloseq objects)
  |     |
  |     +---> [05_LDA_BINARY.R] --> 05_taxa_driving_groups/LDA/binary/
  |     |     (Binary LDA: Captive vs Wild, excluding seized)
  |     |
  |     +---> [06_ANCOMBC_BINARY.R] --> 05_taxa_driving_groups/ANCOMBC/binary/
  |     |     (Binary ANCOMBC: Captive vs Wild, excluding seized)
  |     |
  |     +---> [07_Combine_LDA_ANCOMBC_BINARY.R] --> 05_taxa_driving_groups/combined_lda_ancombc/binary/
  |           (Combines binary LDA and ANCOMBC results, saves CSVs)
  |
  |
[paper_plots.Rmd]
  |
  |-- Requires:
  |   * Phyloseq objects from 01.1_qc_checks/
  |   * ANCOMBC results from 05_taxa_driving_groups/
  |   * LDA results from 05_taxa_driving_groups/LDA/
  |   * utils.R (utility functions)
  |
  |-- Outputs: paper_plots.html (all paper figures)
```

---

## Main Upstream Scripts for the Figures Script 

To generate all figures for the paper, run these scripts in order:

### Step 1: Data Processing & QC
1. **`01.2_EY_name_qc_process.R`**
   - Processes raw sequence data from `00_sequence_data/` and `00_metadata/`
   - Removes problematic samples (S082895, S082806)
   - Filters host DNA, chloroplasts, and mitochondria
   - Performs quality control filtering with multiple thresholds
   - Creates coverage analysis plots
   - Outputs the 4 required phyloseq RDS files to `01.1_qc_checks/`

### Step 2: Differential Abundance Analysis
3. **`05_EY_ANCOMBC.R`**
   - Runs ANCOMBC2 differential abundance analysis
   - Analyzes all 4 datasets at family, genus, and species levels
   - Outputs: `05_taxa_driving_groups/ancom_16m_16n_18m_18m_trimmed.RData`

4. **`05_EY_LDA.R`**
   - Runs LDA (Linear Discriminant Analysis) differential analysis
   - Analyzes all 4 datasets at family, genus, and species levels
   - Outputs: `05_taxa_driving_groups/LDA/lda4_16m_16n_18m_18n.RData`

### Step 3: Generate Paper Plots
5. **`paper_plots.Rmd`**
   - Creates all figures for the paper
   - Requires all files from Steps 1-2
   - Uses `utils.R` for helper functions

---

## Binary Analysis Pipeline (Captive vs Wild)

For binary comparisons excluding "Wild, seized from traffickers" samples, use the scripts in `binary_scripts/`:

### Step 1: Load Phyloseq Objects
Load your phyloseq objects (from `01.1_qc_checks/`) into a named list. See `*_EXAMPLE.R` files for usage.

### Step 2: Run Binary Differential Analysis
1. **`binary_scripts/05_LDA_BINARY.R`**
   - Performs LDA analysis for binary comparison
   - Filters out "Wild, seized from traffickers" samples
   - Analyzes at Family, Genus, and Species levels
   - Outputs: `05_taxa_driving_groups/LDA/binary/lda_binary_results.RData`

2. **`binary_scripts/06_ANCOMBC_BINARY.R`**
   - Performs ANCOMBC2 analysis for binary comparison
   - Filters out "Wild, seized from traffickers" samples
   - Analyzes at Family, Genus, and Species levels
   - Outputs: `05_taxa_driving_groups/ANCOMBC/binary/ancombc_binary_results.RData`

### Step 3: Combine and Visualize Results
3. **`binary_scripts/07_Combine_LDA_ANCOMBC_BINARY.R`**
   - Combines LDA and ANCOMBC binary results
   - Processes all available dataset/taxonomic level combinations
   - Creates combined visualizations
   - Saves input data as CSVs with descriptive filenames:
     * `LDA_ANCOMBC_[dataset]_[tax_level]_LDA_input_data.csv`
     * `LDA_ANCOMBC_[dataset]_[tax_level]_ANCOMBC_input_data.csv`
     * `LDA_ANCOMBC_[dataset]_[tax_level]_combined_LDA_ANCOMBC.csv`
     * `LDA_ANCOMBC_[dataset]_[tax_level]_feature_overlap.csv` (only overlapping features)
   - Outputs plots: `*_combined_plot.png` and `*_separate_plots_p1_p2_p3.png`

---

## How to Use

1. Clone this repository:
   ```bash
   git clone https://github.com/emilynyeo/parrot_project.git
   ```

2. Navigate to the project directory:
   ```bash
   cd parrot_project
   ```

3. Set working directory in R:
   ```r
   setwd("MelissaAnalysis")
   ```

4. Run scripts in order for the full analysis pipeline, or use individual scripts as needed.

---

## Script Dependencies

### Core Dependencies
- **`utils.R`** - Required by `paper_plots.Rmd` and `05_EY_Combine_LDA_ANCOMBC.R`
- **Phyloseq objects** - Created by `01.2_EY_name_qc_process.R`, used by all analysis scripts

### Analysis Dependencies
- **Alpha/Beta Diversity** scripts depend on filtered phyloseq objects
- **ANCOMBC and LDA** scripts depend on filtered phyloseq objects
- **Prediction models** depend on filtered phyloseq objects
- **Paper plots** depend on ANCOMBC and LDA results

---

## Notes

- This project is under active development; folder and script organization may evolve.
- Old scripts have been archived in `old_scripts/` and `melissa_scripts/` folders for reference.
- All active scripts are Emily's updated versions (indicated by `EY` in filename or `.2` suffix). 