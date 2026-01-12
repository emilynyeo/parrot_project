# 🦜 Parrot Project 

Welcome to the **Parrot Project** repository! This project contains scripts and data for the analysis of parrot microbiome samples, focusing on data processing, diversity analysis, and predictive modeling.

We are trying to see if we can identify parrots who have been smuggled from the wild into captivity.

---

## Repository Structure

All analysis scripts are located in the `MelissaAnalysis/` folder. The repository is organized as follows:

### Active Scripts (Current Pipeline)

| Script | Description | Author |
|--------|-------------|--------|
| **01.1_EY_process_and_clean_data.R** | Processes and cleans raw sequence data | Emily |
| **01.2_EY_name_qc_process.R** | Quality control filtering, creates filtered phyloseq objects | Emily |
| **02.2_EY_split_18_classes.R** | Splits 18S data into different eukaryotic classes | Emily |
| **03.2_EY_alpha_diversity.R** | Analyzes alpha diversity between groups | Emily |
| **04.2_EY_beta_diversity.R** | Analyzes beta diversity between groups | Emily |
| **05_EY_ANCOMBC.R** | ANCOMBC2 differential abundance analysis | Emily |
| **05_EY_LDA.R** | LDA differential analysis | Emily |
| **05_EY_Combine_LDA_ANCOMBC.R** | Combines LDA and ANCOMBC results | Emily |
| **05.2_pcoa_levels.R** | PCoA analysis at different taxonomic levels | Emily |
| **06_EY_grp_prediction_models.R** | Prediction models for group classification | Emily |
| **06_EY_grp_prediction_models_18s.R** | 18S-specific prediction models | Emily |
| **06.2_qc_group_pred_models.R** | Quality control for prediction models | Emily |
| **paper_plots.Rmd** | Generates all figures for the paper | Emily |
| **utils.R** | Utility functions used across scripts | Emily |

### Archived Scripts

- **`melissa_scripts/`** - Original scripts by Mel (archived for reference)
- **`old_scripts/`** - Older versions and tutorial scripts (archived)

---

## Analysis Pipeline Flow

The following diagram shows how scripts and data files are connected in the main analysis pipeline:

```
Raw Data
  |
  |-- 00_metadata/          (Sample metadata)
  |-- 00_sequence_data/     (Raw sequence files: 16S, 18S)
  |
  |
[01.1_EY_process_and_clean_data.R]
  |
  |-- Creates initial phyloseq objects
  |-- Outputs to: 01_process_and_clean_data/
  |
  |
[01.2_EY_name_qc_process.R]
  |
  |-- Quality control filtering
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

## Main Pipeline for Figures Script 

To generate all figures for the paper, run these scripts in order:

### Step 1: Data Processing & QC
1. **`01.1_EY_process_and_clean_data.R`**
   - Processes raw sequence data
   - Creates initial phyloseq objects
   - Outputs to `01_process_and_clean_data/`

2. **`01.2_EY_name_qc_process.R`**
   - Performs quality control filtering
   - Creates filtered phyloseq objects with various thresholds
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

For detailed information about dependencies, see `MelissaAnalysis/paper_plots_dependencies.md`.

---

## How to Use

1. Clone this repository:
   ```bash
   git clone <repo-url>
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

## References

- See the `papers/` directory for literature
  - *coda4microbiome: compositional data analysis for microbiome cross‑sectional and longitudinal studies.pdf*

---

## Notes

- This project is under active development; folder and script organization may evolve.
- Old scripts have been archived in `old_scripts/` and `melissa_scripts/` folders for reference.
- All active scripts are Emily's updated versions (indicated by `EY` in filename or `.2` suffix). 