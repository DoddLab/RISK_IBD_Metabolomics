# RISK metabolomics project

This repository contains the R scripts used to process RISK cohort metabolomics data and related analyses to generate the main and extended-data figures reported in the accompanying manuscript.

The code is organized into three components:

1. Untargeted metabolomics data processing, quality control, normalization, and annotation;
2. 16S rRNA sequence processing and microbiome analysis; and
3. Statistical analyses and figure generation for the manuscript.


## Repository contents

| Script | Description |
|---|---|
| `Metabolomics_raw_data_processing.R` | Processes raw LC-MS data for C18 and HILIC chromatography in positive- and negative-ion modes. The workflow uses `xcms` and `DoddLabRawMS` for peak detection and related raw-data processing. |
| `Metabolomics_data_cleaning_normalization.R` | Creates `tidymass` objects, removes low-quality features, imputes missing values, performs SERRF normalization, reviews outliers, and constructs the final merged metabolomics dataset. |
| `Metabolomics_metabolite_annotation.R` | Annotates metabolites by matching accurate mass, retention time, and/or MS/MS spectra against in-house and public spectral libraries. It also merges, dereplicates, and manually curates annotation results. |
| `Metabolomics_lipid_annotation.R` | Performs lipid-focused annotation using MS-DIAL LipidBlast resources, class/adduct constraints, spectral matching, dereplication, and construction of the final lipidomics object. |
| `16S_rRNA_Seq_Data_Processing.R` | Processes 16S rRNA sequencing data with DADA2, assigns taxonomy with SILVA, constructs phyloseq objects and phylogenetic trees, calculates diversity metrics and ordinations, tests differential taxa with MaAsLin2, and evaluates a microbial dysbiosis index. |
| `Source_Code_Fig1_ext_Fig1.R` | Generates Figure 1 and Extended Data Figure 1, including cohort characteristics, data availability summaries, PCA overview plots, and clinical-variable associations. |
| `Source_Code_Fig2_ext_Fig2_5.R` | Generates Figure 2 and Extended Data Figures 2–5, including baseline differential-metabolite analyses, volcano plots, metabolite-class enrichment, pathway-focused visualizations, microbial-metabolite summaries, and selected clinical correlations. |
| `Source_Code_Fig3_ext_Fig6_7.R` | Generates Figure 3 and Extended Data Figures 6–7, including serology–metabolite associations, network analyses, predictive modeling, repeated cross-validation, LASSO feature selection, and comparisons among serology, metabolite, and 16S-based models. |
| `Source_Code_Fig4_ext_Fig8.R` | Generates Figure 4 and Extended Data Figure 8, including cumulative complication analyses, differential metabolite and lipid analyses by Crohn's disease behavior, subclass enrichment, lipid molecular networking, representative phosphatidylethanolamine plots, and the PE index. |
| `Source_Code_Fig5_ext_Fig9.R` | Generates Figure 5 and Extended Data Figure 9, including longitudinal disease-severity models, metabolite-subclass enrichment, tryptophan-metabolite analyses, and development and evaluation of the Microbial Indole Balance Score. |

## Metabolomics processing

### Raw-data processing

`Metabolomics_raw_data_processing.R` processes four LC-MS acquisition modes:

- C18 positive ionization;
- C18 negative ionization;
- HILIC positive ionization; and
- HILIC negative ionization.

The original analysis used `DoddLabRawMS` together with `xcms`. Processing parameters are initialized separately for C18 and HILIC data, and HILIC peak detection uses a 30-ppm setting in the supplied script. The script was originally run on the Stanford Sherlock/SCG computing environment and sets the number of parallel workers to 20.

### Data cleaning and normalization

`Metabolomics_data_cleaning_normalization.R` performs the following major steps for each acquisition mode:

1. merges the peak-area table with the corresponding worklist;
2. constructs a `tidymass` mass-dataset object;
3. removes low-quality or noisy features;
4. imputes missing values;
5. applies SERRF normalization;
6. evaluates sample and feature quality;
7. identifies and manually reviews outliers; and
8. removes designated outliers and creates the final combined analysis object.

Most data-cleaning steps were run locally, whereas SERRF normalization was run on Stanford Sherlock/SCG computing environment.

### Metabolite annotation

`Metabolomics_metabolite_annotation.R` assigns metabolite confidence levels based on combinations of accurate mass, retention time, and MS/MS matching:

- **Level 1:** MS1, retention time, and MS/MS match to the in-house library;
- **Level 2.1:** MS1 and retention-time match to the in-house library; and
- **Level 2.2:** MS1 and MS/MS match to a public library.

The workflow uses the in-house Dodd library and public resources including MS-DIAL, GNPS bile-acid, GNPS acyl-amide, GNPS acyl-ester, and Fiehn peptide libraries. Candidate annotations are merged, prioritized by confidence, dereplicated, and manually reviewed.

### Lipid annotation

`Metabolomics_lipid_annotation.R` performs lipid annotation using MS-DIAL LipidBlast resources and polarity-specific class/adduct combinations. The workflow includes spectral matching, extraction of lipid abbreviations, merging and dereplication of candidates, manual review, and creation of the final lipidomics object used in the disease-behavior analyses.

### Statistics analysis

`Source_Code_Fig_XXXX.R` Statistical analyses were performed in R using methods tailored to study design and data type. Across figures, group-level comparisons primarily used non-parametric tests (for example, Wilcoxon rank-sum/signed-rank) and contingency-table tests where appropriate. Feature-wise differential analyses were commonly fit with generalized linear models (typically adjusting for covariates such as age, sex, race, and antibiotic use). Correlation and trend analyses used Pearson/Spearman correlation and linear models. For prediction analyses, scripts applied LASSO-based feature selection (`cv.glmnet`) and covariate-adjusted logistic regression, with repeated 5-fold cross-validation (100 repeats) and ROC/AUC-based evaluation (`pROC`, including confidence intervals and threshold-derived performance metrics). For complication outcomes, cumulative incidence/competing-risk analyses were performed with `tidycmprsk`/`ggsurvfit`.


## Software requirements

The scripts were written in R and use packages from CRAN, Bioconductor, GitHub, and local laboratory repositories. Major dependencies include:

- data handling and visualization: `tidyverse`, `ggplot2`, `cowplot`, `ggpubr`, `rstatix`, `sjmisc`, `ComplexHeatmap`, and `circlize`;
- metabolomics: `tidymass`, `xcms`, `DoddLabRawMS`, `DoddLabQC`, `DoddLabSERRF`, `DoddLabTool`, `DoddLabMetID`, and `DoddLabDatabase`;
- microbiome analysis: `dada2`, `phyloseq`, `Biostrings`, `DECIPHER`, `phangorn`, `vegan`, and `Maaslin2`;
- statistical modeling and prediction: `lme4`, `lmerTest`, `broom.mixed`, `glmnet`, and `pROC`;
- network and tree visualization: `igraph`, `tidygraph`, `ggraph`, and `ggtree`;
- survival and competing-risk analysis: `tidycmprsk`, `ggsurvfit`, and `gtsummary`.

Custimized in-house Packages that used in this study:

| Package | Version recorded in code |
|---|---:|
| `DoddLabRawMS` | 0.1.16 |
| `DoddLabQC` | 0.1.2 |
| `DoddLabTool` | 0.1.5 |
| `DoddLabMetID` | 0.1.19 |
| `DoddLabDatabase` | 0.2.6 |

These package are available at our group [GitHub](https://github.com/DoddLab) 

## Data availability
- Raw LC-MS/MS data generated in this study have been deposited in [XXXX] under accession number [XXXX]. 
- Processed metabolomics feature tables and MS/MS spectra can be accessed through Zendo [XXX]. 
- Deidentified participant metadata, annotated metabolite and lipid abundance matrices, authentic metabolite reference libraries, experimental MS/MS spectral libraries, detailed LC-MS/MS acquisition and data processing parameters, and complete statistical analysis results are provided in **Supplementary Data Files 1–11** with this paper. **Source data** are provided with this paper.
- **Note:** Previously generated 16S rRNA sequencing data and other paired clinical data from the RISK cohort are available through the **Crohn’s & Colitis Foundation IBD Plexus program** upon approval of a data access request. 

## Citation

When using this code, please cite the accompanying manuscript:

- Zhiwei Zhou, Dylan Dodd*, Plasma metabolomics in a pediatric Crohn’s disease inception cohort links host-microbial tryptophan metabolism to disease activity, Submitted, 2026


## Contact
For questions or issues, please contact **Zhiwei Zhou** (zhouzw@stanford.edu).

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
