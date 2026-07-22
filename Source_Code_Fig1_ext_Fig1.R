setwd('~/Project/00_IBD_project/Data/20241018_enrollment_sample_info/')

################################################################################
# Ext Figure 1

# Enrollment sample information ------------------------------------------------

load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_demographic_241003.RData')

temp_patient_id <- unique(patient_meta_info_demographic$patient_id)

data_patient_enrollment <- lapply(temp_patient_id, function(x){
  temp_data <- patient_meta_info_demographic %>% 
    dplyr::filter(patient_id == x) %>% 
    dplyr::arrange(visit_month)
  
  use_antibiotics <- temp_data$use_antibiotics[1]
  enrollment_age <- temp_data$age[1]
  race <- temp_data$race[1]
  gender <- temp_data$gender[1]
  phenotype_group1 <- temp_data$phenotype_group1[1]
  phenotype_group2 <- temp_data$phenotype_group2[1]
  
  visit_month_0 <- any(temp_data$visit_month == 0)
  visit_month_12 <- any(temp_data$visit_month == 12)
  visit_month_24 <- any(temp_data$visit_month == 24)
  visit_month_36 <- any(temp_data$visit_month == 36)
  
  result <- data.frame(patient_id = x, enrollment_age, race, gender, use_antibiotics, phenotype_group1, phenotype_group2,
                       visit_month_0, visit_month_12, visit_month_24, visit_month_36, stringsAsFactors = FALSE)
}) %>% 
  dplyr::bind_rows()


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(sjmisc)


temp_data <- data_patient_enrollment %>% 
  dplyr::arrange(match(phenotype_group1, c('non_IBD', 'CD')), gender, race, patient_id)

temp_plot_data <- temp_data  %>% 
  select(patient_id, visit_month_0:visit_month_36) %>% 
  column_to_rownames('patient_id') %>%
  mutate(visit_month_0 = ifelse(visit_month_0 == TRUE, 1, 0),
         visit_month_12 = ifelse(visit_month_12 == TRUE, 1, 0),
         visit_month_24 = ifelse(visit_month_24 == TRUE, 1, 0),
         visit_month_36 = ifelse(visit_month_36 == TRUE, 1, 0)) %>% 
  rotate_df() %>% 
  as.matrix()


# define colors
sample_colors <- sequential_hcl(2, palette = 'Blue-Yellow', h = c(250, 60), c = 100, l = c(30, 90))
col_fun = colorRamp2(c(0, 1), c("#F4FAFE", "#06AED5"))

race_color <- c('Unknown' = "#ED90A4",
                'Hispanic/Latino' = "#D8A06A",
                'Asian' = "#ABB150",
                'Native Hawaiian or Other Pacific Islander' = "#62BE79",
                'American Indian or Alaska Native' = "#00C1B2",
                "Black or African American" = "#48B8DE", 
                'Caucasian' = "#ACA2EC", 
                'Others' = "#E190D6")

gender_color <-  c("#90DBF4", "#FDE4CF")
names(gender_color) <- c('Male', 'Female')


test = HeatmapAnnotation(
  phenotype = temp_data$phenotype_group1,
  gender = temp_data$gender,
  race = temp_data$race,
  enrollment_age = anno_points(temp_data$enrollment_age, extend = 0.1, labels_rot = TRUE),
  col = list(enrollment_age = 'black',
             gender = gender_color,
             race = race_color,
             phenotype = c("CD" = "tomato", "non_IBD" = "dodgerblue")
  )
)



# draw(test)
temp_plot <- Heatmap(temp_plot_data, 
                     name = 'Enrollment', 
                     col = col_fun, 
                     # rect_gp = gpar(col = "white", lwd = 0.3),
                     show_row_names = TRUE, 
                     show_column_names = FALSE, 
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE, 
                     top_annotation = test)
dir.create('~/Project/00_IBD_project/Figure/250326', showWarnings = FALSE)
pdf('~/Project/00_IBD_project/Figure/250326/20250326_heatmap_enrollment_demographic.pdf', 
    width = 20, height = 5)
draw(temp_plot)
dev.off()

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_data, file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1a.csv')



################################################################################
# age boxplot ------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

set.seed(250326)

temp_data <- patient_meta_info_demographic %>% 
  filter(visit_number == 'enrollment') %>% 
  mutate(phenotype_group1 = recode_factor(phenotype_group1, 'non_IBD' = 'non_IBD', 'CD' = 'CD')) %>% 
  as.data.frame()

temp_data %>% group_by(phenotype_group1) %>% summarise(age = mean(age, na.rm = TRUE))

group_comparisons <- list(c("CD", "non_IBD"))

temp_plot <- ggplot(temp_data,
                    aes(x = phenotype_group1, y = age)) +
  geom_boxplot(aes(fill = phenotype_group1), 
               na.rm = TRUE, outliers = FALSE) +
  geom_jitter(width = 0.2, alpha = 1) +
  stat_compare_means(comparisons = group_comparisons, 
                     method = 'wilcox.test') +
  scale_fill_manual(values = c('CD' = 'tomato', 'non_IBD' = 'dodgerblue')) +
  scale_x_discrete(labels = c('CD' = 'CD', 'non_IBD' = 'Non-IBD')) +
  ylab('Age (years)') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = 'none')

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure1/20250326_boxplot_enrollment_age.pdf', 
       width = 3, height = 6)

# temp_data <- temp_data %>% dplyr::select(sample_id, age, phenotype_group1)
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_data, file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_age.csv')


################################################################################
# state the data availability ---------------------------------------------------------
library(ggpubr)

load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_raw_241002.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_all_241002.RData')

# 16s ------------------------------------------------------------------------
temp_plot_data_16s <- patient_meta_info_all %>% 
  dplyr::filter(visit_number == 'enrollment') %>%
  dplyr::mutate('data_16s' = case_when(
    !is.na(`16s`) ~ '16s',
    TRUE ~ 'unavailable'
  )) %>%
  count(data_16s)

label <- paste0(temp_plot_data_16s$data_16s, ": ", round(temp_plot_data_16s$n/(sum(temp_plot_data_16s$n))*100, digits = 1), '%') %>% rev()

temp_plot <- ggdonutchart(temp_plot_data_16s, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "data_16s", 
                          color = "white",
                          palette = c('16s' = "#8cd3c7", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20241028_other_omics_data_statistics/16s_donut_241028.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_plot_data_16s, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_16s.csv')

# serology ---------------------------------------------------------------------
temp_plot_data_serology <- patient_meta_info_all %>% 
  dplyr::filter(sample_id != 'B032_S34') %>% 
  dplyr::filter(visit_number == 'enrollment') %>%
  dplyr::mutate('data_serology' = case_when(
    !is.na(iga_asca) ~ 'serology',
    TRUE ~ 'unavailable'
  )) %>%
  count(data_serology)


label <- paste0(temp_plot_data_serology$data_serology, ": ", round(temp_plot_data_serology$n/(sum(temp_plot_data_serology$n))*100, digits = 1), '%') %>% rev()

temp_plot <- ggdonutchart(temp_plot_data_serology, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "data_serology", 
                          color = "white",
                          palette = c('serology' = "#b3de68", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20241028_other_omics_data_statistics/Serology_donut_241028.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# write_csv(temp_plot_data_serology, 
#            file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_serology.csv')


# Endoscopy availability -------------------------------------------------------

load('~/Project/00_IBD_project/Data/20260126_severity_other_criteria/endoscopy_result_table_with_index_cdeis_260127.RData')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')



# assign the severity_cdeis: remission (CDEIS < 3), non_remission (CDEIS >= 3)
endoscopy_summary <- endoscopy_result_table_with_index_cdeis %>% 
  dplyr::left_join(
    object_final %>%
      extract_sample_info() %>%
      dplyr::select(sample_id, phenotype_group1, visit_number, severity_class),
    by = c('sample_id')
  ) %>%
  dplyr::mutate(severity_cdeis = ifelse(CDEIS < 3, 'remission', 'non_remission')) %>% 
  dplyr::rename(severity_wpcdai = severity_class) %>%
  dplyr::select(sample_id, deidentified_master_patient_id, visit_encounter_id, CDEIS, severity_cdeis, severity_wpcdai, phenotype_group1, visit_number)

# Enrollment samples
temp_plot_data_endoscopy <- endoscopy_summary %>% 
  dplyr::filter(visit_number == 'enrollment') %>%
  dplyr::mutate('data_endoscopy' = case_when(
    !is.na(CDEIS) ~ 'Endoscopy',
    TRUE ~ 'unavailable'
  )) %>%
  count(data_endoscopy)

label <- paste0(temp_plot_data_endoscopy$data_endoscopy, ": ", round(temp_plot_data_endoscopy$n/(sum(temp_plot_data_endoscopy$n))*100, digits = 1), '%') %>% rev()

library(ggpubr)
temp_plot <- ggdonutchart(temp_plot_data_endoscopy, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "data_endoscopy", 
                          color = "white",
                          palette = c('Endoscopy' = "#fdb462", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260318_update_figure/endoscopy_donut_260318.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE))
# readr::write_csv(temp_plot_data_endoscopy, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_endoscopy.csv')


# Blood test availability -------------------------------------------------------

temp_plot_data_blood_test <- patient_meta_info_all %>% 
  dplyr::filter(sample_id != 'B032_S34') %>% 
  dplyr::mutate('blood_test' = case_when(
    (!is.na(c_reactive_protein) & !is.na(erythrocyte_sedimentation_rate_esr) & !is.na(albumin)) ~ 'blood_test',
    TRUE ~ 'unavailable'
  )) %>%
  count(blood_test)


label <- paste0(temp_plot_data_blood_test$blood_test, ": ", round(temp_plot_data_blood_test$n/(sum(temp_plot_data_blood_test$n))*100, digits = 1), '%') %>% rev()

temp_plot <- ggdonutchart(temp_plot_data_blood_test, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "blood_test", 
                          color = "white",
                          palette = c('blood_test' = "#fb8072", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20241028_other_omics_data_statistics/blood_test_donut_241028.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_plot_data_blood_test, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_blood_test.csv')


# Anti-TNF availability -------------------------------------------------------

temp_plot_data_anti_TNF <- patient_meta_info_all %>% 
  dplyr::filter(sample_id != 'B032_S34') %>% 
  dplyr::mutate('anti_TNF_alpha' = case_when(
    (!is.na(ongoing_infliximab_unspecified) | !is.na(infliximab_unspecified) | !is.na(ongoing_adalimumab) | !is.na(adalimumab) | !is.na(ongoing_certolizumab_pegol) | !is.na(certolizumab_pegol) | !is.na(ongoing_natalizumab) | !is.na(natalizumab)) ~ 'receive_anti_TNF_alpha',
    TRUE ~ 'unavailable'
  )) %>%
  count(anti_TNF_alpha)


label <- paste0(temp_plot_data_anti_TNF$anti_TNF_alpha, ": ", round(temp_plot_data_anti_TNF$n/(sum(temp_plot_data_anti_TNF$n))*100, digits = 1), '%') %>% rev()

temp_plot <- ggdonutchart(temp_plot_data_anti_TNF, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "anti_TNF_alpha", 
                          color = "white",
                          palette = c('receive_anti_TNF_alpha' = "#7fb1d3", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20241028_other_omics_data_statistics/anti_TNF_alpha_donut_241028.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_plot_data_anti_TNF, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_anti_TNF_alpha.csv')



# Disease severity ---------------------------------------------------------------
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_demographic_241003.RData')
temp_plot_data_severity <- patient_meta_info_demographic %>% 
  dplyr::mutate('severity' = case_when(
    severity_class != 'unavailable' ~ 'severity_evaluation',
    TRUE ~ 'unavailable'
  )) %>%
  count(severity)


label <- paste0(temp_plot_data_severity$severity, ": ", round(temp_plot_data_severity$n/(sum(temp_plot_data_severity$n))*100, digits = 1), '%') %>% rev()

temp_plot <- ggdonutchart(temp_plot_data_severity, "n", 
                          label = label,
                          lab.pos = "in",
                          fill = "severity", 
                          color = "white",
                          palette = c('severity_evaluation' = "#bebada", 'unavailable' = "#c4c7c9"))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20241028_other_omics_data_statistics/severity_test_donut_241028.pdf', 
       width = 6, height = 6)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_plot_data_severity,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_extended_fig1_severity.csv')



################################################################################
# PCA overview analysis ------------------------------------------------------
load('~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/object_total_merge_241024.RData')
# load('~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/pca_merge_total_PCA_241024.RData')

library(sjmisc)
library(mixOmics)
temp_object <- object_total_merge %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(p_value_adjust <= 0.05)
temp_data_trans <- temp_object@expression_data %>% 
  rotate_df()
pca_merge <- pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$phenotype_group2
pca_var_id <- temp_object@sample_info$sample_id
save(pca_merge, 
     file = '~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/pca_merge_total_PCA_241024.RData')

# PCA view
temp_plot_data <- pca_merge$variates$X %>% 
  as.data.frame()

temp_plot_data <- temp_plot_data %>% 
  rownames_to_column(var = 'sample_id') %>% 
  left_join(temp_object@sample_info, by = 'sample_id')

pca_data <- temp_plot_data


# visualization colored by IBD visit_phenotype -------------------------------
temp_plot_data <- temp_plot_data %>% 
  mutate(visit_phenotype = paste(visit_number, phenotype_group1, sep = '_'))
temp_plot <- ggplot(temp_plot_data, aes(x = PC1, y = PC2, colour = visit_phenotype)) +
  geom_point(shape = 16, size = 3, alpha = 0.7) + 
  # stat_ellipse(type = "norm", linetype=1, size = 1) + 
  labs(x=paste("PC1 (", format(pca_merge$prop_expl_var$X[1]*100, digits=4), "%)", sep=""),
       y=paste("PC2 (", format(pca_merge$prop_expl_var$X[2]*100, digits=4), "%)", sep="")) +
  scale_colour_manual(values = c('enrollment_non_IBD' = '#8dd2c7',
                                 'enrollment_CD' = '#ffed70',
                                 '12_month_CD' = '#bebada',
                                 '24_month_CD' = '#fb8071',
                                 '36_month_CD' = '#80b1d3'),
                      name = 'Group') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.1, 0.1), 
        axis.text.y.left = element_text(hjust = 0.5, vjust = 0.5), 
        axis.ticks.length = unit(2, "mm"))


library(cowplot)
temp_plot2 <- axis_canvas(temp_plot, axis = "x") +
  geom_density(data = temp_plot_data, aes(x = PC1, fill = visit_phenotype), alpha = 0.5) +
  scale_fill_manual(values = c('enrollment_non_IBD' = '#8dd2c7',
                               'enrollment_CD' = '#ffed70',
                               '12_month_CD' = '#bebada',
                               '24_month_CD' = '#fb8071',
                               '36_month_CD' = '#80b1d3'),
                    name = 'Group') +
  ZZWtool::ZZWTheme(type = 'classic') 


temp_plot3 <- temp_plot %>%
  insert_xaxis_grob(temp_plot2, grid::unit(1, "in"), position = "top") %>%
  ggdraw()

ggsave(plot = temp_plot3, 
       filename = '~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/global_pca_plot_with_density_and_visit_number_241024.pdf', 
       width = 6, height = 6.5)

temp_data <- temp_plot_data %>% 
  dplyr::select(sample_id, visit_number, phenotype_group1, PC1, PC2)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE)
# readr::write_csv(temp_data, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_fig1_pca.csv')


# PC1 boxplots for different visits --------------------------------------------

temp_plot_data$visit_phenotype <- recode_factor(temp_plot_data$visit_phenotype, 
                                                'enrollment_non_IBD' = 'enrollment_non_IBD',
                                                'enrollment_CD' = 'enrollment_CD',
                                                '12_month_CD' = '12_month_CD',
                                                '24_month_CD' = '24_month_CD',
                                                '36_month_CD' = '36_month_CD')

list_comparsion <- list(c('enrollment_CD', 'enrollment_non_IBD'),
                        c('12_month_CD', 'enrollment_non_IBD'),
                        c('12_month_CD', 'enrollment_CD'),
                        c('12_month_CD', '24_month_CD'),
                        c('24_month_CD', '36_month_CD'))

temp_plot <- ggplot(temp_plot_data,
                    aes(x = visit_phenotype, y = PC1, fill = visit_phenotype)) +
  geom_jitter(position = position_jitter(height = 0.05),
              shape = 21, color = 'black', stroke = 0,
              size = 2, alpha = 1) +
  geom_boxplot(outliers = FALSE, alpha = 0.7) +
  stat_compare_means(comparisons = list_comparsion, 
                     method = 'wilcox.test', 
                     label.y = c(105, 120, 105, 105, 105)) +
  scale_fill_manual(values = c('enrollment_non_IBD' = '#8dd2c7',
                               'enrollment_CD' = '#ffed70',
                               '12_month_CD' = '#bebada',
                               '24_month_CD' = '#fb8071',
                               '36_month_CD' = '#80b1d3'),
                    name = 'Group') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = 'none')

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure1/PC1_boxplot_with_visit_number_2_250326.pdf', 
       width = 6, height = 6)


################################################################################
load('~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/object_total_merge_241024.RData')
load('~/Project/00_IBD_project/Data/20241024_enrollment_total_PCA/pca_merge_total_PCA_241024.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/01_table_lab_test/table_lab_tests_230621.RData')

library(tidyverse)
library(mixOmics)
library(sjmisc)

temp_object <- object_total_merge %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(p_value_adjust <= 0.05)
# temp_data_trans <- temp_object@expression_data %>% 
#   rotate_df()
# pca_merge <- pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$phenotype_group1
pca_var_id <- temp_object@sample_info$sample_id

# PCA view
temp_plot_data <- pca_merge$variates$X %>% 
  as.data.frame()

temp_plot_data <- temp_plot_data %>% 
  rownames_to_column(var = 'sample_id') %>% 
  left_join(temp_object@sample_info, by = 'sample_id')

pca_data <- temp_plot_data



pca_data <- pca_data %>% 
  left_join(table_lab_tests, 
            by = c('subject.id' = 'deidentified_master_patient_id',
                   'visit_encounter_id' = 'visit_encounter_id'))

pca_data <- pca_data %>% 
  dplyr::rename('esr' = 'Erythrocyte Sedimentation Rate (ESR)',
                'platelet_count' = 'Platelet Count',
                'lymphocytes' = 'Lymphocytes',
                'albumin' = 'Albumin',
                'urea' = 'Urea',
                'alt' = 'Alanine Aminotransferase (ALT)',
                'white_blood_cell_count' = 'White Blood Cell Count',
                'hematocrit' = 'Hematocrit',
                'eosinophil' = 'Eosinophil',
                'crp' = 'C Reactive Protein',
                'ggt' = 'Gamma-Glutamyl Transferase (GGT)',
                'creatinine' = 'Creatinine',
                'ast' = 'Aspartate Aminotransferase (AST)',
                'hemoglobin' = 'Hemoglobin',
                'neutrophil' = 'Neutrophil',
                'alp' = 'Alkaline Phosphatase (ALP)',
                'cbir_fla' = 'CBir Fla',
                'i2' = 'I2',
                'iga_asca' = 'IgA ASCA',
                'gm_csf' = 'GM-CSF',
                'igg_asca' = 'IgG ASCA',
                'ompc' = 'OmpC',
                'anca' = 'ANCA')

# CRP linear plot ------------------------------------------------------------------

temp_mean <- mean(pca_data$crp, na.rm = TRUE)
temp_sd <- sd(pca_data$crp, na.rm = TRUE)

pca_data <- pca_data %>%  
  mutate(visit_phenotype = paste(visit_number, phenotype_group1, sep = '_')) %>% 
  mutate(crp_normalized = crp/temp_sd)

temp_lm <- lm(log2(crp+1)~PC1, data = pca_data)
temp_lm2 <- summary(temp_lm)

temp_data <- pca_data %>% 
  dplyr::filter(!is.na(crp)) %>% 
  dplyr::filter(!is.na(PC1))

cor(log2(temp_data$crp+1), temp_data$PC1, method = 'pearson')


temp_plot <- ggplot(pca_data) +
  geom_point(aes(x = PC1, y = log2(crp + 1), color = visit_phenotype), 
             size = 2, alpha = 1) +
  geom_smooth(aes(x = PC1, y = log2(crp + 1)), method = 'lm', color = 'black') +
  scale_color_manual(values = c('enrollment_non_IBD' = '#8dd2c7',
                                'enrollment_CD' = '#ffed70',
                                '12_month_CD' = '#bebada',
                                '24_month_CD' = '#fb8071',
                                '36_month_CD' = '#80b1d3'), 
                     na.value = 'gray',
                     name = 'Group') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = 'none')

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure1/lm_pc1_crp_250326.pdf', 
       width = 6, height = 6)

temp_data <- temp_data %>% 
  dplyr::select(sample_id, visit_number, phenotype_group1, PC1, crp)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data', showWarnings = FALSE))
# write_csv(temp_data, 
#           file = '~/Project/00_IBD_project/Data/20260718_source_data/20260718_source_data_ext_fig1_pc1_crp.csv')
