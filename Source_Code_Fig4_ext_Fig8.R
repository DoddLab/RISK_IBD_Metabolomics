################################################################################
# Cumulative complication plot -------------------------------------------------

library(sjmisc)
library(tidyverse)

load('~/Project/00_IBD_project/Data/20241205_CD_complicates/behaviour_time_table_241205.RData')

# cumulative number of patients with complication accoding to the visit types
behaviour_change_time_table_CD_patient <- behaviour_time_table_update %>%
  mutate(complication = case_when(
    complication == 'Stricturing;Penetrating' ~ 'Stricturing;Penetrating',
    complication == 'Penetrating;Stricturing' ~ 'Stricturing;Penetrating',
    TRUE ~ complication
  ))

behaviour_change_time_table_CD_patient <- behaviour_change_time_table_CD_patient %>%
  mutate(encounter_month_change_behavior = as.numeric(stringr::str_extract(encounter_type_change_behavior, '\\d+')))

behaviour_change_time_table_summary <- behaviour_change_time_table_CD_patient %>%
  mutate(encounter_group = case_when(
    encounter_month_change_behavior <= 12 ~ '0-12',
    encounter_month_change_behavior <= 24 ~ '13-24',
    encounter_month_change_behavior <= 36 ~ '25-36',
    encounter_month_change_behavior > 36 ~ '>36',
    TRUE ~ 'none'
  )) %>%
  group_by(complication, encounter_group) %>%
  count() %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = encounter_group, values_from = n)


# col 1: 0-12, col 2: 13-24, col 3: 25-36, col 4: >36
# Col 5: 0-12 + 13-24, Col 6: 0-12 + 13-24 + 25-36, Col 7: 0-12 + 13-24 + 25-36 + >36
# add a row to sum each column

behaviour_change_time_table_summary <- behaviour_change_time_table_summary %>%
  tidyr::replace_na(list(`0-12` = 0, `13-24` = 0, `25-36` = 0, `>36` = 0)) %>%
  slice(1:3) %>%
  # select(-none) %>%
  mutate(cummu_year0 = 0,
         cummu_year1 = `0-12`,
         cummu_year2 = `0-12` + `13-24`,
         cummu_year3 = `0-12` + `13-24` + `25-36`,
         cummu_year4 = `0-12` + `13-24` + `25-36` + `>36`) %>%
  column_to_rownames('complication')

temp_sum <- apply(behaviour_change_time_table_summary, 2, sum) %>% as.data.frame() %>% rotate_df()
temp_non_complication <- nrow(behaviour_change_time_table_CD_patient) - temp_sum
rownames(temp_non_complication) <- 'none_complication'

behaviour_change_time_cumulative <- behaviour_change_time_table_summary %>%
  bind_rows(temp_non_complication) %>%
  select(cummu_year0:cummu_year4)


temp_plot_data <- behaviour_change_time_cumulative %>%
  rownames_to_column('complication') %>%
  gather('year', 'count', -complication) %>%
  mutate(year = paste0('year', as.numeric(stringr::str_extract(year, '\\d+')))) %>%
  # mutate(year = as.numeric(stringr::str_extract(year, '\\d+'))) %>%
  mutate(complication = factor(complication, levels = c('none_complication', 'Stricturing', 'Stricturing;Penetrating', 'Penetrating')))


color <-c( "#c4c7c9","#7fb1d3", "#fcb463", "#fa8072")

# ggalluvial plot
library(ggplot2)
library(reshape2)
library(ggalluvial)
library(vegan)
library(ggbreak)

temp_plot <- ggplot(temp_plot_data, aes(x = year, y = count, fill = complication, stratum = complication, alluvium = complication)) +
  geom_stratum(width = 0.5,alpha=1) +
  geom_flow(alpha = 0.5) +
  geom_text(stat = "stratum", aes(label = count), size = 3, color = "black") +
  scale_fill_manual(values=color,
                    labels = c('none_complication' = 'B1',
                               'Stricturing' = 'B2',
                               'Stricturing;Penetrating' = 'B2+B3',
                               'Penetrating' = 'B3'))+
  scale_x_discrete(limits = c('year0', 'year1', 'year2', 'year3', 'year4'),
                   labels = c('Enrollment', '0-12 month', '12-24 month', '23-36 month', '>36 month')) +
  labs(x = '', y = 'Cummulative number', fill="CD compliation") +
  theme_bw() +
  theme(axis.text=element_text(colour='black',size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add break
temp_plot <- ggplot(temp_plot_data, aes(x = year, y = count, fill = complication, stratum = complication, alluvium = complication)) +
  geom_stratum(width = 0.5,alpha=1) +
  geom_flow(alpha = 0.5) +
  geom_text(stat = "stratum", aes(label = count), size = 3, color = "black") +
  scale_fill_manual(values=color,
                    labels = c('none_complication' = 'B1',
                               'Stricturing' = 'B2',
                               'Stricturing;Penetrating' = 'B2+B3',
                               'Penetrating' = 'B3'))+
  scale_x_discrete(limits = c('year0', 'year1', 'year2', 'year3', 'year4'),
                   labels = c('Enrollment', '0-12 month', '12-24 month', '23-36 month', '>36 month')) +
  labs(x = '', y = 'Cummulative number', fill="CD compliation") +
  theme_bw() +
  theme(axis.text=element_text(colour='black',size=9), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_break(c(100, 300), scales = 1.5, space = 0.15)

ggplot2::ggsave(temp_plot,
                filename = '~/Project/00_IBD_project/Figure/250326/Figure1/cumulative_complication_plot_2_260709.pdf',
                width = 9.5, height = 6)


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
readr::write_csv(temp_plot_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig_4_cumulative_complication_260721.csv')



################################################################################
# Differential analyses for CD compliation -------------------------------------
library(tidyverse)
library(tidymass)
library(sjmisc)

dir.create('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis', showWarnings = FALSE, recursive = TRUE)

load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')

# B2 vs B1 ---------------------------------------------------------------------

id_B1 <- object_stat %>% extract_sample_info() %>% filter(phenotype_group2 == 'B1') %>% pull(sample_id)
id_B2 <- object_stat %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B2')) %>% pull(sample_id)
id_B3 <- object_stat %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B3')) %>% pull(sample_id)

object_stat_B1_B2 <- object_stat %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B2'))

fc <- apply(object_stat_B1_B2@expression_data, 1, function(x){
  idx_control <- match(id_B1, names(x))
  idx_case <- match(id_B2, names(x))
  
  mean_control <- mean(x[idx_control])
  mean_case <- mean(x[idx_case])
  
  fc <- mean_case/mean_control
})


# calculate P-value
object_stat_temp <- object_stat %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B2')) %>% 
  scale_data(center = TRUE, method = 'auto')


# adjusted p-values with gender, age, race, use_antibiotics
temp_data_stat <- object_stat_temp@expression_data %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_stat_temp@sample_info, by = 'sample_id') %>% 
  mutate(phenotype_group3 = case_when(phenotype_group2 == 'B1' ~ 'B1',
                                      phenotype_group2 == 'B2' ~ 'B2'))

variable_name <- object_stat_temp@annotation_table$variable_id

fix_vs <- c('phenotype_group3', 'gender', 'age', 'race', 'use_antibiotics')

adjusted_result_B1_B2 <- pbapply::pblapply(variable_name, function(x){
  form <- as.formula(paste0(paste0(x, " ~ "), paste(fix_vs, collapse=' + ' )))
  model <- glm(form, data = temp_data_stat)
  
  stat_result <- summary(model)
  p_value <- stat_result$coefficients[, 4]
  coefficient <- stat_result$coefficients[,1]
  
  temp_result <- c(p_value, coefficient)
  names(temp_result) <- c(paste0('p_value_', names(p_value)), paste0('coefficient_', names(coefficient)))
  result <- temp_result
  
  return(result)
})

adjusted_result_B1_B2 <- adjusted_result_B1_B2 %>% 
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  dplyr::mutate(variable_name = variable_name) %>% 
  dplyr::select(variable_name, everything())

adjusted_result_B1_B2 <- adjusted_result_B1_B2 %>% 
  mutate(p_adjusted_intercept = p.adjust(`p_value_(Intercept)`, method = 'BH'),
         p_adjusted_phenotype_group3B2 = p.adjust(p_value_phenotype_group3B2, method = 'BH'),
         p_adjusted_genderMale = p.adjust(p_value_genderMale, method = 'BH'),
         p_adjusted_age = p.adjust(p_value_age, method = 'BH')) %>% 
  mutate(fc = fc)

temp_annot_table <- object_stat %>% 
  extract_annotation_table() %>% 
  select(variable_id, Compound.name:rt, id:adduct, confidence_level:metabolon_subclass)

adjusted_result_B1_B2 <- adjusted_result_B1_B2 %>% 
  bind_cols(temp_annot_table)

save(adjusted_result_B1_B2, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B2_250401.RData')



# B3 vs B1 ---------------------------------------------------------------------

object_stat_B1_B3 <- object_stat %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B3'))

object_stat_B1_B3@sample_info <- object_stat_B1_B3@sample_info %>% 
  mutate(phenotype_group3 = case_when(phenotype_group2 == 'B1' ~ 'B1',
                                      phenotype_group2 == 'B3' ~ 'B3'))

object_stat_B1_B3@sample_info_note <- data.frame(sample_id = colnames(object_stat_B1_B3@sample_info),
                                                 note = colnames(object_stat_B1_B3@sample_info), 
                                                 stringsAsFactors = FALSE)

fc <- apply(object_stat_B1_B3@expression_data, 1, function(x){
  idx_control <- match(id_B1, names(x))
  idx_case <- match(id_B3, names(x))
  
  mean_control <- mean(x[idx_control])
  mean_case <- mean(x[idx_case])
  
  fc <- mean_case/mean_control
})


# calculate P-value
object_stat_temp <- object_stat %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B3')) %>% 
  scale_data(center = TRUE, method = 'auto')


# adjusted p-values with gender, age, race, use_antibiotics
temp_data_stat <- object_stat_temp@expression_data %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_stat_temp@sample_info, by = 'sample_id') %>% 
  mutate(phenotype_group3 = case_when(phenotype_group2 == 'B1' ~ 'B1',
                                      phenotype_group2 == 'B3' ~ 'B3'))

variable_name <- object_stat_temp@annotation_table$variable_id

fix_vs <- c('phenotype_group3', 'gender', 'age', 'race', 'use_antibiotics')

adjusted_result_B1_B3 <- pbapply::pblapply(variable_name, function(x){
  form <- as.formula(paste0(paste0(x, " ~ "), paste(fix_vs, collapse=' + ' )))
  model <- glm(form, data = temp_data_stat)
  
  stat_result <- summary(model)
  p_value <- stat_result$coefficients[, 4]
  coefficient <- stat_result$coefficients[,1]
  
  temp_result <- c(p_value, coefficient)
  names(temp_result) <- c(paste0('p_value_', names(p_value)), paste0('coefficient_', names(coefficient)))
  result <- temp_result
  
  return(result)
})

adjusted_result_B1_B3 <- adjusted_result_B1_B3 %>% 
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  dplyr::mutate(variable_name = variable_name) %>% 
  dplyr::select(variable_name, everything())

adjusted_result_B1_B3 <- adjusted_result_B1_B3 %>% 
  mutate(p_adjusted_intercept = p.adjust(`p_value_(Intercept)`, method = 'BH'),
         p_adjusted_phenotype_group3B3 = p.adjust(p_value_phenotype_group3B3, method = 'BH'),
         p_adjusted_genderMale = p.adjust(p_value_genderMale, method = 'BH'),
         p_adjusted_age = p.adjust(p_value_age, method = 'BH')) %>% 
  mutate(fc = fc)

temp_annot_table <- object_stat %>% 
  extract_annotation_table() %>% 
  select(variable_id, Compound.name:rt, id:adduct, confidence_level:metabolon_subclass)

adjusted_result_B1_B3 <- adjusted_result_B1_B3 %>% 
  bind_cols(temp_annot_table)

save(adjusted_result_B1_B3, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B3_250401.RData')



library(ZZWtool)
temp_plot <- ZZWVolcanoPlot2(object = adjusted_result_B1_B3,
                             fc_column_name = 'fc',
                             p_value_column_name = 'p_adjusted_phenotype_group3B3',
                             fc_up_cutoff = 1,
                             fc_down_cutoff = 1,
                             up_color = 'tomato',
                             down_color = 'dodgerblue',
                             add_text = TRUE,
                             text_for = 'marker',
                             text_from = 'Compound.name') +
  theme(legend.position = c(0.1,0.8))



# Visualization ----------------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B2_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B3_250401.RData')

adjusted_result_B1_B2 <- adjusted_result_B1_B2 %>% 
  select(variable_id, p_adjusted_phenotype_group3B2, fc, Compound.name, id:metabolon_subclass) %>% 
  rename(p_adjusted = p_adjusted_phenotype_group3B2) %>% 
  mutate(p_adjusted_log10 = -log10(p_adjusted)) %>%
  mutate(fc2 = ifelse(fc > 1, fc, 1/fc),
         p_adjusted_log10_2 = ifelse(fc > 1, p_adjusted_log10, -p_adjusted_log10)) %>%
  mutate(phenotype_group = 'B1B2') %>% 
  select(variable_id:fc, fc2, p_adjusted_log10, p_adjusted_log10_2, phenotype_group, everything())

adjusted_result_B1_B3 <- adjusted_result_B1_B3 %>% 
  select(variable_id, p_adjusted_phenotype_group3B3, fc, Compound.name, id:metabolon_subclass) %>% 
  rename(p_adjusted = p_adjusted_phenotype_group3B3) %>% 
  mutate(p_adjusted_log10 = -log10(p_adjusted)) %>%
  mutate(fc2 = ifelse(fc > 1, fc, 1/fc),
         p_adjusted_log10_2 = ifelse(fc > 1, p_adjusted_log10, -p_adjusted_log10)) %>%
  mutate(phenotype_group = 'B1B3') %>% 
  select(variable_id:fc, fc2, p_adjusted_log10, p_adjusted_log10_2, phenotype_group, everything())

adjusted_result <- bind_rows(adjusted_result_B1_B2, adjusted_result_B1_B3)

fc_cutoff <- 1.5
p_adjusted_cutoff <- 0.05

temp_plot_data <- adjusted_result %>%
  mutate(
    comparison = recode(phenotype_group,
                        'B1B2' = 'B2 vs. B1',
                        'B1B3' = 'B3 vs. B1'),
    log2_fc = log2(fc),
    significant = p_adjusted < p_adjusted_cutoff &
      abs(log2_fc) > log2(fc_cutoff),
    plot_group = case_when(
      !significant ~ 'N.S.',
      metabolon_subclass == 'Drug' ~ 'Medications',
      metabolon_subclass == 'Phosphatidylethanolamines' ~ 'PE Lipids',
      metabolon_class == 'Lipids' ~ 'Non-PE Lipids',
      TRUE ~ 'Other metabolites'
    ),
    plot_group = factor(
      plot_group,
      levels = c('N.S.', 'Other metabolites', 'Non-PE Lipids',
                 'PE Lipids', 'Medications')
    )
  ) %>%
  arrange(plot_group != 'N.S.')

temp_plot_data <- temp_plot_data %>%
  select(variable_id, Compound.name, log2_fc, p_adjusted, p_adjusted_log10,
         significant, plot_group, comparison)

text_data <- temp_plot_data %>% 
  dplyr::filter(plot_group %in% c('PE Lipids', 'Medications'))

panel_annotation <- tibble(
  comparison = c('B2 vs. B1', 'B3 vs. B1'),
  higher_label = c('Higher in B2', 'Higher in B3'),
  fc_label = '|FC| > 1.5',
  p_label = 'italic(P)[adj.] < 0.05'
)

temp_plot <- ggplot(temp_plot_data, aes(x = log2_fc, y = p_adjusted_log10)) +
  geom_hline(yintercept = -log10(p_adjusted_cutoff),
             color = '#cfcfd1', linewidth = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
             color = '#cfcfd1', linewidth = 0.6, linetype = 'dashed') +
  geom_point(aes(fill = plot_group, color = plot_group), shape = 21,
             size = 3.2, stroke = 0.45, alpha = 0.95) +
  geom_text(
    data = panel_annotation,
    aes(x = 1.05, y = 7.75, label = higher_label),
    hjust = 0, size = 4.2, inherit.aes = FALSE
  ) +
  geom_text(
    data = panel_annotation,
    aes(x = 0, y = 8.25, label = fc_label),
    size = 4.2, inherit.aes = FALSE
  ) +
  geom_text(
    data = panel_annotation,
    aes(x = -5.75, y = 0.7, label = p_label),
    hjust = 0, size = 4.2, parse = TRUE, inherit.aes = FALSE
  ) +
  ggrepel::geom_text_repel(aes(log2_fc, p_adjusted_log10,
                               label = Compound.name),
                           data = text_data) +
  facet_wrap(~comparison, nrow = 1) +
  scale_fill_manual(
    values = c(
      'N.S.' = '#d1d2d5',
      'Other metabolites' = '#a7a7aa',
      'Non-PE Lipids' = '#f9c43a',
      'PE Lipids' = '#dc4943',
      'Medications' = '#5d82c7'
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      'N.S.' = '#d1d2d5',
      'Other metabolites' = '#a7a7aa',
      'Non-PE Lipids' = '#dfa91f',
      'PE Lipids' = '#bd3935',
      'Medications' = '#3f5da5'
    ),
    drop = FALSE,
    guide = 'none'
  ) +
  scale_x_continuous(
    breaks = seq(-6, 6, 2),
    limits = c(-6, 6),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 8, 2),
    limits = c(0, 8.4),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = expression(log[2](fold-change)),
    y = expression(-Log[10](italic(P)[adjusted])),
    fill = NULL
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic(base_size = 17) +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.16, 'cm'),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15, color = 'black'),
    strip.background = element_blank(),
    strip.text = element_text(size = 24, face = 'bold', margin = margin(b = 26)),
    panel.spacing.x = unit(3.0, 'cm'),
    legend.position = 'top',
    legend.justification = 'center',
    legend.text = element_text(size = 16),
    legend.key.width = unit(0.65, 'cm'),
    plot.margin = margin(12, 34, 10, 10)
  )

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure5/complication_disease_result_250401.pdf', 
       width = 14, 
       height = 6.8)


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
readr::write_csv(temp_plot_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig4b_complication_disease_result_260721.csv')



################################################################################
# Lipid class enrichment analysis ----------------------------------------------
# Part 1: Enrichment analysis of chemical classes ------------------------------
# Load the enrichment data

load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B2_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_adjusted_result_B1_B3_250401.RData')


adjusted_result_B1_B2 <- adjusted_result_B1_B2 %>% 
  select(variable_id, p_value_phenotype_group3B2, p_adjusted_phenotype_group3B2, fc, Compound.name, id:metabolon_subclass) %>% 
  rename(p = p_value_phenotype_group3B2, 
         p_adjusted = p_adjusted_phenotype_group3B2) %>% 
  mutate(p_adjusted_log10 = -log10(p_adjusted)) %>%
  mutate(fc2 = ifelse(fc > 1, fc, 1/fc),
         p_adjusted_log10_2 = ifelse(fc > 1, p_adjusted_log10, -p_adjusted_log10)) %>%
  mutate(phenotype_group = 'B1B2') %>% 
  select(variable_id:fc, p_adjusted_log10, phenotype_group, everything())

adjusted_result_B1_B3 <- adjusted_result_B1_B3 %>% 
  select(variable_id, p_value_phenotype_group3B3, p_adjusted_phenotype_group3B3, fc, Compound.name, id:metabolon_subclass) %>% 
  rename(p = p_value_phenotype_group3B3,
         p_adjusted = p_adjusted_phenotype_group3B3) %>% 
  mutate(p_adjusted_log10 = -log10(p_adjusted)) %>%
  mutate(fc2 = ifelse(fc > 1, fc, 1/fc),
         p_adjusted_log10_2 = ifelse(fc > 1, p_adjusted_log10, -p_adjusted_log10)) %>%
  mutate(phenotype_group = 'B1B3') %>% 
  select(variable_id:fc, p_adjusted_log10, phenotype_group, everything())

adjusted_result <- bind_rows(adjusted_result_B1_B2, adjusted_result_B1_B3)

# save(adjusted_result, 
#      file = '~/Project/00_IBD_project/Data/20260422_PE_lipids_update/adjusted_result_combined_20260422.RData')


# Part 2: Enrichment on metabolon_subclass -------------------------------------

# enrichment using fgsea 

fgsea_subclass <- function(df,
                           class_cutoff = 2,
                           use_raw_p = FALSE,
                           ignore_direction = FALSE) {
  if (use_raw_p && !"p" %in% colnames(df)) {
    stop("Column 'p' is required when use_raw_p = TRUE.")
  }
  
  p_col <- if (use_raw_p) "p" else "p_adjusted"
  
  dfx <- df %>%
    dplyr::filter(
      !is.na(variable_id),
      !is.na(metabolon_subclass), metabolon_subclass != "",
      !is.na(fc), fc > 0,
      !is.na(.data[[p_col]])
    ) %>%
    dplyr::mutate(
      p_for_stat = .data[[p_col]],
      direction_weight = if (ignore_direction) 1 else sign(log2(fc)),
      stat = direction_weight * (-log10(pmax(p_for_stat, 1e-300)))
    )
  
  dfx <- dfx %>% dplyr::filter(!is.na(p_for_stat), !is.na(stat))
  
  # Ensure one statistic per metabolite id.
  rank_df <- dfx %>%
    dplyr::group_by(variable_id) %>%
    dplyr::summarise(stat = mean(stat, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(stat))
  
  stats <- rank_df$stat
  names(stats) <- rank_df$variable_id
  
  pathways <- dfx %>%
    dplyr::group_by(metabolon_subclass) %>%
    dplyr::summarise(members = list(unique(variable_id)), .groups = "drop") %>%
    dplyr::filter(lengths(members) >= class_cutoff)
  
  pathways_list <- pathways$members
  names(pathways_list) <- pathways$metabolon_subclass
  
  fgsea::fgseaMultilevel(
    pathways = pathways_list,
    stats = stats,
    minSize = class_cutoff,
    eps = 1e-10
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      stat_p_source = dplyr::if_else(use_raw_p, "raw_p", "adjusted_p"),
      direction_mode = dplyr::if_else(ignore_direction, "absolute", "signed")
    ) %>%
    dplyr::rename(metabolon_subclass = pathway) %>%
    dplyr::arrange(padj, pval)
}

fgsea_subclass_result <- adjusted_result %>%
  dplyr::group_by(phenotype_group) %>%
  dplyr::group_modify(~ fgsea_subclass(.x, 
                                       class_cutoff = 10,
                                       use_raw_p = TRUE,
                                       ignore_direction = TRUE)) %>%
  dplyr::ungroup()

save(fgsea_subclass_result, 
     file = '~/Project/00_IBD_project/Data/20260422_PE_lipids_update/fgsea_subclass_enrichment_result_20260422.RData')


# fgsea_subclass_sig <- fgsea_subclass_result %>%
#   dplyr::arrange(desc(NES))

rm(list = ls());gc()


# Part 3: Visualization of subclass enrichment to introduce the PE -------------
load('~/Project/00_IBD_project/Data/20260422_PE_lipids_update/adjusted_result_combined_20260422.RData')
load('~/Project/00_IBD_project/Data/20260422_PE_lipids_update/fgsea_subclass_enrichment_result_20260422.RData')

library(cowplot)

# Top 10 classes ranked by max |NES| across both phenotype groups
top10_class_fgsea <- fgsea_subclass_result %>%
  group_by(metabolon_subclass) %>%
  summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = 'drop') %>%
  slice_max(order_by = max_abs_NES, n = 10) %>%
  arrange(desc(max_abs_NES)) %>%       # ascending → bottom-to-top after coord_flip
  pull(metabolon_subclass)

class_levels <- rev(top10_class_fgsea)

# Shared add-on theme: suppress y-axis for non-label panels
no_y_theme <- theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# Panel 1: NES lollipop — one dot per subclass (max |NES| across groups) 
fgsea_top10_max <- fgsea_subclass_result %>%
  filter(metabolon_subclass %in% top10_class_fgsea) %>%
  mutate(metabolon_subclass = factor(metabolon_subclass, levels = class_levels)) %>%
  group_by(metabolon_subclass) %>%
  slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    sig_alpha = ifelse(padj < 0.05, 1, 0.4)
  )

p_nes <- ggplot(fgsea_top10_max,
                aes(x = 1, y = metabolon_subclass,
                    size = abs(NES))) +
  geom_point(color = "#8cd4c8") +
  scale_size_continuous(name = 'Normalized\nenrichment score (NES)', range = c(2, 10)) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  xlab('NES') + ylab('') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(
    axis.text.y  = element_text(hjust = 1, angle = 0),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    panel.border = element_rect(color = 'grey80', fill = NA, linewidth = 0.5),
    plot.margin  = margin(5, 8, 5, 12)
  )

# Panel 2: Total metabolite count (lollipop)
count_top10 <- adjusted_result %>%
  filter(metabolon_subclass %in% top10_class_fgsea, phenotype_group == 'B1B2') %>%
  group_by(metabolon_subclass) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(metabolon_subclass = factor(metabolon_subclass, levels = class_levels))

p_count <- ggplot(count_top10, aes(x = metabolon_subclass, y = n)) +
  geom_segment(aes(xend = metabolon_subclass, y = 0, yend = n),
               color = "#7fb1d3", linewidth = 0.8) +
  geom_point(size = 3, color = "#7fb1d3") +
  coord_flip() +
  xlab('') + ylab('Number of detected metabolites') +
  ZZWtool::ZZWTheme(type = 'classic') +
  no_y_theme

# P-value heatmap (metabolites sorted by P within each subclass)
# fgsea_df: one row per subclass for a single phenotype group
make_pval_heatmap <- function(fgsea_df, title_label, show_legend = TRUE) {
  p_breaks <- c('P < 0.05', 'P > 0.05')
  p_colors <- c('#e8cf79', '#d0d0d0')
  
  plot_df <- tibble::tibble(metabolon_subclass = class_levels) %>%
    left_join(
      fgsea_df %>%
        filter(metabolon_subclass %in% top10_class_fgsea) %>%
        select(metabolon_subclass, pval),
      by = 'metabolon_subclass'
    ) %>%
    mutate(
      metabolon_subclass = factor(metabolon_subclass, levels = class_levels),
      p_cat = dplyr::case_when(
        !is.na(pval) & pval < 0.05 ~ 'P < 0.05',
        TRUE                       ~ 'P > 0.05'
      ),
      p_cat = factor(p_cat, levels = p_breaks)
    )
  
  plot_df %>%
    ggplot(aes(x = 1, y = metabolon_subclass, fill = p_cat)) +
    geom_tile(width = 0.9, height = 0.95, color = '#2b2b2b', linewidth = 0.3) +
    scale_fill_manual(values = setNames(p_colors, p_breaks), name = 'Class enrichment',
                      drop = FALSE) +
    scale_x_continuous(limits = c(0.55, 1.45), expand = c(0, 0)) +
    ggtitle(title_label) +
    xlab('P-value') +
    ZZWtool::ZZWTheme(type = 'classic') +
    no_y_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      axis.line.y  = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8),
      plot.title   = element_text(size = 9, hjust = 0.5),
      plot.margin  = margin(5, 2, 5, 2),
      legend.position = if (show_legend) 'right' else 'none'
    )
}

# Stacked bar (proportion significant by FDR-adjusted P) 
make_stacked_bar <- function(df, p_cutoff = 0.05, show_legend = FALSE) {
  df %>%
    filter(metabolon_subclass %in% top10_class_fgsea) %>%
    mutate(
      metabolon_subclass = factor(metabolon_subclass, levels = class_levels),
      sig = ifelse(p_adjusted < p_cutoff, 'Significant', 'Insignificant')
    ) %>%
    count(metabolon_subclass, sig) %>%
    ggplot(aes(x = metabolon_subclass, y = n, fill = sig)) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(
      values = c('Significant' = '#fb7f72', 'Insignificant' = '#c5c7c9'),
      labels = c('Significant' = 'P < 0.05', 'Insignificant' = 'P > 0.05'),
      name = 'Changed metabolites'
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_flip() +
    xlab('') + ylab('Percentage (%)') +
    ZZWtool::ZZWTheme(type = 'classic') +
    no_y_theme +
    theme(
      legend.position = if (show_legend) 'right' else 'none',
      axis.line.x     = element_blank(),
      axis.line.y     = element_blank(),
      panel.border    = element_rect(color = 'grey80', fill = NA, linewidth = 0.5)
    )
}

# Build per-group panels 
df_B1B2 <- adjusted_result %>% filter(phenotype_group == 'B1B2') %>% select(-c('fc2', 'p_adjusted_log10_2'))
df_B1B3 <- adjusted_result %>% filter(phenotype_group == 'B1B3') %>% select(-c('fc2', 'p_adjusted_log10_2'))

p_heatmap_B1B2 <- make_pval_heatmap(fgsea_subclass_result %>% filter(phenotype_group == 'B1B2'), 'B2 vs B1', show_legend = FALSE)
p_stacked_B1B2 <- make_stacked_bar(df_B1B2, show_legend = FALSE)
p_heatmap_B1B3 <- make_pval_heatmap(fgsea_subclass_result %>% filter(phenotype_group == 'B1B3'), 'B3 vs B1', show_legend = FALSE)
p_stacked_B1B3 <- make_stacked_bar(df_B1B3, show_legend = FALSE)

# Combine: [NES lollipop | Count | B1B2 heatmap | B1B2 bar | B1B3 heatmap | B1B3 bar] 
fgsea_enrich_plot <- plot_grid(
  p_nes, p_count,
  p_heatmap_B1B2, p_stacked_B1B2,
  p_heatmap_B1B3, p_stacked_B1B3,
  ncol = 6,
  align = 'h',
  axis = 'tb',
  rel_widths = c(1.05, 0.45, 0.2, 0.45, 0.2, 0.45)
)

ggsave(fgsea_enrich_plot,
       filename = '~/Project/00_IBD_project/Data/20260422_PE_lipids_update/fgsea_subclass_enrichment_3_20260422.pdf',
       width = 24, height = 6.5)


# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# 
# temp_data1 <- p_nes$data %>% arrange(desc(abs(NES))) 
# temp_data2 <- p_heatmap_B1B2$data %>% rename(pval_B1B2 = pval, p_cat_B1B2 = p_cat)
# temp_data3 <- p_heatmap_B1B3$data %>% rename(pval_B1B3 = pval, p_cat_B1B3 = p_cat)
# temp_data4 <- p_stacked_B1B2$data %>% 
#   rename(n_B1B2 = n, sig_B1B2 = sig) %>% 
#   pivot_wider(names_from = sig_B1B2, values_from = n_B1B2) %>% 
#   tidyr::replace_na(list(Significant = 0, Insignificant = 0)) %>% 
#   rename(n_B1B2_Significant = Significant, n_B1B2_Insignificant = Insignificant)
# temp_data5 <- p_stacked_B1B3$data %>% 
#   rename(n_B1B3 = n, sig_B1B3 = sig) %>% 
#   pivot_wider(names_from = sig_B1B3, values_from = n_B1B3) %>% 
#   tidyr::replace_na(list(Significant = 0, Insignificant = 0)) %>% 
#   rename(n_B1B3_Significant = Significant, n_B1B3_Insignificant = Insignificant)
# 
# temp_data <- temp_data1 %>% 
#   select(-leadingEdge) %>% 
#   left_join(temp_data2, by = 'metabolon_subclass') %>% 
#   left_join(temp_data3, by = 'metabolon_subclass') %>% 
#   left_join(temp_data4, by = 'metabolon_subclass') %>% 
#   left_join(temp_data5, by = 'metabolon_subclass')
# 
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# readr::write_csv(temp_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig8_fgsea_subclass_enrichment_260721.csv')
# 
# rm(list = ls());gc()



################################################################################
# Differential analyses for lipids in CD complication --------------------------
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_250401.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_demographic_241003.RData')

object_lipid_enrollment <- object_lipid %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(visit_number == 'enrollment') %>% 
  filter(sample_id %in% patient_meta_info_demographic$sample_id)

object_lipid_enrollment %>% 
  extract_annotation_table() %>% 
  count(lipid_class)


object_lipid_enrollment@sample_info <- object_lipid_enrollment@sample_info %>% 
  select(sample_id) %>% 
  left_join(patient_meta_info_demographic, by = 'sample_id')

object_lipid_enrollment@sample_info_note <- data.frame(sample_id = colnames(object_lipid_enrollment@sample_info),
                                                       note = colnames(object_lipid_enrollment@sample_info), 
                                                       stringsAsFactors = FALSE)

save(object_lipid_enrollment, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_enrollment_250401.RData')

rm(list = ls());gc()

# B2 vs B1 
# Exclude B2+B3
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_enrollment_250401.RData')

id_B1 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 == 'B1') %>% pull(sample_id)
id_B2 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B2')) %>% pull(sample_id)
id_B3 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B3')) %>% pull(sample_id)

object_lipid_enrollment_B1_B2 <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B2'))

fc <- apply(object_lipid_enrollment_B1_B2@expression_data, 1, function(x){
  idx_control <- match(id_B1, names(x))
  idx_case <- match(id_B2, names(x))
  
  mean_control <- mean(x[idx_control])
  mean_case <- mean(x[idx_case])
  
  fc <- mean_case/mean_control
})


# calculate P-value
object_lipid_enrollment_temp <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B2')) %>% 
  scale_data(center = TRUE, method = 'auto')


# adjusted p-values with gender, age, race, use_antibiotics
temp_data_stat <- object_lipid_enrollment_temp@expression_data %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_lipid_enrollment_temp@sample_info, by = 'sample_id') %>% 
  mutate(phenotype_group3 = case_when(phenotype_group2 == 'B1' ~ 'B1',
                                      phenotype_group2 == 'B2' ~ 'B2'))

variable_name <- object_lipid_enrollment_temp@annotation_table$variable_id

fix_vs <- c('phenotype_group3', 'gender', 'age', 'race', 'use_antibiotics')

lipid_adjusted_result_B1_B2 <- pbapply::pblapply(variable_name, function(x){
  form <- as.formula(paste0(paste0(x, " ~ "), paste(fix_vs, collapse=' + ' )))
  model <- glm(form, data = temp_data_stat)
  stat_result <- summary(model)
  
  p_value <- stat_result$coefficients[, 4]
  coefficient <- stat_result$coefficients[,1]
  
  temp_result <- c(p_value, coefficient)
  names(temp_result) <- c(paste0('p_value_', names(p_value)), paste0('coefficient_', names(coefficient)))
  result <- temp_result
  
  return(result)
})

lipid_adjusted_result_B1_B2 <- lipid_adjusted_result_B1_B2 %>% 
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  dplyr::mutate(variable_name = variable_name) %>% 
  dplyr::select(variable_name, everything())

lipid_adjusted_result_B1_B2 <- lipid_adjusted_result_B1_B2 %>% 
  mutate(p_adjusted_intercept = p.adjust(`p_value_(Intercept)`, method = 'BH'),
         p_adjusted_phenotype_group3B2 = p.adjust(p_value_phenotype_group3B2, method = 'BH'),
         p_adjusted_genderMale = p.adjust(p_value_genderMale, method = 'BH'),
         p_adjusted_age = p.adjust(p_value_age, method = 'BH')) %>% 
  mutate(fc = fc)

temp_annot_table <- object_lipid_enrollment %>% 
  extract_annotation_table()

lipid_adjusted_result_B1_B2 <- lipid_adjusted_result_B1_B2 %>% 
  bind_cols(temp_annot_table)

save(lipid_adjusted_result_B1_B2, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_lipid_adjusted_result_B1_B2_250401.RData')


lipid_adjusted_result_B1_B2 %>% 
  select(variable_name, p_value_phenotype_group3B2, p_adjusted_phenotype_group3B2, fc, everything()) %>% 
  filter(p_adjusted_phenotype_group3B2 <= 0.2) %>% 
  arrange(p_value_phenotype_group3B2) %>% 
  count(lipid_class)

rm(list = ls());gc()


# B3 vs B1 
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_enrollment_250401.RData')

id_B1 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 == 'B1') %>% pull(sample_id)
id_B2 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B2')) %>% pull(sample_id)
id_B3 <- object_lipid_enrollment %>% extract_sample_info() %>% filter(phenotype_group2 %in% c('B3')) %>% pull(sample_id)

object_lipid_enrollment_B1_B3 <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B3'))

fc <- apply(object_lipid_enrollment_B1_B3@expression_data, 1, function(x){
  idx_control <- match(id_B1, names(x))
  idx_case <- match(id_B3, names(x))
  
  mean_control <- mean(x[idx_control])
  mean_case <- mean(x[idx_case])
  
  fc <- mean_case/mean_control
})


# calculate P-value
object_lipid_enrollment_temp <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B3')) %>% 
  scale_data(center = TRUE, method = 'auto')


# adjusted p-values with gender, age, race, use_antibiotics
temp_data_stat <- object_lipid_enrollment_temp@expression_data %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_lipid_enrollment_temp@sample_info, by = 'sample_id') %>% 
  mutate(phenotype_group3 = case_when(phenotype_group2 == 'B1' ~ 'B1',
                                      phenotype_group2 == 'B3' ~ 'B3'))

variable_name <- object_lipid_enrollment_temp@annotation_table$variable_id

fix_vs <- c('phenotype_group3', 'gender', 'age', 'race', 'use_antibiotics')

lipid_adjusted_result_B1_B3 <- pbapply::pblapply(variable_name, function(x){
  form <- as.formula(paste0(paste0(x, " ~ "), paste(fix_vs, collapse=' + ' )))
  model <- glm(form, data = temp_data_stat)
  stat_result <- summary(model)
  
  p_value <- stat_result$coefficients[, 4]
  coefficient <- stat_result$coefficients[,1]
  
  temp_result <- c(p_value, coefficient)
  names(temp_result) <- c(paste0('p_value_', names(p_value)), paste0('coefficient_', names(coefficient)))
  result <- temp_result
  
  return(result)
})

lipid_adjusted_result_B1_B3 <- lipid_adjusted_result_B1_B3 %>% 
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  dplyr::mutate(variable_name = variable_name) %>% 
  dplyr::select(variable_name, everything())

lipid_adjusted_result_B1_B3 <- lipid_adjusted_result_B1_B3 %>% 
  mutate(p_adjusted_intercept = p.adjust(`p_value_(Intercept)`, method = 'BH'),
         p_adjusted_phenotype_group3B3 = p.adjust(p_value_phenotype_group3B3, method = 'BH'),
         p_adjusted_genderMale = p.adjust(p_value_genderMale, method = 'BH'),
         p_adjusted_age = p.adjust(p_value_age, method = 'BH')) %>% 
  mutate(fc = fc)

temp_annot_table <- object_lipid_enrollment %>% 
  extract_annotation_table()

lipid_adjusted_result_B1_B3 <- lipid_adjusted_result_B1_B3 %>% 
  bind_cols(temp_annot_table)

save(lipid_adjusted_result_B1_B3, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_lipid_adjusted_result_B1_B3_250401.RData')

lipid_adjusted_result_B1_B3 %>% 
  select(variable_name, p_value_phenotype_group3B3, p_adjusted_phenotype_group3B3, fc, everything()) %>% 
  filter(p_adjusted_phenotype_group3B3 <= 0.2) %>% 
  arrange(p_value_phenotype_group3B3) %>% 
  count(lipid_class)


rm(list = ls());gc()


# merge results to one table
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_dereplication_table_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_lipid_adjusted_result_B1_B2_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/linear_regression_lipid_adjusted_result_B1_B3_250401.RData')

temp_result_B2 <- lipid_adjusted_result_B1_B2 %>% 
  select(variable_name, p_value_phenotype_group3B2, p_adjusted_phenotype_group3B2, fc) %>% 
  rename(variable_id = variable_name,
         p_value_B2 = p_value_phenotype_group3B2,
         p_value_adj_B2 = p_adjusted_phenotype_group3B2,
         fc_B2 = fc)

temp_result_B3 <- lipid_adjusted_result_B1_B3 %>% 
  select(variable_name, p_value_phenotype_group3B3, p_adjusted_phenotype_group3B3, fc) %>% 
  rename(variable_id = variable_name,
         p_value_B3 = p_value_phenotype_group3B3,
         p_value_adj_B3 = p_adjusted_phenotype_group3B3,
         fc_B3 = fc)

lipid_annot_result <- lipid_annot_table %>% 
  left_join(temp_result_B2, by = 'variable_id') %>%
  left_join(temp_result_B3, by = 'variable_id')

save(lipid_annot_result, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_annot_result_merge_250401.RData')


# lipid analysis visualization -------------------------------------------------
lipid_class_enrichment_table <- lipid_annot_result %>% 
  group_by(lipid_class) %>% 
  summarise(n = n(),
            n_B3_increase = sum(p_value_adj_B3 <= 0.2 & fc_B3 >= 6/5),
            n_B3_decrease = sum(p_value_adj_B3 <= 0.2 & fc_B3 <= 5/6),
            n_B3_no_change = sum(p_value_adj_B3 > 0.2)) %>% 
  mutate(percentage_B3_increase = n_B3_increase/n,
         percentage_B3_decrease = n_B3_decrease/n,
         percentage_B3_no_change = n_B3_no_change/n) %>% arrange(desc(n)) %>% 
  mutate(category = case_when(lipid_class %in% c('FA', 'OxFA', 'Car', "FAHFA") ~ 'Fatty_Acyls',
                              lipid_class %in% c('MG', 'DG', 'TG') ~ 'Glycerolipids',
                              lipid_class %in% c("Cer", "HexCer", "SM", "Sph", 'PhytoSph', 'SphP') ~ 'Sphingolipids',
                              lipid_class %in% c('LPC', 'EtherLPC', 'LPE', 'LPG', 'LPS', 'PC', 'EtherPC', 'PE', 'EtherPE', 'PG', 'PI', 'PS', 'LNAPE', 'OxPE') ~ 'Glycerophospholipids',
                              lipid_class %in% c('Vitamin D', 'Vitamin E') ~ 'Sterol Lipids'))



temp_data1 <- lipid_class_enrichment_table %>%
  filter(n > 2) %>%
  dplyr::arrange(desc(percentage_B3_increase)) %>% 
  mutate(y = match(lipid_class, rev(unique(lipid_class))))

temp_plot1 <- ggplot(temp_data1) +
  geom_bar(aes(x = y, y = n, fill = category), stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c('Fatty_Acyls' = '#b5de68', 
                               'Glycerolipids' = '#fdb562', 
                               'Sphingolipids' = '#ffed6f', 
                               'Glycerophospholipids' = '#8cd4c8', 
                               'Sterol Lipids' = '#beb9da')) +
  scale_x_continuous(breaks = seq(16, 1, -1), labels = unique(temp_data1$lipid_class)) +
  geom_text(aes(x = y, y = n, label = n), vjust = 0.5, hjust = -0.5, size = 3) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab('Lipid class') +
  coord_flip() +
  theme(rect = element_blank(), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, hjust = 0.5), 
        # axis.ticks.x = element_blank(), 
        legend.position = c(0.8, 0.2))

temp_data2 <- lipid_class_enrichment_table %>%
  filter(n > 2) %>% 
  dplyr::arrange(desc(percentage_B3_increase)) %>%
  select(lipid_class, percentage_B3_increase, percentage_B3_decrease, percentage_B3_no_change) %>%
  pivot_longer(cols = c(percentage_B3_increase, percentage_B3_decrease, percentage_B3_no_change), names_to = 'pcg', values_to = 'percentage') %>% 
  mutate(x = case_when(pcg == 'percentage_B3_increase' ~ 1,
                       pcg == 'percentage_B3_no_change' ~ 2,
                       pcg == 'percentage_B3_decrease' ~ 3),
         y = match(lipid_class, rev(unique(lipid_class))))

# library(ggnewscale)
# temp_plot2 <- ggplot(mapping = aes(x = x, y = y, fill = percentage)) +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_increase'),
#             color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#fc8070'), name = 'Increase') +
#   new_scale_fill() +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_no_change'),
#             aes(fill = percentage), color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#b5b5b6'), name = 'No change') +
#   new_scale_fill() +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_decrease'),
#             aes(fill = percentage), color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#81b2d3'), name = 'Decrease') +
#   scale_x_continuous(breaks = c(1, 2, 3), labels = c('Increase', 'No change', 'Decrease')) +
#   scale_y_continuous(breaks = seq(16, 1, -1), labels = unique(temp_data2$lipid_class)) +
#   coord_cartesian(xlim = c(0.5, 3.5)) +
#   theme(rect = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 10),
#         axis.ticks.x = element_blank(),
#         legend.position = 'none')



# tile plot with numerical values
library(ggnewscale)
temp_plot2 <- ggplot(mapping = aes(x = x, y = y, fill = percentage)) +
  geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_increase'),
            color = 'black', size = 0.4) +
  geom_text(data = filter(temp_data2, pcg == 'percentage_B3_increase'),
            aes(label = sprintf('%.0f%%', percentage * 100)), size = 3) +
  scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#fc8070'), name = 'Increase') +
  new_scale_fill() +
  geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_no_change'),
            aes(fill = percentage), color = 'black', size = 0.4) +
  geom_text(data = filter(temp_data2, pcg == 'percentage_B3_no_change'),
            aes(label = sprintf('%.0f%%', percentage * 100)), size = 3) +
  scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#b5b5b6'), name = 'No change') +
  new_scale_fill() +
  geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_decrease'),
            aes(fill = percentage), color = 'black', size = 0.4) +
  geom_text(data = filter(temp_data2, pcg == 'percentage_B3_decrease'),
            aes(label = sprintf('%.0f%%', percentage * 100)), size = 3) +
  scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#81b2d3'), name = 'Decrease') +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c('Increase', 'No change', 'Decrease')) +
  scale_y_continuous(breaks = seq(16, 1, -1), labels = unique(temp_data2$lipid_class)) +
  coord_cartesian(xlim = c(0.5, 3.5)) +
  theme(rect = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 10),
        axis.ticks.x = element_blank(),
        legend.position = 'none')


library(cowplot)
library(ggtree)

pp <- list(temp_plot1, temp_plot2)
temp_plot_merge <- plot_grid(plotlist=pp, ncol=2, align ='h', rel_widths = c(0.7, 0.3))

ggsave(temp_plot_merge, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure4/lipid_class_heatmap_260708.pdf', 
       width = 6, height = 10)



# 
# library(ggnewscale)
# temp_plot3 <- ggplot(mapping = aes(x = x, y = y, fill = percentage)) +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_increase'),
#             color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#fc8070'), name = 'Increase') +
#   new_scale_fill() +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_no_change'),
#             aes(fill = percentage), color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#b5b5b6'), name = 'No change') +
#   new_scale_fill() +
#   geom_tile(data = filter(temp_data2, pcg == 'percentage_B3_decrease'),
#             aes(fill = percentage), color = 'black', size = 0.4) +
#   scale_fill_gradientn(limits = c(0, 1), colors = c('white', '#81b2d3'), name = 'Decrease') +
#   scale_x_continuous(breaks = c(1, 2, 3), labels = c('Increase', 'No change', 'Decrease')) +
#   scale_y_continuous(breaks = seq(16, 1, -1), labels = unique(temp_data2$lipid_class)) +
#   coord_cartesian(xlim = c(0.5, 3.5)) +
#   theme(rect = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 10),
#         axis.ticks.x = element_blank(),
#         legend.position = 'top')
# 
# 
# 
# ggsave(temp_plot3, 
#        filename = '~/Project/00_IBD_project/Figure/250326/Figure4/lipid_class_heatmap_legend_260708.pdf', 
#        width = 6, height = 10)

rm(list = ls());gc()


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
readr::write_csv(lipid_class_enrichment_table, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/figure4_lipid_class_enrichment_260721.csv')



################################################################################
# Lipid molecular networks -----------------------------------------------------
# convert_spectra_data -------------------------------------------------------

#' @title convert_spectra_data
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

convert_spectra_data2 <- function(ms2_data) {
  options(readr.num_columns = 0)
  ms2_data$info <- ms2_data$info %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = 'variable') %>% 
    tidyr::pivot_wider(names_from = 'variable', values_from = 'V1') %>% 
    dplyr::mutate(PRECURSORMZ = as.numeric(PRECURSORMZ),
                  PRECURSORRT = as.numeric(PRECURSORRT))
  
  temp_info <- ms2_data$info %>%
    dplyr::rename(name = NAME,
                  mz = PRECURSORMZ) %>%
    dplyr::select(name:mz) %>%
    readr::type_convert()
  
  temp_ms2_data <- ms2_data$spec
  
  result <- new('SpectraData',
                info = temp_info,
                spectra = list(temp_ms2_data))
  
  return(result)
}



# tidy_network ---------------------------------------------------------------
tidy_network <- function(object_graph, score_cutoff = 0.5, topK = 10) {
  # browser()
  object_graph <- object_graph %>% 
    tidygraph::activate(edges) %>%
    tidygraph::filter(score >= score_cutoff) %>%
    tidygraph::activate(nodes) %>% 
    tidygraph::mutate(index_id = seq(n())) %>% 
    tidygraph::mutate(degree = centrality_degree())
  
  # if the edge table is empty, directly return the score_cutoff filtered graph
  temp_edge <- object_graph %>% activate(edges) %>% as_tibble()
  if (nrow(temp_edge) == 0) {
    return(object_graph)
  }
  
  # extract subgraphs for each node
  #   1. extract according to node id
  #   2. order by score
  #   3. calculate the rank according to the rank: used in top 10
  
  # extract the subgraph list for each node
  temp_nodes <- object_graph %>% as_tibble()
  
  list_pc_subnetworks <- lapply(temp_nodes$index_id, function(x){
    cat(x, ' ')
    temp_graph <- object_graph %>% 
      activate(edges) %>% 
      filter(from == x | to == x)
    
    temp <- temp_graph %>% activate(edges) %>% as_tibble()
    if (nrow(temp) == 0) {
      return(NULL)
    }
    
    temp <- temp_graph %>% 
      activate(edges) %>% 
      mutate(score_2 = n_frag_total*score) %>% 
      arrange(desc(score_2)) %>% 
      mutate(rank = seq(n()))
    return(temp)
  })
  
  names(list_pc_subnetworks) <- temp_nodes$index_id
  
  # add the rank to the edge table
  temp <- object_graph %>% 
    activate(edges) %>% 
    as_tibble()
  
  rank_result <- mapply(function(id_from, id_to){
    rank1 <- list_pc_subnetworks[[id_from]] %>% 
      activate(edges) %>% 
      as_tibble() %>% 
      filter(from == id_from & to == id_to) %>% 
      pull(rank)
    
    rank2 <- list_pc_subnetworks[[id_to]] %>% 
      activate(edges) %>% 
      as_tibble() %>% 
      filter(from == id_from & to == id_to) %>% 
      pull(rank)
    
    result <- tibble(rank_from = rank1, rank_to = rank2)
    
    return(result)
  },
  id_from = temp$from,
  id_to = temp$to, 
  SIMPLIFY = FALSE) %>% 
    dplyr::bind_rows()
  
  
  # merge to the graph
  object_graph <- object_graph %>% 
    activate(edges) %>% 
    mutate(rank_from = rank_result$rank_from, rank_to = rank_result$rank_to)
  
  result_object_graph <- object_graph %>% 
    activate(edges) %>% 
    filter(rank_from <= topK & rank_to <= topK)
  
  return(result_object_graph)
}

library(DoddLabMetID)
library(tidygraph)

load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_annot_result_merge_250401.RData')
load('~/Project/00_IBD_project/Data/20241120_CD_complicates/ms2_lipids_241125.RData')

unique_lipid_class <- lipid_annot_result %>% 
  count(lipid_class) %>% 
  arrange(desc(n)) %>% 
  pull(lipid_class)


list_lipid_networks <- pbapply::pblapply(seq_along(unique_lipid_class), function(i){
  cat('Class ', i, '\n')
  temp_lipid_class <- unique_lipid_class[i]
  feature_lipids <- lipid_annot_result %>% 
    filter(lipid_class == temp_lipid_class) %>% 
    pull(variable_id)
  
  if (length(feature_lipids) == 1) {
    return(NULL)
  }
  
  # lipid_class molecular network 
  lipid_class_pairs <- combn(feature_lipids, m = 2) %>% 
    t() %>% 
    as_tibble() %>% 
    rename('Spec1' = 'V1',
           'Spec2' = 'V2')
  
  lipid_class_networks <- lapply(seq_along(lipid_class_pairs$Spec1), function(i){
    temp_feature1 <- lipid_class_pairs$Spec1[i]
    temp_feature2 <- lipid_class_pairs$Spec2[i]
    temp_spec1 <- convert_spectra_data2(ms2_data = ms2_lipids[[temp_feature1]])
    temp_spec2 <- convert_spectra_data2(ms2_data = ms2_lipids[[temp_feature2]])
    result_ms2 <- runSpecMatch(obj_ms2_cpd1 = temp_spec1, obj_ms2_cpd2 = temp_spec2, mz_tol_ms2 = 35, scoring_approach = 'gnps')
    result <- result_ms2@info %>% dplyr::select(score:n_frag_total)
    rownames(result) <- NULL
    return(result)
  }) %>% bind_rows()
  
  lipid_class_networks <- lipid_class_pairs %>% 
    bind_cols(lipid_class_networks)
  
  # filter the networks
  edge_table <- lipid_class_networks %>% 
    # filter(score >= 0.5) %>% 
    rename('from' = 'Spec1',
           'to' = 'Spec2')
  
  node_table <- lipid_annot_result %>% 
    filter(lipid_class == temp_lipid_class) %>% 
    mutate(direction_B2 = case_when(p_value_adj_B2 <= 0.2 & fc_B2 > 6/5 ~ 'enriched_B2',
                                    p_value_adj_B2 <= 0.2 & fc_B2 < 5/6 ~ 'depleted_B2',
                                    p_value_adj_B2 > 0.2 ~ 'no_change_B2'),
           direction_B3 = case_when(p_value_adj_B3 <= 0.2 & fc_B3 > 6/5 ~ 'enriched_B3',
                                    p_value_adj_B3 <= 0.2 & fc_B3 < 5/6 ~ 'depleted_B3',
                                    p_value_adj_B3 > 0.2 ~ 'no_change_B3')) %>%
    dplyr::rename('name' = 'variable_id')
  
  # Visualization
  lipid_class_graph <- as_tbl_graph(
    data.frame(
      from = edge_table$from,
      to = edge_table$to,
      weight = edge_table$score
    )
  )
  
  lipid_class_graph <- lipid_class_graph %>% 
    activate(nodes) %>% 
    left_join(node_table, by = 'name') %>% 
    activate(edges) %>% 
    mutate(score = edge_table$score,
           n_frag_match = edge_table$n_frag_match,
           n_frag_nl = edge_table$n_frag_nl,
           n_frag_total = edge_table$n_frag_total)
  
  lipid_class_graph <- tidy_network(object_graph = lipid_class_graph, score_cutoff = 0.5, topK = 10)
  
  final_node_table <- lipid_class_graph %>% 
    activate(nodes) %>% 
    as_tibble()
  
  final_edge_table <- lipid_class_graph %>% 
    activate(edges) %>% 
    as_tibble()
  
  result <- list(node_table = final_node_table,
                 edge_table = final_edge_table,
                 graph = lipid_class_graph)
  
  return(result)
})

names(list_lipid_networks) <- unique_lipid_class


save(list_lipid_networks, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/list_lipid_molecular_network_260708.RData')

# Visualize the whole network --------------------------------------------------

temp_graph_list <- lapply(list_lipid_networks, function(x){
  x$graph
})

temp_null_idx <- sapply(temp_graph_list, function(x){
  is.null(x)
}) %>% which()

temp_graph_list <- temp_graph_list[-temp_null_idx]

combined_graph <- as_tbl_graph(temp_graph_list[[1]])

for(i in 2:length(temp_graph_list)) {
  tbl_graph <- as_tbl_graph(temp_graph_list[[i]])
  combined_graph <- bind_graphs(combined_graph, tbl_graph)
}

set.seed(241125)
temp_plot <- ggraph(combined_graph, layout = 'fr') +
  geom_edge_link(aes(size = score), color = 'gray') +
  geom_node_point(aes(size = centrality_pagerank(), color = direction_B3)) +
  geom_node_text(aes(label = abbr_name), repel = TRUE) +
  scale_color_manual(values = c('enriched_B3' = '#fc8070', 'depleted_B3' = '#7fb1d3', 'no_change_B3' = 'gray')) +
  scale_size_continuous(range = c(1, 10)) +
  theme_void() +
  theme(legend.position = 'top',
        family = 'Arial')

# ggsave(temp_plot, 
#        filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure5/lipid_molecular_network_250401.pdf', 
#        width = 14, height = 8)

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure5/lipid_molecular_network_260708.pdf', 
       width = 12, height = 8)

rm(list = ls());gc()


as.data.frame(combined_graph)
combined_graph %>% activate(nodes) %>% as_tibble() 


################################################################################
# Boxplots of representative PE lipids  ----------------------------------------
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_enrollment_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_annot_result_merge_250401.RData')

library(ggpubr)
library(rstatix)
library(sjmisc)
library(tidymass)
library(introdataviz)
library(tidyverse)

# extract significanltly changed PEs
lipid_annot_result_PE <- lipid_annot_result %>% 
  filter(lipid_class == 'PE')

changed_pes <- lipid_annot_result_PE %>% 
  select(variable_id, abbr_name, compound_name, adduct, confidence_level,
         p_value_B2, p_value_adj_B2, p_value_B3, p_value_adj_B3) %>%
  dplyr::filter(p_value_B3 < 0.05) %>% 
  dplyr::filter(abbr_name %in% c("PE 34:1", "PE 34:2", "PE 36:3", "PE 36:4", "PE 38:4", "PE 38:5", "PE 42:10"))

lipid_pe <- changed_pes$abbr_name
temp_feature <- changed_pes$variable_id

# B1 vs B3 comparison ----------------------------------------------------------
temp_data <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% temp_feature) %>%
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B3')) %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_lipid_enrollment@sample_info) %>% 
  mutate(phenotype_complicate = phenotype_group2)

temp_data <- temp_data %>% 
  pivot_longer(cols = temp_feature, 
               names_to = 'variable_id', 
               values_to = 'value') %>% 
  left_join(lipid_annot_result_PE, by = 'variable_id') %>% 
  select(sample_id:visit_number, phenotype_complicate, variable_id:value, id, abbr_name, lipid_class:source)

stat_wilcox <- temp_data %>%
  group_by(abbr_name) %>%
  wilcox_test(value ~ phenotype_complicate) %>%
  ungroup() %>%
  mutate(
    p = signif(p, digits = 2),
    p.label = format(p, digits = 2, scientific = TRUE, trim = TRUE)
  )

y_pos <- temp_data %>%
  group_by(abbr_name) %>%
  summarise(
    y_min = min(value, na.rm = TRUE),
    y_max = max(value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(y_span = pmax(y_max - y_min, 0.4))

stat_wilcox <- stat_wilcox %>%
  left_join(y_pos, by = 'abbr_name') %>%
  mutate(
    x_id = as.numeric(factor(abbr_name, levels = unique(temp_data$abbr_name))),
    xmin = x_id - 0.2,
    xmax = x_id + 0.2,
    y.position = 2.5
  )

# replace p-values with adjusted p-values for B3 comparison
stat_wilcox  <- stat_wilcox %>% 
  mutate(p = changed_pes$p_value_adj_B3) %>% 
  mutate(p.label = format(p, digits = 2, scientific = TRUE, trim = TRUE))

temp_plot <- ggplot(temp_data, aes(x = abbr_name, y = value)) +
  geom_boxplot(aes(fill = phenotype_complicate), outliers = FALSE, position = position_dodge(width = 0.7), width = 0.65) +
  stat_pvalue_manual(
    stat_wilcox,
    label = 'p.label',
    xmin = 'xmin',
    xmax = 'xmax',
    y.position = 'y.position',
    inherit.aes = FALSE,
    tip.length = 0,
    bracket.size = 0.35,
    size = 3
  ) +
  scale_fill_manual(name = 'Group', values = c('B1' = '#c5c7c9',
                                               'B3' = '#fc8070')) +
  xlab('PE lipids') +
  ylab('Z-Score') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.9, 0.7))


ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260422_PE_lipids_update/PE_boxplot_with_points_2_20260422.pdf', 
       width = 12, height = 6)


# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# readr::write_csv(temp_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/figure4e_PE_lipid_boxplot_B1B3_260721.csv')
# 

# B1 vs B2 comparison ----------------------------------------------------------
temp_data <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% temp_feature) %>%
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(phenotype_group2 %in% c('B1', 'B2')) %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_lipid_enrollment@sample_info) %>% 
  mutate(phenotype_complicate = phenotype_group2)

temp_data <- temp_data %>% 
  pivot_longer(cols = temp_feature, 
               names_to = 'variable_id', 
               values_to = 'value') %>% 
  left_join(lipid_annot_result_PE, by = 'variable_id') %>% 
  select(sample_id:visit_number, phenotype_complicate, variable_id:value, id, abbr_name, lipid_class:source)


library(introdataviz)
library(tidyverse)
library(ggpubr)
library(rstatix)

stat_wilcox <- temp_data %>%
  group_by(abbr_name) %>%
  wilcox_test(value ~ phenotype_complicate) %>%
  ungroup() %>%
  mutate(
    p = signif(p, digits = 2),
    p.label = format(p, digits = 2, scientific = TRUE, trim = TRUE)
  )

y_pos <- temp_data %>%
  group_by(abbr_name) %>%
  summarise(
    y_min = min(value, na.rm = TRUE),
    y_max = max(value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(y_span = pmax(y_max - y_min, 0.4))

stat_wilcox <- stat_wilcox %>%
  left_join(y_pos, by = 'abbr_name') %>%
  mutate(
    x_id = as.numeric(factor(abbr_name, levels = unique(temp_data$abbr_name))),
    xmin = x_id - 0.2,
    xmax = x_id + 0.2,
    y.position = 2.5
  )

# replace p-values with adjusted p-values for B2 comparison
stat_wilcox  <- stat_wilcox %>% 
  mutate(p = changed_pes$p_value_adj_B2) %>% 
  mutate(p.label = format(p, digits = 2, scientific = TRUE, trim = TRUE))

temp_plot <- ggplot(temp_data, aes(x = abbr_name, y = value)) +
  # introdataviz::geom_split_violin(alpha = .4, trim = FALSE, scale = 'width') +
  geom_boxplot(aes(fill = phenotype_complicate), outliers = FALSE, position = position_dodge(width = 0.7), width = 0.65) +
  stat_pvalue_manual(
    stat_wilcox,
    label = 'p.label',
    xmin = 'xmin',
    xmax = 'xmax',
    y.position = 'y.position',
    inherit.aes = FALSE,
    tip.length = 0,
    bracket.size = 0.35,
    size = 3
  ) +
  scale_fill_manual(name = 'Group', values = c('B1' = '#c5c7c9',
                                               'B2' = '#7fb1d3')) +
  xlab('PE lipids') +
  ylab('Z-Score') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.9, 0.7))


ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260422_PE_lipids_update/PE_boxplot_with_points_2_20260422.pdf', 
       width = 12, height = 6)

# 
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# readr::write_csv(temp_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig8b_PE_lipid_boxplot_B1B2_260721.csv')



################################################################################
# PE index ---------------------------------------------------------------------

# calculate the PE index for each sample ---------------------------------------
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_enrollment_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_annot_result_merge_250401.RData')
load('~/Project/00_IBD_project/Data/20241205_CD_complicates/behaviour_time_table_241205.RData')

# Confident PE IDs (level 1):
# PE 34:1: M717T730_c18_neg
# PE 36:4: M739T701_c18_neg
# PE 38:4: M767T732_c18_neg
# PE 38:6: M763T694_1_c18_neg

feature_id_pe <- lipid_annot_result %>% 
  filter(lipid_class == 'PE') %>%
  arrange(p_value_adj_B3) %>% 
  filter(confidence_level == 'Level1') %>% 
  select(variable_id:msms_matched_frag, p_value_adj_B3, everything()) %>%
  pull(variable_id)


# confident PC IDs: 
# PC 34:1: M761T185_hilic_pos
# PC 36:4: M783T685_c18_pos
# PC 38:4: M811T716_c18_pos
# PC 38:6: M807T701_1_c18_pos

feature_id_pc <- lipid_annot_result %>% 
  filter(lipid_class == 'PC') %>%
  arrange(p_value_adj_B3) %>% 
  filter(confidence_level == 'Level1') %>% 
  select(variable_id:msms_matched_frag, p_value_adj_B3, everything()) %>%
  arrange(match(abbr_name, c('PC 34:1', 'PC 36:4', 'PC 38:4', 'PC 38:6'))) %>%
  filter(abbr_name %in% c('PC 34:1', 'PC 36:4', 'PC 38:4', 'PC 38:6')) %>% 
  pull(variable_id)


behaviour_time_table_update <- behaviour_time_table_update %>% 
  filter(status %in% c('B1', 'B2', 'B3'))

temp_update_lipid_result <- object_lipid_enrollment %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(sample_id %in% behaviour_time_table_update$sample_id) %>%
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(variable_id %in% c(feature_id_pc, feature_id_pe)) %>% 
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column('sample_id') %>%
  as_tibble()


pe_index1 <- lapply(2:5, function(x){
  temp_update_lipid_result[[x]]/temp_update_lipid_result[[x+4]]
}) %>% 
  do.call(cbind, .) %>% 
  apply(., 1, mean)


pe_index1 <- tibble(sample_id = temp_update_lipid_result$sample_id, 
                    pe_index1 = pe_index1)


pe_index_table <- behaviour_time_table_update %>% 
  left_join(pe_index1, by = 'sample_id') %>% 
  arrange(desc(pe_index1)) %>% 
  mutate(idx = seq(n())) %>% 
  arrange(status) %>% 
  mutate(status = as.factor(status))

save(pe_index_table, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/pe_index_table_250401.RData')

# visualize the PE index distribution ------------------------------------------
temp_data <- behaviour_time_table_update %>% 
  left_join(pe_index1, by = 'sample_id') %>% 
  arrange(desc(pe_index1)) %>% 
  mutate(idx = seq(n())) %>% 
  arrange(status) %>% 
  mutate(status = as.factor(status))


temp_plot <- ggplot(temp_data) +
  geom_point(aes(x = idx, y = -pe_index1, color = status, size = status)) +
  scale_color_manual(values = c('B3' = '#fb8172', 
                                'B2' = '#c5c7c9', 
                                'B1' = '#c5c7c9')) +
  scale_size_manual(values = c('B3' = 2,
                               'B2' = 1,
                               'B1' = 1)) +
  geom_hline(yintercept = -300) +
  xlab('Sample index') +
  ylab('-PE index') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.9, 0.2))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure4/plot_pe_index_250401.pdf', 
       width = 6, height = 6)


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_source_data <- temp_data %>% 
  select(sample_id, patient_id, status, pe_index1, idx)

readr::write_csv(temp_source_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig8c_PE_index_260721.csv')

################################################################################
# Normalization approach + new label + Old lipids ------------------------------
# competing risk model ---------------------------------------------------------

library(tidycmprsk)
library(ggsurvfit)
library(gtsummary)

behaviour_time_table2 <- temp_data %>% 
  mutate(pe_index1_label = case_when(pe_index1 >= 300 ~ 'Low PE',
                                     pe_index1 < 300 ~ 'High PE')) %>% 
  mutate(status_id = case_when(
    status == 'B1' ~ '0',
    status == 'B2' ~ '1',
    status == 'B3' ~ '2')) %>%
  mutate(status_id = recode_factor(status_id, '0' = '0', '1' = '1', '2' = '2'))

# contain B2 & B3
cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
  ggcuminc(outcome = c("1", "2")) +
  labs(
    x = "Days"
  ) +
  geom_vline(xintercept = 365) +
  geom_vline(xintercept = 1095) +
  add_confidence_interval() +
  add_risktable(risktable_group = 'risktable_stats')

# only B2
year_day <- 365
temp_plot <- cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
  ggcuminc(outcome = c("1")) +
  scale_x_continuous(breaks = c(0, seq(7)*year_day), labels = c(0, seq(7)*year_day), limits = c(0, 2900)) +
  labs(
    x = "Days"
  ) +
  # geom_vline(xintercept = 365) +
  # geom_vline(xintercept = 1095) +
  add_confidence_interval() +
  ZZWtool::ZZWTheme() + 
  add_risktable(risktable_group = 'risktable_stats') +
  theme(legend.position = c(0.1, 0.9))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure5/plot_pe_index_B2_250401.pdf', 
       width = 6, height = 6)

cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
  tbl_cuminc(times = c(seq(7)*year_day),
             outcomes = '1') %>% 
  add_p()

# 
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- behaviour_time_table2 %>%
#   dplyr::filter(status %in% c('B1', 'B2')) %>%
#   select(sample_id, patient_id, status, status_id, change_time, pe_index1, pe_index1_label)
# 
# readr::write_csv(temp_source_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig8d_PE_index_B2_260721.csv')

# only B3

year_day <- 365
temp_plot <- cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
  ggcuminc(outcome = c("2")) +
  scale_x_continuous(breaks = c(0, seq(7)*year_day), labels = c(0, seq(7)*year_day), limits = c(0, 2900)) +
  # xlim(0, 2800) +
  # xlim(0, 7*year_day) +
  labs(
    x = "Days"
  ) +
  add_confidence_interval() +
  ZZWtool::ZZWTheme() +
  add_risktable(risktable_group = 'risktable_stats') + 
  theme(legend.position = c(0.1, 0.9))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure4/plot_pe_index_B3_250401.pdf', 
       width = 6, height = 6)


cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
  tbl_cuminc(times = c(seq(7)*year_day),
             outcomes = '2') %>% 
  add_p()


# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- behaviour_time_table2 %>% 
#   dplyr::filter(status %in% c('B1', 'B3')) %>% 
#   select(sample_id, patient_id, status, status_id, change_time, pe_index1, pe_index1_label)
# 
# readr::write_csv(temp_source_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/fig4f_PE_index_B3_260721.csv')

# # show B2 and B3 in the one table
# cuminc(Surv(change_time, status_id) ~ pe_index1_label, data = behaviour_time_table2) %>% 
#   tbl_cuminc(times = c(seq(7)*year_day),
#              outcomes = c('1', '2')) %>% 
#   add_p()
# 
# crr(Surv(change_time, status_id) ~ pe_index1_label, 
#     data = behaviour_time_table2,
#     failcode = c('1'))
# 
# crr(Surv(change_time, status_id) ~ pe_index1_label, 
#     data = behaviour_time_table2,
#     failcode = c('2'))