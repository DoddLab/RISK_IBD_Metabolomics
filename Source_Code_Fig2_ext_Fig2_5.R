setwd('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD')


################################################################################
# CD onset analyses ------------------------------------------------------------

dir.create('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD', showWarnings = FALSE, recursive = TRUE)

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')


# Differential analysis  -------------------------------------------------------
library(tidyverse)
library(tidymass)
library(sjmisc)

object_final2 <- object_final %>% scale_data(center = TRUE, method = 'auto')

temp_data_stat <- object_final2@expression_data %>% 
  rotate_df() %>% 
  tibble::rownames_to_column(var = 'sample_id') %>% 
  left_join(object_final2@sample_info, by = 'sample_id') %>% 
  filter(visit_number == 'enrollment')

variable_name <- object_final2@annotation_table$variable_id

fix_vs <- c('phenotype_group1', 'gender', 'age', 'use_antibiotics', 'race')

adjusted_result <- pbapply::pblapply(variable_name, function(x){
  form <- as.formula(paste0(paste0(x, " ~ "), paste(fix_vs, collapse=' + ' )))
  model <- glm(form, data = temp_data_stat)
  stat_result <- summary(model)
  
  p_value <- stat_result$coefficients[2, 4]
  coefficient <- stat_result$coefficients[2,1]
  
  result <- tibble::tibble(variable_id = x, 
                           coefficient = coefficient,
                           p_value = p_value)
  
  return(result)
})

adjusted_result <- adjusted_result %>% 
  dplyr::bind_rows()
adjust_p <- adjusted_result$p_value %>% p.adjust(method = 'BH')
adjusted_result <- adjusted_result %>% mutate(p_value_adjusted = adjust_p)
sum(adjusted_result$p_value_adjusted <= 0.05)

dir.create('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/', showWarnings = FALSE, recursive = TRUE)
save(adjusted_result, 
     file = '~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/adjusted_result_250328.RData')

rm(list = ls())


# merge into the object_stat ----------------------------------------------------------
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/adjusted_result_250328.RData')
adjusted_result %>% count(p_value_adjusted <= 0.05)
rm(adjusted_result);gc()

object_enrollment <- object_final %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(visit_number == 'enrollment')

id_non_ibd <- object_enrollment@sample_info %>% dplyr::filter(phenotype_group1 == 'non_IBD') %>% pull(sample_id)
id_cd <- object_enrollment@sample_info %>% dplyr::filter(phenotype_group1 == 'CD') %>% pull(sample_id)

object_stat <-
  mutate_p_value(
    object = object_enrollment,
    control_sample_id = id_non_ibd,
    case_sample_id = id_cd,
    method = "wilcox.test",
    p_adjust_methods = "BH"
  )

object_stat <-
  mutate_fc(object = object_stat, 
            control_sample_id = id_non_ibd, 
            case_sample_id = id_cd, 
            mean_median = "mean")

# replace P-values with generalized linear model (age and gender adjusted)
# remove NA statistics
object_stat@variable_info <- object_stat@variable_info %>% select(variable_id:rt, p_value:fc)
object_stat@variable_info$p_value <- adjusted_result$p_value
object_stat@variable_info$p_value_adjust <- adjusted_result$p_value_adjusted
object_stat@variable_info <- object_stat@variable_info %>% 
  mutate(coefficient = adjusted_result$coefficient)

save(object_stat, 
     file = '~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')


################################################################################
# Vocano plot Fig 2a -----------------------------------------------------------
# Highlight with representative metabolite classes -----------------------------

load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')

# colored with highted metabolite classes

result_annot_adjust <- object_stat %>% 
  extract_variable_info() %>% 
  dplyr::rename('mz' = 'mz.x',
                'rt' = 'rt.x' ) %>% 
  dplyr::select(-c('mz.y', 'rt.y'))

# tryptophan metabolites
tryptophan_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Tryptophan Metabolism') %>% 
  pull('variable_id')


# Phenylalanine and Tyrosine Metabolism, **Phenyl/PAGln**, Phenylacetylglutamine
phenyl_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Phenylalanine and Tyrosine Metabolism' | Compound.name %in% c('2-Phenylpropionic acid', '4-Hydroxyhippuric acid', 'Hippuric acid')) %>% 
  pull('variable_id')

# Bile acids
bile_acid_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Steroids') %>% 
  pull('variable_id')


# LCFA.VLCFA
fatty_acid_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass %in% c('Long-chain fatty acids', 'Very long-chain fatty acids')) %>% 
  pull('variable_id')


# Hemo metabolites
heme_feature_id <- result_annot_adjust %>%
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Hemoglobin and Porphyrin Metabolism') %>% 
  pull('variable_id')



# assign different colors to different classes of metabolites
temp_data <- result_annot_adjust %>% 
  mutate(class_color = case_when(
    variable_id %in% tryptophan_feature_id ~ 'Tryptophan',
    variable_id %in% phenyl_feature_id ~ 'Phenyl/PAGln',
    variable_id %in% bile_acid_feature_id ~ 'Bile acids',
    variable_id %in% fatty_acid_feature_id ~ 'LCFA/VLCFA',
    variable_id %in% heme_feature_id ~ 'Heme',
    TRUE ~ 'Other'
  )) %>% 
  arrange(match(class_color, c('Other', 'Tryptophan', 'Phenyl/PAGln', 'Bile acids', 'LCFA/VLCFA', 'Heme')))


text_data <- temp_data %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(class_color != 'Other') %>% 
  mutate(log2_fc = log2(fc),
         log10_p = -log10(p_value_adjust)) %>% 
  select(variable_id, Compound.name, log2_fc, log10_p)


library(tidyverse)
temp_plot <- ggplot(temp_data) +
  geom_point(aes(x = log2(fc), y = -log10(p_value_adjust), color = class_color, size = -log10(p_value_adjust)), shape = 19, alpha = 0.8) +
  scale_size_continuous(range = c(3, 5),
                        name = 'Log10(P)') +
  scale_colour_manual(values = c('Tryptophan' = "#fc8070",
                                 'Phenyl/PAGln' = '#fdb562',
                                 'Bile acids' = '#ffed6f',
                                 'LCFA/VLCFA' = "#7fb1d3",
                                 'Heme' = "#8cd4c8",
                                 'Other' = "#c5c7c9")) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  xlim(-6, 6) +
  ylab('-Log10(P-adjusted)') +
  xlab('Log2(Fold change)') +
  ZZWtool::ZZWTheme() +
  ggrepel::geom_text_repel(aes(log2_fc, log10_p,
                               label = Compound.name),
                           data = text_data) +
  theme(legend.position = c(0.9, 0.7))

dir.create('~/Project/00_IBD_project/Figure/260302/260325_update/', recursive = TRUE, showWarnings = FALSE)
ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/260302/260325_update/volcano_plot_annotated_260325.pdf', 
       width = 12, 
       height = 6)



temp_data <- temp_data %>% 
  dplyr::select(variable_id, Compound.name, fc, p_value_adjust, class_color) %>% 
  dplyr::arrange(p_value_adjust)

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)

readr::write_csv(temp_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2_volcano_plot_260722.csv')


################################################################################
# Vocano plot Extended Fig 3a ------------------------------------------------------------------
load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')

library(tidyverse)

result_annot_adjust <- object_stat %>% 
  extract_variable_info() %>% 
  dplyr::rename('mz' = 'mz.x',
                'rt' = 'rt.x' ) %>% 
  dplyr::select(-c('mz.y', 'rt.y'))


temp_data <- result_annot_adjust %>% 
  mutate(class_color = case_when(
    (fc > 1 & p_value_adjust <= 0.05) ~ 'increased',
    (fc < 1 & p_value_adjust <= 0.05) ~ 'decreased',
    p_value_adjust > 0.05 ~ 'no_change'
  ))

text_data <- temp_data %>% 
  dplyr::filter(p_value_adjust <= 0.05) %>% 
  # filter(class_color != 'no_change') %>% 
  mutate(log2_fc = log2(fc),
         log10_p = -log10(p_value_adjust)) %>% 
  dplyr::select(variable_id, Compound.name, log2_fc, log10_p)


library(tidyverse)
temp_plot <- ggplot(temp_data) +
  geom_point(aes(x = log2(fc), y = -log10(p_value_adjust), color = class_color, size = -log10(p_value_adjust)), shape = 19, alpha = 0.8) +
  scale_size_continuous(range = c(1, 3),
                        name = 'Log10(P)') +
  scale_colour_manual(values = c('increased' = "#fc8070",
                                 'decreased' = '#7fb1d3',
                                 'no_change' = '#c5c7c9')) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  xlim(-6, 6) +
  ylab('-Log10(P-adjusted)') +
  xlab('Log2(Fold change)') +
  ZZWtool::ZZWTheme() +
  ggrepel::geom_text_repel(aes(log2_fc, log10_p,
                               label = Compound.name),
                           data = text_data, size = 2,
                           max.overlaps = 20) +
  theme(legend.position = c(0.9, 0.7))

dir.create('~/Project/00_IBD_project/Figure/260710/', recursive = TRUE, showWarnings = FALSE)

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/260710/volcano_plot_annotated_260710.pdf', 
       width = 8, 
       height = 8)

temp_data <- temp_data %>% 
  dplyr::select(variable_id, Compound.name, fc, p_value_adjust, class_color) %>% 
  dplyr::arrange(p_value_adjust)

readr::write_csv(temp_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/volcano_plot_extended_fig3a_260710.csv')


################################################################################
# tryptophan metabolites -------------------------------------------------------
# tryptophan metabolites
tryptophan_features <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Tryptophan Metabolism') %>% 
  arrange(p_value_adjust)

tryptophan_feature_id <- tryptophan_features %>% pull('variable_id') %>% head(5)
tryptophan_feature_name <- tryptophan_features %>% pull('Compound.name') %>% head(5)


plot_data <- object_stat %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1) %>% 
  left_join(object_stat %>% 
              activate_mass_dataset(what = 'variable_info') %>%
              filter(variable_id %in% tryptophan_feature_id) %>% 
              extract_expression_data() %>% 
              rownames_to_column(var = 'variable_id') %>%
              arrange(match(variable_id, tryptophan_feature_id)) %>% 
              pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'z_score'),
            by = 'sample_id'
  ) %>% 
  group_by(variable_id) %>%
  mutate(z_score = as.numeric(scale(z_score)),
         phenotype_group1 = factor(phenotype_group1, levels = c('non_IBD', 'CD'))) %>% 
  ungroup() %>% 
  mutate(compound_name = match(variable_id, tryptophan_feature_id) %>% tryptophan_feature_name[.],
         compound_name = factor(compound_name, levels = tryptophan_feature_name))


library(ggplot2)
library(dplyr)
library(ggbeeswarm) 
library(ggpubr)
library(ggsignif)
library(rstatix)


stat.test <- plot_data %>%
  group_by(compound_name) %>%
  wilcox_test(z_score ~ phenotype_group1) %>%
  add_significance() 

stat.test <- stat.test %>%
  add_xy_position(
    x = "compound_name", 
    dodge = 0.7, 
    fun = "max",     
    step.increase = 0.1
  )

stat.test$y.position <- 5.5

temp_plot <- ggplot(plot_data, aes(x = compound_name, y = z_score)) +
  geom_point(aes(color = phenotype_group1), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), 
             shape = 19, size = 1.5, alpha = 0.7) +
  geom_boxplot(aes(fill = phenotype_group1), outlier.shape = NA, alpha = 0.6, width = 0.7, color = "black") +
  # scale_fill_manual(values = c("CD" = "#A9CCE3", "non_IBD" = "#F5CBA7")) +
  # scale_color_manual(values = c("CD" = "#5DADE2", "non_IBD" = "#EB984E")) +
  scale_fill_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  scale_color_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  xlab("Compound name") +
  ylab("Z-score") +
  ZZWtool::ZZWTheme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_blank(),     
    legend.position = c(0.8, 0.9),            
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  ) +
  ylim(-2, 6) +
  stat_pvalue_manual(
    stat.test,
    label = "p", 
    tip.length = 0.02,
    hide.ns = FALSE 
  )

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/260302/260325_update/tryptophan_boxplot_260325.pdf', 
       width = 10, 
       height = 6)


readr::write_csv(plot_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2b_tryptophan_boxplot_260325.csv')

################################################################################
# Phenyl metabolites -----------------------------------------------------------

phenyl_features <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Phenylalanine and Tyrosine Metabolism' | Compound.name %in% c('2-Phenylpropionic acid', '4-Hydroxyhippuric acid', 'Hippuric acid')) %>% 
  arrange(p_value_adjust)

phenyl_feature_id <- phenyl_features %>% pull('variable_id') %>% head(9)
phenyl_feature_name <- phenyl_features %>% pull('Compound.name') %>% head(9)

plot_data <- object_stat %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1) %>% 
  left_join(object_stat %>% 
              activate_mass_dataset(what = 'variable_info') %>%
              filter(variable_id %in% phenyl_feature_id) %>% 
              extract_expression_data() %>% 
              rownames_to_column(var = 'variable_id') %>%
              arrange(match(variable_id, phenyl_feature_id)) %>% 
              pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'z_score'),
            by = 'sample_id'
  ) %>% 
  group_by(variable_id) %>%
  mutate(z_score = as.numeric(scale(z_score)),
         phenotype_group1 = factor(phenotype_group1, levels = c('non_IBD', 'CD'))) %>% 
  ungroup() %>% 
  mutate(compound_name = match(variable_id, phenyl_feature_id) %>% phenyl_feature_name[.],
         compound_name = factor(compound_name, levels = phenyl_feature_name))


library(ggplot2)
library(dplyr)
library(ggbeeswarm) 
library(ggpubr)
library(ggsignif)
library(rstatix)


stat.test <- plot_data %>%
  group_by(compound_name) %>%
  wilcox_test(z_score ~ phenotype_group1) %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(
    x = "compound_name", 
    dodge = 0.7, 
    fun = "max",     
    step.increase = 0.1
  )

stat.test$y.position <- 7
# replace the P with orginal P values
stat.test$p <- match(stat.test$compound_name, phenyl_feature_name) %>% phenyl_features$p_value_adjust[.] %>% signif(3)

temp_plot <- ggplot(plot_data, aes(x = compound_name, y = z_score)) +
  geom_point(aes(color = phenotype_group1), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), 
             shape = 19, size = 1.5, alpha = 0.7) +
  geom_boxplot(aes(fill = phenotype_group1), outlier.shape = NA, alpha = 0.6, width = 0.7, color = "black") +
  scale_fill_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  scale_color_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  xlab("Compound name") +
  ylab("Z-score") +
  ZZWtool::ZZWTheme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_blank(),     
    legend.position = c(0.8, 0.9),            
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  ) +
  ylim(-2, 8) +
  stat_pvalue_manual(
    stat.test,
    label = "p", 
    tip.length = 0.02,
    hide.ns = FALSE 
  )

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/260302/260325_update/phenyl_boxplot_260325.pdf', 
       width = 12, 
       height = 6)

readr::write_csv(plot_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2c_phenyl_boxplot_260720.csv')


################################################################################
# Heme metabolites -----------------------------------------------------------

heme_features <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Hemoglobin and Porphyrin Metabolism') %>% 
  arrange(p_value_adjust)

heme_feature_id <- heme_features %>% pull('variable_id') %>% head(3)
heme_feature_name <- heme_features %>% pull('Compound.name') %>% head(3)


plot_data <- object_stat %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1) %>% 
  left_join(object_stat %>% 
              activate_mass_dataset(what = 'variable_info') %>%
              filter(variable_id %in% heme_feature_id) %>% 
              extract_expression_data() %>% 
              rownames_to_column(var = 'variable_id') %>%
              arrange(match(variable_id, heme_feature_id)) %>% 
              pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'z_score'),
            by = 'sample_id'
  ) %>% 
  group_by(variable_id) %>%
  mutate(z_score = as.numeric(scale(z_score)),
         phenotype_group1 = factor(phenotype_group1, levels = c('non_IBD', 'CD'))) %>% 
  ungroup() %>% 
  mutate(compound_name = match(variable_id, heme_feature_id) %>% heme_feature_name[.],
         compound_name = factor(compound_name, levels = heme_feature_name))


library(ggplot2)
library(dplyr)
library(ggbeeswarm) 
library(ggpubr)
library(ggsignif)
library(rstatix)

stat.test <- plot_data %>%
  group_by(compound_name) %>%
  wilcox_test(z_score ~ phenotype_group1) %>%
  add_significance()

stat.test <- stat.test %>%
  add_xy_position(
    x = "compound_name", 
    dodge = 0.7, 
    fun = "max",     
    step.increase = 0.1
  )

stat.test$y.position <- 7
stat.test$p <- match(stat.test$compound_name, heme_feature_name) %>% heme_features$p_value_adjust[.] %>% signif(3)

temp_plot <- ggplot(plot_data, aes(x = compound_name, y = z_score)) +
  geom_point(aes(color = phenotype_group1), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), 
             shape = 19, size = 1.5, alpha = 0.7) +
  geom_boxplot(aes(fill = phenotype_group1), outlier.shape = NA, alpha = 0.6, width = 0.7, color = "black") +
  scale_fill_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  scale_color_manual(values = c("CD" = "#fc8070", "non_IBD" = "#7fb1d3")) +
  xlab("Compound name") +
  ylab("Z-score") +
  ZZWtool::ZZWTheme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_blank(),     
    legend.position = c(0.8, 0.9),            
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  ) +
  ylim(-2, 8) +
  stat_pvalue_manual(
    stat.test,
    label = "p",
    tip.length = 0.02,
    hide.ns = FALSE 
  )

ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/260302/260325_update/heme_boxplot_260325.pdf', 
       width = 8, 
       height = 6)

readr::write_csv(plot_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2d_heme_boxplot_260720.csv')


################################################################################
# Bile acid metabolites -----------------------------------------------------------

library(ggjoy)

bile_acid_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass == 'Steroids') %>% 
  pull('variable_id')

bile_acid_result <- object_stat %>%
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(variable_id %in% bile_acid_feature_id) %>% 
  dplyr::arrange(match(variable_id, bile_acid_feature_id)) %>%
  extract_variable_info()



temp_id <- bile_acid_result$variable_id
temp_compound_name <- bile_acid_result$Compound.name
temp_variable_id <- bile_acid_result$variable_id
temp_p_value <- bile_acid_result$p_value_adjust
temp_fc <- bile_acid_result$fc

plot_list <- lapply(seq_along(temp_id), function(i){
  temp_p <- ggplot_mass_dataset(object = object_stat, 
                                direction = 'variable',
                                variable_id = temp_id[i])
  
  non_ibd_median <- temp_p$data %>% 
    dplyr::filter(phenotype_group1 == 'non_IBD') %>% 
    dplyr::pull(value) %>% 
    median()
  
  
  temp_data <- temp_p$data %>% 
    select(sample_id, phenotype_group1, value) %>% 
    mutate(log_value = log2(value/non_ibd_median),
           compound_name = temp_compound_name[i],
           variable_id = temp_variable_id[i])
  
  temp_label <- paste0(temp_compound_name[i], '\n', 'p-value: ', signif(temp_p_value[i], 4), '\n', 'FC: ', signif(temp_fc[i], 2))
  
  temp_plot <- ggplot(temp_data, aes(x = log_value, y = phenotype_group1, fill = phenotype_group1, height = ..density..)) +
    geom_density_ridges2(scale = 30, stat = "density") +
    scale_y_discrete(expand = c(0.01, 0), name = 'Density') +
    scale_x_continuous(expand = c(0.01, 0), name = 'Log2 (Abundance relative to non-IBD median)') +
    scale_fill_manual(values = c('CD' = '#ffed6f',
                                 'non_IBD' = '#8cd4c7'),
                      label = c('CD' = 'CD',
                                'non_IBD' = 'Non-IBD'),
                      name = 'Group') +
    xlim(-7, 7)+
    ggtitle(temp_label) +
    xlab('Log2 (Abundance relative to non-IBD median)') +
    # ZZWtool::ZZW_annotate_text2(label = temp_label, x = 0.1, y = 0.9) +
    theme_ridges(center_axis_labels = TRUE, grid = FALSE) + 
    theme(legend.position = c(0.8, 0.9))
  
  return(temp_plot)
})


names(plot_list) <- temp_compound_name

dir.create('~/Project/00_IBD_project/Figure/260302/260325_update/bile_acid_distribution', 
           showWarnings = FALSE, recursive = TRUE)  
walk(seq_along(plot_list), function(i){
  temp_plot1 <- plot_list[[i]]
  ggsave(plot = temp_plot1, 
         filename = file.path('~/Project/00_IBD_project/Figure/260302/260325_update/bile_acid_distribution', 
                              paste0(temp_compound_name[i], '.pdf')), 
         width = 6, height = 2)
})

# 
# temp_idx <- c("Cholic acid", "Deoxycholic acid", "cholic acid-C10:0", "lithocholic acid-C12:0", "Ursodeoxycholic acid 3-sulfate", "Glycocholic acid", "Taurocholic acid", "Glycodeoxycholic acid", "Taurodeoxycholic acid", "Glycochenodeoxycholic acid")
# temp_idx <- match(temp_idx, names(plot_list))
# 
# temp_data <- lapply(temp_idx, function(i){
#   plot_list[[i]]$data
# }) %>% bind_rows()
# 
# readr::write_csv(temp_data, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2e_bile_acid_distribution_260720.csv')

################################################################################
# Long chain/Very long chain fatty acid metabolites ---------------------------------------

fatty_acid_feature_id <- result_annot_adjust %>% 
  filter(p_value_adjust <= 0.05) %>% 
  filter(metabolon_subclass %in% c('Long-chain fatty acids', 'Very long-chain fatty acids')) %>% 
  pull('variable_id')

fatty_acid_result <- object_stat %>%
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(variable_id %in% fatty_acid_feature_id) %>% 
  dplyr::arrange(match(variable_id, bile_acid_feature_id)) %>%
  extract_variable_info()

fatty_acid_result_manual_check <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241104_metabolite_class_distribution_CD_nonIBD/annot_fatty_acids_manual_2_241104.xlsx')
temp_data <- fatty_acid_result_manual_check %>% select(variable_id, C_number:COOH)
fatty_acid_result <- fatty_acid_result %>% 
  left_join(temp_data, by = 'variable_id')


temp_data_2 <- fatty_acid_result

temp_id <- temp_data_2 %>% pull(variable_id)
temp_compound_name <- temp_data_2 %>% pull(Compound.name)
temp_c_number <- temp_data_2 %>% pull(C_number)
temp_db_number <- temp_data_2 %>% pull(Double_bound_number)
temp_oh_number <- temp_data_2 %>% pull(OH)
temp_cooh_number <- temp_data_2 %>% pull(COOH)
temp_p_value <- temp_data_2 %>% pull(p_value_adjust)
temp_fc <- temp_data_2 %>% pull(fc)

# plot all joys 
library(ggjoy)

plot_list <- lapply(seq_along(temp_id), function(i){
  temp_p <- ggplot_mass_dataset(object = object_stat, 
                                direction = 'variable',
                                variable_id = temp_id[i])
  
  non_ibd_median <- temp_p$data %>% 
    dplyr::filter(phenotype_group2 == 'non_IBD') %>% 
    dplyr::pull(value) %>% 
    median()
  
  
  temp_data <- temp_p$data %>% 
    select(sample_id, phenotype_group1, value) %>% 
    mutate(log_value = log2(value/non_ibd_median),
           compound_name = temp_compound_name[i],
           variable_id = temp_id[i])
  
  
  temp_label <- paste0('FA ', temp_c_number[i], ":", temp_db_number[i])
  if (temp_oh_number[i] > 0) {
    if (temp_oh_number[i] > 1) {
      temp_label <- paste0(temp_label, ';O', temp_oh_number[i])
    } else {
      temp_label <- paste0(temp_label, ';O')
    }
  }
  
  if (temp_cooh_number[i] > 0) {
    if (temp_cooh_number[i] > 1) {
      temp_label <- paste0(temp_label, ';O', temp_cooh_number[i])
    } else {
      temp_label <- paste0(temp_label, ';O')
    }
  }
  
  temp_label <- paste0(temp_compound_name[i], '\n', temp_label, '\n', 'p-value: ', signif(temp_p_value[i], 4), '\n', 'FC: ', signif(temp_fc[i], 2))
  
  temp_plot <- ggplot(temp_data, aes(x = log_value, y = phenotype_group1, fill = phenotype_group1, height = ..density..)) +
    geom_density_ridges2(scale = 30, stat = "density") +
    scale_y_discrete(expand = c(0.01, 0), name = 'Density') +
    scale_x_continuous(expand = c(0.01, 0), name = 'Log2 (Abundance relative to non-IBD median)') +
    # scale_fill_manual(values = c('CD' = '#FCC900',
    #                              'non_IBD' = '#6293CF'),
    #                   label = c('CD' = 'CD',
    #                             'non_IBD' = 'Non-IBD'),
    #                   name = 'Group') +
    scale_fill_manual(values = c('CD' = '#ffed6f',
                                 'non_IBD' = '#8cd4c7'),
                      label = c('CD' = 'CD',
                                'non_IBD' = 'Non-IBD'),
                      name = 'Group') +
    xlim(-5, 5)+
    ggtitle(temp_label) +
    xlab('Log2 (Abundance relative to non-IBD median)') +
    # ZZWtool::ZZW_annotate_text2(label = temp_label, x = 0.1, y = 0.9) +
    theme_ridges(center_axis_labels = TRUE, grid = FALSE) + 
    theme(legend.position = c(0.8, 0.9))
  
  return(temp_plot)
})



names(plot_list) <- temp_compound_name

dir.create('~/Project/00_IBD_project/Figure/260302/260325_update/LCFA_VLCFA_distribution', 
           showWarnings = FALSE, recursive = TRUE)
walk(seq_along(plot_list), function(i){
  temp_plot1 <- plot_list[[i]]
  ggsave(plot = temp_plot1, 
         filename = file.path('~/Project/00_IBD_project/Figure/260302/260325_update/LCFA_VLCFA_distribution', 
                              paste0(temp_compound_name[i], '.pdf')), 
         width = 6, height = 2)
})



temp_idx <- c("Docosatrienoic acid", "Adrenic Acid", "Docosapentaenoic acid", "Docosahexaenoic acid", "Hexacosanoic acid", "Palmitoleic Acid", "Heptadecanoate", "Elaidic Acid", "Arachidic Acid", "Arachidonic Acid")
temp_idx <- match(temp_idx, names(plot_list))

temp_data <- lapply(temp_idx, function(i){
  plot_list[[i]]$data
}) %>% bind_rows()

readr::write_csv(temp_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2f_LCFA_VLCFA_distribution_260720.csv')





# consider all related FA, like Hydroxy fatty acids, Branched fatty acids ------
fatty_acid_feature_id <- result_annot_adjust %>% 
  filter(metabolon_subclass %in% c("Medium-chain fatty acids",
                                   "Long-chain fatty acids",
                                   "Very long-chain fatty acids",
                                   'Lineolic Acid Metabolism')) %>% 
  pull('variable_id')


fatty_acid_result <- object_stat %>%
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(variable_id %in% fatty_acid_feature_id) %>% 
  dplyr::arrange(match(variable_id, bile_acid_feature_id)) %>%
  extract_variable_info()

fatty_acid_result_manual_check <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241104_metabolite_class_distribution_CD_nonIBD/annot_fatty_acids_manual_2_241104.xlsx')
temp_data <- fatty_acid_result_manual_check %>% select(variable_id, C_number:COOH)
fatty_acid_result <- fatty_acid_result %>% 
  left_join(temp_data, by = 'variable_id')

temp_data_2 <- fatty_acid_result

temp_fa <- fatty_acid_result %>%
  dplyr::filter(metabolon_subclass %in%
                  c("Medium-chain fatty acids",
                    "Long-chain fatty acids",
                    "Very long-chain fatty acids",
                    'Lineolic Acid Metabolism'))

# view statistics
temp_fa %>%
  group_by(C_number > 12) %>%
  count(p_value_adjust <= 0.05)

# contingency table
temp_table <- matrix(c(19,9,1,13), ncol=2,
                     dimnames = list(c('Sig','Unsig'),
                                     c('Long-chain','Medium-chain')))
chisq.test(temp_table)

temp_plot <- temp_fa %>%
  group_by(C_number > 12) %>%
  count(p_value_adjust <= 0.05) %>% 
  mutate(percentage = n/sum(n)) %>% 
  group_by(`C_number > 12`) %>% 
  mutate(pcg = n/sum(n)) %>% 
  ggplot(aes(x = `C_number > 12`, y = pcg, fill = `p_value_adjust <= 0.05`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0) +
  xlab('FA chain length') +
  ylab('Percentage (%)') +
  scale_fill_manual(values = c('TRUE' = '#257a3a',
                               'FALSE' = '#aba8a8'),
                    label = c('TRUE' = 'P < 0.05',
                              'FALSE' = 'P > 0.05'),
                    name = 'P-value') +
  scale_x_discrete(labels = c('TRUE' = 'C > 12',
                              'FALSE' = 'C <= 12')) +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(legend.position = 'none')

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure2/FA_p_value_distribution_250328.pdf', 
       width = 3, height = 6)

temp_data <- temp_fa %>%
  group_by(C_number > 12) %>%
  count(p_value_adjust <= 0.05) %>% 
  mutate(percentage = n/sum(n)) %>% 
  group_by(`C_number > 12`) %>% 
  mutate(pcg = n/sum(n))

readr::write_csv(temp_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/fig2f_right_FA_p_value_distribution_260720.csv')



################################################################################
# Heatmap of the significant metabolites ---------------------------------------


# enrichClass2 -----------------------------------------------------------------

#' @title enrichClass2
#' @author Zhiwei Zhou
#' @description this function was reported by ChemRich (DOI: xxxx)
#' @param table_class_file 'table_class_enrichment.csv'. a csv file, 1-5th colnames need to be named as "compound_name", "order", "pvalue", "foldchange", "class"
#' @param path '.'
#' @export

# path <- '/home/zhouzw/Data_processing/20210522_debug/test/04_biology_intepretation/'
# enrichment_result <- enrichClass(table_class_file = 'table_class_enrichment.csv',
#                                  path = '/home/zhouzw/Data_processing/20210522_debug/test/04_biology_intepretation/')

setGeneric(name = 'enrichClass2',
           def = function(
    # table_class_file = 'table_class_enrichment.csv',
             table_compound,
             path = '.',
             class_cutoff = 3,
             p_cutoff = 0.05
           ){
             # browser()
             # xdf <- readr::read_csv(file.path(path, table_class_file))
             xdf <- table_compound
             clusterids <- names(which(table(xdf$class) >= class_cutoff))
             clusterids <- clusterids[which(clusterids!="")]
             cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
               cl.member <- which(xdf$class==x)
               if( length(which(xdf$pvalue[cl.member] < p_cutoff)) >1 ){
                 pval.cl.member <- xdf$pvalue[cl.member]
                 p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
                 p.test.results$p.value
               } else {
                 1
               }
             })
             cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
             clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)
             
             clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
               dfx <- xdf[which(xdf$class==x),]
               dfx$compound_name[which.min(dfx$pvalue)]
             })
             altrat <- sapply(clusterdf$name, function (k) {
               length(which(xdf$class==k & xdf$pvalue < p_cutoff))/length(which(xdf$class==k))
             })
             uprat <-sapply(clusterdf$name, function (k) {
               length(which(xdf$class==k & xdf$pvalue < p_cutoff & xdf$foldchange > 1.00000000))/length(which(xdf$class==k & xdf$pvalue < p_cutoff))
             })
             clust_s_vec <- sapply(clusterdf$name, function (k) {
               length(which(xdf$class==k))
             })
             clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(xdf$class==k & xdf$pvalue < p_cutoff))})
             clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(xdf$class==k & xdf$pvalue < p_cutoff & xdf$foldchange > 1.00000000))})
             clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(xdf$class==k & xdf$pvalue < p_cutoff & xdf$foldchange < 1.00000000))})
             clusterdf$upratio <- uprat
             clusterdf$altratio <- altrat
             clusterdf$csize <- clust_s_vec
             clusterdf <- clusterdf[which(clusterdf$csize>class_cutoff),]
             clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
             clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(xdf$order[which(xdf$class==x)]) }))
             clusterdf
             clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
               dfx <- xdf[which(xdf$class==x),]
               paste(dfx$compound_name, collapse="<br>")
             }) ## this one is the label on the tooltip of the ggplotly plot.
             clustdf <- clusterdf[which(clusterdf$pvalues!=1),]
             
             clustdf <- tibble::as_tibble(clustdf)
             return(clustdf)
           })


# class heatmap between IBD and Non-IBD ----------------------------------------
library(circlize)
library(scales)
library(ComplexHeatmap)
library(tidymass)

load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')


# calculate average values for manual classes and patient group ----------------
# 1. scale each metabolite to z-score
# 2. sum the all metabolite within the same class
# 3. use the average metabolite class value to represent the patient group, then calculate the z-score again

# calculate z-score for each metabolite
temp_object <- object_stat %>% scale_data(method = 'auto')

# for each sample, sum the same class metabolites (that are significant)
temp_object <- temp_object %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(p_value_adjust <= 0.05)

# top 25 classes maxtrix
top_classes <- temp_object@annotation_table %>% count(metabolon_subclass) %>% arrange(desc(n)) %>% slice(1:25) %>% as.data.frame()
variable_id_class <- lapply(top_classes$metabolon_subclass, function(x){
  temp_object@annotation_table %>% 
    dplyr::filter(metabolon_subclass == x) %>% 
    dplyr::pull(variable_id)
})
names(variable_id_class) <- top_classes$metabolon_subclass


matrix_class <- lapply(seq_along(top_classes$metabolon_subclass), function(i){
  temp_variable <- variable_id_class[[i]]
  temp_object %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(variable_id %in% temp_variable) %>% 
    activate_mass_dataset(what = 'expression_data') %>%
    colSums()
}) %>% do.call(rbind, .)

rownames(matrix_class) <- top_classes$metabolon_subclass

# calculate average value for each patient group ----------------------------

sample_id_ibd <- temp_object@sample_info %>% filter(phenotype_group1 == 'CD') %>% pull(sample_id)
sample_id_nonibd <- temp_object@sample_info %>% filter(phenotype_group1 == 'non_IBD') %>% pull(sample_id)

idx_ibd <- match(sample_id_ibd, colnames(matrix_class))
idx_nonibd <- match(sample_id_nonibd, colnames(matrix_class))

matrix_class_group <- apply(matrix_class, 1, function(x){
  mean_nonibd <- x[idx_nonibd] %>% mean()
  mean_ibd <- x[idx_ibd] %>% mean()
  c(mean_nonibd, mean_ibd)
}) %>% t()

colnames(matrix_class_group) <- c('NonIBD', 'CD')
matrix_class_group2 <- scale(matrix_class_group)


rm(idx_ibd, idx_nonibd, sample_id_ibd, sample_id_nonibd, matrix_class_group)

# for the side annotation, we can use the number of metabolites in each class, and use different size of points to represent the number of metabolites in each class

sum_counts <- temp_object@annotation_table %>% count(metabolon_subclass) %>% arrange(desc(n)) %>% slice(1:25) %>% pull(n)

direction_counts <- lapply(seq_along(variable_id_class), function(i){
  temp_variable_id <- variable_id_class[[i]]
  temp_label <- temp_object %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(variable_id %in% temp_variable_id) %>%
    extract_variable_info() %>%
    mutate(increase_label = case_when(fc > 1 ~ 'increase',
                                      fc < 1 ~ 'decrease',
                                      TRUE ~ 'no_change'))
  
  increase_count <- temp_label %>% filter(increase_label == 'increase') %>% nrow()
  decrease_count <- temp_label %>% filter(increase_label == 'decrease') %>% nrow()
  
  result <- data.frame(increase_count = increase_count,
                       decrease_count = decrease_count, 
                       stringsAsFactors = FALSE)
  
  return(result)
}) %>% bind_rows()

increase_counts <- direction_counts$increase_count
decrease_counts <- direction_counts$decrease_count

rm(direction_counts)

# chemical class enrichment analysis function -------------------------------------------------------

met_table <- object_stat %>%
  extract_variable_info() %>%
  select(variable_id, Compound.name, p_value_adjust, fc, metabolon_subclass) %>%
  rename(pvalue = p_value_adjust, 
         foldchange = fc,
         class = metabolon_subclass) %>% 
  arrange(pvalue) %>% 
  mutate(order = seq(n()))

enrichment_result <- enrichClass2(table_compound = met_table, 
                                  path = '.', 
                                  class_cutoff = 2, 
                                  p_cutoff = 0.2)


temp_idx <- match(top_classes$metabolon_subclass, enrichment_result$name)
pvals <- enrichment_result$adjustedpvalue[temp_idx]
names(pvals) <- enrichment_result$name[temp_idx]
stars <- ifelse(pvals < 0.001, "***",
                ifelse(pvals < 0.01, "**",
                       ifelse(pvals < 0.05, "*", "")))

rm(temp_idx, enrichClass2, met_table);gc()

# Construct custom plotting functions (bar plot & circles) ---------------------

# 1. custom_diverging_bar: a custom function to draw diverging bar plot for increase and decrease counts
custom_diverging_bar <- AnnotationFunction(
  fun = function(index) {  
    n <- length(index)
    inc <- increase_counts[index]
    dec <- decrease_counts[index]
    
    pushViewport(viewport(xscale = c(-15, 15), yscale = c(0, 1)))
    
    y_pos <- (n - seq_len(n) + 0.5) / n
    bar_height <- 0.8 / n
    
    grid.rect(x = 0, y = y_pos, width = dec, height = bar_height,
              just = c("right", "center"), default.units = "native",
              gp = gpar(fill = "#7fb1d3", col = NA))
    
    grid.rect(x = 0, y = y_pos, width = inc, height = bar_height,
              just = c("left", "center"), default.units = "native",
              gp = gpar(fill =  "#fc8070", col = NA))
    
    grid.lines(x = 0, y = c(0, 1), default.units = "native", 
               gp = gpar(col = "black", lty = 2, lwd = 1.5))
    
    popViewport()
  },
  which = "row",        
  width = unit(5, "cm") 
)


# 2. custom_circles: a custom function to draw circles with size representing the number of metabolites in each class

custom_circles <- AnnotationFunction(
  fun = function(index) {
    n <- length(index)
    sizes <- sum_counts[index]
    max_s <- max(sum_counts)
    
    y_pos <- (n - seq_len(n) + 0.5) / n
    
    grid.points(x = rep(0.5, n), y = y_pos, pch = 16,
                size = unit(sizes / max_s * 6, "mm"),
                gp = gpar(col = "#b5de68"))
  },
  which = "row",
  width = unit(2, "cm")
)



# complex heatmap of the chemical class ----------------------------------------


# col_fun <- colorRamp2(c(-0.5, 0, 0.5), c('#1982c4', "white", "#ff595e"), space = "RGB")

col_fun <- colorRamp2(c(-1, 0, 1),
                      c(
                        viridis::viridis(n = 3)[1],
                        viridis::viridis(n = 3)[2],
                        viridis::viridis(n = 3)[3]
                      ))

patient_group  <-  HeatmapAnnotation(patient_group = colnames(matrix_class_group2), 
                                     col = list(
                                       patient_group = c('NonIBD' = '#d0ccd0', 
                                                         'CD' = '#274156')
                                     ), 
                                     annotation_label = 'Patient',
                                     show_legend = FALSE)

right_anno <- rowAnnotation(
  BarChart = custom_diverging_bar,
  SumCircles = custom_circles, 
  Pvalue = anno_text(
    stars,
    gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
    just = "center", location = 0.5, width = unit(1, "cm")
  ),
  annotation_name_side = "top" 
)

met_name <- rownames(matrix_class_group2)

temp_heatmap <- Heatmap(matrix_class_group2, 
                        name = "Z-score", 
                        col = col_fun,
                        cluster_columns = FALSE,
                        cluster_rows = TRUE,
                        clustering_method_rows = 'ward.D',
                        # clustering_method_rows = 'complete',
                        top_annotation = patient_group,
                        row_labels = met_name, 
                        right_annotation = right_anno,
                        # border = TRUE,  
                        rect_gp = gpar(col = "white", lwd = 1),
                        row_names_gp = gpar(fontsize = 7),
                        row_names_side = 'left',
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        heatmap_legend_param = list(direction = "vertical"), 
                        # show_row_names = FALSE, 
                        show_column_names = FALSE)

lgd_bubble <- Legend(
  labels = c("5", "10", "15"), title = "Sum of Counts",
  type = "points", pch = 16,
  size = unit(c(5, 10, 15) / max(sum_counts) * 6, "mm"),
  legend_gp = gpar(col ="#b5de68")
)

lgd_bar <- Legend(
  labels = c("Increase", "Decrease"), title = "Regulation",
  legend_gp = gpar(fill = c("#fc8070",  "#7fb1d3"))
)

grid.newpage() 
temp_plot <- draw(temp_heatmap, 
                  annotation_legend_list = list(lgd_bar, lgd_bubble), 
                  ht_gap = unit(5, "mm"),
                  padding = unit(c(2, 2, 12, 2), "mm"))


# plot the diverging bar plot first to ensure the correct order of rows, then add decorations to the bar plot
decorate_annotation("BarChart", {
  grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  pushViewport(viewport(xscale = c(-15, 15), yscale = c(0, 1), clip = "off"))
  grid.xaxis(at = c(-15, -10, -5, 0, 5, 10, 15), 
             label = c("15", "10", '5', "0", '5', '10', "15"),
             gp = gpar(fontsize = 7))
  popViewport()
})

# Plot the heatmap first to ensure the correct order of rows, then add decorations to circles and P-value annotations
decorate_annotation("SumCircles", { grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1)) })
decorate_annotation("Pvalue", { grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1)) })

# temp_data <- matrix_class_group2 %>% as.data.frame()
# 
# readr::write_csv(temp_data, file = '~/Project/00_IBD_project/Data/20260718_source_data/extended_fig3b_heatmap_class_zscore_260720.csv')



################################################################################
# Differential abundance of microbial metabolites between CD and non-IBD ---------------------------------
setwd('~/Project/00_IBD_project/Data/20251107_boxplot_microbial_met/')

library(tidyverse)
library(tidymass)
library(sjmisc)

load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')

met_info_manual <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20251104_correlation_metabolite_with_microbial_dysbiosis_inex/annotation_table_251104.xlsx')

microbial_met <- met_info_manual %>% 
  dplyr::select(variable_id, `Manual label`) %>% 
  rename(manual_label = `Manual label`)

# calculate the CV, number of l
met_data <- object_stat %>% 
  extract_expression_data()

met_value <- object_stat %>% 
  extract_variable_info()

sample_info <- object_stat %>% 
  extract_sample_info()

idx_CD <- which(sample_info$phenotype_group1 == 'CD')
idx_non_IBD <- which(sample_info$phenotype_group1 == 'non_IBD')

# calculate CV for metabolite according to group
cv_CD <- apply(met_data[, idx_CD], 1, function(x) sd(x) / mean(x))
cv_non_IBD <- apply(met_data[, idx_non_IBD], 1, function(x) sd(x) / mean(x))

# count the sample number with value <= 1000 for each metabolite according to group
num_low_CD <- apply(met_data[, idx_CD], 1, function(x) sum(x <= 1000))
num_low_non_IBD <- apply(met_data[, idx_non_IBD], 1, function(x) sum(x <= 1000))

# calculate the low value ratio
ratio_low_CD <- num_low_CD / length(idx_CD)
ratio_low_non_IBD <- num_low_non_IBD / length(idx_non_IBD)

met_stat <- met_value %>% 
  dplyr::mutate(cv_CD = cv_CD[variable_id],
                cv_non_IBD = cv_non_IBD[variable_id],
                num_low_CD = num_low_CD[variable_id],
                num_low_non_IBD = num_low_non_IBD[variable_id],
                ratio_low_CD = ratio_low_CD[variable_id],
                ratio_low_non_IBD = ratio_low_non_IBD[variable_id]) %>% 
  select(variable_id:Compound.name, confidence_level, hmdb_id, cv_CD:ratio_low_non_IBD) %>% 
  rename(mz = mz.x,
         rt = rt.x) %>% 
  left_join(microbial_met, by = 'variable_id')

# plot the boxplot for selected microbial metabolites
selected_metabolites <- microbial_met %>% dplyr::filter(!is.na(manual_label)) %>% pull(variable_id)


# convert to microbial data ----------------------------------------------------
library(sjmisc)

# extract the data for selected metabolites
object_stat_microbial_met <- object_stat %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(variable_id %in% selected_metabolites)

sample_info <- object_stat_microbial_met %>% 
  extract_sample_info()

idx_CD <- which(sample_info$phenotype_group1 == 'CD')
idx_non_IBD <- which(sample_info$phenotype_group1 == 'non_IBD')

# convert to z-score
z_score_microbial_met <- object_stat_microbial_met %>% 
  extract_expression_data() %>% 
  apply(1, function(x) (x - mean(x)) / sd(x)) %>% 
  as.data.frame() %>% 
  rotate_df()

z_score_mean_microbial_met <- z_score_microbial_met %>% 
  dplyr::mutate(mean_CD = rowMeans(select(., all_of(colnames(z_score_microbial_met)[idx_CD]))),
                mean_non_IBD = rowMeans(select(., all_of(colnames(z_score_microbial_met)[idx_non_IBD])))) %>%
  rownames_to_column(var = 'variable_id') %>%
  select(variable_id, mean_CD, mean_non_IBD)


fold_change_microbial_met <- object_stat_microbial_met %>% 
  extract_expression_data() %>% 
  apply(1, function(x) mean(x[idx_CD]) / mean(x[idx_non_IBD])) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'variable_id') %>%
  rename(fold_change = '.')


microbial_met_stat <- z_score_mean_microbial_met %>% 
  left_join(fold_change_microbial_met, by = 'variable_id') %>% 
  left_join(met_stat, by = 'variable_id') %>% 
  arrange(desc(abs(mean_CD - mean_non_IBD))) %>% 
  select(variable_id:p_value_adjust, Compound.name, confidence_level, hmdb_id, manual_label)

temp <- microbial_met_stat %>% 
  select(variable_id:mean_non_IBD) %>% 
  pivot_longer(cols = mean_CD:mean_non_IBD, names_to = 'group', values_to = 'z_score_mean') %>% 
  mutate(scale_z_score_mean = as.numeric(scale(z_score_mean))) %>% 
  select(-z_score_mean) %>%
  pivot_wider(names_from = group, values_from = scale_z_score_mean) %>% 
  rename(scaled_z_score_mean_CD = mean_CD, scaled_z_score_mean_non_IBD = mean_non_IBD)

microbial_met_stat <- microbial_met_stat %>% 
  left_join(temp, by = 'variable_id') %>% 
  select(variable_id:mean_non_IBD, scaled_z_score_mean_CD, scaled_z_score_mean_non_IBD, everything()) %>% 
  arrange(match(manual_label, c('Trp_met', 'Bile_acids', 'Phenyls/Phenols', 'Other_microMet')))

dir.create('~/Project/00_IBD_project/Data/20260114_microbial_met_summary', recursive = TRUE, showWarnings = FALSE)
save(microbial_met_stat, file = '~/Project/00_IBD_project/Data/20260114_microbial_met_summary/microbial_met_stat_260114.RData')

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE))
readr::write_csv(microbial_met_stat, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig4_microbial_met_stat_260720.csv')

# plot 4 layer circular heatmap with significance and category ----------------------

library(tidyverse)
library(circlize)
library(ComplexHeatmap)

# 1. Data Preparation 

met_matrix <- microbial_met_stat %>% 
  select(variable_id, scaled_z_score_mean_non_IBD, scaled_z_score_mean_CD) %>% 
  column_to_rownames(var = 'variable_id')
rownames(met_matrix) <- microbial_met_stat$Compound.name

split <- factor(microbial_met_stat$manual_label, 
                levels = c('Trp_met', 'Bile_acids', 'Phenyls/Phenols', 'Other_microMet'))

category_colors <- structure(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), 
                             names = levels(split))
col_fun1 = colorRamp2(c(-2, 0, 2), c("dodgerblue", "white", "tomato"))


# 2. Plotting Starts
circos.clear()

# Set track margins
circos.par(track.margin = c(0.01, 0.01)) 

# === Track 1: Heatmap ===
circos.heatmap(met_matrix, 
               split = split, 
               col = col_fun1, 
               cell.border = 'black',
               track.height = 0.15, 
               dend.side = "outside")

fc <- log2(microbial_met_stat$fold_change)

# === Track 2: Fold Change ===
circos.track(ylim = range(fc), track.height = 0.1, panel.fun = function(x, y) {
  y = fc[CELL_META$subset][CELL_META$row_order]
  circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
  circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue"), cex = 0.6)
}, cell.padding = c(0.01, 0, 0.01, 0))

# === Track 3: Category + Significance (Merged Layer) ===
# Slightly increase the track height (e.g., from 0.05 to 0.06) to accommodate asterisks
circos.track(ylim = c(0, 1), track.height = 0.06, bg.border = NA, panel.fun = function(x, y) {
  
  # 1. Prepare data
  # Category data
  current_labels <- microbial_met_stat$manual_label[CELL_META$subset][CELL_META$row_order]
  current_cols <- category_colors[as.character(current_labels)]
  # P-value data
  pvals = microbial_met_stat$p_value_adjust[CELL_META$subset][CELL_META$row_order]
  sig_text = ifelse(pvals < 0.001, '***', ifelse(pvals < 0.01, "**", ifelse(pvals < 0.05, "*", "")))
  
  n = length(current_labels)
  
  # 2. Draw color blocks first (Background)
  circos.rect(xleft = 0:(n-1), ybottom = rep(0, n), xright = 1:n, ytop = rep(1, n), 
              col = current_cols, border = NA) 
  
  # 3. Draw asterisks (Foreground)
  # Note: To make asterisks visible on dark backgrounds, use col = "white" or "black"
  circos.text(x = seq_along(pvals) - 0.5, 
              y = rep(0.5, length(pvals)), # Display in center
              labels = sig_text, 
              facing = "downward", 
              adj = c(0.5, 0.6), # Fine-tune vertical alignment to center asterisks in the block
              cex = 0.8, 
              col = "white") # White asterisks are usually clearer on category color blocks
})

# === Track 4: Labels (Innermost Text) ===

# Adjust margins: Leave space at the top to separate text from the color blocks above
# (Note: You might want to change the second value to 0.05 if you want a larger gap)
circos.par(track.margin = c(0.01, 0.01)) 

circos.track(ylim = c(0, 1), track.height = 0.1, bg.border = NA, panel.fun = function(x, y) {
  labels = rownames(met_matrix)[CELL_META$subset][CELL_META$row_order]
  
  circos.text(x = seq_along(labels) - 0.5, 
              y = rep(1, length(labels)),  
              labels = labels,
              facing = "clockwise", 
              niceFacing = TRUE, 
              adj = c(1, 0.5), 
              cex = 0.5)
}, cell.padding = c(0, 0, 0, 0))

circos.clear()

# --- Draw Legends ---

# 1. Create legend objects
lgd_cat = Legend(title = "Class", 
                 at = names(category_colors), 
                 type = "grid", 
                 legend_gp = gpar(fill = category_colors))

lgd_hm = Legend(title = "Z-score", 
                col_fun = col_fun1)

# 2. Pack legends
lgd_list = packLegend(lgd_hm, lgd_cat)

# 3. Draw to the top-right corner
# x = unit(1, "npc"): Horizontal position at the far right of the canvas (100%)
# y = unit(1, "npc"): Vertical position at the very top of the canvas (100%)
# just = c("right", "top"): Align using the "top-right" corner of the legend itself as the anchor
# gap = unit(2, "mm"): If it's too close to the edge, subtract a small offset

# Move 5mm down-left to leave some margin
draw(lgd_list, 
     x = unit(1, "npc") - unit(2, "mm"), 
     y = unit(1, "npc") - unit(2, "mm"), 
     just = c("right", "top"))



################################################################################
# eGFR adjusted indole-sulfate comparsion between CD and non-IBD ---------------
dir.create('~/Project/00_IBD_project/Data/20251106_IS_adjustment', showWarnings = F, recursive = T)

setwd('~/Project/00_IBD_project/Data/20251106_IS_adjustment')

library(tidyverse)
library(tidymass)
library(IbdData)
library(sjmisc)


# eGFR function ----------------------------------------------------------------

#' Calculate estimated Glomerular Filtration Rate (eGFR) using the CKD-EPI equation


calculate_eGFR_2021 <- function(con_creatinine, age, gender) {
  if (any(c(is.na(con_creatinine), is.na(gender), is.na(age)))) {
    return(NA)
  }
  k_factor <- ifelse(gender == 'Male', 0.9, 0.7)
  alpha_factor <- ifelse(gender == 'Male', -0.302, -0.241)
  min_conc <- min(con_creatinine/k_factor, 1)
  max_conc <- max(con_creatinine/k_factor, 1)
  adj_factor <- ifelse(gender == 'Male', 1, 1.012)
  
  result <- 142*min_conc^alpha_factor*max_conc^-1.200*0.9938^age*adj_factor
  return(result)
}


# calculate eGFR for all samples ----------------------------------------------

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/20251106_IS_adjustment/patient_height_251106.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_lab_test_241002.RData')

sample_table_info <- object_final %>% extract_sample_info()

height_table <- patient_height %>% 
  select(sample_id, height_cm)
creatine_table <- patient_meta_info_lab_test %>% 
  select(sample_id, creatinine)

sample_table_info %>% 
  left_join(height_table, 
            by = c('sample_id')) %>% 
  left_join(creatine_table, 
            by = c('sample_id'))



table_eGFR <- sample_table_info %>% 
  left_join(height_table, 
            by = c('sample_id')) %>% 
  left_join(creatine_table, 
            by = c('sample_id')) %>% 
  select(sample_id, patient_id, visit_encounter_id, gender, age, race, creatinine, height_cm, phenotype_group1, visit_number)


# calculate_eGFR(con_creatinine = 0.6, age = 12, race = 'Caucasian', gender = 'Male')

eGFR_values <- sapply(seq_along(table_eGFR$sample_id), function(i){
  cat(i, ' ')
  temp_creatinine <- as.numeric(table_eGFR$creatinine[i])
  temp_age <- table_eGFR$age[i]
  temp_race <- table_eGFR$race[i]
  temp_gender <- table_eGFR$gender[i]
  
  result <- calculate_eGFR_2021(con_creatinine = temp_creatinine,
                                age = temp_age,
                                gender = temp_gender)
  
  return(result)
})

table_eGFR <- table_eGFR %>% mutate(eGFR = eGFR_values) 
save(table_eGFR, 
     file = '~/Project/00_IBD_project/Data/20251106_IS_adjustment/table_eGFR_260118.RData')



# adjust the IS values by eGFR -------------------------------------------------

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/20251106_IS_adjustment/table_eGFR_260118.RData')

IS_concentration <- object_final %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(variable_id == 'M212T285_c18_neg') %>% 
  extract_expression_data() %>% 
  as.data.frame() %>%
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id') %>% 
  as_tibble()

table_eGFR <- table_eGFR %>% 
  left_join(IS_concentration, 
            by = c('sample_id'))

table_eGFR_plot <- table_eGFR %>% 
  dplyr::filter(visit_number == 'enrollment') %>% 
  dplyr::filter(!is.na(M212T285_c18_neg) & !is.na(eGFR)) %>%
  dplyr::mutate(zscore_IS = scale(M212T285_c18_neg))

# adjust IS by eGFR using linear regression
temp_data <- table_eGFR_plot %>% 
  dplyr::filter(!is.na(M212T285_c18_neg) & !is.na(eGFR))

lm_model <- lm(M212T285_c18_neg ~ eGFR, data = temp_data)
temp_data <- temp_data %>% 
  mutate(adjusted_IS = lm_model$residuals) %>% 
  mutate(zscore_adjusted_IS = scale(adjusted_IS))

temp_data <- temp_data %>% 
  select(sample_id, phenotype_group1, zscore_IS, zscore_adjusted_IS) %>% 
  tidyr::pivot_longer(cols = c(zscore_IS, zscore_adjusted_IS), 
                      names_to = 'adjust_type', 
                      values_to = 'zscore_value') %>%
  mutate(adjust_type = recode_factor(adjust_type,
                                     'zscore_IS' = 'Before eGFR Adjustment',
                                     'zscore_adjusted_IS' = 'After eGFR Adjustment'))

temp_data$zscore_value <- as.numeric(temp_data$zscore_value)

temp_plot <- temp_data %>% 
  ggplot(aes(x = adjust_type, y = zscore_value, colour = phenotype_group1)) +
  geom_quasirandom(
    method = "quasirandom",
    alpha = 0.5, 
    size = 2,
    width = 0.15,
    dodge.width = 0.8
  ) +
  stat_summary( # 2. plot quartile range
    aes(group = phenotype_group1),
    fun.min = function(x) quantile(x, 0.25),
    fun.max = function(x) quantile(x, 0.75),
    geom = "errorbar", # errorbar
    width = 0.1,       # control the width of the error bars
    color = "black",
    linewidth = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  # 3. median bar
  stat_summary(
    aes(group = phenotype_group1),
    fun = median,
    geom = "crossbar",
    width = 0.2,
    color = "black",
    linewidth = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  stat_compare_means(
    aes(group = phenotype_group1),
    # comparisons = list(c("CD", "non_IBD")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(temp_data$zscore_value, na.rm = TRUE) + 0.5
  ) +
  # stat_compare_means(comparisons = list(c("CD", "non_IBD")),
  #                    method = "wilcox.test",
  #                    label = "p.format") +
  scale_colour_manual(values = c('CD' = 'tomato',
                                 'non_IBD' = 'dodgerblue')) +
  xlab('Phenotype Group') +
  ylab('Z-score') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.8, 0.8))


ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20251106_IS_adjustment/IS_phenotype_compare_before_after_correction_260115.pdf', 
       width = 6, height = 6)



dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)
readr::write_csv(temp_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig5_IS_phenotype_compare_before_after_correction_260720.csv')



################################################################################
# Compare the IS value with Albumin --------------------------------------------

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_lab_test_241002.RData')

value_albumin <- patient_meta_info_lab_test %>% 
  select(sample_id, albumin) %>% 
  mutate(albumin = as.numeric(albumin))


IS_concentration <- object_final %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(variable_id == 'M212T285_c18_neg') %>% 
  extract_expression_data() %>% 
  as.data.frame() %>%
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id') %>% 
  as_tibble()


sample_info_table <- object_final %>% 
  extract_sample_info() %>% 
  select(sample_id:phenotype_group1)


table_IS_albumin <- sample_info_table %>%
  left_join(IS_concentration, 
            by = c('sample_id')) %>% 
  left_join(value_albumin, 
            by = c('sample_id'))


# add quartiles of albumin
temp_IS_albumin <- table_IS_albumin %>%
  filter(!is.na(albumin)) %>% 
  as_tibble() %>% 
  mutate(zscore_IS = as.numeric(scale(M212T285_c18_neg)))


# enrollment visit - CD --------------------------------------------------------
# add t-test between albumin quartiles

library(ggpubr)
IS_albumin_CD_enrollment <- temp_IS_albumin %>%
  filter(visit_number == 'enrollment' & phenotype_group1 == 'CD') %>%
  mutate(albumin_quartile = ntile(albumin, 4)) %>% 
  mutate(albumin_quartile = factor(albumin_quartile,
                                   levels = c(1, 2, 3, 4),
                                   labels = c('Q1', 'Q2', 'Q3', 'Q4')))


temp_plot <- ggplot(IS_albumin_CD_enrollment, aes(x = as.factor(albumin_quartile), y = zscore_IS, colour = albumin_quartile)) +
  # ggplot(aes(x = adjust_type, y = zscore_value, colour = phenotype_group1)) +
  geom_quasirandom(
    method = "quasirandom",
    alpha = 0.5, 
    size = 2,
    width = 0.15,
    dodge.width = 0.8
  ) +
  stat_summary( # 2. plot quartile range
    aes(group = phenotype_group1),
    fun.min = function(x) quantile(x, 0.25),
    fun.max = function(x) quantile(x, 0.75),
    geom = "errorbar", # errorbar
    width = 0.1,       # control the width of the error bars
    color = "black",
    linewidth = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  # 3. median bar
  stat_summary(
    aes(group = phenotype_group1),
    fun = median,
    geom = "crossbar",
    width = 0.2,
    color = "black",
    linewidth = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  stat_compare_means(
    # aes(group = albumin_quartile),
    comparisons = list(c("Q4", "Q1"), c("Q4", "Q2"), c("Q4", "Q3")),
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.9
  ) +
  xlab('Albumin (Quartile)') +
  ylab('3-Indoxyl sulfate (Z-score)') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.1, 0.9))


dir.create('~/Project/00_IBD_project/Data/20260116_albumin_adjusted_indoxyl_sulfate',
           recursive = TRUE,
           showWarnings = FALSE)
ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260116_albumin_adjusted_indoxyl_sulfate/enrollment_CD_albumin_3_indoxyl_sulfate_260116.pdf',
       width = 6, height = 6)


# enrollment visit - non-IBD ---------------------------------------------------

IS_albumin_non_IBD_enrollment <- temp_IS_albumin %>%
  filter(visit_number == 'enrollment' & phenotype_group1 == 'non_IBD') %>%
  mutate(albumin_quartile = ntile(albumin, 4)) %>% 
  mutate(albumin_quartile = factor(albumin_quartile,
                                   levels = c(1, 2, 3, 4),
                                   labels = c('Q1', 'Q2', 'Q3', 'Q4')))


temp_plot <- ggplot(IS_albumin_non_IBD_enrollment, aes(x = as.factor(albumin_quartile), y = zscore_IS, colour = albumin_quartile)) +
  # ggplot(aes(x = adjust_type, y = zscore_value, colour = phenotype_group1)) +
  geom_quasirandom(
    method = "quasirandom",
    alpha = 0.5, 
    size = 2,
    width = 0.15,
    dodge.width = 0.8
  ) +
  stat_summary( # 2. plot quartile range
    aes(group = phenotype_group1),
    fun.min = function(x) quantile(x, 0.25),
    fun.max = function(x) quantile(x, 0.75),
    geom = "errorbar", # errorbar
    width = 0.1,       # control the width of the error bars
    color = "black",
    linewidth = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  # 3. median bar
  stat_summary(
    aes(group = phenotype_group1),
    fun = median,
    geom = "crossbar",
    width = 0.2,
    color = "black",
    linewidth = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  stat_compare_means(
    # aes(group = albumin_quartile),
    comparisons = list(c("Q4", "Q1"), c("Q4", "Q2"), c("Q4", "Q3")),
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.9
  ) +
  xlab('Albumin (Quartile)') +
  ylab('3-Indoxyl sulfate (Z-score)') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.1, 0.9))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260116_albumin_adjusted_indoxyl_sulfate/enrollment_non_IBD_albumin_3_indoxyl_sulfate_260116.pdf',
       width = 6, height = 6)

temp_data <- IS_albumin_non_IBD_enrollment %>% 
  dplyr::select(sample_id, albumin, zscore_IS, albumin_quartile) %>% 
  dplyr::mutate(label = 'non_IBD_enrollment') %>% 
  bind_rows(IS_albumin_CD_enrollment %>% 
              dplyr::select(sample_id, albumin, zscore_IS, albumin_quartile) %>% 
              dplyr::mutate(label = 'CD_enrollment'))


readr::write_csv(temp_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig5_IS_albumin_correlation_260720.csv')



################################################################################
# Correlation analysis between hemoglobin and its related metabolites ----------
dir.create('~/Project/00_IBD_project/Data/20260401_correlation_hb_hemo_met', showWarnings = FALSE, recursive = TRUE)
setwd('~/Project/00_IBD_project/Data/20260401_correlation_hb_hemo_met/')

library(tidyverse)
library(tidymass)
library(sjmisc)

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_lab_test_241002.RData')


sample_info <- object_stat %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  extract_sample_info() %>% 
  select(sample_id, patient_id, phenotype_group1)


# hemoglobin
temp_hb <- patient_meta_info_lab_test %>% 
  select(sample_id, patient_id, hemoglobin) %>%
  filter(!is.na(hemoglobin)) %>% 
  mutate(hemoglobin = as.numeric(hemoglobin)) %>%
  mutate(relative_abundance = as.numeric(scale(hemoglobin))) %>% 
  select(sample_id, relative_abundance) %>%
  rename(relative_abundance_hb = relative_abundance)


# Bilirubin, M585T45_hilic_pos
temp_data_bilirubin <- object_stat %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(variable_id == 'M585T45_hilic_pos') %>% 
  extract_expression_data() %>% 
  rownames_to_column('variable_id') %>%
  pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'abundance') %>% 
  mutate(relative_abundance = as.numeric(scale(abundance))) %>% 
  select(sample_id, relative_abundance) %>% 
  rename(relative_abundance_bilirubin = relative_abundance)


# Biliverdin, M583T482_c18_pos
temp_data_biliverdin <- object_stat %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(variable_id == 'M583T482_c18_pos') %>% 
  extract_expression_data() %>% 
  rownames_to_column('variable_id') %>%
  pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'abundance') %>% 
  mutate(relative_abundance = as.numeric(scale(abundance))) %>% 
  select(sample_id, relative_abundance) %>% 
  rename(relative_abundance_biliverdin = relative_abundance)


# Urobilin, M591T385_c18_pos
temp_data_urobilin <- object_stat %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(variable_id == 'M591T385_c18_pos') %>% 
  extract_expression_data() %>% 
  rownames_to_column('variable_id') %>%
  pivot_longer(cols = -variable_id, names_to = 'sample_id', values_to = 'abundance') %>% 
  mutate(relative_abundance = as.numeric(scale(abundance))) %>% 
  select(sample_id, relative_abundance) %>% 
  rename(relative_abundance_urobilin = relative_abundance)

# merge data
temp_data_merge <- sample_info %>%
  left_join(temp_hb, by = 'sample_id') %>%
  left_join(temp_data_bilirubin, by = 'sample_id') %>% 
  left_join(temp_data_biliverdin, by = 'sample_id') %>% 
  left_join(temp_data_urobilin, by = 'sample_id') %>% 
  filter(!is.na(relative_abundance_hb))


# plot the linear correlation between hb and each metabolite
# CD samples, x: hb, y: each metabolite, each metabolite assign one color
temp_data_merge_cd <- temp_data_merge %>% 
  filter(phenotype_group1 == 'CD')

# non-IBD samples, x: hb, y: each metabolite, each metabolite assign one color
temp_data_merge_nonibd <- temp_data_merge %>% 
  filter(phenotype_group1 == 'non_IBD')

# print the correlation results in a table, include 3 groups: all samples, CD samples, non-IBD samples; include correlation coefficient and p value
# column: group, bilirubin_cor, bilirubin_p, biliverdin_cor, biliverdin_p, urobilin_cor, urobilin_p

cor_results <- data.frame(
  group = c('all_samples', 'CD_samples', 'nonIBD_samples'),
  bilirubin_cor = c(cor(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_bilirubin, method = 'pearson'),
                    cor(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_bilirubin, method = 'pearson'),
                    cor(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_bilirubin, method = 'pearson')),
  bilirubin_p = c(cor.test(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_bilirubin, method = 'pearson')$p.value,
                  cor.test(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_bilirubin, method = 'pearson')$p.value,
                  cor.test(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_bilirubin, method = 'pearson')$p.value),
  biliverdin_cor = c(cor(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_biliverdin, method = 'pearson'),
                     cor(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_biliverdin, method = 'pearson'),
                     cor(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_biliverdin, method = 'pearson')),
  biliverdin_p = c(cor.test(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_biliverdin, method = 'pearson')$p.value,
                   cor.test(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_biliverdin, method = 'pearson')$p.value,
                   cor.test(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_biliverdin, method = 'pearson')$p.value),
  urobilin_cor = c(cor(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_urobilin, method = 'pearson'),
                   cor(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_urobilin, method = 'pearson'),
                   cor(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_urobilin, method = 'pearson')),
  urobilin_p = c(cor.test(temp_data_merge$relative_abundance_hb, temp_data_merge$relative_abundance_urobilin, method = 'pearson')$p.value,
                 cor.test(temp_data_merge_cd$relative_abundance_hb, temp_data_merge_cd$relative_abundance_urobilin, method = 'pearson')$p.value,
                 cor.test(temp_data_merge_nonibd$relative_abundance_hb, temp_data_merge_nonibd$relative_abundance_urobilin, method = 'pearson')$p.value),
  stringsAsFactors = FALSE
)



# plot the linear correlation between hb and each metabolite

ZZWtool::ZZWcolors()

# all samples, x: hb, y: each metabolite
# phenotype_group1 as color
# add 3 lines, each line represents one group, CD, non-IBD, all samples
# add p values, R square values to the plot

temp_color <- c('CD' = "#fc8070", 'non_IBD' =  "#7fb1d3")


walk(c('bilirubin', 'biliverdin', 'urobilin'), function(x){
  cat('Processing correlation between hb and ', x, '...\n')
  
  cor_value <- cor(temp_data_merge[[paste0('relative_abundance_hb')]], temp_data_merge[[paste0('relative_abundance_', x)]], method = 'pearson')
  p_value <- cor.test(temp_data_merge[[paste0('relative_abundance_hb')]], temp_data_merge[[paste0('relative_abundance_', x)]], method = 'pearson')$p.value
  r_square <- cor_value^2
  
  p <- ggplot(temp_data_merge, 
              aes_string(x = "relative_abundance_hb", 
                         y = paste0("relative_abundance_", x))) +
    geom_point(aes(color = phenotype_group1), size = 2, shape = 16, 
               alpha = 0.8) +
    geom_smooth(method = "lm", data = temp_data_merge, color = "black") +
    # geom_smooth(method = "lm", 
    #             data = temp_data_merge %>% 
    #               filter(phenotype_group1 == "CD"), 
    #             color = "red") +
    # geom_smooth(method = "lm", 
    #             data = temp_data_merge %>% 
    #               filter(phenotype_group1 == "non_IBD"), 
    #             color = "blue") +
    scale_color_manual(values = temp_color) +
    labs(x = "Relative Abundance of Hemoglobin",
         y = paste0("Relative Abundance of ", str_to_title(x))) +
    ZZWtool::ZZWTheme() +
    annotate("text", 
             x = min(temp_data_merge$relative_abundance_hb), 
             y = max(temp_data_merge[[paste0("relative_abundance_", x)]]), 
             label = paste0("R = ", round(cor_value, 2), "\n", 
                            "R² = ", round(r_square, 2), "\n", 
                            "p = ", signif(p_value, 2)), 
             hjust = 0, vjust = 1, size = 4) +
    theme(legend.position = "none", 
          axis.text.y = element_text(angle = 0))
  
  
  ggsave(file.path('~/Project/00_IBD_project/Data/20260401_correlation_hb_hemo_met/', paste0('all_enrolled_samples_correlation_hb_', x, '_with_group_260401.pdf')), 
         plot = p, width = 6, height = 6)
  
  print(p)
  
  
})


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)
readr::write_csv(temp_data_merge, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig5_hb_hemo_met_correlation_260720.csv')
