


################################################################################
# serology marker between CD vs non-IBD ----------------------------------------

dir.create('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model', recursive = FALSE, showWarnings = TRUE)
setwd('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model')

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_lab_test_241002.RData')


temp_info <- object_final@sample_info %>% 
  select(sample_id, phenotype_group1:severity_class)

patient_meta_info_lab_test_enrollment <- patient_meta_info_lab_test %>% 
  filter(sample_id != 'B032_S34') %>% 
  filter(visit_number == 'enrollment') %>% 
  select(sample_id, patient_id, visit_encounter_id, visit_number, gm_csf:anca) %>%
  mutate(gm_csf = as.numeric(gm_csf), 
         iga_asca = as.numeric(iga_asca),
         igg_asca = as.numeric(igg_asca), 
         i2 = as.numeric(i2), 
         ompc = as.numeric(ompc), 
         cbir_fla = as.numeric(cbir_fla),
         anca = as.numeric(anca)) %>% 
  mutate(gm_csf = scale(gm_csf)[,1], 
         iga_asca = scale(iga_asca)[,1], 
         igg_asca = scale(igg_asca)[,1], 
         i2 = scale(i2)[,1], 
         ompc = scale(ompc)[,1], 
         cbir_fla = scale(cbir_fla)[,1],
         anca = scale(anca)[,1])


temp_data <- patient_meta_info_lab_test_enrollment %>% 
  left_join(temp_info, by = 'sample_id') %>% 
  select(sample_id, patient_id, phenotype_group1:severity_class, gm_csf:anca) %>% 
  mutate(phenotype_group1 = as.factor(phenotype_group1),
         severity_class = as.factor(severity_class))


temp_data <- temp_data %>% 
  tidyr::pivot_longer(cols = gm_csf:anca, names_to = 'marker', values_to = 'value') %>% 
  filter(!is.na(value))


temp_data$marker <- temp_data$marker %>% 
  recode_factor('gm_csf' = 'GM-CSF',
                'iga_asca' = 'IgA ASCA',
                'igg_asca' = 'IgG ASCA',
                'i2' = 'I2',
                'ompc' = 'OmpC',
                'cbir_fla' = 'CBir Fla',
                'anca' = 'ANCA')

library(introdataviz)
library(tidyverse)
library(ggpubr)
library(rstatix)


stat_test <- temp_data %>%
  group_by(marker) %>%
  wilcox_test(value ~ phenotype_group1) %>% 
  adjust_pvalue() %>%
  add_significance()


temp_data$phenotype_group1 <- temp_data$phenotype_group1 %>% 
  recode_factor('non_IBD' = 'non_IBD',
                'CD' = 'CD')

temp_plot <- ggplot(temp_data, aes(x = marker, y = value, fill = phenotype_group1)) +
  # introdataviz::geom_split_violin(alpha = .4, trim = FALSE, scale = 'width') +
  geom_boxplot(outliers = FALSE) +
  geom_pwc(
    aes(group = phenotype_group1),
    y.position = 1.8,
    tip.length = 0,
    method = "wilcox_test",
    label = "p"
  ) +
  scale_fill_manual(values = c('CD' = '#fb8172',
                               'non_IBD' = '#7fb2d4'),
                    label = c('CD' = 'CD',
                              'non_IBD' = 'Non-IBD'), 
                    name = 'Phenotype group') +
  xlab('Serology marker') +
  ylab('Z-Score') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.8, 0.8))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Figure3/serology_marker_compare_CD_nonIBD_250402.pdf', 
       width = 5.5, height = 4)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
# readr::write_csv(temp_data, '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure6_serology_marker_compare_CD_nonIBD_260720.csv')




# serology marker correlation network ------------------------------------------
# linear association between serology marker and metabolite --------------------

library(tidyverse)
rm(list = ls())

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/00_meta_data/00_meta_info/20241003/patient_meta_info_lab_test_241002.RData')

object_enrollment <- object_final %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(visit_number == 'enrollment')

serology_marker <- patient_meta_info_lab_test %>% 
  filter(visit_number == 'enrollment') %>% 
  select(sample_id, gm_csf:anca) %>%
  mutate(gm_csf = as.numeric(gm_csf), 
         iga_asca = as.numeric(iga_asca),
         igg_asca = as.numeric(igg_asca), 
         i2 = as.numeric(i2), 
         ompc = as.numeric(ompc), 
         cbir_fla = as.numeric(cbir_fla),
         anca = as.numeric(anca)) %>% 
  mutate(gm_csf = scale(gm_csf)[,1], 
         iga_asca = scale(iga_asca)[,1], 
         igg_asca = scale(igg_asca)[,1], 
         i2 = scale(i2)[,1], 
         ompc = scale(ompc)[,1], 
         cbir_fla = scale(cbir_fla)[,1],
         anca = scale(anca)[,1])



object_enrollment@sample_info <- object_enrollment@sample_info %>% 
  left_join(serology_marker, by = 'sample_id')

save(object_enrollment, 
     file = '~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')



variable_interest <- c('gm_csf', 'iga_asca', 'igg_asca', 'i2', 'ompc', 'cbir_fla', 'anca')

library(Maaslin2)
library(tidymass)
# dir.create('association_bmi', showWarnings = FALSE, recursive = TRUE)

walk(seq_along(variable_interest), function(i){
  cat(i, '\n')
  
  # scaling z-score for metabolite data
  temp_object <- object_enrollment %>%
    scale_data(center = TRUE, method = "auto") %>%
    activate_mass_dataset(what = 'annotation_table') %>%
    filter(!is.na(variable_id))
  
  input_annot_table <- temp_object %>%
    extract_expression_data() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'ID')
  
  temp_variable <- variable_interest[i]
  
  # modify meta info
  meta_info <- temp_object %>%
    extract_sample_info() %>% 
    dplyr::select(sample_id, temp_variable, gender, age, race, use_antibiotics)
  
  # remove na & outliers (mean +/- 3sd)
  temp_idx_na <- meta_info[[2]] %>% is.na() %>% which()
  meta_info <- meta_info[-temp_idx_na,]
  
  # remove NA and outlier samples in the input_data
  input_annot_table <- input_annot_table %>%
    dplyr::filter(ID %in% meta_info$sample_id) %>%
    dplyr::arrange(match(meta_info$sample_id, ID))
  
  # write these files
  temp_path <- file.path('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/', temp_variable)
  dir.create(temp_path, showWarnings = FALSE, recursive = TRUE)
  write_tsv(input_annot_table, file = file.path(temp_path, 'input_data_metabolite.tsv'))
  write_tsv(meta_info, file = file.path(temp_path, 'input_data_meta.tsv'))
  
  fit_data <- Maaslin2(
    input_data = file.path(temp_path, 'input_data_metabolite.tsv'),
    input_metadata = file.path(temp_path, 'input_data_meta.tsv'),
    output = temp_path,
    fixed_effects = c(temp_variable, 'gender', 'age', 'race', 'use_antibiotics'),
    normalization = 'NONE',
    transform = 'NONE',
    # random_effects = c('site', 'subject'),
    reference = "race,Caucasian",
    standardize = FALSE,
    save_models = TRUE)
  
  cat('\n\n\n')
})

rm(list = ls());gc()


# read association result ------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')

variable_interest <- c('gm_csf', 'iga_asca', 'igg_asca', 'i2', 'ompc', 'cbir_fla', 'anca')
metabolite_table <- object_enrollment@annotation_table

list_association_result <- lapply(variable_interest, function(x){
  temp_path <- file.path('.', x)
  result <- readr::read_tsv(file.path(temp_path, 'significant_results.tsv'))
  
  cpd_name <- match(result$feature, metabolite_table$variable_id) %>% 
    metabolite_table$Compound.name[.]
  result <- result %>% 
    dplyr::mutate(metabolite = cpd_name,
                  marker_name = x) %>% 
    dplyr::select(metabolite, marker_name, everything())
  
  return(result)
})

names(list_association_result) <- variable_interest
save(list_association_result, 
     file = './list_metabolite_serology_marker_association_result_250331.RData')


table_association_met_serology <- list_association_result %>% 
  bind_rows() %>% 
  filter(value %in% variable_interest)

table_association_met_serology %>% 
  count(qval <= 0.05)


# network analysis -------------------------------------------------------------

library(igraph)
library(tidygraph)

# edge table
edge_table <- table_association_met_serology %>% 
  dplyr::filter(qval <= 0.05) %>% 
  dplyr::mutate(from = feature, to = marker_name) %>% 
  dplyr::select(from, to, everything())

readr::write_tsv(edge_table, 
                 file = '~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/edge_table_250331.tsv')

# node table
node_table <- object_enrollment@annotation_table %>% 
  select(variable_id, Compound.name, metabolon_class, metabolon_subclass) %>% 
  mutate(class = 'metabolite') %>% 
  rename('attribute' = 'metabolon_class',
         'attribute2' = 'metabolon_subclass') %>% 
  add_row(variable_id = variable_interest, 
          Compound.name = variable_interest, 
          class = rep('serology_marker', length(variable_interest)),
          attribute = rep('serology_marker', length(variable_interest)), 
          attribute2 = rep('serology_marker', length(variable_interest))) %>% 
  dplyr::filter(variable_id %in% unique(c(edge_table$from, edge_table$to)))

readr::write_tsv(node_table, 
                 file = '~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/node_table_250331.tsv')

network_association <- graph_from_data_frame(d = edge_table, directed = FALSE, vertices = node_table)

save(network_association, 
     file = './network_association_met_serology_250331.RData')

# network visualization --------------------------------------------------------

library(tidygraph)
network_all <- network_association %>% as_tbl_graph()

# extract serology marker edges 

temp_serology_marker <- network_all %>% 
  activate('edges') %>% 
  # filter(edge_class == 'ibd_serology_marker') %>% 
  mutate(cor_direction = case_when(coef > 0 ~ 'pos',
                                   coef < 0 ~ 'neg')) %>%
  mutate(cor_index_absolute = abs(coef)) %>% 
  activate('nodes') %>% 
  mutate(neighbour_count = centrality_degree()) %>% 
  mutate(node_size = case_when(class == 'metabolite' ~ 1,
                               class != 'metabolite' ~ 3)) %>% 
  filter(neighbour_count > 0)

set.seed(20241106)
temp_plot <- ggraph(temp_serology_marker, layout = 'fr') + 
  geom_edge_fan(aes(width = cor_index_absolute, colour = cor_direction)) + 
  scale_edge_colour_manual(values = c('pos' = '#ff595e',
                                      'neg' = '#1982c4'),
                           label = c('pos' = 'Positive',
                                     'neg' = 'Negative'),
                           name = 'Correlation direction') +
  scale_edge_width_continuous(range = c(0.5, 2.5),
                              name = 'Spearman correlation coefficient') +
  geom_node_point(aes(size = node_size, fill = class), shape = 21) +
  scale_size_continuous(range = c(3.5, 9)) +
  scale_fill_manual(values = c('serology_marker' = '#eddea4',
                               'metabolite' = '#cce3de'), 
                    label = c('serology_marker' = 'Serology marker',
                              'metabolite' = 'Metabolite'),
                    name = 'Node class') +
  geom_node_text(aes(label = Compound.name), repel = TRUE, label.r = 0) +
  theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(), 
        legend.key = ggplot2::element_blank())


ggsave(plot = temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/network_association_met_serology_250331.pdf', 
       width = 14, height = 8)


# temp_source_data <- edge_table
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
# readr::write_csv(temp_source_data, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure6_network_association_met_serology_260720.csv')

# stack bar plot of nodes and edges --------------------------------------------
temp_data_order <- temp_serology_marker %>% 
  activate('edges') %>% 
  as_tibble() %>% 
  count(marker_name) %>% 
  arrange(desc(n))

temp_data <- temp_serology_marker %>% 
  activate('edges') %>% 
  as_tibble() %>% 
  count(marker_name, cor_direction) %>% 
  arrange(match(marker_name, rev(temp_data_order$marker_name)))


temp_data$marker_name <- recode_factor(temp_data$marker_name, !!!setNames(rev(temp_data_order$marker_name), rev(temp_data_order$marker_name)))

temp_plot <- ggplot(temp_data, aes(x = marker_name, y = n, fill = cor_direction)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('pos' = '#ff595e',
                               'neg' = '#1982c4'),
                    label = c('pos' = 'Positive',
                              'neg' = 'Negative'),
                    name = 'Correlation direction') +
  geom_text(aes(label = n), color = "white", size = 5, 
            position = position_stack(vjust = 0.5), 
            show.legend = FALSE) +
  coord_flip() +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(legend.position = c(0.9, 0.2))


ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/serology_marker_association_stat_250331.pdf', 
       width = 3, height = 6)

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
temp_source_data <- temp_data
readr::write_csv(temp_source_data, 
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure6_serology_marker_association_stat_260720.csv')


# circular association plot ----------------------------------------------------
# Chord Diagram ----------------------------------------------------------------


edge_table_serology <- temp_serology_marker %>% 
  activate('edges') %>% 
  as_tibble()

node_table_serology <- temp_serology_marker %>% 
  activate('nodes') %>% 
  as_tibble()

idx <- match(edge_table_serology$feature, node_table_serology$name)
edge_table_serology <- edge_table_serology %>% 
  mutate(attribute = node_table_serology$attribute[idx],
         attribute2 = node_table_serology$attribute2[idx])

temp_data <- edge_table_serology %>% 
  count(metadata, attribute, cor_direction) %>% 
  # replace_na(value = 'Other') %>% 
  mutate(metadata = case_when(metadata == 'anca' ~ 'ANCA',
                              metadata == 'cbir_fla' ~ 'CBir-Fla',
                              metadata == 'gm_csf' ~ 'GM-SCF',
                              metadata == 'i2' ~ 'I2',
                              metadata == 'iga_asca' ~ 'IgA ASCA',
                              metadata == 'igg_asca' ~ 'IgG ASCA',
                              metadata == 'ompc' ~ 'OmpC'))

# temp_data <- edge_table_serology %>% 
#   count(metadata, class, cor_direction) %>% 
#   # replace_na(value = 'Other') %>% 
#   mutate(metadata = case_when(metadata == 'anca' ~ 'ANCA',
#                               metadata == 'cbir_fla' ~ 'CBir-Fla',
#                               metadata == 'gm_scf' ~ 'GM-SCF',
#                               metadata == 'i2' ~ 'I2',
#                               metadata == 'lga_asca' ~ 'IgA ASCA',
#                               metadata == 'lgg_asca' ~ 'IgG ASCA',
#                               metadata == 'ompc' ~ 'OmpC'))

library(circlize)
# set.seed(12345678)


grid.col = c('ANCA' = "#8ed3c6", 
             'CBir-Fla' = "#ffed6f", 
             'GM-SCF' = "#bebada",
             'I2' = "#fdb461", 
             'IgA ASCA' = "#fb8072", 
             'IgG ASCA' = "#7fb1d3", 
             'OmpC' = "#b5de68")

set.seed(20241108)
temp_plot <- chordDiagram(temp_data, grid.col = grid.col)



pdf('~/Project/00_IBD_project/Figure/250326/Figure3/circular_association_plot_250331.pdf', width = 10, height = 10)
grid.col = c('ANCA' = "#8ed3c6", 
             'CBir-Fla' = "#ffed6f", 
             'GM-SCF' = "#bebada",
             'I2' = "#fdb461", 
             'IgA ASCA' = "#fb8072", 
             'IgG ASCA' = "#7fb1d3", 
             'OmpC' = "#b5de68")

set.seed(20241108)
temp_plot <- chordDiagram(temp_data, grid.col = grid.col)
dev.off()

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
# temp_source_data <- temp_data
# readr::write_csv(temp_source_data, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure6_circular_association_plot_260720.csv')

rm(list = ls());gc()





################################################################################
# Serology marker and metabolite association example ----------------------------
# 2-oxindole association -------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/network_association_met_serology_250331.RData')


library(tidygraph)
library(sjmisc)
network_all <- network_association %>% as_tbl_graph()

temp_2_oxindole <- object_enrollment %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(Compound.name == '2-Oxindole') %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id')


temp_serology_marker <- object_enrollment %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1, iga_asca)


temp_data <- temp_2_oxindole %>% 
  left_join(temp_serology_marker, by = 'sample_id') %>% 
  filter(!is.na(iga_asca))

temp_plot <- ggplot(data = temp_data) +
  geom_point(aes(x = iga_asca, y = M134T46_hilic_pos, color = phenotype_group1), size = 2, data = temp_data) +
  geom_smooth(aes(x = iga_asca, y = M134T46_hilic_pos), method = "lm", se = TRUE) +
  scale_color_manual(values = c('CD' = '#fb8172',
                                'non_IBD' = '#7fb2d4'),
                     label = c('CD' = 'CD',
                               'non_IBD' = 'Non-IBD'), 
                     name = 'Phenotype group') +
  xlab('IgA ASCA (z-score)') +
  ylab('2-Oxindole (z-score)') + 
  ZZWtool::ZZWTheme() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text = element_text(hjust = 0.5, vjust = 0.5))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/2_oxindole_iga_asca_association_250321.pdf', 
       width = 6, height = 6)

temp_lm <- lm(M134T46_hilic_pos ~ iga_asca, data = temp_data)
summary(temp_lm)

cor(temp_data$M134T46_hilic_pos, temp_data$iga_asca, method = 'pearson')

temp_2_oxindole_source_data <- temp_data

# PC(34:0) association -------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/network_association_met_serology_250331.RData')


library(tidygraph)
library(sjmisc)
network_all <- network_association %>% as_tbl_graph()

temp_pc_34_0 <- object_enrollment %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(Compound.name == 'PC(34:0)') %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id')


temp_serology_marker <- object_enrollment %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1, igg_asca)


temp_data <- temp_pc_34_0 %>% 
  left_join(temp_serology_marker, by = 'sample_id') %>% 
  filter(!is.na(igg_asca))

temp_plot <- ggplot(data = temp_data) +
  geom_point(aes(x = igg_asca, y = M763T746_c18_pos, color = phenotype_group1), size = 2, data = temp_data) +
  geom_smooth(aes(x = igg_asca, y = M763T746_c18_pos), method = "lm", se = TRUE) +
  scale_color_manual(values = c('CD' = '#fb8172',
                                'non_IBD' = '#7fb2d4'),
                     label = c('CD' = 'CD',
                               'non_IBD' = 'Non-IBD'), 
                     name = 'Phenotype group') +
  xlab('IgG ASCA (z-score)') +
  ylab('PC(34:0) (z-score)') + 
  ZZWtool::ZZWTheme() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text = element_text(hjust = 0.5, vjust = 0.5))

ggsave(temp_plot,  
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/pc_34_0_iga_asca_association_250331.pdf', 
       width = 6, height = 6)

temp_lm <- lm(M763T746_c18_pos ~ igg_asca, data = temp_data)
summary(temp_lm)

cor(temp_data$M763T746_c18_pos, temp_data$igg_asca, method = 'pearson')
temp_pc_34_0_source_data <- temp_data


# Arg-C20:2 --------------------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/network_association_met_serology_250331.RData')


library(tidygraph)
library(sjmisc)
network_all <- network_association %>% as_tbl_graph()

temp_arg_c20_2 <- object_enrollment %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(Compound.name == 'Arg-C20:2') %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id')


temp_serology_marker <- object_enrollment %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1, igg_asca)


temp_data <- temp_arg_c20_2 %>% 
  left_join(temp_serology_marker, by = 'sample_id') %>% 
  filter(!is.na(igg_asca))

temp_plot <- ggplot(data = temp_data) +
  geom_point(aes(x = igg_asca, y = M465T614_c18_pos, color = phenotype_group1), size = 2, data = temp_data) +
  geom_smooth(aes(x = igg_asca, y = M465T614_c18_pos), method = "lm", se = TRUE) +
  scale_color_manual(values = c('CD' = '#fb8172',
                                'non_IBD' = '#7fb2d4'),
                     label = c('CD' = 'CD',
                               'non_IBD' = 'Non-IBD'), 
                     name = 'Phenotype group') +
  xlab('IgG ASCA (z-score)') +
  ylab('Arg-C20:2 (z-score)') + 
  ZZWtool::ZZWTheme() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text = element_text(hjust = 0.5, vjust = 0.5))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/arg_c20_2_iga_asca_association_250331.pdf', 
       width = 6, height = 6)

temp_lm <- lm(M465T614_c18_pos ~ igg_asca, data = temp_data)
summary(temp_lm)

cor(temp_data$M465T614_c18_pos, temp_data$igg_asca, method = 'pearson')
network_association %>% 
  as_tibble()

temp_arg_c20_2_source_data <- temp_data

# Bilirubin --------------------------------------------------------------------
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/network_association_met_serology_250331.RData')


library(tidygraph)
library(sjmisc)
network_all <- network_association %>% as_tbl_graph()

temp_bilirubin <- object_enrollment %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  filter(Compound.name == 'Bilirubin') %>% 
  scale_data(center = TRUE, method = 'auto') %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'sample_id')


temp_serology_marker <- object_enrollment %>% 
  extract_sample_info() %>% 
  select(sample_id, phenotype_group1, iga_asca)


temp_data <- temp_bilirubin %>% 
  left_join(temp_serology_marker, by = 'sample_id') %>% 
  filter(!is.na(iga_asca))

temp_plot <- ggplot(data = temp_data) +
  geom_point(aes(x = iga_asca, y = M585T45_hilic_pos, color = phenotype_group1), size = 2, data = temp_data) +
  geom_smooth(aes(x = iga_asca, y = M585T45_hilic_pos), method = "lm", se = TRUE) +
  scale_color_manual(values = c('CD' = '#fb8172',
                                'non_IBD' = '#7fb2d4'),
                     label = c('CD' = 'CD',
                               'non_IBD' = 'Non-IBD'), 
                     name = 'Phenotype group') +
  xlab('IgA ASCA (z-score)') +
  ylab('Bilirubin (z-score)') + 
  ZZWtool::ZZWTheme() + 
  theme(legend.position = c(0.9,0.9), 
        axis.text = element_text(hjust = 0.5, vjust = 0.5))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure4/bilirubin_iga_asca_association_250331.pdf', 
       width = 6, height = 6)

temp_lm <- lm(M585T45_hilic_pos ~ iga_asca, data = temp_data)
summary(temp_lm)
cor(temp_data$M585T45_hilic_pos, temp_data$iga_asca, method = 'pearson')

network_association %>% 
  as_tibble()

temp_bilirubin_source_data <- temp_data

# 
# temp_2_oxindole_source_data <- temp_2_oxindole_source_data %>% 
#   select(sample_id, iga_asca, M134T46_hilic_pos) %>% 
#   rename('marker_value' = 'iga_asca', 'metabolite_value' = 'M134T46_hilic_pos') %>% 
#   mutate(metabolite = '2-Oxindole', marker_name = 'iga_asca')
# 
# temp_pc_34_0_source_data <- temp_pc_34_0_source_data %>%
#   select(sample_id, igg_asca, M763T746_c18_pos) %>% 
#   rename('marker_value' = 'igg_asca', 'metabolite_value' = 'M763T746_c18_pos') %>% 
#   mutate(metabolite = 'PC(34:0)', marker_name = 'igg_asca')
# 
# temp_arg_c20_2_source_data <- temp_arg_c20_2_source_data %>%
#   select(sample_id, igg_asca, M465T614_c18_pos) %>% 
#   rename('marker_value' = 'igg_asca', 'metabolite_value' = 'M465T614_c18_pos') %>% 
#   mutate(metabolite = 'Arg-C20:2', marker_name = 'igg_asca')
# 
# temp_bilirubin_source_data <- temp_bilirubin_source_data %>%
#   select(sample_id, iga_asca, M585T45_hilic_pos) %>%
#   rename('marker_value' = 'iga_asca', 'metabolite_value' = 'M585T45_hilic_pos') %>%
#   mutate(metabolite = 'Bilirubin', marker_name = 'iga_asca')
#   
# 
# temp_data <- bind_rows(temp_2_oxindole_source_data, 
#                        temp_pc_34_0_source_data, 
#                        temp_arg_c20_2_source_data, 
#                        temp_bilirubin_source_data)
# 
# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
# readr::write_csv(temp_data,
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure6_metabolite_serology_marker_example_260720.csv')



################################################################################
# Individual serology marker predictive model ----------------------------------
# Part 0: set up environment and functions -------------------------------------

# summarize the odds results
summarize_odds <- function(logistic_model) {
  CI_lower <- coefficients(logistic_model) - 1.96*summary(logistic_model)$coefficients[,2]
  CI_upper <- coefficients(logistic_model) + 1.96*summary(logistic_model)$coefficients[,2]
  odds <- exp(coef(logistic_model))
  odds_ci_min <- exp(CI_lower)
  odds_ci_max <- exp(CI_upper)
  # odds_ci <- paste(exp(CI_lower), exp(CI_upper), sep = ';')
  result_table <- tibble::tibble(variable = names(coef(logistic_model)),
                                 odds = odds, 
                                 odds_ci_min = odds_ci_min, 
                                 odds_ci_max = odds_ci_max)
}

# Part 1: predictive model for serology markers --------------------------------

# load data and set up environment 
setwd('~/Project/00_IBD_project/20260420_serology_marker_evaluation_Fig3/')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(ggraph)
library(glmnet)
library(pROC)

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')

sample_info_enrollemnt_serology <- object_enrollment@sample_info %>% 
  dplyr::select(sample_id, phenotype_group1, age, gender, race, use_antibiotics, gm_csf:anca) %>% 
  dplyr::mutate(disease = case_when(phenotype_group1 == 'CD' ~ 1,
                                    phenotype_group1 == 'non_IBD' ~ 0),
                gender = as.factor(gender),
                
                # Group rare races to prevent cross-validation train/test split errors
                race = case_when(
                  race %in% c('Caucasian', 'Black or African American', 'Asian', 'Hispanic/Latino') ~ race,
                  TRUE ~ 'Others/Unknown' # Groups Native Hawaiian, American Indian, Unknown, and Others
                ),
                race = as.factor(race),
                
                use_antibiotics = as.factor(use_antibiotics)) %>% 
  dplyr::select(sample_id, disease, age, gender, race, use_antibiotics, gm_csf:anca)

save(sample_info_enrollemnt_serology, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_serology_260420.RData')

# Convert to dataframe with row names (avoiding as.matrix to preserve factor types for glm)
sample_info_model_data <- sample_info_enrollemnt_serology %>% 
  tibble::column_to_rownames(var = 'sample_id')

save(sample_info_model_data, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_info_enrollemnt_serology_scaled_260420.RData')


# 10-fold cross validation for serology markers --------------------------------

sample_ids <- rownames(sample_info_model_data)

# generate seed list to reproduce the results
set.seed(20241111)
index <- sample(10000000, 100)

list_serology_10_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index[z]
  
  set.seed(seed_id)
  
  # Stratified Cross-Validation
  # obtain the patient IDs for CD (disease = 1) and non_IBD (disease = 0)
  id_case <- rownames(sample_info_model_data)[sample_info_model_data$disease == 1]
  id_control <- rownames(sample_info_model_data)[sample_info_model_data$disease == 0]
  
  # randomly shuffle the patient IDs within each group to ensure randomness in fold assignment
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  # divide the patient IDs into 10 folds for each group separately
  groups_case <- rep(1:10, length.out = length(id_case))
  groups_control <- rep(1:10, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  # combine the folds from both groups to create the final list of patient IDs for each fold
  cv_splits <- lapply(1:10, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_10_fold_cv <- lapply(1:10, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    # split data sets into training and testing sets directly
    data_train <- sample_info_model_data[ids_training, ]
    data_test <- sample_info_model_data[ids_test, ]
    
    # Data Leakage Fix: Compute scaling paraers on the training set ONLY
    # Scale age
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    serology_names <- colnames(data_train)[-c(1:5)] # cols 1-5: disease, age, gender, race, use_antibiotics
    
    # Scale serology markers (based on training set)
    for(marker in serology_names){
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    list_roc_serology <- lapply(seq_along(serology_names), function(i){
      x <- serology_names[i]
      
      # Robust column renaming
      temp_data_train <- data_train %>% 
        as_tibble() %>% 
        dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(x)) %>% 
        dplyr::rename(x = dplyr::all_of(x)) %>% 
        dplyr::filter(!is.na(x))
      
      temp_data_test <- data_test %>% 
        as_tibble() %>% 
        dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(x)) %>% 
        dplyr::rename(x = dplyr::all_of(x)) %>% 
        dplyr::filter(!is.na(x))
      
      # logistic regression - training set
      logistic_model <- glm(disease ~ x + gender + age + race + use_antibiotics, data = temp_data_train, family = binomial(link = 'logit'))
      
      CI_lower <- coefficients(logistic_model)[2] - 1.96*summary(logistic_model)$coefficients[2,2]
      CI_upper <- coefficients(logistic_model)[2] + 1.96*summary(logistic_model)$coefficients[2,2]
      odds <- exp(coef(logistic_model)[2])
      odds_ci <- paste(c(exp(CI_lower), exp(CI_upper)), collapse = ';')
      
      # ROC curve --- training set
      possibility_train <- predict(logistic_model, newdata = temp_data_train, type = "response")
      roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
      auc_train <- roc_obj_train$auc[[1]]
      auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
      
      # ROC curve --- testing set
      possibility_test <- predict(logistic_model, newdata = temp_data_test, type = "response")
      roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
      auc_test <- roc_obj_test$auc[[1]]
      auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
      
      stat_values <- tibble::tibble(variable_id = x,
                                    odds = odds,
                                    odds_ci = odds_ci,
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
      
      result <- list(stat_result = stat_values, 
                     data_train = temp_data_train,
                     data_test = temp_data_test,
                     logistic_model = logistic_model,
                     roc_train = roc_obj_train,
                     roc_test = roc_obj_test)
      
      return(result)
    })
    
    names(list_roc_serology) <- serology_names
    
    return(list_roc_serology)
  })
  
  names(result_list_10_fold_cv) <- paste0('fold_', 1:10)
  
  return(result_list_10_fold_cv)
})

names(list_serology_10_fold_cv) <- paste0('repeat_', 1:100)

save(list_serology_10_fold_cv, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/list_serology_10_fold_cv_260420.RData')

result_serology_10_fold_cv <- lapply(list_serology_10_fold_cv, function(z){
  z <- lapply(z, function(x){
    x <- lapply(x, function(y){
      y$stat_result
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows() %>% 
    mutate(fold = rep(paste0('fold_', 1:10), each = 7))
}) %>% 
  bind_rows() %>% 
  mutate(times = rep(paste0('repeat_', 1:100), each = 70))

# split odds_ci, auc_ci_train, and auc_ci_test into two columns
result_serology_10_fold_cv <- result_serology_10_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq(n())) %>% 
  dplyr::rename('odds_ratio' = 'odds') %>% 
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_serology_10_fold_cv, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/result_serology_10_fold_cv_260420.RData')

list_serology_10_fold_cv$repeat_1$fold_1$gm_csf$roc_test$original.response
list_serology_10_fold_cv$repeat_1$fold_1$gm_csf$roc_test$original.predictor

# calculate the mean odds ratio, odds_ci_min, odds_ci_max, auc_train, auc_test, auc_ci_train_min, auc_ci_train_max, auc_ci_test_min, auc_ci_test_max
temp_serology_table <- result_serology_10_fold_cv %>% 
  group_by(variable_id) %>% 
  summarise(mean_odds_ratio = mean(odds_ratio),
            mean_odds_ci_min = mean(odds_ci_min),
            mean_odds_ci_max = mean(odds_ci_max),
            mean_auc_train = mean(auc_train),
            mean_auc_test = mean(auc_test),
            mean_auc_ci_train_min = mean(auc_ci_train_min),
            mean_auc_ci_train_max = mean(auc_ci_train_max),
            mean_auc_ci_test_min = mean(auc_ci_test_min),
            mean_auc_ci_test_max = mean(auc_ci_test_max))


# plot serology marker 100 times 10-fold cross validation ----------------------

result_serology_marker <- temp_serology_table %>%
  dplyr::mutate(idx = seq(n())) %>%
  dplyr::rename('odds_ratio' = 'mean_odds_ratio',
                'odds_ci_min' = 'mean_odds_ci_min',
                'odds_ci_max' = 'mean_odds_ci_max') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

temp_plot <- ggplot(result_serology_marker, aes(x = rev(idx), y = log10(odds_ratio), color = color)) +
  geom_pointrange(aes(ymin = log10(odds_ci_min), ymax = log10(odds_ci_max)), size = 1.3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(result_serology_marker$variable_id)), 
                     label = c('OmpC', 'IgG ASCA', 'IgA ASCA', 'I2', 'GM-CSF', 'CBir Fla', 'ANCA')
  ) +
  scale_colour_manual(values = c('protective factor' = '#7fb1d3',
                                 'risk factor' = '#fc8070'),
                      label = c('protective factor' = 'Protective factor',
                                'risk factor' = 'Risk factor')) +
  scale_y_continuous(breaks = c(-1, 0, 1), label = c(0.1, 1, 10), limits = c(-1.3, 2.1)) +
  xlab('Serology marker') +
  ylab('Odds Ratio (OR)') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

ggsave(plot = temp_plot,
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/serology_odds_ratio_100cv_260420.pdf',
       width = 9, height = 6)

# Plot ROC curve ---------------------------------------------------------------

result_serology_marker_10_fold_cv <- lapply(seq_along(list_serology_10_fold_cv), function(i){
  z <- list_serology_10_fold_cv[[i]]
  z <- lapply(seq_along(z), function(n){
    x <- z[[n]]
    x <- lapply(seq_along(x), function(m){
      y <- x[[m]]
      temp_result <- tibble::tibble(variable_id = y$stat_result$variable_id,
                                    disease = y$data_test$disease, 
                                    possibility = y$roc_test$predictor)
    }) %>% 
      bind_rows() %>% 
      mutate(times = paste0('repeat_', i),
             fold = paste0('fold_', n))
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

temp_serology_marker <- unique(result_serology_marker_10_fold_cv$variable_id)

list_roc_test <- lapply(temp_serology_marker, function(x){
  temp_data <- result_serology_marker_10_fold_cv %>% 
    filter(variable_id == x)
  
  roc_test <- roc(temp_data$disease ~ temp_data$possibility, 
                  plot = FALSE, 
                  print.auc = TRUE, 
                  ci = TRUE)
  
  return(roc_test)
})

names(list_roc_test) <- temp_serology_marker
save(list_roc_test, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/list_roc_test_serology_260420.RData')

temp_plot <- ggroc(list_roc_test, size = 1) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed", size = 1) +
  scale_color_discrete(label = temp_serology_marker) +
  xlab('Specificity') +
  ylab('Sensitivity') +
  ZZWtool::ZZWTheme() +
  theme(legend.position = c(0.8, 0.3))

ggsave(plot = temp_plot,
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/serology_roc_10cv_260420.pdf',
       width = 6,
       height = 6)

table_auc_serology <- lapply(list_roc_test, function(x){
  tibble::tibble(auc = as.numeric(x$auc), auc_ci_min = x$ci[1], auc_ci_max = x$ci[3])
}) %>% 
  dplyr::bind_rows() %>%
  mutate(serology_name = temp_serology_marker) %>%
  select(serology_name, everything()) %>%
  arrange(desc(serology_name)) %>% 
  dplyr::mutate(idx = seq(n()))

save(table_auc_serology, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_serology_260420.RData')

mean_auc <- table_auc_serology %>% 
  summarise(mean_auc = mean(auc),
            mean_auc_ci_min = mean(auc_ci_min),
            mean_auc_ci_max = mean(auc_ci_max)) %>% 
  pull(mean_auc)

temp_plot <- ggplot(table_auc_serology, aes(x = idx, y = auc)) +
  geom_pointrange(aes(ymin = auc_ci_min, ymax = auc_ci_max), colour = 'black', size = 1.3) +
  geom_hline(yintercept = 0.7376, linetype = 'dashed', size = 1) +
  scale_x_continuous(breaks = seq(length(table_auc_serology$serology_name)), 
                     label = table_auc_serology$serology_name) +
  ylim(0.6, 0.9) +
  xlab('Serological marker') +
  ylab('AUC') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

ggsave(plot = temp_plot,
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/serology_auc_10cv_260420.pdf',
       width = 9,
       height = 6)


rm(list = ls());gc()

temp_data <- temp_serology_table %>% 
  left_join(table_auc_serology, by = c('variable_id' = 'serology_name')) %>% 
  dplyr::select(-idx)

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
readr::write_csv(temp_data, 
                 path = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure7_individual_serology_marker_predictive_model_260720.csv')









################################################################################
# Individual Metabolite marker predictive model --------------------------------
load('~/Project/00_IBD_project/Data/20250328_metabolic_signiture_enrollment_CD_nonIBD/object_stat_250328.RData')

# use metabolites with p-value < 0.5
met_sig <- object_stat %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(p_value < 0.5) %>%
  extract_expression_data() %>% 
  sjmisc::rotate_df() %>% 
  rownames_to_column(var = 'sample_id')

# extract the cofounder
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')

sample_info_enrollemnt_cofounder <- object_enrollment@sample_info %>% 
  dplyr::select(sample_id, phenotype_group1, age, gender, race, use_antibiotics) %>% 
  dplyr::mutate(disease = case_when(phenotype_group1 == 'CD' ~ 1,
                                    phenotype_group1 == 'non_IBD' ~ 0),
                gender = as.factor(gender),
                race = case_when(
                  race %in% c('Caucasian', 'Black or African American', 'Asian', 'Hispanic/Latino') ~ race,
                  TRUE ~ 'Others/Unknown'
                ),
                race = as.factor(race),
                use_antibiotics = as.factor(use_antibiotics)) %>% 
  dplyr::select(sample_id, disease, age, gender:use_antibiotics)


sample_info_enrollment_metabolite <- sample_info_enrollemnt_cofounder %>% 
  left_join(met_sig, by = 'sample_id')

save(sample_info_enrollment_metabolite, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_metabolite_260420.RData')


# 10-fold cross validation for metabolite markers ------------------------------
sample_ids <- sample_info_enrollment_metabolite$sample_id

set.seed(20250106)
sample_ids_discovery <- sample(sample_ids, floor(length(sample_ids) * 0.8))
sample_ids_validation <- sample_ids[!sample_ids %in% sample_ids_discovery]


# generate seed list to reproduce the results:
# Note: we can set a smaller number of randomization (e.g., 10) for testing the code, and then set a larger number (e.g., 100) for the final result.


# set.seed(20241111)
# # index <- sample(10000000, 100)
# index <- sample(10000000, 10)

# list_metabolite_10_fold_cv <- lapply(seq_along(index), function(z){
#   cat('Perform the ', z, 'th randomization\n')
#   seed_id <- index[z]

#   set.seed(seed_id)
#   # stratified cross-validation within discovery set
#   id_case <- sample_info_enrollment_metabolite %>% 
#     dplyr::filter(sample_id %in% sample_ids_discovery, disease == 1) %>% 
#     dplyr::pull(sample_id)
#   id_control <- sample_info_enrollment_metabolite %>% 
#     dplyr::filter(sample_id %in% sample_ids_discovery, disease == 0) %>% 
#     dplyr::pull(sample_id)

#   id_case <- sample(id_case)
#   id_control <- sample(id_control)

#   groups_case <- rep(1:10, length.out = length(id_case))
#   groups_control <- rep(1:10, length.out = length(id_control))

#   split_case <- split(id_case, groups_case)
#   split_control <- split(id_control, groups_control)

#   cv_splits <- lapply(1:10, function(f) {
#     c(split_case[[f]], split_control[[f]])
#   })


#   result_list_10_fold_cv <- lapply(1:10, function(j){
#     cat('Perform the ', j, 'th cross validation\n')
#     ids_test <- unlist(cv_splits[[j]])
#     ids_training <- unlist(cv_splits[-j])

#     data_train <- sample_info_enrollment_metabolite %>% 
#       dplyr::filter(sample_id %in% ids_training)
#     data_test <- sample_info_enrollment_metabolite %>% 
#       dplyr::filter(sample_id %in% ids_test)

#     # scale age based on training set only
#     age_mean <- mean(data_train$age, na.rm = TRUE)
#     age_sd <- sd(data_train$age, na.rm = TRUE)
#     data_train$age <- (data_train$age - age_mean) / age_sd
#     data_test$age <- (data_test$age - age_mean) / age_sd

#     # train logistic model
#     metabolite_names <- colnames(data_train)[-c(1:6)]

#     list_roc_metabolite <- lapply(seq_along(metabolite_names), function(i){
#       # cat(i, ' ')
#       x <- metabolite_names[i]
#       temp_data_train <- data_train %>% 
#         as_tibble() %>% 
#         dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(x)) %>% 
#         dplyr::rename(x = dplyr::all_of(x))
#       temp_data_train <- temp_data_train %>% 
#         dplyr::filter(!is.na(x))

#       temp_data_test <- data_test %>% 
#         as_tibble() %>% 
#         dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(x)) %>% 
#         dplyr::rename(x = dplyr::all_of(x))
#       temp_data_test <- temp_data_test %>% 
#         dplyr::filter(!is.na(x))

#       # scale metabolite based on training set only
#       x_mean <- mean(temp_data_train$x, na.rm = TRUE)
#       x_sd <- sd(temp_data_train$x, na.rm = TRUE)
#       if (!is.na(x_sd) && x_sd > 0) {
#         temp_data_train$x <- (temp_data_train$x - x_mean) / x_sd
#         temp_data_test$x <- (temp_data_test$x - x_mean) / x_sd
#       } else {
#         temp_data_train$x <- 0
#         temp_data_test$x <- 0
#       }


#       # logistic regression - training set
#       logistic_model <- glm(disease ~ x + gender + age + race + use_antibiotics, data = temp_data_train, family = binomial(link = 'logit'))
#       CI_lower <- coefficients(logistic_model)[2] - 1.96*summary(logistic_model)$coefficients[2,2]
#       CI_upper <- coefficients(logistic_model)[2] + 1.96*summary(logistic_model)$coefficients[2,2]
#       odds <- exp(coef(logistic_model)[2])
#       odds_ci <- paste(c(exp(CI_lower), exp(CI_upper)), collapse = ';')

#       # ROC curve --- training set
#       possibility <- predict(logistic_model, newdata = temp_data_train, type = "response")
#       roc_obj_train <- roc(temp_data_train$disease ~ possibility, plot = FALSE, print.auc = TRUE, ci = TRUE)
#       auc_train <- roc_obj_train$auc[[1]]
#       auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')

#       # ROC curve --- testing set
#       possibility <- predict(logistic_model, newdata = temp_data_test, type = "response")
#       roc_obj_test <- roc(temp_data_test$disease ~ possibility, plot = FALSE, print.auc = TRUE, ci = TRUE)
#       auc_test <- roc_obj_test$auc[[1]]
#       auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')

#       stat_values <- tibble::tibble(variable_id = x,
#                                     odds = odds,
#                                     odds_ci = odds_ci,
#                                     auc_train = auc_train,
#                                     auc_ci_train = auc_ci_train,
#                                     auc_test = auc_test,
#                                     auc_ci_test = auc_ci_test)

#       result <- list(stat_result = stat_values, 
#                      data_train = temp_data_train,
#                      data_test = temp_data_test,
#                      logistic_model = logistic_model,
#                      roc_train = roc_obj_train,
#                      roc_test = roc_obj_test)

#       return(result)
#     })

#     names(list_roc_metabolite) <- metabolite_names

#     return(list_roc_metabolite)

#   })

#   names(result_list_10_fold_cv) <- paste0('fold_', 1:10)

#   return(result_list_10_fold_cv)
# })


# names(list_metabolite_10_fold_cv) <- paste0('repeat_', 1:10)

# save(list_metabolite_10_fold_cv, 
#      file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/list_metabolite_10_fold_cv_10_times_260420.RData')


# I ran the 100 randomization with 10-fold CV on SCG cluster, so I loaded the results below
# The script: /Users/zhouzw/Project/00_IBD_project/Data/20260420_single_met_logistic_evaluation/demo_code.r

load('~/Project/00_IBD_project/Data/20260420_single_met_logistic_evaluation/list_metabolite_10_fold_cv_100_repeat_260421.RData')

# split odds_ci, auc_ci_train, and auc_ci_test into two columns
result_metabolite_10_fold_cv <- list_metabolite_10_fold_cv$stat_values %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq(n())) %>% 
  dplyr::rename('odds_ratio' = 'odds') %>% 
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))


save(result_metabolite_10_fold_cv, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/result_metabolite_10_fold_cv_100_times_260423.RData')


metabolite_info <- object_stat@annotation_table %>% 
  select(variable_id, Compound.name)

# calculate the mean odds ratio, odds_ci_min, odds_ci_max, auc_train, auc_test, auc_ci_train_min, auc_ci_train_max, auc_ci_test_min, auc_ci_test_max
temp_metabolite_table <- result_metabolite_10_fold_cv %>% 
  group_by(variable_id) %>% 
  summarise(mean_odds_ratio = mean(odds_ratio),
            mean_odds_ci_min = mean(odds_ci_min),
            mean_odds_ci_max = mean(odds_ci_max),
            mean_auc_train = mean(auc_train),
            mean_auc_test = mean(auc_test),
            mean_auc_ci_train_min = mean(auc_ci_train_min),
            mean_auc_ci_train_max = mean(auc_ci_train_max),
            mean_auc_ci_test_min = mean(auc_ci_test_min),
            mean_auc_ci_test_max = mean(auc_ci_test_max)) %>% 
  arrange(desc(mean_auc_test)) %>% 
  left_join(metabolite_info, by = 'variable_id') %>% 
  select(variable_id, Compound.name, everything())


# plot metabolite marker 100 times 10-fold cross validation ----------------------


result_metabolite_marker <- temp_metabolite_table %>%
  dplyr::mutate(idx = seq(n())) %>%
  dplyr::rename('odds_ratio' = 'mean_odds_ratio',
                'odds_ci_min' = 'mean_odds_ci_min',
                'odds_ci_max' = 'mean_odds_ci_max') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))


# Plot ROC curve ---------------------------------------------------------------
result_metabolite_marker_10_fold_cv <- list_metabolite_10_fold_cv$prediction_result_test

temp_metabolite_marker <- result_metabolite_marker$variable_id
temp_metabolite_name <- result_metabolite_marker$Compound.name

list_roc_test <- lapply(temp_metabolite_marker, function(x){
  temp_data <- result_metabolite_marker_10_fold_cv %>% 
    filter(variable_id == x)
  
  roc_test <- roc(temp_data$disease ~ temp_data$possibility, 
                  plot = FALSE, 
                  print.auc = TRUE, 
                  ci = TRUE)
  
  return(roc_test)
})


names(list_roc_test) <- temp_metabolite_marker


table_auc <- lapply(list_roc_test, function(x){
  tibble::tibble(auc = as.numeric(x$auc), auc_ci_min = x$ci[1], auc_ci_max = x$ci[3])
}) %>% 
  dplyr::bind_rows() %>%
  mutate(variable_id = temp_metabolite_marker,
         compound_name = temp_metabolite_name) %>%
  select(compound_name, variable_id, everything())

table_auc_sig <- table_auc %>% 
  arrange(desc(auc)) %>% 
  filter(auc >= 0.7376) %>%
  dplyr::mutate(idx = seq(n()))

save(table_auc_sig, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_sig_260423.RData')

temp_plot <- ggplot(table_auc_sig, aes(x = idx, y = auc)) +
  geom_pointrange(aes(ymin = auc_ci_min, ymax = auc_ci_max), colour = 'black') +
  geom_hline(yintercept = 0.7376, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(table_auc_sig$variable_id)), 
                     label = table_auc_sig$compound_name
                     # label = c('OmpC', 'IgG ASCA', 'IgA ASCA', 'I2', 'GM-CSF', 'CBir Fla', 'ANCA')
  ) +
  ylim(0.6, 0.9) +
  xlab('Metabolite marker') +
  ylab('AUC') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

dir.create('~/Project/00_IBD_project/Figure/250326/Figure3/',
           showWarnings = FALSE, recursive = TRUE)
ggsave(plot = temp_plot,
       filename = '~/Project/00_IBD_project/Figure/250326/Figure3/metabolite_roc_10cv_250403.pdf',
       width = 16,
       height = 5)



# Odds ratio  ------------------------------------------------------------------
result_metabolite_marker_sig <- result_metabolite_marker %>% 
  filter(variable_id %in% table_auc_sig$variable_id) %>%
  arrange(match(variable_id, table_auc_sig$variable_id)) %>% 
  mutate(idx = seq(n()))

writexl::write_xlsx(result_metabolite_marker_sig, 
                    path = '~/Project/00_IBD_project/Data/20260420_single_met_logistic_evaluation/result_metabolite_marker_sig_260423.xlsx')

save(result_metabolite_marker_sig, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/result_metabolite_marker_sig_260423.RData')

table_auc_sig
writexl::write_xlsx(table_auc_sig, 
                    path = '~/Project/00_IBD_project/Data/20260420_single_met_logistic_evaluation/result_table_auc_sig_260423.xlsx')

temp_plot <- ggplot(result_metabolite_marker_sig, aes(x = idx, y = log10(odds_ratio))) +
  geom_pointrange(aes(ymin = log10(odds_ci_min), ymax = log10(odds_ci_max), colour = color)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(result_metabolite_marker_sig$variable_id)), 
                     label = result_metabolite_marker_sig$Compound.name
  ) +
  # ylim(-2.5, 2.5) +
  scale_y_continuous(breaks = c(-1, 0, 1), label = c(0.1, 1, 10), limits = c(-1.3, 2.1)) +
  # coord_flip() +
  scale_colour_manual(values = c('protective factor' = '#7fb1d3',
                                 'risk factor' = '#fc8070'),
                      label = c('protective factor' = 'Protective factor',
                                'risk factor' = 'Risk factor')) +
  xlab('metabolite marker') +
  ylab('Odds Ratio (OR)') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

ggsave(plot = temp_plot,
       filename = '~/Project/00_IBD_project/Figure/250326/Figure3/metabolite_odds_250403.pdf',
       width = 16, height = 5)

met_table <- table_auc_sig %>% 
  left_join(temp_metabolite_table) %>% 
  select(compound_name, variable_id, mean_odds_ratio:mean_odds_ci_max, auc:auc_ci_max)

# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = FALSE, showWarnings = TRUE)
# readr::write_csv(met_table,
#                  path = '~/Project/00_IBD_project/Data/20260718_source_data/ext_figure7_individual_metabolite_marker_predictive_model_260720.csv')

rm(list = ls());gc()






################################################################################
# alignment plot ---------------------------------------------------------------

library(pROC)
library(tidyverse)


# serology data
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/result_serology_10_fold_cv_260420.RData')

temp_serology_table <- result_serology_10_fold_cv %>% 
  group_by(variable_id) %>% 
  summarise(mean_odds_ratio = mean(odds_ratio),
            mean_odds_ci_min = mean(odds_ci_min),
            mean_odds_ci_max = mean(odds_ci_max),
            mean_auc_train = mean(auc_train),
            mean_auc_test = mean(auc_test),
            mean_auc_ci_train_min = mean(auc_ci_train_min),
            mean_auc_ci_train_max = mean(auc_ci_train_max),
            mean_auc_ci_test_min = mean(auc_ci_test_min),
            mean_auc_ci_test_max = mean(auc_ci_test_max))

result_serology_marker <- temp_serology_table %>%
  dplyr::mutate(idx = seq(n())) %>%
  dplyr::rename('odds_ratio' = 'mean_odds_ratio',
                'odds_ci_min' = 'mean_odds_ci_min',
                'odds_ci_max' = 'mean_odds_ci_max') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_serology_260420.RData')

mean_auc <- table_auc_serology %>% 
  summarise(mean_auc = mean(auc),
            mean_auc_ci_min = mean(auc_ci_min),
            mean_auc_ci_max = mean(auc_ci_max)) %>% 
  pull(mean_auc)

# metabolite data
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/result_metabolite_marker_sig_260423.RData')

load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_sig_260423.RData')

# merge AUC figures

merge_point_size <- 1.3 * 0.6

temp_plot_auc_serology_merge <- ggplot(table_auc_serology, aes(x = idx, y = auc)) +
  geom_pointrange(aes(ymin = auc_ci_min, ymax = auc_ci_max), colour = 'black', size = merge_point_size) +
  geom_hline(yintercept = 0.7376, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(table_auc_serology$serology_name)),
                     label = table_auc_serology$serology_name) +
  scale_y_continuous(limits = c(0.6, 0.9), breaks = c(0.6, 0.7, 0.8, 0.9)) +
  xlab('Serological marker') +
  ylab('AUC') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

temp_plot_auc_met_merge <- ggplot(table_auc_sig, aes(x = idx, y = auc)) +
  geom_pointrange(aes(ymin = auc_ci_min, ymax = auc_ci_max), colour = 'black', size = merge_point_size) +
  geom_hline(yintercept = 0.7376, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(table_auc_sig$variable_id)),
                     label = table_auc_sig$compound_name,
                     limits = c(0.5, length(result_metabolite_marker_sig$variable_id) + 0.5)) +
  scale_y_continuous(limits = c(0.6, 0.9), breaks = c(0.6, 0.7, 0.8, 0.9)) +
  xlab('Metabolite marker') +
  ylab('AUC') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))

temp_plot_auc_merge <- cowplot::plot_grid(
  temp_plot_auc_serology_merge,
  temp_plot_auc_met_merge,
  nrow = 1,
  align = 'hv',
  axis = 'tblr',
  rel_widths = c(0.2, 0.8)
)

ggsave(plot = temp_plot_auc_merge,
       filename = '~/Project/00_IBD_project/Data/20260417_AUC_Odds_serology_met_individual/serology_metabolite_auc_merge_20_80_260423.pdf',
       width = 20,
       height = 6)



# Odds ratio merge

merge_point_size_odds <- 1.3 * 0.6

temp_plot_odds_serology_merge <- ggplot(result_serology_marker, aes(x = rev(idx), y = log10(odds_ratio), color = color)) +
  geom_pointrange(aes(ymin = log10(odds_ci_min), ymax = log10(odds_ci_max)), size = merge_point_size_odds) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(result_serology_marker$variable_id)),
                     label = c('OmpC', 'IgG ASCA', 'IgA ASCA', 'I2', 'GM-CSF', 'CBir Fla', 'ANCA')) +
  scale_colour_manual(values = c('protective factor' = '#7fb1d3',
                                 'risk factor' = '#fc8070'),
                      label = c('protective factor' = 'Protective factor',
                                'risk factor' = 'Risk factor')) +
  scale_y_continuous(breaks = c(-1, 0, 1), label = c(0.1, 1, 10), limits = c(-1.4, 2.3)) +
  xlab('Serology marker') +
  ylab('Odds Ratio (OR)') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = 'none')

temp_plot_odds_met_merge <- ggplot(result_metabolite_marker_sig, aes(x = idx, y = log10(odds_ratio))) +
  geom_pointrange(aes(ymin = log10(odds_ci_min), ymax = log10(odds_ci_max), colour = color), size = merge_point_size_odds) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(length(result_metabolite_marker_sig$variable_id)),
                     label = result_metabolite_marker_sig$Compound.name,
                     limits = c(0.5, length(result_metabolite_marker_sig$variable_id) + 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1), label = c(0.1, 1, 10), limits = c(-1.4, 2.3)) +
  scale_colour_manual(values = c('protective factor' = '#7fb1d3',
                                 'risk factor' = '#fc8070'),
                      label = c('protective factor' = 'Protective factor',
                                'risk factor' = 'Risk factor')) +
  xlab('metabolite marker') +
  ylab('Odds Ratio (OR)') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = 'none')

temp_plot_odds_merge <- cowplot::plot_grid(
  temp_plot_odds_serology_merge,
  temp_plot_odds_met_merge,
  nrow = 1,
  align = 'hv',
  axis = 'tblr',
  rel_widths = c(0.2, 0.8)
)

ggsave(plot = temp_plot_odds_merge,
       filename = '~/Project/00_IBD_project/Data/20260417_AUC_Odds_serology_met_individual/serology_metabolite_odds_merge_20_80_260423.pdf',
       width = 20,
       height = 6)







################################################################################
# LASSO model to select the most important variables for prediction model construction --------
library(tidyverse)
library(tidymass)

load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_sig_260423.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_metabolite_260420.RData')


sig_metabolite_table <- table_auc_sig

sample_ids <- object_enrollment %>% extract_sample_info() %>% pull(sample_id)

set.seed(20250106)
sample_ids_discovery <- sample(sample_ids, length(sample_ids)*0.8)
sample_ids_validation <- sample_ids[!sample_ids %in% sample_ids_discovery]


# create train data set for metabolite-only LASSO model ------------------------
#   Keep the metabolites that have AUC larger than the mean AUC of serology markers

metabolite_data <- sample_info_enrollment_metabolite %>% 
  select(sample_id:disease, sig_metabolite_table$variable_id)

metabolite_data_discovery <- metabolite_data %>% 
  filter(sample_id %in% sample_ids_discovery) %>% 
  column_to_rownames(var = 'sample_id')

metabolite_data_validation <- metabolite_data %>%
  filter(sample_id %in% sample_ids_validation) %>%
  column_to_rownames(var = 'sample_id')

# calculate and savle the mean, sd for each metabolite in the discovery set, which will be used for scaling the validation set
metabolite_mean_sd <- metabolite_data_discovery %>% 
  dplyr::select(dplyr::any_of(sig_metabolite_table$variable_id)) %>% 
  summarise(across(everything(), list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))) %>% 
  pivot_longer(cols = everything(), names_to = c('variable_id', 'stat'), names_pattern = '^(.*)_(mean|sd)$') %>% 
  pivot_wider(names_from = stat, values_from = value)

save(metabolite_mean_sd, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_mean_sd_260423.RData')

# scale the discovery set based on the mean and sd of the discovery set
mean_vec <- setNames(metabolite_mean_sd$mean, metabolite_mean_sd$variable_id)
sd_vec   <- setNames(metabolite_mean_sd$sd, metabolite_mean_sd$variable_id)

metabolite_data_discovery_scaled <- metabolite_data_discovery %>%
  mutate(across(any_of(sig_metabolite_table$variable_id), ~{
    m <- mean_vec[cur_column()]
    s <- sd_vec[cur_column()]
    if (is.na(s) || s == 0) 0 else (. - m) / s
  }))

# scale the validation set using the mean and sd from the discovery set
metabolite_data_validation_scaled <- metabolite_data_validation %>%
  mutate(across(any_of(sig_metabolite_table$variable_id), ~{
    m <- mean_vec[cur_column()]
    s <- sd_vec[cur_column()]
    if (is.na(s) || s == 0) 0 else (. - m) / s
  }))

save(metabolite_data_discovery_scaled, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_data_discovery_scaled_260423.RData')

save(metabolite_data_validation_scaled, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_data_validation_scaled_260423.RData')

# LASSO approach first ----------------------------------------------------
# use the all features (q_value <= 0.05)

x_train <- metabolite_data_discovery_scaled %>% select(-disease) %>% as.matrix()
y_train <- metabolite_data_discovery_scaled %>% pull(disease)
names(y_train) <- rownames(metabolite_data_discovery_scaled)

# use all significant variables to train the model 
library(glmnet)

set.seed(241113)
seed_idx <- sample(10000000, 100)


lasso_table <- pbapply::pblapply(seq_along(seed_idx), function(i){
  x <- seed_idx[i]
  set.seed(x)
  lasso_model <- cv.glmnet(x = x_train, 
                           y = y_train, 
                           family="binomial", 
                           type.measure="class", 
                           alpha = 1, 
                           nfolds = 10)
  # plot(lasso_model)
  
  # Extract coefficients for selected variables
  optimal_lambda <- lasso_model$lambda.min
  selected_coeffs <- coef(lasso_model, s = optimal_lambda)
  ids_selected <- selected_coeffs@Dimnames[[1]][selected_coeffs@i + 1]
  annot_table <- object_enrollment %>% 
    extract_annotation_table()
  
  selected_met_id <- ids_selected[stringr::str_detect(ids_selected, 'M\\d+T\\d+')]
  selected_other_var <- ids_selected[!stringr::str_detect(ids_selected, 'M\\d+T\\d+')]
  selected_met <- match(selected_met_id, annot_table$variable_id)  %>% 
    annot_table$Compound.name[.]
  selected_met_confidence <- match(selected_met_id, annot_table$variable_id)  %>% 
    annot_table$confidence_level[.]
  selected_met_metabolon_subclass <- match(selected_met_id, annot_table$variable_id)  %>% 
    annot_table$metabolon_subclass[.]
  selected_met_table <- tibble(met = selected_met, 
                               confidence = selected_met_confidence,
                               metabolon_subclass = selected_met_metabolon_subclass,
                               id = selected_met_id,
                               efficient = selected_coeffs@x[stringr::str_detect(ids_selected, 'M\\d+T\\d+')])  
  
  selected_met_table2 <- selected_met_table %>% 
    filter(confidence == 'Level1') %>% 
    arrange(desc(abs(efficient))) %>% 
    mutate(time = i) %>% 
    mutate(rank = seq(n()))
  
  list_result <- list(selected_met_table2 = selected_met_table2,
                      selected_coeffs = selected_coeffs)
  
  return(list_result)
})

save(lasso_table,
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/lasso_metabolite_table_260423.RData')

annot_table <- object_enrollment %>% 
  extract_annotation_table() %>% 
  select(variable_id, Compound.name, confidence_level, metabolon_subclass)

stat_variable <- lapply(lasso_table, function(x){
  selected_coeffs <- x[[2]]
  selected_coeffs@Dimnames[[1]][selected_coeffs@i + 1]
}) %>% unlist() %>% table()

stat_variable_table <- tibble(variable_id = names(stat_variable),
                              n = as.numeric(stat_variable)) %>% 
  left_join(annot_table, by = 'variable_id') %>% 
  filter(variable_id != '(Intercept)') %>% 
  arrange(desc(n)) %>% 
  select(variable_id, Compound.name, confidence_level, metabolon_subclass, n) %>% 
  filter(confidence_level == 'Level1' | is.na(confidence_level)) %>% 
  mutate(variable_type = case_when(
    is.na(confidence_level) ~ 'demographic_variable',
    confidence_level == 'Level1' ~ 'metabolite_variable'
  ))



# keep variable repeated more than 90 times
stat_variable_table %>% 
  filter(n > 90) %>% 
  pull(variable_id) %>% 
  length()

stat_variable_table <- stat_variable_table %>% 
  mutate(idx = seq(n())) %>% 
  mutate(x_lab_name = case_when(
    is.na(confidence_level) ~ variable_id,
    !is.na(confidence_level) ~ Compound.name
  ))


save(stat_variable_table,
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/lasso_metabolite_table_260423.RData')

optimized_variable_table <- stat_variable_table %>% 
  filter(n > 90)
save(optimized_variable_table,
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/optimized_variable_table_260423.RData')

temp_plot <- ggplot(stat_variable_table, aes(x = idx, y = n)) +
  geom_bar(aes(fill = metabolon_subclass), stat = 'identity', position = 'dodge') +
  geom_hline(yintercept = 90, linetype = 'dashed') +
  scale_fill_manual(values = c("Urea cycle; Arginine and Proline Metabolism" = "#1b9e77",
                               "Methionine, Cysteine, SAM and Taurine Metabolism" = "#d95f02",
                               "Tryptophan Metabolism" = "#7570b3", 
                               "Long-chain fatty acids" = "#e7298a",
                               "Pyrimidine Metabolism, Cytosine containing" = "#66a61e",
                               "Steroids" = "#e6ab02",
                               "Hemoglobin and Porphyrin Metabolism" = "#a6761d",
                               "Phenylalanine and Tyrosine Metabolism" = "#666666",
                               "Very long-chain fatty acids" = "#4575b4",
                               "Bacterial/Fungal" = "#1a9850",
                               "Lineolic Acid Metabolism" = "#4d4d4d")) +
  # scale_fill_manual(values = c('demographic_variable' = '#b5de68',
  #                              'metabolite_variable' = '#beb9da'),
  #                   label = c('demographic_variable' = 'Demographics',
  #                             'metabolite_variable' = 'Metabolite')) +
  scale_x_continuous(breaks = stat_variable_table$idx, 
                     label = stat_variable_table$x_lab_name) +
  ylab('Frequency') +
  xlab('Variables') +
  ZZWtool::ZZWTheme() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        legend.position = c(0.9, 0.8))

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/lasso_variable_selection_260423.pdf',
       width = 8, height = 6)

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)

readr::write_csv(stat_variable_table, 
                 path = '~/Project/00_IBD_project/Data/ext_figure7_lasso_variable_selection_260720.csv')



rm(list = ls());gc()



################################################################################
# Multi-variate logistic regression model for the optimized metabolites --------------------------------------

# summarize the odds results
summarize_odds <- function(logistic_model) {
  CI_lower <- coefficients(logistic_model) - 1.96*summary(logistic_model)$coefficients[,2]
  CI_upper <- coefficients(logistic_model) + 1.96*summary(logistic_model)$coefficients[,2]
  odds <- exp(coef(logistic_model))
  odds_ci_min <- exp(CI_lower)
  odds_ci_max <- exp(CI_upper)
  # odds_ci <- paste(exp(CI_lower), exp(CI_upper), sep = ';')
  result_table <- tibble::tibble(variable = names(coef(logistic_model)),
                                 odds = odds, 
                                 odds_ci_min = odds_ci_min, 
                                 odds_ci_max = odds_ci_max)
}



load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_sig_260423.RData')
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_metabolite_260420.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/optimized_variable_table_260423.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_mean_sd_260423.RData')


sample_ids <- object_enrollment %>% extract_sample_info() %>% pull(sample_id)

set.seed(20250106)
sample_ids_discovery <- sample(sample_ids, length(sample_ids)*0.8)
sample_ids_validation <- sample_ids[!sample_ids %in% sample_ids_discovery]

# create train data set for metabolite-only LASSO model ------------------------
#   Keep the metabolites that have AUC larger than the mean AUC of serology markers

metabolite_data <- sample_info_enrollment_metabolite %>% 
  select(sample_id, disease, age:use_antibiotics, optimized_variable_table$variable_id)

metabolite_data_discovery <- metabolite_data %>% 
  filter(sample_id %in% sample_ids_discovery) %>% 
  column_to_rownames(var = 'sample_id')

metabolite_data_validation <- metabolite_data %>% 
  filter(sample_id %in% sample_ids_validation) %>% 
  column_to_rownames(var = 'sample_id')

age_mean <- mean(metabolite_data_discovery$age, na.rm = TRUE)
age_sd <- sd(metabolite_data_discovery$age, na.rm = TRUE)

# use the saved metabolite_mean_sd to scale the discovery set and validation set
mean_vec <- setNames(metabolite_mean_sd$mean, metabolite_mean_sd$variable_id)
sd_vec   <- setNames(metabolite_mean_sd$sd, metabolite_mean_sd$variable_id)

# only scale variables present in both optimized list and mean/sd table
scale_vars <- intersect(optimized_variable_table$variable_id, names(mean_vec))

metabolite_data_discovery <- metabolite_data_discovery %>%
  mutate(across(all_of(scale_vars), ~{
    m <- mean_vec[cur_column()]
    s <- sd_vec[cur_column()]
    if (is.na(m) || is.na(s) || s == 0) 0 else (.x - m) / s
  }))

metabolite_data_validation <- metabolite_data_validation %>%
  mutate(across(all_of(scale_vars), ~{
    m <- mean_vec[cur_column()]
    s <- sd_vec[cur_column()]
    if (is.na(m) || is.na(s) || s == 0) 0 else (.x - m) / s
  }))

# scale age as well
if (is.na(age_sd) || age_sd == 0) {
  metabolite_data_discovery <- metabolite_data_discovery %>%
    mutate(age = 0)
  metabolite_data_validation <- metabolite_data_validation %>%
    mutate(age = 0)
} else {
  metabolite_data_discovery <- metabolite_data_discovery %>%
    mutate(age = (age - age_mean) / age_sd)
  metabolite_data_validation <- metabolite_data_validation %>%
    mutate(age = (age - age_mean) / age_sd)
}

save(metabolite_data_discovery, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_data_discovery_scaled_260423.RData')

save(metabolite_data_validation, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/metabolite_data_validation_scaled_260423.RData')


temp_formula <- paste0('disease ~ ', paste(colnames(metabolite_data_discovery)[-1], collapse = '+'))

# logistic regression - training set
logistic_model <- glm(as.formula(temp_formula), data = metabolite_data_discovery, family = binomial(link = 'logit'))

logistic_result <- summarize_odds(logistic_model)

save(logistic_result, 
     file = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/logistic_result_260423.RData')


possibility_discovery <- predict(logistic_model, newdata = metabolite_data_discovery, type = "response")

possibility_validation <- predict(logistic_model, newdata = metabolite_data_validation, type = "response")

roc_obj_discovery <- roc(metabolite_data_discovery$disease ~ possibility_discovery, plot = FALSE, print.auc = TRUE, ci = TRUE)

# best Coordinates
coords(roc_obj_discovery, "best", 
       ret = c("threshold", "sensitivity", "specificity"), 
       best.method = "youden")


roc_obj_validation <- roc(metabolite_data_validation$disease ~ possibility_validation, plot = FALSE, print.auc = TRUE, ci = TRUE)

coords(roc_obj_validation, "best", 
       ret = c("threshold", "sensitivity", "specificity"), 
       best.method = "youden")

# visualize the logistic model -------------------------------------------------

annot_table <- object_enrollment %>% 
  extract_annotation_table()

logistic_result_met <- logistic_result %>% 
  left_join(annot_table, by = c('variable' = 'variable_id')) %>% 
  select(variable:odds_ci_max, Compound.name, confidence_level, metabolon_class, metabolon_subclass) %>% 
  filter(!is.na(confidence_level)) %>%
  arrange(metabolon_subclass) %>% 
  dplyr::mutate(idx = seq(n())) %>%
  dplyr::mutate(color = case_when(odds > 1 ~ 'risk factor',
                                  odds <= 1 ~ 'protective factor')) %>% 
  dplyr::mutate(y = 0.1)


temp_plot1 <- ggplot(logistic_result_met, aes(x = idx, y = log10(odds))) +
  geom_pointrange(aes(ymin = log10(odds_ci_min), ymax = log10(odds_ci_max), colour = color)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(breaks = logistic_result_met$idx, 
                     label = logistic_result_met$Compound.name
  ) +
  # ylim(-2.5, 2.5) +
  scale_y_continuous(breaks = c(-1, 0, 1), label = c(0.1, 1, 10), limits = c(-1, 1)) +
  coord_flip() +
  scale_colour_manual(values = c('protective factor' = '#7fb1d3',
                                 'risk factor' = '#fc8070'),
                      label = c('protective factor' = 'Protective factor',
                                'risk factor' = 'Risk factor')) +
  xlab('metabolite marker') +
  ylab('Odds Ratio (OR)') +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 0),
        axis.ticks.length.x.bottom = unit(1.5, 'mm'),
        axis.ticks.length.y.left = unit(1.5, 'mm'),
        legend.position = c(0.85, 0.85))


temp_plot2 <- ggplot(logistic_result_met, aes(x = y, y = idx-0.5, fill = metabolon_subclass)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_manual(values = c("Urea cycle; Arginine and Proline Metabolism" = '#8cd4c8',
                               "Methionine, Cysteine, SAM and Taurine Metabolism" = '#fdb562',
                               "Tryptophan Metabolism" = '#ffed6f',
                               "Long-chain fatty acids" = '#b5de68',
                               "Pyrimidine Metabolism, Cytosine containing" = '#beb9da',
                               "Very long-chain fatty acids" = '#fbcee5',
                               "Steroids" = '#fc8070',
                               "Hemoglobin and Porphyrin Metabolism" = '#bd80bd',
                               "Phenylalanine and Tyrosine Metabolism" = '#7fb1d3')) +
  coord_cartesian(xlim = c(0,13)) +
  scale_y_continuous(expand = c(0,0))+
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none') 


temp_plot3 <- ggplot(logistic_result_met, aes(x = y, y = idx-0.5, fill = metabolon_class)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_manual(values = c("Cofactors and Vitamins" = '#8cd4c8',
                               "Lipids" = '#ffed6f',
                               "Amino acid and derivatives" = '#fbcee5',
                               "Nucleotide" = '#beb9da')) +
  coord_cartesian(xlim = c(0,13)) +
  scale_y_continuous(expand = c(0,0))+
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none') 


library(cowplot)
library(ggtree)

# pp <- list(temp_plot1, temp_plot2)
# temp_plot_merge <- plot_grid(plotlist=pp, ncol=2, align ='h', rel_widths = c(0.8, 0.1))
# 
# ggsave(temp_plot_merge, 
#        filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/logistic_model_odds_met_260423.pdf', 
#        width = 6, height = 6)


pp2 <- list(temp_plot1, temp_plot3)
temp_plot_merge <- plot_grid(plotlist=pp2, ncol=2, align ='h', rel_widths = c(0.8, 0.1))

ggsave(temp_plot_merge, 
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/logistic_model_odds_met_2_260423.pdf', 
       width = 6, height = 6)



# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)
# readr::write_csv(logistic_result_met,
#                  path = '~/Project/00_IBD_project/Data/20260718_source_data/figure3_logistic_model_odds_metabolite_260720.csv')



# plot ROC curve ---------------------------------------------------------------


# discovery set
best_threshold <- coords(roc_obj_discovery, 'best')

ciobj <- ci.se(roc_obj_discovery, specificities=seq(0, 1, l=25))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

temp_plot <- ggroc(roc_obj_discovery, size = 1) +
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper),
              stat = "identity", fill = "#7fb1d3", alpha= 0.3) +
  coord_equal() +
  geom_point(aes(x = specificity, y = sensitivity), data = best_threshold, size = 3) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed", size = 1) +
  xlab('Specificity') +
  ylab('Sensitivity') +
  ZZWtool::ZZW_annotate_text2(x = 0.2, y = 0.8, label = paste0('AUC = ', round(roc_obj_discovery$auc[[1]], 4), 
                                                               '(', round(roc_obj_discovery$ci[[1]], 4), ' - ', round(roc_obj_discovery$ci[[3]], 4), ')')) +
  ZZWtool::ZZWTheme() +
  theme(legend.position = 'right')

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/roc_discovery_metabolite_260423.pdf',
       width = 6, height = 6)


# validation set
best_threshold <- coords(roc_obj_validation, 'best')

ciobj <- ci.se(roc_obj_validation, specificities=seq(0, 1, l=25))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

temp_plot <- ggroc(roc_obj_validation, size = 1) +
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper),
              stat = "identity", fill = "#fc8070", alpha= 0.3) +
  coord_equal() +
  geom_point(aes(x = specificity, y = sensitivity), data = best_threshold, size = 3) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed", size = 1) +
  xlab('Specificity') +
  ylab('Sensitivity') +
  ZZWtool::ZZW_annotate_text2(x = 0.2, y = 0.8, label = paste0('AUC = ', round(roc_obj_validation$auc[[1]], 4), 
                                                               '(', round(roc_obj_validation$ci[[1]], 4), ' - ', round(roc_obj_validation$ci[[3]], 4), ')')) +
  ZZWtool::ZZWTheme() +
  theme(legend.position = 'right')

ggsave(temp_plot, 
       filename = '~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/roc_validation_metabolite_260423.pdf',
       width = 6, height = 6)



# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', recursive = TRUE, showWarnings = FALSE)
# temp_data_discovery <- data.frame(sample_id = rownames(metabolite_data_discovery),
#                                   disease = metabolite_data_discovery$disease,
#                                   possibility = possibility_discovery, 
#                                   label = "discovery_cohort")
# 
# temp_data_validation <- data.frame(sample_id = rownames(metabolite_data_validation),
#                                    disease = metabolite_data_validation$disease,
#                                    possibility = possibility_validation, 
#                                    label = "test_cohort")
# 
# readr::write_csv(temp_data_discovery, 
#                  path = '~/Project/00_IBD_project/Data/20260718_source_data/figure3_roc_discovery_metabolite_260720.csv')
# 
# readr::write_csv(temp_data_validation,
#                  path = '~/Project/00_IBD_project/Data/20260718_source_data/figure3_roc_validation_metabolite_260720.csv')





################################################################################
# Compare the prediction performance of serology and metabolite models ---------
# 100 times cross validation for serology, metabolite and 16s data

library(tidyverse)
library(tidymass)
library(sjmisc)
library(ggraph)
library(glmnet)
library(pROC)
library(ggpubr)

setwd('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/')


# Part 1: Serology markers - 100 times 10-fold cross validation ----------------

# Imputation for missing values in serology variables before modeling
library(VIM)

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_info_enrollemnt_serology_scaled_260420.RData')

# KNN imputation for missing values in serology variables before modeling
serology_variables <- colnames(sample_info_model_data)[-c(1:5)]

sample_info_model_data_imputed <- sample_info_model_data %>%
  tibble::rownames_to_column(var = 'sample_id')

if (anyNA(sample_info_model_data_imputed[, serology_variables])) {
  set.seed(20260423) # for reproducibility of KNN imputation
  sample_info_model_data_imputed <- VIM::kNN(sample_info_model_data_imputed,
                                             variable = serology_variables,
                                             k = 5,
                                             imp_var = FALSE)
}

sample_info_model_data_imputed <- sample_info_model_data_imputed %>%
  tibble::column_to_rownames(var = 'sample_id')

save(sample_info_model_data_imputed, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/sample_info_model_data_imputed_260423.RData')

# 10-fold cross validation for serology markers

sample_ids <- rownames(sample_info_model_data_imputed)
serology_variables <- colnames(sample_info_model_data_imputed)[-c(1:5)]

# generate seed list to reproduce the results
set.seed(20260423)
index <- sample(10000000, 100)

list_multi_serology_10_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index[z]
  
  set.seed(seed_id)
  
  # Stratified Cross-Validation
  # obtain the patient IDs for CD (disease = 1) and non_IBD (disease = 0)
  id_case <- rownames(sample_info_model_data_imputed)[sample_info_model_data_imputed$disease == 1]
  id_control <- rownames(sample_info_model_data_imputed)[sample_info_model_data_imputed$disease == 0]
  
  # randomly shuffle the patient IDs within each group to ensure randomness in fold assignment
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  # divide the patient IDs into 10 folds for each group separately
  groups_case <- rep(1:10, length.out = length(id_case))
  groups_control <- rep(1:10, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  # combine the folds from both groups to create the final list of patient IDs for each fold
  cv_splits <- lapply(1:10, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_10_fold_cv <- lapply(1:10, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    # split data sets into training and testing sets directly
    data_train <- sample_info_model_data_imputed[ids_training, ]
    data_test <- sample_info_model_data_imputed[ids_test, ]
    
    # Data Leakage Fix: Compute scaling parameters on the training set ONLY
    # Scale age
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    serology_names <- colnames(data_train)[-c(1:5)] # cols 1-5: disease, age, gender, race, use_antibiotics
    
    # Scale serology markers (based on training set)
    for(marker in serology_names){
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>% 
      as_tibble() %>% 
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(serology_names))
    
    temp_data_test <- data_test %>% 
      as_tibble() %>% 
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(serology_names))
    
    # logistic regression with all serology markers - training set
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', serology_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(serology_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    # ROC curve --- training set
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    # ROC curve --- testing set
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = serology_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_serology_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    result <- list(stat_result = stat_values,
                   data_train = temp_data_train,
                   data_test = temp_data_test,
                   logistic_model = logistic_model,
                   roc_train = roc_obj_train,
                   roc_test = roc_obj_test)
    
    return(result)
  })
  
  names(result_list_10_fold_cv) <- paste0('fold_', 1:10)
  
  return(result_list_10_fold_cv)
})

names(list_multi_serology_10_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_serology_10_fold_cv, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_serology_10_fold_cv_260423.RData')

result_multi_serology_10_fold_cv <- lapply(list_multi_serology_10_fold_cv, function(z){
  z <- lapply(z, function(x){
    x$stat_result
  }) %>% 
    bind_rows(.id = 'fold')
}) %>% 
  bind_rows(.id = 'times')

# split odds_ci, auc_ci_train, and auc_ci_test into two columns
result_multi_serology_10_fold_cv <- result_multi_serology_10_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>% 
  dplyr::rename('odds_ratio' = 'odds') %>% 
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_serology_10_fold_cv, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_serology_10_fold_cv_260423.RData')


# ROC analysis for 100-times CV (test set)

cv_test_predictions_multi_serology <- lapply(names(list_multi_serology_10_fold_cv), function(rep_id) {
  fold_list <- list_multi_serology_10_fold_cv[[rep_id]]
  
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_serology <- cv_test_predictions_multi_serology %>% 
  dplyr::group_by(times) %>% 
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% 
  dplyr::ungroup()

roc_overall_100cv_multi_serology <- pROC::roc(response = cv_test_predictions_multi_serology$disease,
                                              predictor = cv_test_predictions_multi_serology$prediction,
                                              levels = c(0, 1),
                                              direction = '<',
                                              ci = TRUE,
                                              quiet = TRUE)

roc_overall_summary_100cv_multi_serology <- tibble::tibble(
  auc = as.numeric(roc_overall_100cv_multi_serology$auc),
  auc_ci_min = as.numeric(roc_overall_100cv_multi_serology$ci[[1]]),
  auc_ci_max = as.numeric(roc_overall_100cv_multi_serology$ci[[3]])
)

# Interpolate repeat-level ROC to a common FPR grid and compute mean ROC
fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_serology <- roc_by_repeat_multi_serology %>% 
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x,
                                 x = 'all',
                                 input = 'threshold',
                                 ret = c('specificity', 'sensitivity'),
                                 transpose = FALSE) %>% 
      as_tibble()
    
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    
    tpr_interp <- approx(x = fpr,
                         y = tpr,
                         xout = fpr_grid,
                         method = 'linear',
                         ties = 'ordered',
                         rule = 2)$y
    
    tibble::tibble(fpr = fpr_grid,
                   tpr = tpr_interp)
  })) %>% 
  dplyr::select(times, curve) %>% 
  tidyr::unnest(curve) %>% 
  dplyr::group_by(fpr) %>% 
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_serology <- ggplot(roc_curve_mean_100cv_multi_serology,
                                             aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max),
              fill = '#2C7FB8', alpha = 0.25) +
  geom_line(color = '#08519C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_serology$auc, 3),
                         ' (95% CI: ',
                         round(roc_overall_summary_100cv_multi_serology$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_serology$auc_ci_max, 3), ')'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')

save(cv_test_predictions_multi_serology,
     roc_by_repeat_multi_serology,
     roc_overall_100cv_multi_serology,
     roc_overall_summary_100cv_multi_serology,
     roc_curve_mean_100cv_multi_serology,
     plot_roc_mean_100cv_multi_serology,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_serology_100cv_260423.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_serology_100cv_260423.pdf',
       plot = plot_roc_mean_100cv_multi_serology,
       width = 6,
       height = 5)




# Part 2: Metabolite markers - 100 times 10-fold cross validation --------------

load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/table_auc_sig_260423.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_metabolite_260420.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/optimized_variable_table_260423.RData')

load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/optimized_variable_table_260720.RData')

metabolite_variables <- unique(optimized_variable_table$variable_id)
metabolite_variables <- metabolite_variables[metabolite_variables %in% colnames(sample_info_enrollment_metabolite)]

if (length(metabolite_variables) == 0) {
  stop('No optimized metabolite variables found in sample_info_enrollment_metabolite.')
}

metabolite_data <- sample_info_enrollment_metabolite %>% 
  dplyr::select(sample_id, disease, age, gender, race, use_antibiotics, dplyr::all_of(metabolite_variables))

# KNN imputation for missing values in metabolite variables before modeling
metabolite_data_imputed <- metabolite_data

if (anyNA(metabolite_data_imputed[, metabolite_variables])) {
  set.seed(20260423)
  metabolite_data_imputed <- VIM::kNN(metabolite_data_imputed,
                                      variable = metabolite_variables,
                                      k = 5,
                                      imp_var = FALSE)
}

save(metabolite_data_imputed,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/metabolite_data_imputed_260423.RData')

# 10-fold cross validation for optimized metabolite markers 

sample_ids_metabolite <- metabolite_data_imputed$sample_id

# generate seed list to reproduce the results
set.seed(20260423)
index <- sample(10000000, 100)

list_multi_metabolite_10_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization for metabolite model\n')
  seed_id <- index[z]
  
  set.seed(seed_id)
  
  # Stratified Cross-Validation by disease status
  id_case <- metabolite_data_imputed$sample_id[metabolite_data_imputed$disease == 1]
  id_control <- metabolite_data_imputed$sample_id[metabolite_data_imputed$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:10, length.out = length(id_case))
  groups_control <- rep(1:10, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:10, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_10_fold_cv <- lapply(1:10, function(j){
    cat('Perform the ', j, 'th cross validation for metabolite model\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- metabolite_data_imputed %>% 
      dplyr::filter(sample_id %in% ids_training)
    data_test <- metabolite_data_imputed %>% 
      dplyr::filter(sample_id %in% ids_test)
    
    # Data Leakage Fix: scale features using training-set parameters only
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for(marker in metabolite_variables){
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>% 
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(metabolite_variables))
    
    temp_data_test <- data_test %>% 
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(metabolite_variables))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', metabolite_variables),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(metabolite_variables, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    
    temp_roc_train <- tibble::tibble(disease = temp_data_train$disease,
                                     prediction = as.numeric(possibility_train)) %>% 
      dplyr::filter(!is.na(disease), !is.na(prediction))
    temp_roc_test <- tibble::tibble(disease = temp_data_test$disease,
                                    prediction = as.numeric(possibility_test)) %>% 
      dplyr::filter(!is.na(disease), !is.na(prediction))
    
    roc_obj_train <- pROC::roc(response = temp_roc_train$disease,
                               predictor = temp_roc_train$prediction,
                               levels = c(0, 1),
                               direction = '<',
                               ci = TRUE,
                               quiet = TRUE)
    auc_train <- as.numeric(roc_obj_train$auc)
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    roc_obj_test <- pROC::roc(response = temp_roc_test$disease,
                              predictor = temp_roc_test$prediction,
                              levels = c(0, 1),
                              direction = '<',
                              ci = TRUE,
                              quiet = TRUE)
    auc_test <- as.numeric(roc_obj_test$auc)
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = metabolite_variables[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_optimized_metabolite_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    result <- list(stat_result = stat_values,
                   data_train = temp_data_train,
                   data_test = temp_data_test,
                   logistic_model = logistic_model,
                   roc_train = roc_obj_train,
                   roc_test = roc_obj_test)
    
    return(result)
  })
  
  names(result_list_10_fold_cv) <- paste0('fold_', 1:10)
  return(result_list_10_fold_cv)
})

names(list_multi_metabolite_10_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_metabolite_10_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_metabolite_10_fold_cv_260423.RData')

result_multi_metabolite_10_fold_cv <- lapply(list_multi_metabolite_10_fold_cv, function(z){
  z <- lapply(z, function(x){
    x$stat_result
  }) %>% 
    bind_rows(.id = 'fold')
}) %>% 
  bind_rows(.id = 'times')

result_multi_metabolite_10_fold_cv <- result_multi_metabolite_10_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>% 
  dplyr::rename('odds_ratio' = 'odds') %>% 
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_metabolite_10_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_metabolite_10_fold_cv_260423.RData')

# ROC analysis for 100-times CV (test set) 

cv_test_predictions_multi_metabolite <- lapply(names(list_multi_metabolite_10_fold_cv), function(rep_id) {
  fold_list <- list_multi_metabolite_10_fold_cv[[rep_id]]
  
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows() %>% 
  dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_metabolite <- cv_test_predictions_multi_metabolite %>% 
  dplyr::group_by(times) %>% 
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% 
  dplyr::ungroup()

roc_overall_100cv_multi_metabolite <- pROC::roc(response = cv_test_predictions_multi_metabolite$disease,
                                                predictor = cv_test_predictions_multi_metabolite$prediction,
                                                levels = c(0, 1),
                                                direction = '<',
                                                ci = TRUE,
                                                quiet = TRUE)

roc_overall_summary_100cv_multi_metabolite <- tibble::tibble(
  auc = as.numeric(roc_overall_100cv_multi_metabolite$auc),
  auc_ci_min = as.numeric(roc_overall_100cv_multi_metabolite$ci[[1]]),
  auc_ci_max = as.numeric(roc_overall_100cv_multi_metabolite$ci[[3]])
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_metabolite <- roc_by_repeat_multi_metabolite %>% 
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x,
                                 x = 'all',
                                 input = 'threshold',
                                 ret = c('specificity', 'sensitivity'),
                                 transpose = FALSE) %>% 
      as_tibble()
    
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    
    tpr_interp <- approx(x = fpr,
                         y = tpr,
                         xout = fpr_grid,
                         method = 'linear',
                         ties = 'ordered',
                         rule = 2)$y
    
    tibble::tibble(fpr = fpr_grid,
                   tpr = tpr_interp)
  })) %>% 
  dplyr::select(times, curve) %>% 
  tidyr::unnest(curve) %>% 
  dplyr::group_by(fpr) %>% 
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_metabolite <- ggplot(roc_curve_mean_100cv_multi_metabolite,
                                               aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max),
              fill = '#31A354', alpha = 0.25) +
  geom_line(color = '#006D2C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Metabolite Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_metabolite$auc, 3),
                         ' (95% CI: ',
                         round(roc_overall_summary_100cv_multi_metabolite$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_metabolite$auc_ci_max, 3), ')'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')

save(cv_test_predictions_multi_metabolite,
     roc_by_repeat_multi_metabolite,
     roc_overall_100cv_multi_metabolite,
     roc_overall_summary_100cv_multi_metabolite,
     roc_curve_mean_100cv_multi_metabolite,
     plot_roc_mean_100cv_multi_metabolite,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_metabolite_100cv_260423.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_metabolite_100cv_260423.pdf',
       plot = plot_roc_mean_100cv_multi_metabolite,
       width = 6,
       height = 5)


# Part 3: Compare overall ROC between serology and metabolite ------------------

calc_metrics_at_best_threshold <- function(roc_obj, prediction_df, model_id) {
  best_cut <- pROC::coords(roc_obj,
                           x = 'best',
                           input = 'threshold',
                           best.method = 'youden',
                           ret = c('threshold', 'sensitivity', 'specificity'),
                           transpose = FALSE)
  
  threshold <- as.numeric(best_cut$threshold)
  sensitivity <- as.numeric(best_cut$sensitivity)
  specificity <- as.numeric(best_cut$specificity)
  
  pred_class <- ifelse(prediction_df$prediction >= threshold, 1, 0)
  true_class <- as.numeric(prediction_df$disease)
  
  tp <- sum(pred_class == 1 & true_class == 1, na.rm = TRUE)
  tn <- sum(pred_class == 0 & true_class == 0, na.rm = TRUE)
  fp <- sum(pred_class == 1 & true_class == 0, na.rm = TRUE)
  fn <- sum(pred_class == 0 & true_class == 1, na.rm = TRUE)
  
  n <- tp + tn + fp + fn
  accuracy <- ifelse(n > 0, (tp + tn) / n, NA_real_)
  ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
  f1_score <- ifelse((2 * tp + fp + fn) > 0, 2 * tp / (2 * tp + fp + fn), NA_real_)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  youden_index <- sensitivity + specificity - 1
  
  tibble::tibble(
    model = model_id,
    auc = as.numeric(roc_obj$auc),
    auc_ci_min = as.numeric(roc_obj$ci[[1]]),
    auc_ci_max = as.numeric(roc_obj$ci[[3]]),
    best_threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy,
    f1_score = f1_score,
    balanced_accuracy = balanced_accuracy,
    youden_index = youden_index,
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn
  )
}

best_threshold_metrics_serology <- calc_metrics_at_best_threshold(
  roc_obj = roc_overall_100cv_multi_serology,
  prediction_df = cv_test_predictions_multi_serology,
  model_id = 'Serology'
)

best_threshold_metrics_metabolite <- calc_metrics_at_best_threshold(
  roc_obj = roc_overall_100cv_multi_metabolite,
  prediction_df = cv_test_predictions_multi_metabolite,
  model_id = 'Metabolite'
)

best_threshold_metrics_compare <- dplyr::bind_rows(best_threshold_metrics_serology,
                                                   best_threshold_metrics_metabolite)

best_threshold_points_plot <- best_threshold_metrics_compare %>% 
  dplyr::mutate(fpr = 1 - specificity,
                tpr = sensitivity,
                threshold_label = paste0(model, ' threshold = ', round(best_threshold, 3)))



# overall ROC objects for combined plotting
roc_combined_overall <- list(
  Serology = roc_overall_100cv_multi_serology,
  Metabolite = roc_overall_100cv_multi_metabolite
)

auc_serology <- as.numeric(roc_overall_100cv_multi_serology$auc)
auc_metabolite <- as.numeric(roc_overall_100cv_multi_metabolite$auc)

best_threshold_metrics_compare %>% 
  column_to_rownames('model') %>%
  sjmisc::rotate_df()


# Merged ROC curve
plot_roc_overall_serology_metabolite <- pROC::ggroc(roc_combined_overall, size = 1.2, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  scale_color_manual(values = c(Serology = '#08519C', Metabolite = '#006D2C')) +
  geom_point(data = best_threshold_points_plot,
             aes(x = fpr, y = tpr, color = model),
             size = 3,
             inherit.aes = FALSE) +
  geom_text(data = best_threshold_points_plot,
            aes(x = fpr, 
                y = tpr, 
                label = paste0('AUC = ', round(auc, 3), ' [', round(auc_ci_min, 3), '-', round(auc_ci_max, 3), ']'),
                color = model),
            nudge_y = 0.04,
            size = 3,
            show.legend = FALSE,
            inherit.aes = FALSE) +
  ZZWtool::ZZWTheme() +
  labs(
    x = 'False Positive Rate',
    y = 'True Positive Rate',
    color = 'Model') +
  theme(legend.position = c(0.8, 0.2))


ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_overall_compare_serology_metabolite_260423.pdf',
       plot = plot_roc_overall_serology_metabolite,
       width = 6,
       height = 6)

save(plot_roc_overall_serology_metabolite,
     best_threshold_metrics_serology,
     best_threshold_metrics_metabolite,
     best_threshold_metrics_compare,
     best_threshold_points_plot,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_overall_compare_serology_metabolite_260423.RData')

readr::write_csv(best_threshold_metrics_compare,
                 file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/best_threshold_metrics_compare_serology_metabolite_260423.csv')


# ROC for combined plotting using mean ROC curve from 100-times CV

roc_curve_combined_mean_100cv <- dplyr::bind_rows(
  roc_curve_mean_100cv_multi_serology %>% 
    dplyr::transmute(model = 'Serology',
                     fpr = fpr,
                     mean_tpr = mean_tpr,
                     tpr_ci_min = tpr_ci_min,
                     tpr_ci_max = tpr_ci_max),
  roc_curve_mean_100cv_multi_metabolite %>% 
    dplyr::transmute(model = 'Metabolite',
                     fpr = fpr,
                     mean_tpr = mean_tpr,
                     tpr_ci_min = tpr_ci_min,
                     tpr_ci_max = tpr_ci_max)
)

best_threshold_metrics_for_mean_plot <- best_threshold_metrics_compare %>% 
  dplyr::mutate(
    fpr = 1 - specificity,
    tpr = sensitivity,
    metrics_label = paste0(
      model,
      '\nAUC = ', round(auc, 3), ' [', round(auc_ci_min, 3), '-', round(auc_ci_max, 3), ']',
      '\nCutoff = ', round(best_threshold, 3),
      '\nSens = ', round(sensitivity, 3), '; Spec = ', round(specificity, 3),
      '\nPPV = ', round(ppv, 3), '; NPV = ', round(npv, 3),
      '\nAcc = ', round(accuracy, 3), '; F1 = ', round(f1_score, 3),
      '\nBal Acc = ', round(balanced_accuracy, 3), '; Youden = ', round(youden_index, 3)
    )
  )

plot_roc_mean_compare_serology_metabolite <- ggplot(roc_curve_combined_mean_100cv,
                                                    aes(x = fpr, y = mean_tpr, color = model, fill = model)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max),
              alpha = 0.2,
              color = NA) +
  geom_line(size = 1.2) +
  geom_point(data = best_threshold_metrics_for_mean_plot,
             aes(x = fpr, y = tpr, color = model),
             size = 2.8,
             inherit.aes = FALSE) +
  geom_text(data = best_threshold_metrics_for_mean_plot,
            aes(x = fpr, y = tpr, label = metrics_label, color = model),
            nudge_y = 0.05,
            size = 3,
            lineheight = 0.95,
            show.legend = FALSE,
            inherit.aes = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  scale_color_manual(values = c(Serology = '#08519C', Metabolite = '#006D2C')) +
  scale_fill_manual(values = c(Serology = '#2C7FB8', Metabolite = '#31A354')) +
  ZZWtool::ZZWTheme() +
  labs(x = 'False Positive Rate',
       y = 'True Positive Rate',
       color = 'Model') +
  theme(legend.position = c(0.8, 0.2))

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_mean_compare_serology_metabolite_95ci_260511.pdf',
       plot = plot_roc_mean_compare_serology_metabolite,
       width = 6,
       height = 6)

readr::write_csv(best_threshold_metrics_for_mean_plot,
                 file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/best_threshold_metrics_for_mean_roc_260511.csv')

save(roc_curve_mean_100cv_multi_serology,
     roc_curve_mean_100cv_multi_metabolite,
     roc_curve_combined_mean_100cv,
     best_threshold_metrics_for_mean_plot,
     plot_roc_mean_compare_serology_metabolite,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_mean_compare_serology_metabolite_95ci_260423.RData')


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_source_data <- roc_curve_mean_100cv_multi_serology %>% 
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename(fdr_serology = fpr,
                tpr_serology = mean_tpr,
                tpr_ci_min_serology = tpr_ci_min,
                tpr_ci_max_serology = tpr_ci_max) %>% 
  left_join(roc_curve_mean_100cv_multi_metabolite %>% 
              dplyr::mutate(idx = seq_len(n())) %>%
              dplyr::rename(fdr_metabolite = fpr,
                            tpr_metabolite = mean_tpr,
                            tpr_ci_min_metabolite = tpr_ci_min,
                            tpr_ci_max_metabolite = tpr_ci_max),
            by = 'idx') %>% 
  dplyr::select(idx, everything())

readr::write_csv(best_threshold_metrics_compare,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig7d_source_data_roc_compare_serology_metabolite_260721.csv')






################################################################################
# Compare the prediction performance of 16S and metabolite models --------------

# Part 01: 16S rRNA data variable selection (LASSO) ---------------------------------------------------------

# Centered Log-Ratio transformation + Z-score scaling

library(tidyverse)
library(glmnet)
library(pROC)
library(sjmisc)

load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/seqtab_nochim_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/taxa_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/sample_table_251027.RData')

source('~/Project/00_IBD_project/Code/20251027_function_16S_rRNA_data_processing.R')

# aggregate to the genus level 

# otu table
temp_sample_name <- rownames(seqtab_nochim)
otu_table <- seqtab_nochim %>% 
  as_tibble() %>% 
  rotate_df()
colnames(otu_table) <- temp_sample_name
otu_table <- otu_table %>% 
  rownames_to_column('sequence')
rm(temp_sample_name);gc()

# taxa table
temp_asv <- rownames(taxa)
taxa_table <- taxa %>%
  as.data.frame() %>%
  rownames_to_column('sequence') %>%
  as_tibble()
rm(temp_asv);gc()

# aggregate to the genus levels
otu_table_genus <- aggregate_taxa(otu_table, taxa_table, level = 'Genus', relative = FALSE, remove_na = TRUE) %>% 
  rename('Taxa' = 'Genus') %>% 
  tidyr::replace_na(list(Taxa = 'Unassigned_genus'))

# apply CLR transformation first, then perform discovery/validation Z-scaling
# create the training and test data sets

sample_table <- sample_table %>% 
  dplyr::select(data_16s:visit_encounter_id, phenotype_type, age_at_encounter, gender, antibiotics, characteristics_bio_material) %>% 
  dplyr::mutate(disease = case_when(phenotype_type == 'CD' ~ 1,
                                    phenotype_type == 'non_IBD' ~ 0),
                gender = as.factor(gender),
                use_antibiotics = as.factor(antibiotics),
                characteristics_bio_material = as.factor(characteristics_bio_material)) %>% 
  rename(age = age_at_encounter) %>% 
  select(data_16s, disease)

temp_data <- otu_table_genus %>% 
  column_to_rownames('Taxa') %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('SampleID')

data_16s_genus <- sample_table %>% 
  left_join(temp_data, by = c('data_16s' = 'SampleID')) %>%
  column_to_rownames('data_16s') %>% 
  select(disease, everything())

# CLR transform on genus abundance matrix (row-wise) with pseudocount for zeros
genus_mat <- data_16s_genus %>%
  dplyr::select(-disease) %>%
  as.matrix()
storage.mode(genus_mat) <- 'numeric'

min_positive <- suppressWarnings(min(genus_mat[genus_mat > 0], na.rm = TRUE))
if (!is.finite(min_positive)) {
  min_positive <- 1e-06
}
pseudocount <- min_positive / 2

genus_mat_clr <- t(apply(genus_mat, 1, function(x) {
  log_x <- log(x + pseudocount)
  log_x - mean(log_x, na.rm = TRUE)
}))

colnames(genus_mat_clr) <- colnames(genus_mat)
rownames(genus_mat_clr) <- rownames(genus_mat)

data_16s_genus_clr <- data_16s_genus %>%
  dplyr::select(disease) %>%
  bind_cols(as.data.frame(genus_mat_clr, check.names = FALSE))
rownames(data_16s_genus_clr) <- rownames(genus_mat_clr)

save(data_16s_genus_clr, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_clr_260424.RData')

rm(temp_data);gc()

# use all samples for scaling and LASSO variable selection
data_16s_genus_all <- data_16s_genus_clr %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric))

genus_var_ids <- colnames(data_16s_genus_all)[colnames(data_16s_genus_all) != 'disease']

genus_mean_sd_clr_scaled <- data_16s_genus_all %>%
  dplyr::select(dplyr::any_of(genus_var_ids)) %>%
  summarise(across(everything(), list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))) %>%
  pivot_longer(cols = everything(), names_to = c('variable_id', 'stat'), names_pattern = '^(.*)_(mean|sd)$') %>%
  pivot_wider(names_from = stat, values_from = value)

mean_vec_16s <- setNames(genus_mean_sd_clr_scaled$mean, genus_mean_sd_clr_scaled$variable_id)
sd_vec_16s <- setNames(genus_mean_sd_clr_scaled$sd, genus_mean_sd_clr_scaled$variable_id)

data_16s_genus_all_clr_scaled <- data_16s_genus_all %>%
  mutate(across(any_of(genus_var_ids), ~{
    m <- mean_vec_16s[cur_column()]
    s <- sd_vec_16s[cur_column()]
    if (is.na(s) || s == 0) 0 else (. - m) / s
  }))

save(genus_mean_sd_clr_scaled,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/genus_mean_sd_clr_scaled_260424.RData')

save(data_16s_genus_all_clr_scaled,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_all_clr_scaled_260424.RData')



# LASSO regression 
x <- data_16s_genus_all_clr_scaled %>% dplyr::select(-disease) %>% as.matrix()
y <- data_16s_genus_all_clr_scaled %>% pull(disease)
names(y) <- rownames(data_16s_genus_all_clr_scaled)

set.seed(251029)
seed_idx_16s <- sample(10000000, 100)

lasso_var_result_genus_clr_scaled <- pbapply::pblapply(seq_along(seed_idx_16s), function(i){
  set.seed(seed_idx_16s[i])
  
  idx_training <- sample(seq_along(y), floor(0.8 * length(y)))
  x_train <- x[idx_training, , drop = FALSE]
  y_train <- y[idx_training]
  
  lasso_model <- cv.glmnet(x = x_train,
                           y = y_train,
                           family = 'binomial',
                           type.measure = 'auc',
                           alpha = 1,
                           nfolds = 10)
  
  optimal_lambda <- lasso_model$lambda.1se
  selected_coeffs <- coef(lasso_model, s = optimal_lambda)
  ids_selected <- selected_coeffs@Dimnames[[1]][selected_coeffs@i + 1]
  
  selected_genus_table <- tibble(
    variable_id = ids_selected,
    coefficient = selected_coeffs@x
  ) %>%
    filter(variable_id != '(Intercept)') %>%
    arrange(desc(abs(coefficient))) %>%
    mutate(time = i,
           rank = seq_len(n()))
  
  list_result <- list(selected_genus_table = selected_genus_table,
                      selected_coeffs = selected_coeffs)
  
  return(list_result)
})


save(lasso_var_result_genus_clr_scaled,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/lasso_var_result_genus_clr_scaled_all_samples_260424.RData')


# summarize the optimized variable results 
stat_var_genus_clr_scaled <- lasso_var_result_genus_clr_scaled %>%
  lapply(function(x){
    x$selected_genus_table
  }) %>%
  bind_rows() %>%
  group_by(variable_id) %>%
  summarise(
    Frequency = n(),
    Mean_Coefficient = mean(coefficient)
  ) %>%
  arrange(desc(Frequency), desc(abs(Mean_Coefficient)))

save(stat_var_genus_clr_scaled,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/stat_var_genus_clr_scaled_all_samples_260424.RData')

stat_var_genus_clr_scaled %>% filter(Frequency >= 95)


rm(list = ls());gc()

# Part 02: select the intersected samples ----------------------------------------------------------------

# 16s rRNA data
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/stat_var_genus_clr_scaled_all_samples_260424.RData')

optimized_variable_16s_clr_scaled <- stat_var_genus_clr_scaled %>% 
  filter(Frequency >= 95) %>% 
  pull(variable_id)

# 16s rRNA data
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_260424.RData')
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_clr_260424.RData')

load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/otu_table_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/sample_table_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/taxa_table_ileum_251103.RData')

load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/otu_table_rectum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/sample_table_rectum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/taxa_table_rectum_251103.RData')

load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/otu_table_stool_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/sample_table_stool_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/taxa_table_stool_251103.RData')


# metabolomics data
load('~/Project/00_IBD_project/Data/20250331_serology_marker_predictive_model/object_enrollment_250331.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/sample_expression_metabolite_260420.RData')
load('~/Project/00_IBD_project/Data/20260420_serology_marker_evaluation_Fig3/optimized_variable_table_260423.RData')

optimized_met_variable <- optimized_variable_table %>% 
  pull(variable_id)


# ileum data
data_met_ileum <- sample_info_enrollment_metabolite %>% 
  filter(sample_id %in% sample_table_ileum$sample_id) %>% 
  arrange(match(sample_id, sample_table_ileum$sample_id)) %>% 
  select(sample_id:disease, age:use_antibiotics, all_of(optimized_met_variable))

save(data_met_ileum, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_ileum_260427.RData')

data_16s_genus_ileum <- data_16s_genus_clr %>% 
  rownames_to_column('data_16s') %>% 
  filter(data_16s %in% sample_table_ileum$data_16s) %>%
  arrange(match(data_16s, sample_table_ileum$data_16s)) %>%
  select(data_16s, disease, all_of(optimized_variable_16s_clr_scaled)) %>% 
  left_join(sample_table_ileum %>% select(data_16s, age, gender, race, use_antibiotics), by = 'data_16s') %>% 
  select(data_16s, disease, age, gender, race, use_antibiotics, all_of(optimized_variable_16s_clr_scaled)) %>% 
  mutate(race = case_when(
    race %in% c('Caucasian', 'Black or African American', 'Asian', 'Hispanic/Latino') ~ race,
    TRUE ~ 'Others/Unknown'
  ),
  race = as.factor(race))

x_variable_id_ileum <- data.frame(variable_name = colnames(data_16s_genus_ileum)[-c(1:6)],
                                  variable_id = paste0('genus_', seq_along(colnames(data_16s_genus_ileum)[-c(1:6)])),
                                  stringsAsFactors = FALSE)
colnames(data_16s_genus_ileum)[-c(1:6)] <- x_variable_id_ileum$variable_id

save(data_16s_genus_ileum, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_ileum_all_samples_260427.RData')

save(x_variable_id_ileum, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/x_variable_id_ileum_genus_260427.RData')

rm(list = c('data_met_ileum', 'data_16s_genus_ileum', 'x_variable_id_ileum'));gc()


# rectum data 
data_met_rectum <- sample_info_enrollment_metabolite %>% 
  filter(sample_id %in% sample_table_rectum$sample_id) %>% 
  arrange(match(sample_id, sample_table_rectum$sample_id)) %>% 
  select(sample_id:disease, age:use_antibiotics, all_of(optimized_met_variable))

save(data_met_rectum, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_rectum_260427.RData')

data_16s_genus_rectum <- data_16s_genus_clr %>%
  rownames_to_column('data_16s') %>% 
  filter(data_16s %in% sample_table_rectum$data_16s) %>%
  arrange(match(data_16s, sample_table_rectum$data_16s)) %>%
  select(data_16s, disease, all_of(optimized_variable_16s_clr_scaled)) %>%
  left_join(sample_table_rectum %>% select(data_16s, age, gender, race, use_antibiotics), by = 'data_16s') %>%
  select(data_16s, disease, age, gender, race, use_antibiotics, all_of(optimized_variable_16s_clr_scaled)) %>% 
  mutate(race = case_when(
    race %in% c('Caucasian', 'Black or African American', 'Asian', 'Hispanic/Latino') ~ race,
    TRUE ~ 'Others/Unknown'
  ),
  race = as.factor(race))

x_variable_id_rectum <- data.frame(variable_name = colnames(data_16s_genus_rectum)[-c(1:6)],
                                   variable_id = paste0('genus_', seq_along(colnames(data_16s_genus_rectum)[-c(1:6)])),
                                   stringsAsFactors = FALSE)

colnames(data_16s_genus_rectum)[-c(1:6)] <- x_variable_id_rectum$variable_id
save(data_16s_genus_rectum, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_rectum_all_samples_260427.RData')

save(x_variable_id_rectum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/x_variable_id_rectum_genus_all_samples_260427.RData')


# stool data 
# Note: for stool samples, we assigned 'Hispanic/Latino' to 'Others/Unknown' due to the small sample size of Hispanic/Latino in stool samples (n=2), while we kept 'Hispanic/Latino' as a separate category in ileum and rectum samples since there are more Hispanic/Latino samples in those groups.

data_met_stool <- sample_info_enrollment_metabolite %>% 
  filter(sample_id %in% sample_table_stool$sample_id) %>% 
  arrange(match(sample_id, sample_table_stool$sample_id)) %>% 
  select(sample_id:disease, age:use_antibiotics, all_of(optimized_met_variable)) %>% 
  mutate(race = case_when(
    race %in% c('Caucasian', 'Black or African American', 'Asian') ~ race,
    TRUE ~ 'Others/Unknown'
  ),
  race = as.factor(race))

save(data_met_stool, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_stool_260427.RData')

data_16s_genus_stool <- data_16s_genus_clr %>%
  rownames_to_column('data_16s') %>% 
  filter(data_16s %in% sample_table_stool$data_16s) %>%
  arrange(match(data_16s, sample_table_stool$data_16s)) %>%
  select(data_16s, disease, all_of(optimized_variable_16s_clr_scaled)) %>%
  left_join(sample_table_stool %>% select(data_16s, age, gender, race, use_antibiotics), by = 'data_16s') %>%
  select(data_16s, disease, age, gender, race, use_antibiotics, all_of(optimized_variable_16s_clr_scaled)) %>% 
  mutate(race = case_when(
    race %in% c('Caucasian', 'Black or African American', 'Asian') ~ race,
    TRUE ~ 'Others/Unknown'
  ),
  race = as.factor(race))

x_variable_id_stool <- data.frame(variable_name = colnames(data_16s_genus_stool)[-c(1:6)],
                                  variable_id = paste0('genus_', seq_along(colnames(data_16s_genus_stool)[-c(1:6)])),
                                  stringsAsFactors = FALSE)

colnames(data_16s_genus_stool)[-c(1:6)] <- x_variable_id_stool$variable_id

save(data_16s_genus_stool, 
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_stool_all_samples_260427.RData')

save(x_variable_id_stool,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/x_variable_id_stool_genus_all_samples_260427.RData')

rm(list = ls());gc()



# Part 03: running 5-fold cross validation with 100 repeats ---------------------------

load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_ileum_260427.RData')
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_ileum_all_samples_260427.RData')


library(tidyverse)

# ileum ----------------------------------------------------------------------
# 5-fold cross validation for ileum 16S markers -------------------------------

sample_ids <- data_16s_genus_ileum$data_16s
covariate_names <- c('age', 'gender', 'race', 'use_antibiotics')
marker_16s_names <- setdiff(colnames(data_16s_genus_ileum), c('data_16s', 'disease', covariate_names))

# generate seed list to reproduce the results
set.seed(20260427)
index <- sample(10000000, 100)

list_multi_16s_ileum_5_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index[z]
  
  set.seed(seed_id)
  
  # Stratified cross-validation by disease label
  id_case <- data_16s_genus_ileum$data_16s[data_16s_genus_ileum$disease == 1]
  id_control <- data_16s_genus_ileum$data_16s[data_16s_genus_ileum$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_16s_genus_ileum %>%
      dplyr::filter(data_16s %in% ids_training)
    data_test <- data_16s_genus_ileum %>%
      dplyr::filter(data_16s %in% ids_test)
    
    # Keep deterministic sample order in each split
    data_train <- data_train %>%
      dplyr::arrange(match(data_16s, ids_training))
    data_test <- data_test %>%
      dplyr::arrange(match(data_16s, ids_test))
    
    # Scale age and 16S markers based on training set only
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_16s_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_16s_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_16s_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_16s_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_ileum_16s_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    result <- list(stat_result = stat_values,
                   data_train = temp_data_train,
                   data_test = temp_data_test,
                   logistic_model = logistic_model,
                   roc_train = roc_obj_train,
                   roc_test = roc_obj_test)
    
    return(result)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  return(result_list_5_fold_cv)
})

names(list_multi_16s_ileum_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_16s_ileum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_16s_ileum_5_fold_cv_all_samples_260428.RData')

result_multi_16s_ileum_5_fold_cv <- lapply(list_multi_16s_ileum_5_fold_cv, function(z){
  z <- lapply(z, function(x){
    x$stat_result
  }) %>%
    bind_rows(.id = 'fold')
}) %>%
  bind_rows(.id = 'times')

result_multi_16s_ileum_5_fold_cv <- result_multi_16s_ileum_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_16s_ileum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_16s_ileum_5_fold_cv_all_samples_260428.RData')

# ROC analysis for 100-times CV 

cv_test_predictions_multi_16s_ileum <- lapply(names(list_multi_16s_ileum_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_16s_ileum_5_fold_cv[[rep_id]]
  
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_16s_ileum <- cv_test_predictions_multi_16s_ileum %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>%
  dplyr::ungroup()

roc_overall_summary_100cv_multi_16s_ileum <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_16s_ileum$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_16s_ileum$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_16s_ileum$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_16s_ileum <- roc_by_repeat_multi_16s_ileum %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x,
                                 x = 'all',
                                 input = 'threshold',
                                 ret = c('specificity', 'sensitivity'),
                                 transpose = FALSE) %>%
      as_tibble()
    
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    
    tpr_interp <- approx(x = fpr,
                         y = tpr,
                         xout = fpr_grid,
                         method = 'linear',
                         ties = 'ordered',
                         rule = 2)$y
    
    tibble::tibble(fpr = fpr_grid,
                   tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_16s_ileum <- ggplot(roc_curve_mean_100cv_multi_16s_ileum,
                                              aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max),
              fill = '#2C7FB8', alpha = 0.25) +
  geom_line(color = '#08519C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Ileum 16S, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_16s_ileum$auc, 3),
                         ' (95% CI: ',
                         round(roc_overall_summary_100cv_multi_16s_ileum$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_16s_ileum$auc_ci_max, 3), ')'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')

save(cv_test_predictions_multi_16s_ileum,
     roc_by_repeat_multi_16s_ileum,
     roc_overall_summary_100cv_multi_16s_ileum,
     roc_curve_mean_100cv_multi_16s_ileum,
     plot_roc_mean_100cv_multi_16s_ileum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_16s_ileum_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_16s_ileum_100cv_260428.pdf',
       plot = plot_roc_mean_100cv_multi_16s_ileum,
       width = 6,
       height = 6)

# 5-fold cross validation for ileum metabolite markers -------------------------------

sample_ids_met <- data_met_ileum$sample_id
covariate_names_met <- c('age', 'gender', 'race', 'use_antibiotics')
marker_met_names <- setdiff(colnames(data_met_ileum), c('sample_id', 'disease', covariate_names_met))

# generate seed list to reproduce the results
set.seed(20260427)
index_met <- sample(10000000, 100)

list_multi_metabolite_ileum_5_fold_cv <- lapply(seq_along(index_met), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index_met[z]
  
  set.seed(seed_id)
  
  # Stratified cross-validation by disease label
  id_case <- data_met_ileum$sample_id[data_met_ileum$disease == 1]
  id_control <- data_met_ileum$sample_id[data_met_ileum$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_met_ileum %>%
      dplyr::filter(sample_id %in% ids_training)
    data_test <- data_met_ileum %>%
      dplyr::filter(sample_id %in% ids_test)
    
    data_train <- data_train %>%
      dplyr::arrange(match(sample_id, ids_training))
    data_test <- data_test %>%
      dplyr::arrange(match(sample_id, ids_test))
    
    # Scale age and metabolite markers based on training set only
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_met_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_met_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_met_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_met_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_ileum_metabolite_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    result <- list(stat_result = stat_values,
                   data_train = temp_data_train,
                   data_test = temp_data_test,
                   logistic_model = logistic_model,
                   roc_train = roc_obj_train,
                   roc_test = roc_obj_test)
    
    return(result)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  return(result_list_5_fold_cv)
})

names(list_multi_metabolite_ileum_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_metabolite_ileum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_metabolite_ileum_5_fold_cv_260428.RData')

result_multi_metabolite_ileum_5_fold_cv <- lapply(list_multi_metabolite_ileum_5_fold_cv, function(z){
  z <- lapply(z, function(x){
    x$stat_result
  }) %>%
    bind_rows(.id = 'fold')
}) %>%
  bind_rows(.id = 'times')

result_multi_metabolite_ileum_5_fold_cv <- result_multi_metabolite_ileum_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_metabolite_ileum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_metabolite_ileum_5_fold_cv_260428.RData')

# ROC analysis for 100-times CV

cv_test_predictions_multi_metabolite_ileum <- lapply(names(list_multi_metabolite_ileum_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_metabolite_ileum_5_fold_cv[[rep_id]]
  
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_metabolite_ileum <- cv_test_predictions_multi_metabolite_ileum %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>%
  dplyr::ungroup()

roc_overall_100cv_multi_metabolite_ileum <- pROC::roc(response = cv_test_predictions_multi_metabolite_ileum$disease,
                                                      predictor = cv_test_predictions_multi_metabolite_ileum$prediction,
                                                      levels = c(0, 1),
                                                      direction = '<',
                                                      ci = TRUE,
                                                      quiet = TRUE)

roc_overall_summary_100cv_multi_metabolite_ileum <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_metabolite_ileum$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_metabolite_ileum$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_metabolite_ileum$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_metabolite_ileum <- roc_by_repeat_multi_metabolite_ileum %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x,
                                 x = 'all',
                                 input = 'threshold',
                                 ret = c('specificity', 'sensitivity'),
                                 transpose = FALSE) %>%
      as_tibble()
    
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    
    tpr_interp <- approx(x = fpr,
                         y = tpr,
                         xout = fpr_grid,
                         method = 'linear',
                         ties = 'ordered',
                         rule = 2)$y
    
    tibble::tibble(fpr = fpr_grid,
                   tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_metabolite_ileum <- ggplot(roc_curve_mean_100cv_multi_metabolite_ileum,
                                                     aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max),
              fill = '#31A354', alpha = 0.25) +
  geom_line(color = '#006D2C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Ileum Metabolite, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_metabolite_ileum$auc, 3),
                         ' (95% CI: ',
                         round(roc_overall_summary_100cv_multi_metabolite_ileum$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_metabolite_ileum$auc_ci_max, 3), ')'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')

save(cv_test_predictions_multi_metabolite_ileum,
     roc_by_repeat_multi_metabolite_ileum,
     roc_overall_100cv_multi_metabolite_ileum,
     roc_overall_summary_100cv_multi_metabolite_ileum,
     roc_curve_mean_100cv_multi_metabolite_ileum,
     plot_roc_mean_100cv_multi_metabolite_ileum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_metabolite_100cv_ileum_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_metabolite_100cv_ileum_260428.pdf',
       plot = plot_roc_mean_100cv_multi_metabolite_ileum,
       width = 6,
       height = 6)

# Merged ROC and key metrics for ileum -----------------------------------------------------

calc_metrics_at_best_threshold <- function(roc_by_repeat_df, prediction_df, model_id) {
  split_pred <- split(prediction_df, prediction_df[['times']])
  
  repeat_metrics <- lapply(names(split_pred), function(repeat_id) {
    x <- split_pred[[repeat_id]]
    roc_obj <- roc_by_repeat_df$roc_obj[[match(repeat_id, roc_by_repeat_df$times)]]
    
    if (is.null(roc_obj)) {
      roc_obj <- pROC::roc(response = x$disease,
                           predictor = x$prediction,
                           levels = c(0, 1),
                           direction = '<',
                           ci = TRUE,
                           quiet = TRUE)
    }
    
    best_cut <- pROC::coords(roc_obj,
                             x = 'best',
                             input = 'threshold',
                             best.method = 'youden',
                             ret = c('threshold', 'sensitivity', 'specificity'),
                             transpose = FALSE)
    
    threshold <- as.numeric(best_cut$threshold)
    sensitivity <- as.numeric(best_cut$sensitivity)
    specificity <- as.numeric(best_cut$specificity)
    
    pred_class <- ifelse(x$prediction >= threshold, 1, 0)
    true_class <- as.numeric(x$disease)
    
    tp <- sum(pred_class == 1 & true_class == 1, na.rm = TRUE)
    tn <- sum(pred_class == 0 & true_class == 0, na.rm = TRUE)
    fp <- sum(pred_class == 1 & true_class == 0, na.rm = TRUE)
    fn <- sum(pred_class == 0 & true_class == 1, na.rm = TRUE)
    
    n <- tp + tn + fp + fn
    accuracy <- ifelse(n > 0, (tp + tn) / n, NA_real_)
    ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
    npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
    precision <- ppv
    recall <- sensitivity
    f1_score <- ifelse((2 * tp + fp + fn) > 0, 2 * tp / (2 * tp + fp + fn), NA_real_)
    balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
    youden_index <- sensitivity + specificity - 1
    fpr <- ifelse((fp + tn) > 0, fp / (fp + tn), NA_real_)
    fnr <- ifelse((fn + tp) > 0, fn / (fn + tp), NA_real_)
    mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- ifelse(mcc_den > 0, ((tp * tn) - (fp * fn)) / mcc_den, NA_real_)
    
    data.frame(
      auc = as.numeric(roc_obj$auc),
      best_threshold = threshold,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      recall = recall,
      ppv = ppv,
      npv = npv,
      accuracy = accuracy,
      f1_score = f1_score,
      balanced_accuracy = balanced_accuracy,
      youden_index = youden_index,
      fpr = fpr,
      fnr = fnr,
      mcc = mcc,
      tp = tp,
      tn = tn,
      fp = fp,
      fn = fn,
      stringsAsFactors = FALSE
    )
  })
  
  repeat_metrics <- do.call(rbind, repeat_metrics)
  rownames(repeat_metrics) <- NULL
  
  metric_cols <- c('auc', 'best_threshold', 'sensitivity', 'specificity', 'precision',
                   'recall', 'ppv', 'npv', 'accuracy', 'f1_score', 'balanced_accuracy',
                   'youden_index', 'fpr', 'fnr', 'mcc', 'tp', 'tn', 'fp', 'fn')
  metric_means <- colMeans(repeat_metrics[, metric_cols, drop = FALSE], na.rm = TRUE)
  
  tibble::tibble(
    model = model_id,
    auc = metric_means[['auc']],
    auc_ci_min = as.numeric(quantile(repeat_metrics$auc, 0.025, na.rm = TRUE)),
    auc_ci_max = as.numeric(quantile(repeat_metrics$auc, 0.975, na.rm = TRUE)),
    best_threshold = metric_means[['best_threshold']],
    sensitivity = metric_means[['sensitivity']],
    specificity = metric_means[['specificity']],
    precision = metric_means[['precision']],
    recall = metric_means[['recall']],
    ppv = metric_means[['ppv']],
    npv = metric_means[['npv']],
    accuracy = metric_means[['accuracy']],
    f1_score = metric_means[['f1_score']],
    balanced_accuracy = metric_means[['balanced_accuracy']],
    youden_index = metric_means[['youden_index']],
    fpr = metric_means[['fpr']],
    fnr = metric_means[['fnr']],
    mcc = metric_means[['mcc']],
    tp = metric_means[['tp']],
    tn = metric_means[['tn']],
    fp = metric_means[['fp']],
    fn = metric_means[['fn']]
  )
}

best_threshold_metrics_16s_ileum <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_16s_ileum,
  prediction_df = cv_test_predictions_multi_16s_ileum,
  model_id = '16S'
)

best_threshold_metrics_metabolite_ileum <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_metabolite_ileum,
  prediction_df = cv_test_predictions_multi_metabolite_ileum,
  model_id = 'Metabolite'
)

best_threshold_metrics_compare_ileum <- dplyr::bind_rows(best_threshold_metrics_16s_ileum,
                                                         best_threshold_metrics_metabolite_ileum)

best_threshold_metrics_compare_ileum %>% 
  tibble::column_to_rownames('model') %>%
  sjmisc::rotate_df()


confusion_matrix_compare_ileum <- best_threshold_metrics_compare_ileum %>%
  dplyr::select(model, tp, fp, fn, tn) %>%
  tidyr::pivot_longer(cols = c(tp, fp, fn, tn), names_to = 'cell', values_to = 'count') %>%
  dplyr::mutate(actual = dplyr::case_when(
    cell %in% c('tp', 'fn') ~ 'Positive',
    TRUE ~ 'Negative'
  ),
  predicted = dplyr::case_when(
    cell %in% c('tp', 'fp') ~ 'Positive',
    TRUE ~ 'Negative'
  )) %>%
  dplyr::select(model, actual, predicted, count)

best_threshold_points_plot_ileum <- best_threshold_metrics_compare_ileum %>%
  dplyr::mutate(fpr_point = 1 - specificity,
                tpr_point = sensitivity)

roc_curve_compare_ileum <- dplyr::bind_rows(
  roc_curve_mean_100cv_multi_16s_ileum %>% dplyr::mutate(model = '16S'),
  roc_curve_mean_100cv_multi_metabolite_ileum %>% dplyr::mutate(model = 'Metabolite')
)

plot_roc_compare_100cv_ileum <- ggplot(roc_curve_compare_ileum,
                                       aes(x = fpr, y = mean_tpr, color = model, fill = model)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  scale_color_manual(values = c('16S' = '#08519C', 'Metabolite' = '#006D2C')) +
  scale_fill_manual(values = c('16S' = '#2C7FB8', 'Metabolite' = '#31A354')) +
  ZZWtool::ZZWTheme() +
  labs(title = 'Merged ROC Comparison: Ileum 16S vs Metabolite (100-times CV)',
       subtitle = paste0(
         '16S AUC = ', round(best_threshold_metrics_16s_ileum$auc, 3),
         ' [', round(best_threshold_metrics_16s_ileum$auc_ci_min, 3), '-', round(best_threshold_metrics_16s_ileum$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_16s_ileum$best_threshold, 3),
         '\nMetabolite AUC = ', round(best_threshold_metrics_metabolite_ileum$auc, 3),
         ' [', round(best_threshold_metrics_metabolite_ileum$auc_ci_min, 3), '-', round(best_threshold_metrics_metabolite_ileum$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_metabolite_ileum$best_threshold, 3)
       ),
       x = 'False Positive Rate',
       y = 'True Positive Rate',
       color = 'Model',
       fill = 'Model') +
  theme(legend.position = c(0.8, 0.2))

save(plot_roc_compare_100cv_ileum,
     best_threshold_metrics_16s_ileum,
     best_threshold_metrics_metabolite_ileum,
     best_threshold_metrics_compare_ileum,
     confusion_matrix_compare_ileum,
     best_threshold_points_plot_ileum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_ileum_100cv_260428.RData')

# readr::write_csv(best_threshold_metrics_compare_ileum,
#                  file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/best_threshold_metrics_compare_16s_vs_metabolite_ileum_100cv_260427.csv')

# readr::write_csv(confusion_matrix_compare_ileum,
#                  file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/confusion_matrix_compare_16s_vs_metabolite_ileum_100cv_260427.csv')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_ileum_100cv_260428.pdf',
       plot = plot_roc_compare_100cv_ileum,
       width = 7,
       height = 6)

rm(list = ls());gc()


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_data <- plot_roc_compare_100cv_ileum$data %>%
  dplyr::mutate(sample = 'ileum')

readr::write_csv(temp_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig7e_roc_compare_ileum.csv')


# Rectum ---------------------------------------------------------------------------------
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_rectum_260427.RData')
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_rectum_all_samples_260427.RData')

# 5-fold cross validation for rectum 16S markers -------------------------------

sample_ids <- data_16s_genus_rectum$data_16s
covariate_names <- c('age', 'gender', 'race', 'use_antibiotics')
marker_16s_names <- setdiff(colnames(data_16s_genus_rectum), c('data_16s', 'disease', covariate_names))

set.seed(20260427)
index <- sample(10000000, 100)

list_multi_16s_rectum_5_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index[z]
  set.seed(seed_id)
  
  id_case <- data_16s_genus_rectum$data_16s[data_16s_genus_rectum$disease == 1]
  id_control <- data_16s_genus_rectum$data_16s[data_16s_genus_rectum$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_16s_genus_rectum %>%
      dplyr::filter(data_16s %in% ids_training) %>%
      dplyr::arrange(match(data_16s, ids_training))
    data_test <- data_16s_genus_rectum %>%
      dplyr::filter(data_16s %in% ids_test) %>%
      dplyr::arrange(match(data_16s, ids_test))
    
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_16s_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_16s_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_16s_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_16s_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_rectum_16s_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    list(stat_result = stat_values,
         data_train = temp_data_train,
         data_test = temp_data_test,
         logistic_model = logistic_model,
         roc_train = roc_obj_train,
         roc_test = roc_obj_test)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  result_list_5_fold_cv
})

names(list_multi_16s_rectum_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_16s_rectum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_16s_rectum_5_fold_cv_260428.RData')

result_multi_16s_rectum_5_fold_cv <- lapply(list_multi_16s_rectum_5_fold_cv, function(z){
  lapply(z, function(x){x$stat_result}) %>% bind_rows(.id = 'fold')
}) %>% bind_rows(.id = 'times')

result_multi_16s_rectum_5_fold_cv <- result_multi_16s_rectum_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_16s_rectum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_16s_rectum_5_fold_cv_260428.RData')

cv_test_predictions_multi_16s_rectum <- lapply(names(list_multi_16s_rectum_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_16s_rectum_5_fold_cv[[rep_id]]
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% bind_rows()
}) %>% bind_rows() %>% dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_16s_rectum <- cv_test_predictions_multi_16s_rectum %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% dplyr::ungroup()

roc_overall_100cv_multi_16s_rectum <- pROC::roc(response = cv_test_predictions_multi_16s_rectum$disease,
                                                predictor = cv_test_predictions_multi_16s_rectum$prediction,
                                                levels = c(0, 1),
                                                direction = '<',
                                                ci = TRUE,
                                                quiet = TRUE)

roc_overall_summary_100cv_multi_16s_rectum <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_16s_rectum$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_16s_rectum$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_16s_rectum$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_16s_rectum <- roc_by_repeat_multi_16s_rectum %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x, x = 'all', input = 'threshold', ret = c('specificity', 'sensitivity'), transpose = FALSE) %>% as_tibble()
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    tpr_interp <- approx(x = fpr, y = tpr, xout = fpr_grid, method = 'linear', ties = 'ordered', rule = 2)$y
    tibble::tibble(fpr = fpr_grid, tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_16s_rectum <- ggplot(roc_curve_mean_100cv_multi_16s_rectum,
                                               aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), fill = '#2C7FB8', alpha = 0.25) +
  geom_line(color = '#08519C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Rectum 16S, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_16s_rectum$auc, 3),
                         ' (95% CI: ', round(roc_overall_summary_100cv_multi_16s_rectum$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_16s_rectum$auc_ci_max, 3), ')'),
       x = 'False Positive Rate', y = 'True Positive Rate')

save(cv_test_predictions_multi_16s_rectum,
     roc_by_repeat_multi_16s_rectum,
     roc_overall_100cv_multi_16s_rectum,
     roc_overall_summary_100cv_multi_16s_rectum,
     roc_curve_mean_100cv_multi_16s_rectum,
     plot_roc_mean_100cv_multi_16s_rectum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_16s_rectum_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_16s_rectum_100cv_260428.pdf',
       plot = plot_roc_mean_100cv_multi_16s_rectum,
       width = 6,
       height = 6)

# 5-fold cross validation for rectum metabolite markers -------------------------------

sample_ids_met <- data_met_rectum$sample_id
covariate_names_met <- c('age', 'gender', 'race', 'use_antibiotics')
marker_met_names <- setdiff(colnames(data_met_rectum), c('sample_id', 'disease', covariate_names_met))

set.seed(20260427)
index_met <- sample(10000000, 100)

list_multi_metabolite_rectum_5_fold_cv <- lapply(seq_along(index_met), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index_met[z]
  set.seed(seed_id)
  
  id_case <- data_met_rectum$sample_id[data_met_rectum$disease == 1]
  id_control <- data_met_rectum$sample_id[data_met_rectum$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_met_rectum %>%
      dplyr::filter(sample_id %in% ids_training) %>%
      dplyr::arrange(match(sample_id, ids_training))
    data_test <- data_met_rectum %>%
      dplyr::filter(sample_id %in% ids_test) %>%
      dplyr::arrange(match(sample_id, ids_test))
    
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_met_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_met_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_met_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_met_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_rectum_metabolite_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    list(stat_result = stat_values,
         data_train = temp_data_train,
         data_test = temp_data_test,
         logistic_model = logistic_model,
         roc_train = roc_obj_train,
         roc_test = roc_obj_test)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  result_list_5_fold_cv
})

names(list_multi_metabolite_rectum_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_metabolite_rectum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_metabolite_rectum_5_fold_cv_260428.RData')

result_multi_metabolite_rectum_5_fold_cv <- lapply(list_multi_metabolite_rectum_5_fold_cv, function(z){
  lapply(z, function(x){x$stat_result}) %>% bind_rows(.id = 'fold')
}) %>% bind_rows(.id = 'times')

result_multi_metabolite_rectum_5_fold_cv <- result_multi_metabolite_rectum_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_metabolite_rectum_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_metabolite_rectum_5_fold_cv_260428.RData')

cv_test_predictions_multi_metabolite_rectum <- lapply(names(list_multi_metabolite_rectum_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_metabolite_rectum_5_fold_cv[[rep_id]]
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% bind_rows()
}) %>% bind_rows() %>% dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_metabolite_rectum <- cv_test_predictions_multi_metabolite_rectum %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% dplyr::ungroup()

roc_overall_100cv_multi_metabolite_rectum <- pROC::roc(response = cv_test_predictions_multi_metabolite_rectum$disease,
                                                       predictor = cv_test_predictions_multi_metabolite_rectum$prediction,
                                                       levels = c(0, 1),
                                                       direction = '<',
                                                       ci = TRUE,
                                                       quiet = TRUE)

roc_overall_summary_100cv_multi_metabolite_rectum <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_metabolite_rectum$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_metabolite_rectum$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_metabolite_rectum$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_metabolite_rectum <- roc_by_repeat_multi_metabolite_rectum %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x, x = 'all', input = 'threshold', ret = c('specificity', 'sensitivity'), transpose = FALSE) %>% as_tibble()
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    tpr_interp <- approx(x = fpr, y = tpr, xout = fpr_grid, method = 'linear', ties = 'ordered', rule = 2)$y
    tibble::tibble(fpr = fpr_grid, tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_metabolite_rectum <- ggplot(roc_curve_mean_100cv_multi_metabolite_rectum,
                                                      aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), fill = '#31A354', alpha = 0.25) +
  geom_line(color = '#006D2C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Rectum Metabolite, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_metabolite_rectum$auc, 3),
                         ' (95% CI: ', round(roc_overall_summary_100cv_multi_metabolite_rectum$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_metabolite_rectum$auc_ci_max, 3), ')'),
       x = 'False Positive Rate', y = 'True Positive Rate')

save(cv_test_predictions_multi_metabolite_rectum,
     roc_by_repeat_multi_metabolite_rectum,
     roc_overall_100cv_multi_metabolite_rectum,
     roc_overall_summary_100cv_multi_metabolite_rectum,
     roc_curve_mean_100cv_multi_metabolite_rectum,
     plot_roc_mean_100cv_multi_metabolite_rectum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_metabolite_rectum_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_metabolite_rectum_100cv_260428.pdf',
       plot = plot_roc_mean_100cv_multi_metabolite_rectum,
       width = 6,
       height = 5)

# Merged ROC and key metrics for rectum -------------------------------------

best_threshold_metrics_16s_rectum <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_16s_rectum,
  prediction_df = cv_test_predictions_multi_16s_rectum,
  model_id = '16S'
)

best_threshold_metrics_metabolite_rectum <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_metabolite_rectum,
  prediction_df = cv_test_predictions_multi_metabolite_rectum,
  model_id = 'Metabolite'
)

best_threshold_metrics_compare_rectum <- dplyr::bind_rows(best_threshold_metrics_16s_rectum,
                                                          best_threshold_metrics_metabolite_rectum)

best_threshold_metrics_compare_rectum %>% 
  tibble::column_to_rownames('model') %>% 
  sjmisc::rotate_df()

confusion_matrix_compare_rectum <- best_threshold_metrics_compare_rectum %>%
  dplyr::select(model, tp, fp, fn, tn) %>%
  tidyr::pivot_longer(cols = c(tp, fp, fn, tn), names_to = 'cell', values_to = 'count') %>%
  dplyr::mutate(actual = dplyr::case_when(cell %in% c('tp', 'fn') ~ 'Positive', TRUE ~ 'Negative'),
                predicted = dplyr::case_when(cell %in% c('tp', 'fp') ~ 'Positive', TRUE ~ 'Negative')) %>%
  dplyr::select(model, actual, predicted, count)

best_threshold_points_plot_rectum <- best_threshold_metrics_compare_rectum %>%
  dplyr::mutate(fpr_point = 1 - specificity, tpr_point = sensitivity)

roc_curve_compare_rectum <- dplyr::bind_rows(
  roc_curve_mean_100cv_multi_16s_rectum %>% dplyr::mutate(model = '16S'),
  roc_curve_mean_100cv_multi_metabolite_rectum %>% dplyr::mutate(model = 'Metabolite')
)

plot_roc_compare_100cv_rectum <- ggplot(roc_curve_compare_rectum,
                                        aes(x = fpr, y = mean_tpr, color = model, fill = model)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  #   geom_point(data = best_threshold_points_plot_rectum,
  #              aes(x = fpr_point, y = tpr_point, color = model),
  #              size = 3,
  #              inherit.aes = FALSE) +
  coord_equal() +
  scale_color_manual(values = c('16S' = '#08519C', 'Metabolite' = '#006D2C')) +
  scale_fill_manual(values = c('16S' = '#2C7FB8', 'Metabolite' = '#31A354')) +
  ZZWtool::ZZWTheme() +
  labs(title = 'Merged ROC Comparison: Rectum 16S vs Metabolite (100-times CV)',
       subtitle = paste0(
         '16S AUC = ', round(best_threshold_metrics_16s_rectum$auc, 3),
         ' [', round(best_threshold_metrics_16s_rectum$auc_ci_min, 3), '-', round(best_threshold_metrics_16s_rectum$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_16s_rectum$best_threshold, 3),
         '\nMetabolite AUC = ', round(best_threshold_metrics_metabolite_rectum$auc, 3),
         ' [', round(best_threshold_metrics_metabolite_rectum$auc_ci_min, 3), '-', round(best_threshold_metrics_metabolite_rectum$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_metabolite_rectum$best_threshold, 3)
       ),
       x = 'False Positive Rate',
       y = 'True Positive Rate',
       color = 'Model',
       fill = 'Model') +
  theme(legend.position = c(0.8, 0.2))

save(plot_roc_compare_100cv_rectum,
     best_threshold_metrics_16s_rectum,
     best_threshold_metrics_metabolite_rectum,
     best_threshold_metrics_compare_rectum,
     confusion_matrix_compare_rectum,
     best_threshold_points_plot_rectum,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_rectum_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_rectum_100cv_260428.pdf',
       plot = plot_roc_compare_100cv_rectum,
       width = 7,
       height = 6)

rm(list = ls());gc()

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_data <- plot_roc_compare_100cv_rectum$data %>%
  dplyr::mutate(sample = 'rectum')

readr::write_csv(temp_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig7f_roc_compare_rectum.csv')


# Stool ---------------------------------------------------------------------------------
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_met_stool_260427.RData')
load('~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/data_16s_genus_stool_all_samples_260427.RData')

# 5-fold cross validation for stool 16S markers -------------------------------

sample_ids <- data_16s_genus_stool$data_16s
covariate_names <- c('age', 'gender', 'race', 'use_antibiotics')
marker_16s_names <- setdiff(colnames(data_16s_genus_stool), c('data_16s', 'disease', covariate_names))

set.seed(20260427)
index <- sample(10000000, 100)

list_multi_16s_stool_5_fold_cv <- lapply(seq_along(index), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index[z]
  set.seed(seed_id)
  
  id_case <- data_16s_genus_stool$data_16s[data_16s_genus_stool$disease == 1]
  id_control <- data_16s_genus_stool$data_16s[data_16s_genus_stool$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_16s_genus_stool %>%
      dplyr::filter(data_16s %in% ids_training) %>%
      dplyr::arrange(match(data_16s, ids_training))
    data_test <- data_16s_genus_stool %>%
      dplyr::filter(data_16s %in% ids_test) %>%
      dplyr::arrange(match(data_16s, ids_test))
    
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_16s_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_16s_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_16s_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_16s_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_16s_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_stool_16s_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    list(stat_result = stat_values,
         data_train = temp_data_train,
         data_test = temp_data_test,
         logistic_model = logistic_model,
         roc_train = roc_obj_train,
         roc_test = roc_obj_test)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  result_list_5_fold_cv
})

names(list_multi_16s_stool_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_16s_stool_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_16s_stool_5_fold_cv_260428.RData')

result_multi_16s_stool_5_fold_cv <- lapply(list_multi_16s_stool_5_fold_cv, function(z){
  lapply(z, function(x){x$stat_result}) %>% bind_rows(.id = 'fold')
}) %>% bind_rows(.id = 'times')

result_multi_16s_stool_5_fold_cv <- result_multi_16s_stool_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_16s_stool_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_16s_stool_5_fold_cv_260428.RData')

cv_test_predictions_multi_16s_stool <- lapply(names(list_multi_16s_stool_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_16s_stool_5_fold_cv[[rep_id]]
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% bind_rows()
}) %>% bind_rows() %>% dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_16s_stool <- cv_test_predictions_multi_16s_stool %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% dplyr::ungroup()

roc_overall_100cv_multi_16s_stool <- pROC::roc(response = cv_test_predictions_multi_16s_stool$disease,
                                               predictor = cv_test_predictions_multi_16s_stool$prediction,
                                               levels = c(0, 1),
                                               direction = '<',
                                               ci = TRUE,
                                               quiet = TRUE)

roc_overall_summary_100cv_multi_16s_stool <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_16s_stool$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_16s_stool$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_16s_stool$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_16s_stool <- roc_by_repeat_multi_16s_stool %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x, x = 'all', input = 'threshold', ret = c('specificity', 'sensitivity'), transpose = FALSE) %>% as_tibble()
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    tpr_interp <- approx(x = fpr, y = tpr, xout = fpr_grid, method = 'linear', ties = 'ordered', rule = 2)$y
    tibble::tibble(fpr = fpr_grid, tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_16s_stool <- ggplot(roc_curve_mean_100cv_multi_16s_stool,
                                              aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), fill = '#2C7FB8', alpha = 0.25) +
  geom_line(color = '#08519C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Stool 16S, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_16s_stool$auc, 3),
                         ' (95% CI: ', round(roc_overall_summary_100cv_multi_16s_stool$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_16s_stool$auc_ci_max, 3), ')'),
       x = 'False Positive Rate', y = 'True Positive Rate')

save(cv_test_predictions_multi_16s_stool,
     roc_by_repeat_multi_16s_stool,
     roc_overall_100cv_multi_16s_stool,
     roc_overall_summary_100cv_multi_16s_stool,
     roc_curve_mean_100cv_multi_16s_stool,
     plot_roc_mean_100cv_multi_16s_stool,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_16s_stool_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_16s_stool_100cv_260428.pdf',
       plot = plot_roc_mean_100cv_multi_16s_stool,
       width = 6,
       height = 6)

# 5-fold cross validation for stool metabolite markers -------------------------------

sample_ids_met <- data_met_stool$sample_id
covariate_names_met <- c('age', 'gender', 'race', 'use_antibiotics')
marker_met_names <- setdiff(colnames(data_met_stool), c('sample_id', 'disease', covariate_names_met))

set.seed(20260427)
index_met <- sample(10000000, 100)

list_multi_metabolite_stool_5_fold_cv <- lapply(seq_along(index_met), function(z){
  cat('Perform the ', z, 'th randomization\n')
  seed_id <- index_met[z]
  set.seed(seed_id)
  
  id_case <- data_met_stool$sample_id[data_met_stool$disease == 1]
  id_control <- data_met_stool$sample_id[data_met_stool$disease == 0]
  
  id_case <- sample(id_case)
  id_control <- sample(id_control)
  
  groups_case <- rep(1:5, length.out = length(id_case))
  groups_control <- rep(1:5, length.out = length(id_control))
  
  split_case <- split(id_case, groups_case)
  split_control <- split(id_control, groups_control)
  
  cv_splits <- lapply(1:5, function(f) {
    c(split_case[[f]], split_control[[f]])
  })
  
  result_list_5_fold_cv <- lapply(1:5, function(j){
    cat('Perform the ', j, 'th cross validation\n')
    
    ids_test <- unlist(cv_splits[[j]])
    ids_training <- unlist(cv_splits[-j])
    
    data_train <- data_met_stool %>%
      dplyr::filter(sample_id %in% ids_training) %>%
      dplyr::arrange(match(sample_id, ids_training))
    data_test <- data_met_stool %>%
      dplyr::filter(sample_id %in% ids_test) %>%
      dplyr::arrange(match(sample_id, ids_test))
    
    age_mean <- mean(data_train$age, na.rm = TRUE)
    age_sd <- sd(data_train$age, na.rm = TRUE)
    if (is.na(age_sd) || age_sd == 0) {
      age_sd <- 1
    }
    data_train$age <- (data_train$age - age_mean) / age_sd
    data_test$age <- (data_test$age - age_mean) / age_sd
    
    for (marker in marker_met_names) {
      m_mean <- mean(data_train[[marker]], na.rm = TRUE)
      m_sd <- sd(data_train[[marker]], na.rm = TRUE)
      if (is.na(m_sd) || m_sd == 0) {
        m_sd <- 1
      }
      data_train[[marker]] <- (data_train[[marker]] - m_mean) / m_sd
      data_test[[marker]] <- (data_test[[marker]] - m_mean) / m_sd
    }
    
    temp_data_train <- data_train %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    temp_data_test <- data_test %>%
      as_tibble() %>%
      dplyr::select(disease, age, gender, race, use_antibiotics, dplyr::all_of(marker_met_names))
    
    logistic_formula <- reformulate(termlabels = c('gender', 'age', 'race', 'use_antibiotics', marker_met_names),
                                    response = 'disease')
    logistic_model <- glm(logistic_formula,
                          data = temp_data_train,
                          family = binomial(link = 'logit'))
    
    coef_table <- summary(logistic_model)$coefficients
    coef_terms_clean <- gsub('`', '', rownames(coef_table), fixed = TRUE)
    marker_idx <- match(marker_met_names, coef_terms_clean)
    valid_marker <- !is.na(marker_idx)
    
    possibility_train <- predict(logistic_model, newdata = temp_data_train, type = 'response')
    roc_obj_train <- roc(temp_data_train$disease ~ possibility_train, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_train <- roc_obj_train$auc[[1]]
    auc_ci_train <- paste(c(roc_obj_train$ci[[1]], roc_obj_train$ci[[3]]), collapse = ';')
    
    possibility_test <- predict(logistic_model, newdata = temp_data_test, type = 'response')
    roc_obj_test <- roc(temp_data_test$disease ~ possibility_test, plot = FALSE, print.auc = TRUE, ci = TRUE)
    auc_test <- roc_obj_test$auc[[1]]
    auc_ci_test <- paste(c(roc_obj_test$ci[[1]], roc_obj_test$ci[[3]]), collapse = ';')
    
    if (any(valid_marker)) {
      est <- coef_table[marker_idx[valid_marker], 'Estimate']
      se <- coef_table[marker_idx[valid_marker], 'Std. Error']
      
      stat_values <- tibble::tibble(variable_id = marker_met_names[valid_marker],
                                    odds = exp(est),
                                    odds_ci = paste(exp(est - 1.96 * se), exp(est + 1.96 * se), sep = ';'),
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    } else {
      stat_values <- tibble::tibble(variable_id = 'all_stool_metabolite_markers',
                                    odds = NA_real_,
                                    odds_ci = 'NA;NA',
                                    auc_train = auc_train,
                                    auc_ci_train = auc_ci_train,
                                    auc_test = auc_test,
                                    auc_ci_test = auc_ci_test)
    }
    
    list(stat_result = stat_values,
         data_train = temp_data_train,
         data_test = temp_data_test,
         logistic_model = logistic_model,
         roc_train = roc_obj_train,
         roc_test = roc_obj_test)
  })
  
  names(result_list_5_fold_cv) <- paste0('fold_', 1:5)
  result_list_5_fold_cv
})

names(list_multi_metabolite_stool_5_fold_cv) <- paste0('repeat_', 1:100)

save(list_multi_metabolite_stool_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/list_multi_metabolite_stool_5_fold_cv_260428.RData')

result_multi_metabolite_stool_5_fold_cv <- lapply(list_multi_metabolite_stool_5_fold_cv, function(z){
  lapply(z, function(x){x$stat_result}) %>% bind_rows(.id = 'fold')
}) %>% bind_rows(.id = 'times')

result_multi_metabolite_stool_5_fold_cv <- result_multi_metabolite_stool_5_fold_cv %>%
  tidyr::separate(col = 'odds_ci', sep = ';', into = c('odds_ci_min', 'odds_ci_max')) %>%
  tidyr::separate(col = 'auc_ci_train', sep = ';', into = c('auc_ci_train_min', 'auc_ci_train_max')) %>%
  tidyr::separate(col = 'auc_ci_test', sep = ';', into = c('auc_ci_test_min', 'auc_ci_test_max')) %>%
  dplyr::mutate(odds_ci_min = as.numeric(odds_ci_min),
                odds_ci_max = as.numeric(odds_ci_max),
                auc_ci_train_min = as.numeric(auc_ci_train_min),
                auc_ci_train_max = as.numeric(auc_ci_train_max),
                auc_ci_test_min = as.numeric(auc_ci_test_min),
                auc_ci_test_max = as.numeric(auc_ci_test_max)) %>%
  dplyr::mutate(idx = seq_len(n())) %>%
  dplyr::rename('odds_ratio' = 'odds') %>%
  dplyr::mutate(color = case_when(odds_ratio > 1 ~ 'risk factor',
                                  odds_ratio <= 1 ~ 'protective factor'))

save(result_multi_metabolite_stool_5_fold_cv,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/result_multi_metabolite_stool_5_fold_cv_260428.RData')

cv_test_predictions_multi_metabolite_stool <- lapply(names(list_multi_metabolite_stool_5_fold_cv), function(rep_id) {
  fold_list <- list_multi_metabolite_stool_5_fold_cv[[rep_id]]
  lapply(names(fold_list), function(fold_id) {
    roc_obj <- fold_list[[fold_id]]$roc_test
    response_vec <- roc_obj$original.response
    predictor_vec <- roc_obj$original.predictor
    if (is.factor(response_vec)) {
      response_vec <- as.numeric(as.character(response_vec))
    }
    tibble::tibble(times = rep_id,
                   fold = fold_id,
                   disease = as.numeric(response_vec),
                   prediction = as.numeric(predictor_vec))
  }) %>% bind_rows()
}) %>% bind_rows() %>% dplyr::filter(!is.na(disease), !is.na(prediction))

roc_by_repeat_multi_metabolite_stool <- cv_test_predictions_multi_metabolite_stool %>%
  dplyr::group_by(times) %>%
  dplyr::group_modify(~{
    roc_obj <- pROC::roc(response = .x$disease,
                         predictor = .x$prediction,
                         levels = c(0, 1),
                         direction = '<',
                         ci = TRUE,
                         quiet = TRUE)
    tibble::tibble(auc = as.numeric(roc_obj$auc),
                   auc_ci_min = as.numeric(roc_obj$ci[[1]]),
                   auc_ci_max = as.numeric(roc_obj$ci[[3]]),
                   roc_obj = list(roc_obj))
  }) %>% dplyr::ungroup()

roc_overall_100cv_multi_metabolite_stool <- pROC::roc(response = cv_test_predictions_multi_metabolite_stool$disease,
                                                      predictor = cv_test_predictions_multi_metabolite_stool$prediction,
                                                      levels = c(0, 1),
                                                      direction = '<',
                                                      ci = TRUE,
                                                      quiet = TRUE)

roc_overall_summary_100cv_multi_metabolite_stool <- tibble::tibble(
  auc = mean(roc_by_repeat_multi_metabolite_stool$auc, na.rm = TRUE),
  auc_ci_min = as.numeric(quantile(roc_by_repeat_multi_metabolite_stool$auc, 0.025, na.rm = TRUE)),
  auc_ci_max = as.numeric(quantile(roc_by_repeat_multi_metabolite_stool$auc, 0.975, na.rm = TRUE))
)

fpr_grid <- seq(0, 1, by = 0.01)

roc_curve_mean_100cv_multi_metabolite_stool <- roc_by_repeat_multi_metabolite_stool %>%
  dplyr::mutate(curve = purrr::map(roc_obj, function(x) {
    curve_points <- pROC::coords(x, x = 'all', input = 'threshold', ret = c('specificity', 'sensitivity'), transpose = FALSE) %>% as_tibble()
    fpr <- 1 - curve_points$specificity
    tpr <- curve_points$sensitivity
    ord <- order(fpr, tpr)
    fpr <- fpr[ord]
    tpr <- tpr[ord]
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    tpr_interp <- approx(x = fpr, y = tpr, xout = fpr_grid, method = 'linear', ties = 'ordered', rule = 2)$y
    tibble::tibble(fpr = fpr_grid, tpr = tpr_interp)
  })) %>%
  dplyr::select(times, curve) %>%
  tidyr::unnest(curve) %>%
  dplyr::group_by(fpr) %>%
  dplyr::summarise(mean_tpr = mean(tpr),
                   tpr_ci_min = quantile(tpr, 0.025),
                   tpr_ci_max = quantile(tpr, 0.975),
                   .groups = 'drop')

plot_roc_mean_100cv_multi_metabolite_stool <- ggplot(roc_curve_mean_100cv_multi_metabolite_stool,
                                                     aes(x = fpr, y = mean_tpr)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), fill = '#31A354', alpha = 0.25) +
  geom_line(color = '#006D2C', size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  coord_equal() +
  ZZWtool::ZZWTheme() +
  labs(title = 'Mean ROC Curve of 100-times CV (Stool Metabolite, Test Sets)',
       subtitle = paste0('Overall AUC = ', round(roc_overall_summary_100cv_multi_metabolite_stool$auc, 3),
                         ' (95% CI: ', round(roc_overall_summary_100cv_multi_metabolite_stool$auc_ci_min, 3), '-',
                         round(roc_overall_summary_100cv_multi_metabolite_stool$auc_ci_max, 3), ')'),
       x = 'False Positive Rate', y = 'True Positive Rate')

save(cv_test_predictions_multi_metabolite_stool,
     roc_by_repeat_multi_metabolite_stool,
     roc_overall_100cv_multi_metabolite_stool,
     roc_overall_summary_100cv_multi_metabolite_stool,
     roc_curve_mean_100cv_multi_metabolite_stool,
     plot_roc_mean_100cv_multi_metabolite_stool,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_analysis_multi_metabolite_stool_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_curve_mean_multi_metabolite_stool_100cv_260428.pdf',
       plot = plot_roc_mean_100cv_multi_metabolite_stool,
       width = 6,
       height = 6)

# Merged ROC and key metrics for stool -------------------------------------

best_threshold_metrics_16s_stool <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_16s_stool,
  prediction_df = cv_test_predictions_multi_16s_stool,
  model_id = '16S'
)

best_threshold_metrics_metabolite_stool <- calc_metrics_at_best_threshold(
  roc_by_repeat_df = roc_by_repeat_multi_metabolite_stool,
  prediction_df = cv_test_predictions_multi_metabolite_stool,
  model_id = 'Metabolite'
)

best_threshold_metrics_compare_stool <- dplyr::bind_rows(best_threshold_metrics_16s_stool,
                                                         best_threshold_metrics_metabolite_stool)

best_threshold_metrics_compare_stool %>% 
  tibble::column_to_rownames(var = 'model') %>%
  sjmisc::rotate_df()

confusion_matrix_compare_stool <- best_threshold_metrics_compare_stool %>%
  dplyr::select(model, tp, fp, fn, tn) %>%
  tidyr::pivot_longer(cols = c(tp, fp, fn, tn), names_to = 'cell', values_to = 'count') %>%
  dplyr::mutate(actual = dplyr::case_when(cell %in% c('tp', 'fn') ~ 'Positive', TRUE ~ 'Negative'),
                predicted = dplyr::case_when(cell %in% c('tp', 'fp') ~ 'Positive', TRUE ~ 'Negative')) %>%
  dplyr::select(model, actual, predicted, count)

best_threshold_points_plot_stool <- best_threshold_metrics_compare_stool %>%
  dplyr::mutate(fpr_point = 1 - specificity, tpr_point = sensitivity)

roc_curve_compare_stool <- dplyr::bind_rows(
  roc_curve_mean_100cv_multi_16s_stool %>% dplyr::mutate(model = '16S'),
  roc_curve_mean_100cv_multi_metabolite_stool %>% dplyr::mutate(model = 'Metabolite')
)

plot_roc_compare_100cv_stool <- ggplot(roc_curve_compare_stool,
                                       aes(x = fpr, y = mean_tpr, color = model, fill = model)) +
  geom_ribbon(aes(ymin = tpr_ci_min, ymax = tpr_ci_max), alpha = 0.20, color = NA) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40') +
  #   geom_point(data = best_threshold_points_plot_stool,
  #              aes(x = fpr_point, y = tpr_point, color = model),
  #              size = 3,
  #              inherit.aes = FALSE) +
  coord_equal() +
  scale_color_manual(values = c('16S' = '#08519C', 'Metabolite' = '#006D2C')) +
  scale_fill_manual(values = c('16S' = '#2C7FB8', 'Metabolite' = '#31A354')) +
  ZZWtool::ZZWTheme() +
  labs(title = 'Merged ROC Comparison: Stool 16S vs Metabolite (100-times CV)',
       subtitle = paste0(
         '16S AUC = ', round(best_threshold_metrics_16s_stool$auc, 3),
         ' [', round(best_threshold_metrics_16s_stool$auc_ci_min, 3), '-', round(best_threshold_metrics_16s_stool$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_16s_stool$best_threshold, 3),
         '\nMetabolite AUC = ', round(best_threshold_metrics_metabolite_stool$auc, 3),
         ' [', round(best_threshold_metrics_metabolite_stool$auc_ci_min, 3), '-', round(best_threshold_metrics_metabolite_stool$auc_ci_max, 3),
         '], Best threshold = ', round(best_threshold_metrics_metabolite_stool$best_threshold, 3)
       ),
       x = 'False Positive Rate',
       y = 'True Positive Rate',
       color = 'Model',
       fill = 'Model') +
  theme(legend.position = c(0.8, 0.2))

save(plot_roc_compare_100cv_stool,
     best_threshold_metrics_16s_stool,
     best_threshold_metrics_metabolite_stool,
     best_threshold_metrics_compare_stool,
     confusion_matrix_compare_stool,
     best_threshold_points_plot_stool,
     file = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_stool_100cv_260428.RData')

ggsave(filename = '~/Project/00_IBD_project/Data/20260423_compare_prediction_model_serology_metabolite_16s/roc_merged_compare_16s_vs_metabolite_stool_100cv_260428.pdf',
       plot = plot_roc_compare_100cv_stool,
       width = 7,
       height = 6)

rm(list = ls());gc()

dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_data <- plot_roc_compare_100cv_stool$data %>%
  dplyr::mutate(sample = 'stool')

readr::write_csv(temp_data,
                 file = '~/Project/00_IBD_project/Data/20260718_source_data/ext_fig7g_roc_compare_stool.csv')






