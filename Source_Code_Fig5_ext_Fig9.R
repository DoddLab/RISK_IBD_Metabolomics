################################################################################
# Differential analysis of longitudinal severity in CD patients ------------------------------------
dir.create('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update', showWarnings = FALSE, recursive = TRUE)
setwd('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update')

# Extract CD follow-up samples -------------------------------------------------

library(tidyverse)
rm(list = ls())

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')

object_CD_followup <- object_final %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(sample_id != 'B032_S34') %>% 
  filter(phenotype_group1 == 'CD') %>% 
  filter(severity_class != 'unavailable')

object_CD_followup@sample_info <- object_CD_followup@sample_info %>% 
  mutate(severity_class2 = case_when(
    severity_class == 'remission' ~ 'remission',
    severity_class == 'mild' ~ 'mild',
    severity_class == 'moderate' ~ 'moderate_severe',
    severity_class == 'severe' ~ 'moderate_severe'
  )) %>% 
  select(sample_id:severity_class, severity_class2, everything())

object_CD_followup@sample_info_note <- data.frame(sample_id = colnames(object_CD_followup@sample_info), 
                                                  note = colnames(object_CD_followup@sample_info), 
                                                  stringsAsFactors = FALSE)

save(object_CD_followup, 
     file = '~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/object_CD_followup_250401.RData')

# Linear mixed effect model ----------------------------------------------------
library(Maaslin2)
library(sjmisc)

# scaling z-score for metabolite data
input_met_table <- object_CD_followup %>%
  extract_expression_data() %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'ID')

# modify meta info
input_meta_table <- object_CD_followup %>%
  extract_sample_info() %>% 
  dplyr::select(sample_id, patient_id, severity_class2, gender, age, race, use_antibiotics)



# write these files
temp_path <- '~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/'
dir.create('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/', showWarnings = FALSE, recursive = TRUE)
write_tsv(input_met_table, file = file.path(temp_path, 'input_data_metabolite.tsv'))
write_tsv(input_meta_table, file = file.path(temp_path, 'input_data_meta.tsv'))

fit_data <- Maaslin2(
  input_data = file.path(temp_path, 'input_data_metabolite.tsv'),
  input_metadata = file.path(temp_path, 'input_data_meta.tsv'),
  output = file.path(temp_path, 'Maaslin2_output'),
  fixed_effects = c('severity_class2', 'gender', 'age', 'race', 'use_antibiotics'),
  random_effects = c('patient_id'),
  normalization = 'NONE',
  transform = 'LOG',
  reference = "severity_class2,remission;race,Caucasian",
  standardize = TRUE,
  save_models = TRUE)

# extract the Maaslin2 output --------------------------------------------------

annot_table <- object_CD_followup %>% 
  extract_annotation_table()

result_maaslin2 <- readr::read_tsv('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/Maaslin2_output/all_results.tsv')

result_maaslin2_severity <- result_maaslin2 %>%
  dplyr::filter(metadata == 'severity_class2') %>% 
  left_join(annot_table, by = c('feature' = 'variable_id')) %>% 
  dplyr::select(feature:qval, Compound.name:rt, id:adduct, confidence_level, metabolon_class:metabolon_subclass) %>% 
  rename('variable_id' = feature,
         'compound_name' = 'Compound.name',
         'severity_class' = 'value')

save(result_maaslin2_severity, 
     file = '~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/result_maaslin2_severity_250401.RData')



# result_maaslin2_severity %>% 
#   group_by(severity_class) %>% 
#   count(qval <= 0.05)

# visualization - barplot comparsion -------------------------------------------
# count significant metabolites per severity class, split by direction
plot_data <- result_maaslin2_severity %>%
  group_by(severity_class) %>%
  count(qval <= 0.05, coef > 0) %>%
  rename(is_sig = `qval <= 0.05`,
         is_positive = `coef > 0`) %>%
  ungroup() %>%
  dplyr::filter(is_sig) %>%
  mutate(severity_class = factor(severity_class,
                                 levels = c('moderate_severe', 'mild'),
                                 labels = c('Moderate/severe\nvs. remission', 'Mild vs. remission')),
         # stack so 'Higher in remission' (green) is at the left, orange to its right
         is_positive = factor(is_positive, levels = c('TRUE', 'FALSE')))

# per-class totals for the 'n = ...' labels and the fold-change annotation
plot_totals <- plot_data %>%
  group_by(severity_class) %>%
  summarise(total = sum(n), .groups = 'drop')

fold_change <- max(plot_totals$total) / min(plot_totals$total)

temp_plot <- plot_data %>%
  ggplot(aes(x = severity_class, y = n, fill = is_positive, label = n)) +
  geom_bar(stat = 'identity', width = 0.6) +
  geom_text(position = position_stack(vjust = 0.5), size = 3.2, colour = 'white') +
  # 'n = ...' total to the right of each bar
  geom_text(data = plot_totals, inherit.aes = FALSE,
            aes(x = severity_class, y = total, label = paste('n =', total)),
            hjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c('TRUE' = '#d95f02', 'FALSE' = '#1b9e77'),
                    labels = c('TRUE' = 'Higher in active disease', 'FALSE' = 'Higher in remission'),
                    breaks = c('FALSE', 'TRUE'),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 200), expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = 'off') +
  # fold-change bracket + label spanning the two bars
  annotate('segment', x = 1, xend = 2, y = 175, yend = 175) +
  annotate('segment', x = 1, xend = 1, y = 170, yend = 175) +
  annotate('segment', x = 2, xend = 2, y = 170, yend = 175) +
  annotate('text', x = 1.5, y = 175, vjust = -0.6, size = 3.4,
           label = paste0(round(fold_change, 1), 'x increase')) +
  labs(title = 'Differential metabolites vs remission',
       y = 'Differential metabolites (n)', x = NULL) +
  ZZWtool::ZZWTheme(type = 'classic') +
  theme(legend.position = 'top',
        legend.justification = 'left',
        plot.title = element_text(face = 'bold'),
        plot.margin = margin(10, 40, 10, 10))

ggsave(temp_plot,
       filename = '~/Project/00_IBD_project/Figure/250326/Supplementary_Figure6/plot_severity_metabolite_increase_decrease_250506.pdf',
       width = 6, height = 3)


# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- temp_plot$data
# readr::write_csv(temp_source_data, 
#                    file = '~/Project/00_IBD_project/Data/20260718_source_data/Fig5b_differential_met_num_260721.csv')





################################################################################
# Metabolite subclass enrichment analysis --------------------------------------

# enrichClass2 ------------------------------------------------------------------
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

# Class enrichment according to increase and decrease --------------------------
load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/result_maaslin2_severity_250401.RData')

# moderate severe increase class enrichment
table_severity_increase <- result_maaslin2_severity %>%
  filter(severity_class == 'moderate_severe') %>%
  select(compound_name, qval, coef, metabolon_subclass) %>% 
  filter(coef > 0) %>% 
  rename(pvalue = qval, 
         fc = coef,
         class = metabolon_subclass) %>% 
  arrange(pvalue) %>% 
  mutate(order = seq(n()))


table_severity_decrease <- result_maaslin2_severity %>%
  filter(severity_class == 'moderate_severe') %>%
  select(compound_name, qval, coef, metabolon_subclass) %>% 
  filter(coef < 0) %>% 
  rename(pvalue = qval, 
         fc = coef,
         class = metabolon_subclass) %>% 
  arrange(pvalue) %>% 
  mutate(order = seq(n()))

class_enrich_severe_increase <- enrichClass2(table_severity_increase,
                                             class_cutoff = 1,
                                             p_cutoff = 1) %>% 
  arrange(adjustedpvalue) %>% 
  mutate(logP = -log(adjustedpvalue)) %>%
  slice(1:10) %>% 
  mutate(idx = seq(n(), 1, -1)) %>% 
  mutate(is_sig = case_when(adjustedpvalue < 0.05 ~ 'sig',
                            TRUE ~ 'ns'))

# temp_plot <- ggplot(class_enrich_severe_increase, aes(x = reorder(name, logP), y = logP, fill = is_sig)) +
#   geom_bar(stat = 'identity') +
#   coord_flip() +
#   geom_hline(yintercept = -log(0.05), linetype = 'dashed', color = 'black') +
#   scale_fill_manual(values = c('sig' = '#fc8070', 'ns' = '#c5c7c9'), 
#                     labels = c('sig' = 'Significant', 'ns' = 'Not Significant'),
#                     name = 'vs. Remission') +
#   # geom_text(aes(label = name), hjust = 1.5, size = 3) +
#   scale_y_continuous(expand = c(0, 0)) +
#   ZZWtool::ZZWTheme(type = 'classic') +
#   theme(legend.position = c(0.8, 0.2),,
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
# 
# ggsave(temp_plot, 
#        filename = '~/Project/00_IBD_project/Figure/250326/Figure5/plot_class_enrichment_severe_increase_250404.pdf', 
#        width = 8, height = 6)


class_enrich_severe_decrease <- enrichClass2(table_severity_decrease,
                                             class_cutoff = 1,
                                             p_cutoff = 1) %>% 
  arrange(adjustedpvalue) %>% 
  mutate(logP = -log(adjustedpvalue)) %>%
  slice(1:10) %>% 
  mutate(idx = seq(n(), 1, -1)) %>% 
  mutate(is_sig = case_when(adjustedpvalue < 0.05 ~ 'sig',
                            TRUE ~ 'ns'))


# temp_plot <- ggplot(class_enrich_severe_decrease, aes(x = reorder(name, logP), y = logP, fill = is_sig)) +
#   geom_bar(stat = 'identity') +
#   coord_flip() +
#   geom_hline(yintercept = -log(0.05), linetype = 'dashed', color = 'black') +
#   scale_fill_manual(values = c('sig' = '#7fb1d3', 'ns' = '#c5c7c9'), 
#                     labels = c('sig' = 'Significant', 'ns' = 'Not Significant'),
#                     name = 'vs. Remission') +
#   # geom_text(aes(label = name), hjust = 1.5, size = 3) +
#   scale_y_continuous(expand = c(0, 0)) +
#   ZZWtool::ZZWTheme(type = 'classic') +
#   theme(legend.position = c(0.8, 0.2),,
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5))
# 
# ggsave(temp_plot, 
#        filename = '~/Project/00_IBD_project/Figure/250326/Figure5/plot_class_enrichment_severe_decrease_250404.pdf', 
#        width = 8, height = 6)


# merged pathway enrichment plot ------------------------------------------------

# The reference shows the top 5 strongest remission pathways and strongest active-disease pathways.
enrich_merge <- bind_rows(
  class_enrich_severe_decrease %>%
    slice_min(adjustedpvalue, n = 5, with_ties = FALSE) %>%
    mutate(direction = 'decrease'),
  class_enrich_severe_increase %>%
    slice_min(adjustedpvalue, n = 5, with_ties = FALSE) %>%
    mutate(direction = 'increase')
) %>%
  mutate(
    # Simplify the names of the pathways for display in the plot.
    logP = -log(adjustedpvalue),
    signed_logP = if_else(direction == 'decrease', -logP, logP),
    display_name = recode(
      name,
      'Tryptophan Metabolism' = 'Tryptophan metabolism',
      'Fructose, Mannose and Galactose Metabolism' = 'Hexose metabolism',
      'TCA Cycle' = 'TCA cycle',
      'Phosphatidylethanolamines' = 'PE lipids',
      'Methionine, Cysteine, SAM and Taurine Metabolism' =
        'Sulfur amino acid metabolism',
      'Pyrimidine Metabolism, Uracil containing' = 'Pyrimidine metabolism',
      .default = name
    )
  )

enrich_merge <- bind_rows(
  enrich_merge %>%
    filter(direction == 'decrease'),
  enrich_merge %>%
    filter(direction == 'increase') %>%
    arrange(logP)
) %>%
  mutate(y = rev(seq_len(n())))

lim <- ceiling(max(enrich_merge$logP))
y_top <- max(enrich_merge$y)
top_pathway <- enrich_merge %>%
  filter(direction == 'decrease') %>%
  slice_max(logP, n = 1, with_ties = FALSE) %>%
  pull(display_name)

enrich_merge <- enrich_merge %>%
  mutate(
    fontface = if_else(display_name == top_pathway, 'bold', 'plain'),
    label_x = if_else(direction == 'decrease', lim, -lim),
    label_hjust = if_else(direction == 'decrease', 1, 0)
  )

temp_plot <- ggplot(enrich_merge, aes(x = signed_logP, y = y)) +
  # Dotted leaders occupy the side opposite each bar.
  geom_segment(
    aes(x = 0, xend = label_x, yend = y),
    linetype = 'dotted', lineend = 'round',
    colour = '#a6a6a6', linewidth = 0.75
  ) +
  geom_col(aes(fill = direction), width = 0.90, orientation = 'y') +
  annotate(
    'segment', x = 0, xend = 0, y = 0.35, yend = y_top + 0.60,
    colour = 'black', linewidth = 0.8
  ) +
  # White label backgrounds mask the leaders beneath the pathway names.
  geom_label(
    aes(x = label_x, label = display_name,
        hjust = label_hjust, fontface = fontface),
    fill = 'white', colour = 'black', label.size = 0,
    label.padding = grid::unit(0.06, 'lines'), size = 5.0
  ) +
  scale_fill_manual(
    values = c('decrease' = '#159f7d', 'increase' = '#df6426'),
    guide = 'none'
  ) +
  scale_x_continuous(
    limits = c(-lim, lim),
    breaks = seq(-12, 12, by = 4),
    labels = function(x) abs(x),
    expand = expansion(mult = 0)
  ) +
  annotate(
    'segment', x = -0.25, xend = -lim, y = y_top + 0.75, yend = y_top + 0.75,
    linewidth = 0.8,
    arrow = arrow(length = grid::unit(0.16, 'inch'), type = 'closed')
  ) +
  annotate(
    'segment', x = 0.25, xend = lim, y = y_top + 0.75, yend = y_top + 0.75,
    linewidth = 0.8,
    arrow = arrow(length = grid::unit(0.16, 'inch'), type = 'closed')
  ) +
  annotate(
    'text', x = -lim / 2, y = y_top + 1.25,
    label = 'Higher in remission', size = 5.2
  ) +
  annotate(
    'text', x = lim / 2, y = y_top + 1.25,
    label = 'Higher in active disease', size = 5.2
  ) +
  labs(
    title = paste(top_pathway, 'is\nenriched in remission'),
    x = expression(-log[10](P[adjusted])),
    y = NULL
  ) +
  coord_cartesian(clip = 'off', ylim = c(0.35, y_top + 1.55)) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      face = 'bold', hjust = 0.5, size = 26,
      lineheight = 1.05, margin = margin(b = 15)
    ),
    axis.title.x = element_text(size = 22, margin = margin(t = 8)),
    axis.text.x = element_text(size = 16, colour = 'black'),
    axis.text.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.8),
    axis.line.y = element_blank(),
    plot.margin = margin(8, 8, 5, 8)
  )

ggsave(temp_plot,
       filename = '~/Project/00_IBD_project/Figure/250326/Figure5/plot_class_enrichment_merged_250404.pdf',
       width = 8, height = 8)


dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE)
temp_data <- enrich_merge %>% 
  dplyr::select(name, direction, adjustedpvalue, logP, signed_logP, display_name) %>% 
  rename('pathway' = 'name',
         'direction_vs_remission' = 'direction',
         'adjusted_pvalue' = 'adjustedpvalue',
         'log10_adjusted_pvalue' = 'logP',
         'signed_log10_adjusted_pvalue' = 'signed_logP',
         'display_name_in_plot' = 'display_name')

readr::write_csv(temp_data, 
                   file = '~/Project/00_IBD_project/Data/20260718_source_data/Fig5c_differential_met_pathway_enrichment_260721.csv')





################################################################################
# Tryptophan metabolites across disease severity ------------------------------

library(tidyverse)
library(cowplot)

load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/result_maaslin2_severity_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/object_CD_followup_250401.RData')


# severe vs. remission: determine the order
result_maaslin2_severe <- result_maaslin2_severity %>% 
  filter(severity_class == 'moderate_severe') %>% 
  filter(metabolon_subclass == 'Tryptophan Metabolism') %>% 
  mutate(direction = case_when(coef > 0 ~ 'increase',
                               coef < 0 ~ 'decrease')) %>%
  arrange(direction) %>% 
  mutate(index = seq(n(), 1, -1)) %>% 
  mutate(shape = 'moderate_severe')

# tryptophan_order <- result_maaslin2_severe %>% 
#   arrange(index) %>% 
#   pull(compound_name)

# Modify the compound names to same style
tryptophan_order <- c(
  'Tryptophan',
  'Indole-3-lactic acid',
  'Indole-3-carboxaldehyde',
  'Indole-3-acetaldehyde',
  '5-Methoxy-3-indoleacetic acid',
  'Indole-3-Acetamide',
  'Skatole',
  'Indole-3-propionic acid',
  'Indole-3-acetic acid',
  'Kynurenine',
  '5-Hydroxytryptophan',
  '5-Hydroxyindole',
  '2-Oxindole',
  'Kynurenic acid',
  '3-Indoxyl sulfate',
  '3-Methyl-2-oxindole'
)

# add display_name, group, direction, significance, y, y_point columns
tryptophan_result <- result_maaslin2_severity %>%
  filter(
    metabolon_subclass == 'Tryptophan Metabolism',
    severity_class %in% c('moderate_severe', 'mild')
  ) %>%
  mutate(
    display_name = recode(
      compound_name,
      'Indole-3-Carboxaldehyde' = 'Indole-3-carboxaldehyde',
      'Indolepropionic acid' = 'Indole-3-propionic acid',
      'Indole-3-Acetic Acid' = 'Indole-3-acetic acid',
      '3-Methyl-2-Oxindole' = '3-Methyl-2-oxindole',
      .default = compound_name
    ),
    display_name = factor(display_name, levels = rev(tryptophan_order)),
    group = recode(
      severity_class,
      'mild' = 'Mild vs remission',
      'moderate_severe' = 'Moderate/severe vs remission'
    ),
    direction = if_else(coef < 0,
                        'Higher in remission',
                        'Higher in active disease'),
    significance = case_when(
      qval >= 0.05 ~ 'ns',
      qval < 0.001 ~ '***',
      qval < 0.01 ~ '**',
      qval < 0.05 ~ '*',
      TRUE ~ 'ns'
    ),
    y = as.numeric(display_name),
    y_point = y + if_else(severity_class == 'mild', 0.10, -0.10)
  )

# Mean abundance in each severity group, expressed relative to remission.
tryptophan_features <- tryptophan_result %>%
  distinct(variable_id, display_name) %>%
  arrange(match(as.character(display_name), tryptophan_order))

expression_matrix <- object_CD_followup@expression_data[
  match(tryptophan_features$variable_id,
        rownames(object_CD_followup@expression_data)),
  , drop = FALSE
]

sample_table <- object_CD_followup@sample_info %>%
  select(sample_id, severity_class2)

group_mean <- function(group_name) {
  sample_ids <- sample_table %>%
    filter(severity_class2 == group_name) %>%
    pull(sample_id)
  rowMeans(expression_matrix[, intersect(colnames(expression_matrix), sample_ids),
                             drop = FALSE], na.rm = TRUE)
}

remission_mean <- group_mean('remission')

abundance_data <- tibble(
  variable_id = tryptophan_features$variable_id,
  display_name = tryptophan_features$display_name,
  moderate_severe = group_mean('moderate_severe') / remission_mean,
  mild = group_mean('mild') / remission_mean
) %>%
  pivot_longer(
    cols = c(moderate_severe, mild),
    names_to = 'severity_class', values_to = 'relative_abundance'
  ) %>%
  mutate(
    severity_class = factor(
      severity_class,
      levels = c('moderate_severe', 'mild'),
      labels = c('Moderate/\nsevere', 'Mild')
    ),
    display_name = factor(display_name, levels = rev(tryptophan_order))
  ) %>%
  left_join(
    tryptophan_result %>%
      transmute(
        variable_id,
        severity_class = factor(
          severity_class,
          levels = c('moderate_severe', 'mild'),
          labels = c('Moderate/\nsevere', 'Mild')
        ),
        significance
      ),
    by = c('variable_id', 'severity_class')
  )

# Left panel: relative abundance heatmap with model significance.
heatmap_plot <- ggplot(
  abundance_data,
  aes(x = severity_class, y = display_name, fill = relative_abundance)
) +
  geom_tile(colour = 'black', linewidth = 0.55) +
  geom_text(aes(label = significance), size = 4.0, fontface = 'bold') +
  scale_fill_gradient2(
    limits = c(0.7, 1.5), midpoint = 1.0,
    low = '#b9dfda', mid = '#ffffff', high = '#f9bd73',
    oob = scales::squish,
    breaks = c(0.7, 1.1, 1.5),
    name = 'Relative abundance\n(relative to Remission)'
  ) +
  scale_x_discrete(position = 'top', expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(
    fill = guide_colorbar(
      title.position = 'top', title.hjust = 0.5,
      barwidth = grid::unit(2.0, 'cm'),
      barheight = grid::unit(0.35, 'cm'),
      ticks = FALSE
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = 12) +
  theme(
    axis.text.x.top = element_text(
      size = 12, colour = 'black', lineheight = 0.95,
      margin = margin(b = 6)
    ),
    axis.text.y = element_text(
      size = 11.5, colour = 'black', hjust = 1,
      margin = margin(r = 5)
    ),
    legend.position = 'bottom',
    legend.title = element_text(size = 10.5, hjust = 0.5, lineheight = 0.95),
    legend.text = element_text(size = 9.5),
    legend.margin = margin(t = 7),
    plot.margin = margin(4, 5, 2, 5)
  )


# Right panel: paired coefficient lollipops for mild and moderate/severe CD.
coefficient_plot <- ggplot(
  tryptophan_result,
  aes(y = y_point, colour = direction, shape = group)
) +
  geom_segment(
    aes(x = 0, xend = coef, yend = y_point),
    linewidth = 1.0, show.legend = TRUE
  ) +
  geom_point(aes(x = coef), size = 3.5, stroke = 0.6) +
  geom_vline(xintercept = 0, colour = 'black', linewidth = 0.6) +
  scale_colour_manual(
    values = c(
      'Higher in remission' = '#119f7d',
      'Higher in active disease' = '#dc6229'
    ),
    breaks = c('Higher in remission', 'Higher in active disease'),
    name = 'Direction'
  ) +
  scale_shape_manual(
    values = c(
      'Mild vs remission' = 18,
      'Moderate/severe vs remission' = 16
    ),
    labels = c(
      'Mild vs remission',
      'Moderate/severe\nvs remission'
    ),
    name = 'Group'
  ) +
  scale_x_continuous(
    limits = c(-0.46, 0.52),
    breaks = seq(-0.4, 0.4, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, 16.5), expand = c(0, 0)) +
  annotate(
    'text', x = 0.025, y = c(10.4, 9.9, 9.4, 8.9),
    label = c(
      'ns~~italic(P) > 0.05',
      "'*'~~italic(P) < 0.05",
      "'**'~~italic(P) < 0.01",
      "'***'~~italic(P) < 0.001"
    ),
    parse = TRUE, hjust = 0, size = 3.5, colour = 'black'
  ) +
  labs(x = 'Coefficient', y = NULL, title = 'Coefficient') +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 16)),
    shape = guide_legend(order = 2, override.aes = list(colour = 'black'))
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 5)),
    axis.title.x = element_text(size = 13, margin = margin(t = 7)),
    axis.text.x = element_text(size = 10.5, colour = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.65),
    legend.position = 'inside',
    legend.position.inside = c(0.73, 0.97),
    legend.justification = c(0.5, 1),
    legend.direction = 'vertical',
    legend.box = 'vertical',
    legend.box.spacing = grid::unit(0.05, 'cm'),
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 10.0, lineheight = 0.95),
    legend.key.height = grid::unit(0.46, 'cm'),
    legend.key.width = grid::unit(0.65, 'cm'),
    plot.margin = margin(4, 5, 45, 5)
  )



# dir.create('~/Project/00_IBD_project/Data/20260718_source_data/', showWarnings = FALSE, recursive = TRUE))
# temp_source_data <- abundance_data %>% 
#   left_join(
#     tryptophan_result %>% 
#       mutate(
#         severity_class = factor(
#           severity_class,
#           levels = c('moderate_severe', 'mild'),
#           labels = c('Moderate/\nsevere', 'Mild')
#         ),
#         display_name = factor(display_name, levels = rev(tryptophan_order))
#       ) %>% 
#       dplyr::select(variable_id, display_name, severity_class, coef, qval) %>% 
#       rename('coefficient' = 'coef',
#              'q_value' = 'qval'),
#     by = c('variable_id', 'display_name', 'severity_class')
#   )
# 
# readr::write_csv(temp_source_data, 
#                  file = '~/Project/00_IBD_project/Data/20260718_source_data/Fig5d_tryptophan_metabolite_abundance_260721.csv')


# Place the color scale and coefficient axis title in a shared footer so the
# plotting regions above have identical heights and every row stays aligned.
heatmap_guides <- cowplot::get_plot_component(
  heatmap_plot + theme(legend.position = 'bottom'),
  'guide-box', return_all = TRUE
)
heatmap_legend <- Filter(
  function(x) !inherits(x, 'zeroGrob'),
  heatmap_guides
)[[1]]

heatmap_panel <- heatmap_plot +
  guides(fill = 'none') +
  theme(plot.margin = margin(4, 5, 0, 5))

coefficient_panel <- coefficient_plot +
  labs(x = NULL) +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(4, 5, 0, 5)
  )

# Combine the heatmap and coefficient panels into a single row with aligned axes.
panel_row <- cowplot::plot_grid(
  heatmap_panel, coefficient_panel,
  nrow = 1, align = 'h', axis = 'tb',
  rel_widths = c(1.05, 1.15)
)

# add a shared x-axis title for the coefficient panel
coefficient_footer <- cowplot::ggdraw() +
  cowplot::draw_label('Coefficient', size = 13, x = 0.5, y = 0.82)

footer_row <- cowplot::plot_grid(
  heatmap_legend, coefficient_footer,
  nrow = 1, rel_widths = c(1.05, 1.15)
)

combined_panels <- cowplot::plot_grid(
  panel_row, footer_row,
  ncol = 1, rel_heights = c(0.89, 0.11)
)

figure_title <- cowplot::ggdraw() +
  theme(plot.background = element_rect(fill = 'white', colour = NA)) +
  cowplot::draw_label(
    'Tryptophan metabolites shift with disease severity',
    fontface = 'bold', size = 18, x = 0.5, hjust = 0.5
  )

# Combine the title and the panels into a single figure with a white background.
temp_plot_merge <- cowplot::plot_grid(
  figure_title, combined_panels,
  ncol = 1, rel_heights = c(0.075, 0.925)
)

temp_plot_merge <- cowplot::ggdraw(temp_plot_merge) +
  theme(plot.background = element_rect(fill = 'white', colour = NA))

ggsave(
  temp_plot_merge,
  filename = '~/Project/00_IBD_project/Figure/250326/Figure5/tryptophan_severity_260721.pdf',
  width = 10, height = 9.5
)




################################################################################
# Microbial Indole Balance Score (MIBS) ----------------------------------------
#   - protective metabolites (microbial AhR-active indoles):
#       Indole-3-lactic acid, Indole-3-acetaldehyde, Indole-3-Carboxaldehyde, Indolepropionic acid
#   - adverse metabolites:
#       5-Hydroxyindole, 2-Oxindole

# Part 0: setup ----------------------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/')

library(tidyverse)
library(tidymass)
suppressPackageStartupMessages({
  library(lme4)        # logistic / Gaussian mixed model
  library(broom.mixed) # tidy mixed-model output
  library(cowplot)     # combine plots
  library(lmerTest)    # Satterthwaite p-values for lmer
  library(vegan)       # adonis2 (PERMANOVA)
})
set.seed(20260623)

rm(list = ls()); gc()

# Part 1: load raw data and build dat / vinfo -----------------------------------------------------

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')
load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/object_CD_followup_250401.RData')
load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/result_maaslin2_severity_250401.RData')

# log2 -> z-score (single source of truth)

log2z <- function(x, c = 1) {
  x2 <- log2(x + c)
  as.numeric((x2 - mean(x2, na.rm = TRUE)) / sd(x2, na.rm = TRUE))
}

# 1.1 tryptophan-metabolism panel: variable_ids + abundances (samples x compounds)
temp_trp_met <- result_maaslin2_severity %>%
  dplyr::filter(metabolon_subclass == 'Tryptophan Metabolism')

trp_all <- temp_trp_met %>% dplyr::pull(variable_id)

met_wide <- object_CD_followup %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% trp_all) %>%
  extract_expression_data() %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  as_tibble()

vinfo <- object_CD_followup %>%
  extract_variable_info() %>%
  dplyr::filter(variable_id %in% trp_all) %>%
  dplyr::select(variable_id, Compound.name)

# 1.2 sample-level outcome + covariates
sample_info_raw <- object_CD_followup %>% extract_sample_info()
patient_id_col <- "patient_id"
sinfo <- sample_info_raw %>%
  dplyr::select(dplyr::any_of(c(
    "sample_id", "severity_class", patient_id_col,
    "gender", "age", "race", "use_antibiotics")))

dat <- sinfo %>%
  dplyr::left_join(met_wide, by = "sample_id") %>%
  dplyr::mutate(
    patient_id = .data[[patient_id_col]],
    remission = dplyr::if_else(severity_class == "remission", 1L, 0L)
  )

save(dat, vinfo, patient_id_col,
     file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_input_20260721.RData")


# Part 2: Develop Microbial Indole Balance Score (MIBS) ------------------------
#   members fixed by biochemical knowledge (protective vs adverse)

# protective indoles
protective_indoles <- c(
  "Indole-3-lactic acid",
  "Indole-3-acetaldehyde",
  "Indole-3-Carboxaldehyde",
  "Indolepropionic acid"
)

# adverse indoles
adverse_indoles_core <- c(
  "5-Hydroxyindole",
  "2-Oxindole"
)


name_to_vid <- function(compound_names) {
  vid <- vinfo$variable_id[match(compound_names, vinfo$Compound.name)]
  missing <- compound_names[is.na(vid)]
  if (length(missing) > 0) {
    warning("Compound.name not found in vinfo, dropped from score: ", paste(missing, collapse = ", "))
  }
  vid[!is.na(vid)]
}

protective_vid <- name_to_vid(protective_indoles)
adverse_vid    <- name_to_vid(adverse_indoles_core)



cat("\nMIBS protective members matched:\n")
print(vinfo$Compound.name[match(protective_vid, vinfo$variable_id)])
cat("\nMIBS adverse members matched:\n")
print(vinfo$Compound.name[match(adverse_vid, vinfo$variable_id)])

# 2.1 z-scored matrix for score building
metabolite_cols <- intersect(trp_all, colnames(dat))

# Save mean/sd used by z-score transformation (based on log2(x + 1)) for later reuse.
zscore_params_all_metabolites <- purrr::map_dfr(metabolite_cols, function(vid) {
  x_log2 <- log2(dat[[vid]] + 1)
  tibble::tibble(
    variable_id = vid,
    Compound.name = vinfo$Compound.name[match(vid, vinfo$variable_id)],
    z_mean = mean(x_log2, na.rm = TRUE),
    z_sd = sd(x_log2, na.rm = TRUE),
    n_non_na = sum(!is.na(x_log2))
  )
})

save(
  zscore_params_all_metabolites,
  file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_zscore_params_all_metabolites_20260721.RData"
)

score_members_vid <- Reduce(union, list(protective_vid, adverse_vid))
score_members_vid <- intersect(score_members_vid, colnames(dat))

dz <- dat %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(metabolite_cols), log2z, .names = "z_{.col}"))

zcol <- function(vids) paste0("z_", intersect(vids, score_members_vid))

# 2.1.1 main: MIBS = mean(protective z) - mean(adverse z)
dz <- dz %>%
  dplyr::mutate(
    protective_score              = rowMeans(dplyr::across(dplyr::all_of(zcol(protective_vid))), na.rm = TRUE),
    adverse_score                 = rowMeans(dplyr::across(dplyr::all_of(zcol(adverse_vid))), na.rm = TRUE),
    microbial_indole_balance_score = protective_score - adverse_score,  # higher = healthier
  )


table_score <- dz %>%
  dplyr::select(sample_id, severity_class, microbial_indole_balance_score,
                protective_score, adverse_score, everything())

save(table_score,
     file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_table_20260721.RData")


dz <- dz %>%
  dplyr::mutate(severity_class = factor(severity_class,
                                        levels = c("remission", "mild", "moderate", "severe"), ordered = TRUE)) %>% 
  dplyr::mutate(severity_class2 = case_when(
    severity_class == "remission" ~ "remission",
    severity_class == "mild" ~ "mild",
    severity_class %in% c("moderate", "severe") ~ "mod_severe",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::mutate(severity_class2 = factor(severity_class2,
                                         levels = c("remission", "mild", "mod_severe"), ordered = TRUE))

# scores to run with the same evaluation workflow
score_vars <- c(
  "microbial_indole_balance_score",
  "protective_score",
  "adverse_score"
)

# same model covariates as the 0615 main analysis
covs <- intersect(c("gender", "age", "use_antibiotics"), colnames(dz))

save(dz, score_vars, covs,
     file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_evaluation_input_20260721.RData")


rm(list = ls()); gc()


# Part 3: Evaluation -----------------------------------------------------------

# 3.1: cross-sectional evaluation: score across severity (expect remission > mild > moderate/severe) ------

load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_evaluation_input_20260721.RData")

# Individual boxplots across severity for each score 
score_labels <- c(
  microbial_indole_balance_score   = "Microbial Indole Balance Score (MIBS, main)",
  protective_score                 = "Protective Score (mean z of protective indoles)",
  adverse_score                    = "Adverse Score (mean z of adverse indoles)"
)

# Add wilcoxon test for comparisons between severity groups
severity_colors <- c(
  remission = "#169e77",
  mild = "#7570b3",
  mod_severe = "#d76327"
)

# load the evaluation results from the previous step (20260623)
load("~/Project/00_IBD_project/Data/20260615_tryptophan_met_severity/mibs_score_evaluation_20260623.RData")
score_summary <- purrr::map_dfr(eval_list, "summary_row")

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", format(signif(x, 2), scientific = TRUE))
}

fmt_mean <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

score_plot_labels <- score_summary %>%
  dplyr::mutate(
    plot_label = paste0(
      "P=", fmt_p(trend_p),
      "\nMean remission=", fmt_mean(mean_remission),
      "\nMean mild=", fmt_mean(mean_mild),
      "\nMean mod_severe=", fmt_mean(mean_mod_severe)
    )
  ) %>%
  dplyr::select(score, plot_label)

plot_one_score_change_order_wilcoxon <- function(score_var) {
  d <- dz %>%
    dplyr::select(severity_class2, value = dplyr::all_of(score_var)) %>%
    dplyr::filter(!is.na(value), !is.na(severity_class2)) %>%
    dplyr::mutate(severity_class2 = factor(severity_class2, levels = c("mod_severe", "mild", "remission"), ordered = TRUE))
  
  sub_lab <- score_plot_labels$plot_label[match(score_var, score_plot_labels$score)]
  if (length(sub_lab) == 0 || is.na(sub_lab)) sub_lab <- "P=NA"
  
  # pairwise Wilcoxon comparisons between the three severity groups
  sev_levels  <- c("mod_severe", "mild", "remission")
  present     <- sev_levels[sev_levels %in% unique(as.character(d$severity_class2))]
  comparisons <- if (length(present) >= 2) utils::combn(present, 2, simplify = FALSE) else list()
  
  p <- ggplot(d, aes(x = severity_class2, y = value, fill = severity_class2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(color = severity_class2), position = position_jitter(width = 0.2), color = "black",
               alpha = 0.5, size = 3, shape = 21) +
    annotate("text", x = Inf, y = Inf, label = sub_lab, hjust = 1.02, vjust = 1.1, size = 3.3) +
    scale_fill_manual(values = severity_colors, drop = FALSE) +
    scale_color_manual(values = severity_colors, drop = FALSE) +
    labs(
      x = "Severity (wPCDAI)",
      y = "Score",
      title = score_labels[[score_var]]
    ) +
    ZZWtool::ZZWTheme() +
    theme(legend.position = "none")
  
  # add pairwise Wilcoxon test brackets with p-values
  if (length(comparisons) > 0) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method      = "wilcox.test",
      label       = "p.format",
      tip.length  = 0.01,
      size        = 3
    )
  }
  
  print(p)
  invisible(p)
}

score_boxplots <- purrr::set_names(score_vars) %>%
  purrr::map(plot_one_score_change_order_wilcoxon)

dir.create("~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/score_boxplots_260710", showWarnings = FALSE, recursive = TRUE)

walk(score_boxplots, function(p) {
  ggsave(
    filename = file.path('~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/score_boxplots_260710/', 
                         paste0("mibs_scores_boxplot_", p$labels$title, ".pdf")),
    plot = p,
    width = 6, height = 6
  )
})


# Explicit requested panel handles
fig5_f  <- score_boxplots[["microbial_indole_balance_score"]]
ed9_b1  <- score_boxplots[["protective_score"]]
ed9_b2  <- score_boxplots[["adverse_score"]]

# dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
# readr::write_csv(fig5_f$data, file = "~/Project/00_IBD_project/Data/20260718_source_data/Fig5f_MIBS_score_260721.csv")
# readr::write_csv(ed9_b1$data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9b_protective_score_260721.csv")
# readr::write_csv(ed9_b2$data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9b_adverse_score_260721.csv")



# 3.2: Longitudinal symptom-activity change vs MIBS change (wPCDAI line) ------------------

rm(list = ls()); gc()

load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_input_20260721.RData")
load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_evaluation_input_20260721.RData")
load('~/Project/00_IBD_project/Data/20250401_longitudinal_severity_analysis_update/object_CD_followup_250401.RData')

# Build visit-level table (patient_id, visit_order, severity_class, wPCDAI) 
sample_info_vis <- object_CD_followup %>% extract_sample_info()

visit_id_col    <- "visit_encounter_id"
visit_order_col <- "visit_month"
wpcdai_col      <- "wPCDAI"
sample_id_col   <-"sample_id"
patient_id_col  <-  "patient_id"


visit_df <- sample_info_vis %>%
  dplyr::mutate(
    patient_id = .data[[patient_id_col]],
    sample_id  = if (is.na(sample_id_col)) sample_id else .data[[sample_id_col]],
    visit_id   = if (is.na(visit_id_col)) {
      if (is.na(sample_id_col)) dplyr::row_number() else .data[[sample_id_col]]
    } else .data[[visit_id_col]],
    visit_order = if (is.na(visit_order_col)) dplyr::row_number() else .data[[visit_order_col]],
    severity_class = as.character(severity_class),
    wPCDAI = if (is.na(wpcdai_col)) NA_real_ else suppressWarnings(as.numeric(.data[[wpcdai_col]]))
  ) %>%
  dplyr::filter(!is.na(patient_id))

# attach MIBS (+ all score versions) from dz, ordered within patient by visit
score_cols_present <- intersect(score_vars, colnames(dz))
cov_cols_present   <- intersect(c("gender", "age", "race", "use_antibiotics"), colnames(dz))

dz_for_join <- dz %>%
  dplyr::select(sample_id, dplyr::any_of(c(score_cols_present, cov_cols_present)))

long_df <- visit_df %>%
  dplyr::select(sample_id, patient_id, visit_id, visit_order, severity_class, severity_class2, wPCDAI) %>%
  dplyr::left_join(dz_for_join, by = "sample_id") %>%
  dplyr::group_by(patient_id) %>%
  dplyr::arrange(visit_order, .by_group = TRUE) %>%
  dplyr::ungroup()

save(long_df,
     file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_longitudinal_20260721.RData")


main_score <- "microbial_indole_balance_score"

# adjacent visit-pair table 
pair_df <- long_df %>%
  dplyr::group_by(patient_id) %>%
  dplyr::arrange(visit_order, .by_group = TRUE) %>%
  dplyr::mutate(
    delta_wPCDAI = wPCDAI - dplyr::lag(wPCDAI),
    dplyr::across(dplyr::all_of(score_cols_present), ~ . - dplyr::lag(.), .names = "delta_{.col}"),
    severity_class_prev = dplyr::lag(severity_class)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(severity_class_curr = severity_class) %>%
  dplyr::filter(!is.na(delta_wPCDAI))

# delta wPCDAI vs delta MIBS scatter, adjacent visit pairs 
run_delta_delta <- function(score_var = main_score) {
  dcol <- paste0("delta_", score_var)
  d <- pair_df %>% dplyr::filter(!is.na(delta_wPCDAI), !is.na(.data[[dcol]]))
  form <- reformulate(c("delta_wPCDAI", "(1 | patient_id)"), response = dcol)
  m <- lmerTest::lmer(form, data = d)
  broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
    dplyr::filter(term == "delta_wPCDAI") %>%
    dplyr::mutate(score = score_var)
}

delta_delta_slope_mibs <- run_delta_delta(main_score)
cat("\n3.2 delta-delta slope (delta_MIBS ~ delta_wPCDAI), expect NEGATIVE:\n")
print(delta_delta_slope_mibs)

delta_delta_ols_mibs <- lm(delta_microbial_indole_balance_score ~ delta_wPCDAI, data = pair_df)
delta_ols_sum <- summary(delta_delta_ols_mibs)

# pearson correlation for delta-delta (scatter) check
delta_cor <- cor(pair_df$delta_wPCDAI, pair_df$delta_microbial_indole_balance_score, use = "complete.obs")
cat("\n3.2 delta-delta correlation (Pearson r): ", round(delta_cor, 3), "\n", sep = "")

delta_intercept <- unname(coef(delta_delta_ols_mibs)[1])
delta_slope <- unname(coef(delta_delta_ols_mibs)[2])
delta_p <- delta_ols_sum$coefficients["delta_wPCDAI", "Pr(>|t|)"]
delta_formula_label <- sprintf(
  "y = %.3f %s %.3f x\nP = %s",
  delta_intercept,
  ifelse(delta_slope < 0, "-", "+"),
  abs(delta_slope),
  ifelse(is.na(delta_p), "NA", format(signif(delta_p, 2), scientific = TRUE))
)

plot_delta_wpcdai_vs_delta_mibs <- pair_df %>%
  dplyr::filter(!is.na(delta_wPCDAI), !is.na(delta_microbial_indole_balance_score)) %>%
  ggplot(aes(delta_wPCDAI, delta_microbial_indole_balance_score)) +
  geom_point(alpha = 0.8, size = 3, color = "#4575b4") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", se = TRUE) +
  annotate("text", x = -Inf, y = -Inf, label = delta_formula_label,
           hjust = -0.05, vjust = -0.3, size = 4) +
  labs(x = "Delta wPCDAI",
       y = "Delta MIBS") +
  ZZWtool::ZZWTheme() +
  theme(legend.position = "none")

ggsave(plot_delta_wpcdai_vs_delta_mibs, 
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_delta_delta_plot_20260630.pdf", 
       width = 6, height = 6)

# dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- plot_delta_wpcdai_vs_delta_mibs$data %>% 
#     dplyr::select(sample_id, patient_id, visit_order, severity_class_prev, severity_class_curr,delta_wPCDAI, delta_microbial_indole_balance_score)

# readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/Fig5g_delta_wPCDAI_vs_delta_MIBS_260721.csv")



# 3.3: within-patient PAIRED transition plots (prev visit -> curr visit), MIBS --------------------
#   Severity collapsed to 3 levels: Remission / Mild / Moderate/Severe (moderate + severe merged).
#   9 directed sub-classes grouped into 3 macro-classes:
#     Improving : Mod/Sev->Remission, Mod/Sev->Mild, Mild->Remission
#     Stable    : Remission->Remission, Mild->Mild, Mod/Sev->Mod/Sev
#     Worsening : Remission->Mild, Remission->Mod/Sev, Mild->Mod/Sev
#   For each sub-class: paired point-line (MIBS at prev vs curr, one line per visit-pair) + boxplot,
#   with a paired Wilcoxon signed-rank test annotated (P, and mean MIBS at each timepoint).

load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_longitudinal_20260721.RData")

sev3_of <- function(x) {
  dplyr::case_when(
    x == "remission"               ~ "Remission",
    x == "mild"                    ~ "Mild",
    x %in% c("moderate", "severe") ~ "Moderate/Severe",
    TRUE                           ~ NA_character_
  )
}
sev3_lvls <- c("Remission", "Mild", "Moderate/Severe")

paired_trans <- long_df %>%
  dplyr::filter(!is.na(microbial_indole_balance_score)) %>%
  dplyr::group_by(patient_id) %>%
  dplyr::arrange(visit_order, .by_group = TRUE) %>%
  dplyr::mutate(
    sev3_curr  = sev3_of(as.character(severity_class)),
    sev3_prev  = dplyr::lag(sev3_curr),
    mibs_curr  = microbial_indole_balance_score,
    mibs_prev  = dplyr::lag(microbial_indole_balance_score),
    visit_order_prev = dplyr::lag(visit_order)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(sev3_prev), !is.na(sev3_curr),
                !is.na(mibs_prev), !is.na(mibs_curr)) %>%
  dplyr::mutate(
    sub_class = paste0(sev3_prev, " \u2192 ", sev3_curr),   # "prev -> curr"
    macro_class = dplyr::case_when(
      (sev3_prev == "Moderate/Severe" & sev3_curr == "Remission") |
        (sev3_prev == "Moderate/Severe" & sev3_curr == "Mild") |
        (sev3_prev == "Mild"            & sev3_curr == "Remission") ~ "Improving",
      sev3_prev == sev3_curr ~ "Stable",
      (sev3_prev == "Remission" & sev3_curr == "Mild") |
        (sev3_prev == "Remission" & sev3_curr == "Moderate/Severe") |
        (sev3_prev == "Mild"      & sev3_curr == "Moderate/Severe") ~ "Worsening",
      TRUE ~ NA_character_
    ),
    macro_class = factor(macro_class, levels = c("Improving", "Stable", "Worsening")),
    pair_id = dplyr::row_number()
  ) %>%
  dplyr::filter(!is.na(macro_class))

# fixed sub-class ordering within each macro-class (for consistent facet order)
sub_levels <- c(
  # Improving
  "Moderate/Severe \u2192 Remission", "Moderate/Severe \u2192 Mild", "Mild \u2192 Remission",
  # Stable
  "Remission \u2192 Remission", "Mild \u2192 Mild", "Moderate/Severe \u2192 Moderate/Severe",
  # Worsening
  "Remission \u2192 Mild", "Remission \u2192 Moderate/Severe", "Mild \u2192 Moderate/Severe"
)
paired_trans <- paired_trans %>%
  dplyr::mutate(sub_class = factor(sub_class, levels = intersect(sub_levels, unique(sub_class))))

cat("\n3.4 paired transition counts (visit-pairs per sub-class):\n")
print(table(paired_trans$macro_class, paired_trans$sub_class))

macro_colors <- c("Improving" = "#91bfdb", "Stable" = "#e6ab02", "Worsening" = "#e7298a")

# per-sub-class paired stats: paired Wilcoxon (prev vs curr) + means
paired_trans_long <- paired_trans %>%
  dplyr::select(pair_id, patient_id, macro_class, sub_class, mibs_prev, mibs_curr) %>%
  tidyr::pivot_longer(c(mibs_prev, mibs_curr),
                      names_to = "timepoint", values_to = "mibs") %>%
  dplyr::mutate(timepoint = factor(dplyr::recode(timepoint,
                                                 mibs_prev = "prev", mibs_curr = "curr"),
                                   levels = c("prev", "curr")))

paired_stats <- paired_trans %>%
  dplyr::group_by(macro_class, sub_class) %>%
  dplyr::group_modify(function(g, key) {
    n_pairs <- nrow(g)
    mean_prev <- mean(g$mibs_prev, na.rm = TRUE)
    mean_curr <- mean(g$mibs_curr, na.rm = TRUE)
    p_wilcox <- if (n_pairs >= 2 && stats::sd(g$mibs_curr - g$mibs_prev, na.rm = TRUE) > 0) {
      tryCatch(stats::wilcox.test(g$mibs_curr, g$mibs_prev, paired = TRUE)$p.value,
               error = function(e) NA_real_)
    } else NA_real_
    tibble::tibble(n_pairs = n_pairs, mean_prev = mean_prev, mean_curr = mean_curr,
                   mean_delta = mean_curr - mean_prev, p_wilcox = p_wilcox)
  }) %>%
  dplyr::ungroup()

cat("\n3.4 paired Wilcoxon (curr vs prev MIBS) per sub-class:\n")
print(paired_stats)

# annotation label per facet: P + means
paired_annot <- paired_stats %>%
  dplyr::mutate(
    label = sprintf("P = %s\nmean: %.2f \u2192 %.2f\n(n = %d)",
                    ifelse(is.na(p_wilcox), "NA", format(signif(p_wilcox, 2), scientific = TRUE)),
                    mean_prev, mean_curr, n_pairs)
  )

# one combined paired point-line + boxplot, faceted by sub-class (grouped by macro) 
build_paired_plot <- function(df_long, df_annot) {
  ggplot(df_long, aes(x = timepoint, y = mibs)) +
    geom_boxplot(aes(fill = macro_class), outlier.shape = NA, alpha = 0.35, width = 0.5) +
    geom_line(aes(group = pair_id), color = "grey55", alpha = 0.5, linewidth = 0.35) +
    geom_point(aes(color = macro_class), size = 1.6, alpha = 0.75) +
    geom_text(data = df_annot,
              aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.1, size = 3, lineheight = 0.95) +
    scale_fill_manual(values = macro_colors, guide = "none") +
    scale_color_manual(values = macro_colors, guide = "none") +
    labs(x = NULL, y = "MIBS") +
    ZZWtool::ZZWTheme()
}

plot_paired_transition <- build_paired_plot(
  paired_trans_long, paired_annot)

plot_paired_transition <- plot_paired_transition +
  facet_wrap(~ macro_class + sub_class, nrow = 1, scales = "free_x")

print(plot_paired_transition)
ggsave(plot_paired_transition,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_paired_transition_20260630.pdf",
       width = 10, height = 5)

# dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- plot_paired_transition$data %>% 
#     dplyr::select(pair_id, patient_id, macro_class, sub_class, timepoint, mibs)
# readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9d_MIBS_paired_transition_260721.csv")


# focused paired boxplot: Mod/Sev→Rem, Mild→Rem, Rem→Rem 
focus_groups <- c(
  "Moderate/Severe → Remission",
  "Mild → Remission",
  "Remission → Remission"
)

paired_trans_long_focus <- paired_trans_long %>%
  dplyr::filter(sub_class %in% focus_groups) %>%
  dplyr::mutate(sub_class = factor(sub_class, levels = focus_groups))

paired_annot_focus <- paired_annot %>%
  dplyr::filter(sub_class %in% focus_groups) %>%
  dplyr::mutate(sub_class = factor(sub_class, levels = focus_groups))

temp_colors <- c("Moderate/Severe → Remission" = "#d76327", "Mild → Remission" = "#7570b3", "Remission → Remission" = "#199e77")


temp_plot <- ggplot(paired_trans_long_focus, aes(x = timepoint, y = mibs)) +
  geom_boxplot(aes(fill = sub_class), outlier.shape = NA, alpha = 0.35, width = 0.5) +
  geom_line(aes(group = pair_id), color = "grey55", alpha = 0.5, linewidth = 0.35) +
  geom_point(aes(color = sub_class), size = 1.6, alpha = 0.75) +
  geom_text(data = paired_annot_focus,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.1, size = 3, lineheight = 0.95) +
  scale_fill_manual(values = temp_colors, guide = "none") +
  scale_color_manual(values = temp_colors, guide = "none") +
  labs(x = NULL, y = "MIBS") +
  ZZWtool::ZZWTheme() +
  facet_wrap(~ sub_class, nrow = 1, scales = "free_x")

save(paired_trans, paired_stats,
     file = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/paired_transition_mibs_20260630.RData")


ggsave(temp_plot,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_paired_transition_focus_20260630.pdf",
       width = 6, height = 6)

dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
temp_source_data <- temp_plot$data
readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/fig5h_MIBS_paired_transition_focus_260721.csv")

rm(list = ls());gc()




# 3.4: evaluate MIBS against CDEIS (endoscopic / mucosal activity) -----------------------------
#   CDEIS = endoscopic severity (mucosal healing). Independent from wPCDAI (symptoms).
#   The MIBS formula and its z-score parameters are FROZEN from the Part 2 wPCDAI-training run
#   (zscore_params_all_metabolites): the CDEIS samples are scored with those constants, never
#   re-centered on the CDEIS distribution -> this is a held-out test, not a re-fit.
#
# Analysis plan:
#     (1) MIBS vs CDEIS_continuous scatter + lm assessment (expect NEGATIVE slope; higher MIBS
#         = healthier = lower CDEIS)
#     (2) logistic models using MIBS as predictor, OR per +1 MIBS, for two CDEIS-severity splits:
#           - remission (CDEIS<3) vs non-remission
#           - moderate/severe (CDEIS>=9) vs rest

rm(list = ls()); gc()

# frozen z-params + score members + source objects + wPCDAI sample set 
load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_zscore_params_all_metabolites_20260721.RData") # zscore_params_all_metabolites
load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_longitudinal_20260721.RData")             # long_df (has sample_id + wPCDAI)
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData') # object_final (full cohort)
load('~/Project/00_IBD_project/Data/20260126_severity_other_criteria/endoscopy_result_table_with_index_cdeis_260127.RData') # endoscopy_result_table_with_index_cdeis

# redeclare the FROZEN MIBS members (same as Part 2 main score)
protective_indoles <- c(
  "Indole-3-lactic acid",
  "Indole-3-acetaldehyde",
  "Indole-3-Carboxaldehyde",
  "Indolepropionic acid"
)
adverse_indoles_core <- c(
  "5-Hydroxyindole",
  "2-Oxindole"
)
score_members <- c(protective_indoles, adverse_indoles_core)

vinfo_final <- object_final %>%
  extract_variable_info() %>%
  dplyr::select(variable_id, Compound.name) %>%
  dplyr::filter(Compound.name %in% score_members)

protective_vid   <- vinfo_final$variable_id[match(protective_indoles, vinfo_final$Compound.name)]
adverse_vid      <- vinfo_final$variable_id[match(adverse_indoles_core, vinfo_final$Compound.name)]
score_member_vid <- c(protective_vid, adverse_vid)

# frozen z params for exactly these members
freeze_params <- zscore_params_all_metabolites %>%
  dplyr::filter(variable_id %in% score_member_vid) %>%
  dplyr::select(variable_id, z_mean, z_sd)

# CDEIS table: CD samples with non-missing CDEIS + remission flag
temp_CDEIS <- endoscopy_result_table_with_index_cdeis %>%
  dplyr::left_join(
    object_final %>% extract_sample_info() %>%
      dplyr::select(sample_id, phenotype_group1),
    by = "sample_id"
  ) %>%
  dplyr::filter(!is.na(CDEIS), phenotype_group1 == "CD") %>%
  dplyr::mutate(severity_cdeis = dplyr::if_else(CDEIS < 3, "remission", "non_remission")) %>%
  dplyr::select(sample_id, visit_encounter_id, CDEIS, severity_cdeis, phenotype_group1, dplyr::everything())

cat("\ntemp_CDEIS: ", nrow(temp_CDEIS), " CD samples with non-missing CDEIS.\n", sep = "")
print(table(temp_CDEIS$severity_cdeis))

# pull CDEIS-sample abundances for the MIBS members, apply the FROZEN z
cdeis_met_wide <- object_final %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% score_member_vid) %>%
  extract_expression_data() %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  as_tibble()

cdeis_covariates <- object_final %>%
  extract_sample_info() %>%
  dplyr::select(dplyr::any_of(c("sample_id", "deidentified_master_patient_id",
                                "gender", "age", "use_antibiotics")))

cdeis_dat <- temp_CDEIS %>%
  dplyr::select(sample_id, visit_encounter_id, CDEIS, severity_cdeis) %>%
  dplyr::rename(CDEIS_continuous = CDEIS) %>%
  dplyr::mutate(severity_cdeis = factor(severity_cdeis, levels = c("non_remission", "remission"))) %>%
  dplyr::left_join(cdeis_met_wide, by = "sample_id") %>%
  dplyr::left_join(cdeis_covariates, by = "sample_id")

# apply the FROZEN (wPCDAI-trained) mean/sd; never re-derived from the CDEIS samples
for (i in seq_len(nrow(freeze_params))) {
  vid <- freeze_params$variable_id[i]
  cdeis_dat[[paste0("z_", vid)]] <-
    (log2(cdeis_dat[[vid]] + 1) - freeze_params$z_mean[i]) / freeze_params$z_sd[i]
}
zcol_frozen <- function(vids) paste0("z_", vids)

cdeis_dat <- cdeis_dat %>%
  dplyr::mutate(
    protective_score_cdeis         = rowMeans(dplyr::across(dplyr::all_of(zcol_frozen(protective_vid))), na.rm = TRUE),
    adverse_score_cdeis            = rowMeans(dplyr::across(dplyr::all_of(zcol_frozen(adverse_vid))), na.rm = TRUE),
    microbial_indole_balance_score = protective_score_cdeis - adverse_score_cdeis  # SAME frozen formula as Part 2
  ) %>%
  dplyr::filter(!is.na(CDEIS_continuous), !is.na(microbial_indole_balance_score))

# CDEIS severity outcomes (pre-specified literature cutoffs: <3 rem, 3-9 mild, 9-12 mod, >12 sev)
cdeis_dat <- cdeis_dat %>%
  dplyr::mutate(
    remission_cdeis     = dplyr::if_else(CDEIS_continuous < 3, 0L, 1L),    # 1 = non-remission (CDEIS>=3)
    cdeis_modsev         = dplyr::if_else(CDEIS_continuous >= 9, 1L, 0L),  # 1 = moderate/severe
    cdeis_modsev_vs_rem = dplyr::case_when(                                # mod/severe (>=9) vs remission (<3); mild (3-9) excluded
      CDEIS_continuous >= 9 ~ 1L,
      CDEIS_continuous < 3  ~ 0L,
      TRUE                  ~ NA_integer_
    )
  )

wpcdai_sample_ids <- long_df %>%
  dplyr::filter(!is.na(wPCDAI)) %>%
  dplyr::pull(sample_id) %>% unique()

cdeis_dat <- cdeis_dat %>%
  dplyr::mutate(has_wpcdai = sample_id %in% wpcdai_sample_ids)

# analysis functions (split into continuous and logistic OR)
covs_cdeis_default <- c("age", "gender", "use_antibiotics")

# continuous: MIBS ~ CDEIS_continuous (lm + Pearson/Spearman + scatter)
analyse_cdeis_continuous <- function(d, set_label, set_tag) {
  d <- d %>% dplyr::filter(!is.na(CDEIS_continuous), !is.na(microbial_indole_balance_score))
  covs_cdeis <- intersect(covs_cdeis_default, colnames(d))
  n_set <- nrow(d)
  cat("\n==================== ", set_label, " (n = ", n_set, ") ====================\n", sep = "")
  
  lm_mibs_cdeis <- lm(reformulate(c("CDEIS_continuous", covs_cdeis),
                                  response = "microbial_indole_balance_score"), data = d)
  tidy_lm <- broom::tidy(lm_mibs_cdeis, conf.int = TRUE) %>%
    dplyr::filter(term == "CDEIS_continuous") %>%
    dplyr::mutate(set = set_label)
  cat("\n(1) MIBS ~ CDEIS_continuous (adjusted), expect NEGATIVE slope:\n")
  print(tidy_lm)
  if (nrow(tidy_lm) > 0 && !is.na(tidy_lm$estimate[1]) && tidy_lm$estimate[1] > 0) {
    warning(set_label, ": MIBS-CDEIS slope is POSITIVE (expected negative). Inspect.")
  }
  
  pear  <- suppressWarnings(cor.test(d$microbial_indole_balance_score, d$CDEIS_continuous, method = "pearson"))
  spear <- suppressWarnings(cor.test(d$microbial_indole_balance_score, d$CDEIS_continuous, method = "spearman"))
  cat(sprintf("    Pearson r = %.3f (p = %.3g);  Spearman rho = %.3f (p = %.3g)\n",
              unname(pear$estimate), pear$p.value, unname(spear$estimate), spear$p.value))
  
  r_lab <- sprintf("Pearson r = %.2f, p = %.3g\nadj. beta = %.3f [%.3f, %.3f]",
                   unname(pear$estimate), pear$p.value,
                   tidy_lm$estimate[1], tidy_lm$conf.low[1], tidy_lm$conf.high[1])
  
  p_scatter <- ggplot(d, aes(CDEIS_continuous, microbial_indole_balance_score)) +
    geom_point(alpha = 0.45, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = r_lab, hjust = 1.05, vjust = 1.3, size = 3.6) +
    labs(x = "CDEIS (endoscopic / mucosal activity)",
         y = "MIBS (frozen, wPCDAI-trained parameters)",
         title = paste0("4.x ", set_label, ": MIBS vs CDEIS")) +
    ZZWtool::ZZWTheme()
  
  print(p_scatter)
  
  list(d = d, n = n_set, lm = tidy_lm, pearson = pear, spearman = spear, r_lab = r_lab, scatter = p_scatter)
}

# logistic OR per +1 MIBS, three CDEIS-severity splits
analyse_cdeis_logit <- function(d, set_label, set_tag) {
  d <- d %>% dplyr::filter(!is.na(CDEIS_continuous), !is.na(microbial_indole_balance_score))
  covs_cdeis <- intersect(covs_cdeis_default, colnames(d))
  n_set <- nrow(d)
  cat("\n==================== ", set_label, " (n = ", n_set, ") — logistic OR ====================\n", sep = "")
  
  # ensure binary outcome columns exist (defensive: create from CDEIS_continuous if absent)
  if (!"remission_cdeis" %in% colnames(d)) {
    d <- d %>% dplyr::mutate(
      remission_cdeis     = dplyr::if_else(CDEIS_continuous < 3,  0L, 1L),   # 1 = non-remission
      cdeis_modsev         = dplyr::if_else(CDEIS_continuous >= 9, 1L, 0L),
      cdeis_modsev_vs_rem = dplyr::case_when(                                # mod/severe (>=9) vs remission (<3); mild excluded
        CDEIS_continuous >= 9 ~ 1L,
        CDEIS_continuous < 3  ~ 0L,
        TRUE                  ~ NA_integer_
      )
    )
  }
  
  run_logit <- function(outcome_var, outcome_label, expected_dir) {
    dd <- d %>% dplyr::filter(!is.na(.data[[outcome_var]]), !is.na(microbial_indole_balance_score))
    n0 <- sum(dd[[outcome_var]] == 0, na.rm = TRUE)
    n1 <- sum(dd[[outcome_var]] == 1, na.rm = TRUE)
    if (n0 < 5 || n1 < 5) {
      cat("\n   [", outcome_label, "] skipped - group too small (n0=", n0, ", n1=", n1, ")\n", sep = "")
      return(tibble::tibble(set = set_label, outcome = outcome_label,
                            estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                            p.value = NA_real_, n0 = n0, n1 = n1))
    }
    m <- glm(reformulate(c("microbial_indole_balance_score", covs_cdeis), response = outcome_var),
             data = dd, family = binomial)
    tout <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term == "microbial_indole_balance_score") %>%
      dplyr::transmute(set = set_label, outcome = outcome_label,
                       estimate, conf.low, conf.high, p.value, n0 = n0, n1 = n1)
    if (nrow(tout) > 0 && !is.na(tout$estimate[1])) {
      if (expected_dir == "OR>1" && tout$estimate[1] < 1)
        warning(set_label, " [", outcome_label, "]: OR < 1 (expected > 1). Inspect.")
      if (expected_dir == "OR<1" && tout$estimate[1] > 1)
        warning(set_label, " [", outcome_label, "]: OR > 1 (expected < 1). Inspect.")
    }
    tout
  }
  
  or_remission     <- run_logit("remission_cdeis",     "non-remission (CDEIS>=3) vs remission",     "OR<1")
  or_modsev        <- run_logit("cdeis_modsev",        "moderate/severe (CDEIS>=9) vs rest",        "OR<1")
  or_modsev_vs_rem <- run_logit("cdeis_modsev_vs_rem", "moderate/severe (CDEIS>=9) vs remission",   "OR<1")
  or_tab <- dplyr::bind_rows(or_remission, or_modsev, or_modsev_vs_rem)
  cat("\n(2) logistic OR per +1 MIBS:\n")
  print(or_tab)
  
  list(n = n_set, or_table = or_tab)
}

res_cdeis_all  <- analyse_cdeis_continuous(cdeis_dat, "all_CDEIS", "all_CDEIS")
res_cdeis_all_logit <- analyse_cdeis_logit(cdeis_dat, "all_CDEIS", "all_CDEIS")


save(temp_CDEIS, freeze_params, cdeis_dat, 
     res_cdeis_all, cdeis_both, cdeis_only, 
     res_cdeis_all_logit, res_cdeis_both_logit, res_cdeis_only_logit, 
     file = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_cdeis_validation_20260625.RData")

# combined OR forest across both datasets + all three thresholds
or_all <- res_cdeis_all_logit$or_table %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::mutate(
    row_label = paste0(set, "\n", outcome, " (n1=", n1, "/", n0 + n1, ")"),
    row_label = factor(row_label, levels = rev(row_label)),
    p_label = paste0("p=", formatC(p.value, format = "f", digits = 3))
  )

p_or_forest <- ggplot(or_all, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "grey35") +
  geom_point(size = 3, color = "steelblue") +
  geom_text(aes(label = paste0("OR=", formatC(estimate, format = "f", digits = 2), 
                               " [", formatC(conf.low, format = "f", digits = 2), 
                               ", ", formatC(conf.high, format = "f", digits = 2), "]",
                               "\n", p_label)),
            nudge_y = 0.28, hjust = 0.5, size = 3.2) +
  scale_x_log10() +
  labs(x = "Odds Ratios (OR)", y = NULL) +
  ZZWtool::ZZWTheme() + 
  theme(axis.text.y = element_text(size = 9, angle = 0, hjust = 1))

print(p_or_forest)
ggsave(p_or_forest,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_cdeis_OR_forest_20260625.pdf",
       width = 8, height = 8)

# all CDEIS scatter plot
p_scatter_a <- res_cdeis_all$d %>% 
  ggplot(aes(CDEIS_continuous, microbial_indole_balance_score)) +
  geom_point(alpha = 0.8, size = 3, color = "#4575b4") +
  geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = res_cdeis_all$r_lab, hjust = 1.05, vjust = 1.3, size = 3.6) +
  geom_vline(xintercept = 3, linetype = 2, size = 1, color = "#1a9850") +
  geom_vline(xintercept = 9, linetype = 2, size = 1, color = "#fc8d59") +
  labs(x = "CDEIS",
       y = "MIBS") +
  ZZWtool::ZZWTheme()

ggsave(p_scatter_a,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_cdeis_all_samples_20260630.pdf",
       width = 6, height = 6)

or_scatter_a <- res_cdeis_all_logit$or_table %>% 
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::mutate(
    row_label = paste0(set, "\n", outcome, " (n1=", n1, "/", n0 + n1, ")"),
    row_label = factor(row_label, levels = rev(row_label)),
    p_label = paste0("p=", formatC(p.value, format = "f", digits = 3))
  )

p_or_forest_a <- ggplot(or_scatter_a, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "grey35") +
  geom_point(size = 3, color = "black") +
  geom_text(aes(label = paste0("OR=", formatC(estimate, format = "f", digits = 2),
                               " [", formatC(conf.low, format = "f", digits = 2),
                               ", ", formatC(conf.high, format = "f", digits = 2), "]",
                               "\n", p_label)),
            nudge_y = 0.28, hjust = 0.5, size = 3.2) +
  scale_x_log10() +
  labs(x = "Odds Ratios (OR)", y = NULL) +
  ZZWtool::ZZWTheme() +
  theme(axis.text.y = element_text(size = 9, angle = 90, hjust = 1))

ggsave(p_or_forest_a,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_cdeis_all_samples_forest_OR_20260630.pdf",
       width = 6, height = 6)

p_cdeis_combined <- cowplot::plot_grid(
  p_scatter_a, p_or_forest_a,
  ncol = 2, rel_widths = c(2, 1),
  align = "h", axis = "tb",
  # labels = c("A", "B"), 
  label_size = 12
)
print(p_cdeis_combined)
ggsave(p_cdeis_combined,
       filename = "~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_cdeis_combined_20260630.pdf",
       width = 9, height = 6)

# dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
# temp_source_data <- res_cdeis_all$scatter$data %>% 
#   dplyr::select(sample_id, visit_encounter_id, CDEIS_continuous, microbial_indole_balance_score)
# 
# readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9e_MIBS_vs_CDEIS_260721.csv")
# 
# 
# temp_source_data <- p_or_forest_a$data
# readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9e_MIBS_vs_CDEIS_logit_260721.csv")

rm(list = ls()); gc()




# 3.5: overall MIBS landscape  -------------------------------------------------
#   Datasets: ALL wPCDAI samples  (dz, from Part 2; severity = severity_class2)
#   Two visualizations per dataset:
#     (1) rank plot - samples ordered by MIBS (ascending) as the x index; y = MIBS; color = severity
#         (remission / mild / moderate_severe)
#     (2) PCA on ALL tryptophan metabolites (z-scored), colored by MIBS (continuous) and by severity

rm(list = ls()); gc()

load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_score_evaluation_input_20260721.RData")      # dz, score_vars, covs
load("~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity/mibs_zscore_params_all_metabolites_20260721.RData") # zscore_params_all_metabolites
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData') # object_final

severity_colors3 <- c(remission = "#169e77", mild = "#7570b3", moderate_severe = "#d76327")

sev3_levels <- c("remission", "mild", "moderate_severe")

# all tryptophan-metabolite z-score columns (frozen, wPCDAI-trained)
trp_met_vids <- zscore_params_all_metabolites$variable_id
zcol_all <- paste0("z_", trp_met_vids)

# dataset 1: wPCDAI samples (dz)
ds_wpcdai <- dz %>%
  dplyr::mutate(
    severity3 = dplyr::recode(as.character(severity_class2),
                              "mod_severe" = "moderate_severe"),
    severity3 = factor(severity3, levels = sev3_levels)
  ) %>%
  dplyr::filter(!is.na(microbial_indole_balance_score))

save(ds_wpcdai, 
     file = "~/Project/00_IBD_project/Data/20260721_MIBS_CD_severity//mibs_overall_landscape_20260630.RData")

# reusable visualizations 
# (1) rank plot: order samples by MIBS ascending -> x = order index, y = MIBS, color = severity
plot_mibs_rank <- function(d, cohort_label, file_tag, n_quantiles = 4) {
  d_rank <- d %>%
    dplyr::filter(!is.na(microbial_indole_balance_score), !is.na(severity3)) %>%
    dplyr::arrange(microbial_indole_balance_score) %>%
    dplyr::mutate(order_index = dplyr::row_number())
  
  p_main <- ggplot(d_rank, aes(order_index, microbial_indole_balance_score, color = severity3)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(values = severity_colors3, drop = FALSE, name = "severity") +
    labs(x = "Sample order index (MIBS ascending)",
         y = "MIBS") +
    ZZWtool::ZZWTheme() +
    theme(legend.position = c(0.82, 0.18))
  
  d_q <- d_rank %>%
    dplyr::mutate(
      quantile_grp = dplyr::ntile(order_index, n_quantiles),
      quantile_grp = factor(quantile_grp, labels = paste0("Q", seq_len(n_quantiles)))
    ) %>%
    dplyr::count(quantile_grp, severity3, .drop = FALSE) %>%
    dplyr::group_by(quantile_grp) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
  
  p_inset <- ggplot(d_q, aes(quantile_grp, prop, fill = severity3)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = severity_colors3, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "MIBS quantile", y = "Proportion") +
    ZZWtool::ZZWTheme() +
    theme(
      legend.position = "none",
      axis.text  = element_text(size = 6),
      axis.title = element_text(size = 7),
      plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3)
    )
  
  p <- cowplot::ggdraw(p_main) +
    cowplot::draw_plot(p_inset, x = 0.1, y = 0.55, width = 0.36, height = 0.38)
  
  print(p)
  p
}

#   ggsave(p, filename = paste0("~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_rank_",
#                               file_tag, "_20260630.pdf"), width = 6, height = 6)

# (2) PCA on ALL tryptophan metabolites (frozen z); two copies colored by MIBS and by severity
plot_trp_pca <- function(d, cohort_label, file_tag) {
  zc <- intersect(zcol_all, colnames(d))
  mat <- d %>% dplyr::select(dplyr::all_of(zc)) %>% as.data.frame()
  
  # keep samples with complete z across the metabolite panel; drop zero-variance columns
  complete_rows <- stats::complete.cases(mat)
  mat_c <- mat[complete_rows, , drop = FALSE]
  keep_col <- vapply(mat_c, function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1))
  mat_c <- mat_c[, keep_col, drop = FALSE]
  cat("\n", cohort_label, " PCA: ", nrow(mat_c), " samples x ", ncol(mat_c), " metabolites\n", sep = "")
  
  pr <- prcomp(mat_c, center = TRUE, scale. = FALSE)
  var_explained <- (pr$sdev^2) / sum(pr$sdev^2)
  
  pca_df <- d[complete_rows, ] %>%
    dplyr::mutate(PC1 = pr$x[, 1], PC2 = pr$x[, 2])
  
  axis_lab <- function(k) sprintf("PC%d (%.1f%%)", k, 100 * var_explained[k])
  
  # PERMANOVA: does severity3 explain variance in the trp-metabolite space (Euclidean dist, same geometry prcomp() decomposes)
  sev_idx <- !is.na(pca_df$severity3)
  set.seed(20260623)
  permanova_res <- vegan::adonis2(
    mat_c[sev_idx, , drop = FALSE] ~ severity3,
    data        = pca_df[sev_idx, , drop = FALSE],
    method      = "euclidean",
    permutations = 999
  )
  permanova_r2 <- permanova_res$R2[1]
  permanova_p  <- permanova_res$`Pr(>F)`[1]
  cat("\n", cohort_label, " PERMANOVA (severity3, n = ", sum(sev_idx), "): R2 = ",
      round(permanova_r2, 3), ", p = ", round(permanova_p, 4), "\n", sep = "")
  permanova_label <- sprintf("PERMANOVA: R² = %.3f, %s",
                             permanova_r2,
                             ifelse(permanova_p < 0.001, "p < 0.001",
                                    sprintf("p = %.3f", permanova_p)))
  
  # colored by MIBS (continuous)
  p_mibs <- ggplot(pca_df, aes(PC1, PC2, color = microbial_indole_balance_score)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_viridis_c(option = "D", name = "MIBS") +
    labs(x = axis_lab(1), y = axis_lab(2)) +
    ZZWtool::ZZWTheme() +
    theme(legend.position = c(0.85, 0.85))
  
  # colored by severity (3 levels); annotate with the PERMANOVA result (severity3 ~ metabolite space)
  p_sev <- ggplot(pca_df %>% dplyr::filter(!is.na(severity3)),
                  aes(PC1, PC2, color = severity3)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(values = severity_colors3, drop = FALSE, name = "severity") +
    labs(x = axis_lab(1), y = axis_lab(2), subtitle = permanova_label) +
    ZZWtool::ZZWTheme() +
    theme(legend.position = c(0.85, 0.85),
          plot.subtitle = element_text(size = 9, hjust = 0))
  
  p_combined <- cowplot::plot_grid(p_sev, p_mibs, ncol = 2, align = "h", axis = "tb")
  print(p_combined)
  #   ggsave(p_combined, filename = paste0("~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_trp_pca_",
  #                                        file_tag, "_20260630.pdf"), width = 10, height = 5)
  
  list(pca = pr, var_explained = var_explained, pca_df = pca_df,
       permanova = permanova_res, permanova_r2 = permanova_r2, permanova_p = permanova_p,
       plot_mibs = p_mibs, plot_sev = p_sev, plot_combined = p_combined)
}

# run both visualizations on both datasets
rank_wpcdai <- plot_mibs_rank(ds_wpcdai, "wPCDAI cohort", "wpcdai") 
pca_wpcdai <- plot_trp_pca(ds_wpcdai, "wPCDAI cohort", "wpcdai")

# ggsave(pca_wpcdai$plot_combined, filename = paste0("~/Project/00_IBD_project/Data/20260623_MIBS_CD_severity/mibs_trp_pca_wpcdai_20260630.pdf"), width = 10, height = 5)

dir.create("~/Project/00_IBD_project/Data/20260718_source_data/", showWarnings = FALSE, recursive = TRUE)
temp_source_data <- pca_wpcdai$plot_sev$data %>% 
  dplyr::select(sample_id, PC1, PC2, severity3)
readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9a_PCA_severity_260721.csv")

temp_source_data <- pca_wpcdai$plot_mibs$data %>% 
  dplyr::select(sample_id, PC1, PC2, microbial_indole_balance_score)
readr::write_csv(temp_source_data, file = "~/Project/00_IBD_project/Data/20260718_source_data/ext_fig9c_PCA_MIBS_260721.csv")