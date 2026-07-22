################################################################################
# Lipidomics workflow ----------------------------------------------------------

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms1_data')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms2_data_combined')
raw_table <- ms1_data$info %>% 
  bind_cols(ms1_data$subject) %>% 
  select(name:rt, contains("poolQC"))

readr::write_csv(raw_table, 
                 file = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms1_c18_pos.csv')

file.copy(from = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/ms2.msp', 
          to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms2_c18_pos.msp', 
          recursive = TRUE, 
          overwrite = TRUE)


load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms1_data')
raw_table <- ms1_data$info %>% 
  bind_cols(ms1_data$subject) %>% 
  select(name:rt, contains("poolQC"))

readr::write_csv(raw_table, file = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms1_c18_neg.csv')

file.copy(from = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/ms2.msp', 
          to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms2_c18_neg.msp', 
          recursive = TRUE, 
          overwrite = TRUE)


load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms1_data')
raw_table <- ms1_data$info %>% 
  bind_cols(ms1_data$subject) %>% 
  select(name:rt, contains("poolQC"))

readr::write_csv(raw_table, 
                 file = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms1_hilic_pos.csv')

file.copy(from = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/ms2.msp', 
          to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms2_hilic_pos.msp', 
          recursive = TRUE, 
          overwrite = TRUE)


load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms1_data')
raw_table <- ms1_data$info %>% 
  bind_cols(ms1_data$subject) %>% 
  select(name:rt, contains("poolQC"))

readr::write_csv(raw_table, file = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms1_hilic_neg.csv')

file.copy(from = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/ms2.msp', 
          to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/ms2_hilic_neg.msp', 
          recursive = TRUE, 
          overwrite = TRUE)


# class_adduct -----------------------------------------------------------------
select_class_adduct_pos <- readxl::read_xlsx('~/Project/04_package/00_Database/MSDIAL_lipidBlast/msdial_lipidblast_lib_class_adduct_241121_pos_selected.xlsx', sheet = 2)
select_class_adduct_pos <- select_class_adduct_pos %>% 
  mutate(class_adduct = paste0(compound_class, '_', adduct)) %>% 
  pull(class_adduct)


select_class_adduct_neg <- readxl::read_xlsx('~/Project/04_package/00_Database/MSDIAL_lipidBlast/msdial_lipidblast_lib_class_adduct_241121_neg_selected.xlsx', sheet = 2)
select_class_adduct_neg <- select_class_adduct_neg %>% 
  mutate(class_adduct = paste0(compound_class, '_', adduct)) %>% 
  pull(class_adduct)

list_class_adduct <- list('pos' = select_class_adduct_pos, 
                          'neg' = select_class_adduct_neg)

save(list_class_adduct, 
     file = '~/Project/00_IBD_project/Data/20241120_CD_complicates/list_class_adduct_241121.RData')

library(DoddLabMetID)
library(DoddLabDatabase)
DoddLabMetID::annotate_metabolite(ms1_file = 'ms1_c18_pos.csv', 
                                  ms2_file = 'ms2_c18_pos.msp',
                                  ms2_type = 'msp',
                                  polarity = 'positive',
                                  path = '~/Project/00_IBD_project/Data/20241120_CD_complicates',
                                  column = 'c18',
                                  ce = '20',
                                  lib =  'msdial_lipid',
                                  class_adduct_list = select_class_adduct_pos, 
                                  mz_tol = 10, 
                                  mz_ppm_thr = 150, 
                                  is_rt_score = FALSE, 
                                  is_include_precursor = FALSE,
                                  matched_frag_cutoff = 2, 
                                  scoring_approach = 'dp')


file.rename(from = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation/', 
            to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_c18_pos/')




library(DoddLabMetID)
library(DoddLabDatabase)
DoddLabMetID::annotate_metabolite(ms1_file = 'ms1_c18_neg.csv', 
                                  ms2_file = 'ms2_c18_neg.msp',
                                  ms2_type = 'msp',
                                  polarity = 'negative',
                                  path = '~/Project/00_IBD_project/Data/20241120_CD_complicates',
                                  column = 'c18',
                                  ce = '20',
                                  lib =  'msdial_lipid',
                                  class_adduct_list = select_class_adduct_neg, 
                                  mz_tol = 10, 
                                  mz_ppm_thr = 150, 
                                  is_rt_score = FALSE, 
                                  is_include_precursor = FALSE,
                                  matched_frag_cutoff = 2, 
                                  scoring_approach = 'dp')


file.rename(from = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation/', 
            to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_c18_neg/')


library(DoddLabMetID)
library(DoddLabDatabase)
DoddLabMetID::annotate_metabolite(ms1_file = 'ms1_hilic_pos.csv', 
                                  ms2_file = 'ms2_hilic_pos.msp',
                                  ms2_type = 'msp',
                                  polarity = 'positive',
                                  path = '~/Project/00_IBD_project/Data/20241120_CD_complicates',
                                  column = 'hilic',
                                  ce = '20',
                                  lib =  'msdial_lipid',
                                  class_adduct_list = select_class_adduct_pos, 
                                  mz_tol = 10, 
                                  mz_ppm_thr = 150, 
                                  is_rt_score = FALSE, 
                                  is_include_precursor = FALSE,
                                  matched_frag_cutoff = 2, 
                                  scoring_approach = 'dp')


file.rename(from = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation/', 
            to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_hilic_pos/')




library(DoddLabMetID)
library(DoddLabDatabase)
DoddLabMetID::annotate_metabolite(ms1_file = 'ms1_hilic_neg.csv', 
                                  ms2_file = 'ms2_hilic_neg.msp',
                                  ms2_type = 'msp',
                                  polarity = 'negative',
                                  path = '~/Project/00_IBD_project/Data/20241120_CD_complicates',
                                  column = 'hilic',
                                  ce = '20',
                                  lib =  'msdial_lipid',
                                  class_adduct_list = select_class_adduct_neg, 
                                  mz_tol = 10, 
                                  mz_ppm_thr = 150, 
                                  is_rt_score = FALSE, 
                                  is_include_precursor = FALSE,
                                  matched_frag_cutoff = 2, 
                                  scoring_approach = 'dp')


file.rename(from = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation/', 
            to = '~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_hilic_neg/')


# merge and dereplication for lipids -------------------------------------------
lipid_c18_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_c18_pos/annotation_summary.xlsx')
lipid_c18_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_c18_neg/annotation_summary.xlsx')
lipid_hilic_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_hilic_pos/annotation_summary.xlsx')
lipid_hilic_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/01_metabolite_annotation_hilic_neg/annotation_summary.xlsx')

extract_lipid_abbr <- function(x) {
  # browser()
  temp_label <- x %>% 
    str_split(' ') %>% 
    do.call(rbind, .)
  
  # extract lipid class
  temp_class <- temp_label[,1]
  
  # sum chain length
  temp_chain <- temp_label[,2] %>% 
    stringr::str_extract_all('\\d+\\:\\d+')
  
  temp_chain_table <- lapply(temp_chain, function(x){
    temp <- length(x)
    if (temp == 0){
      x <- c("0:0", "0:0", "0:0")
    } else if (temp == 1){
      x <- c(x, "0:0", "0:0")
    } else if (temp == 2){
      x <- c(x, "0:0")
    } else {
      x <- x
    }
    return(x)
  }) %>% 
    do.call(rbind, .)
  
  temp_chain_table <- temp_chain_table %>% 
    as.data.frame() %>% 
    separate(V1, into = c("chain_length_1", "db_1"), sep = ":") %>% 
    separate(V2, into = c("chain_length_2", "db_2"), sep = ":") %>%
    separate(V3, into = c("chain_length_3", "db_3"), sep = ":") %>%
    mutate_all(as.numeric) %>% 
    as_tibble() %>% 
    mutate(chain_length = chain_length_1 + chain_length_2 + chain_length_3) %>%
    mutate(db = db_1 + db_2 + db_3) %>% 
    select(chain_length, db)
  
  # extract O-, ;2O, -SN1, -SN2
  temp_subclass <- temp_label[,2] %>% 
    str_extract(pattern = 'O-|;2O|-SN1|-SN2')
  
  
  lipid_table <- temp_chain_table %>% 
    mutate(lipid_class = temp_class,
           lipid_subclass = temp_subclass) %>% 
    select(lipid_class, chain_length, db, lipid_subclass) %>% 
    mutate(abbr = case_when(
      lipid_subclass == 'O-' ~ paste0(lipid_class, ' O-', chain_length, ':', db),
      lipid_subclass == ';2O' ~ paste0(lipid_class, ' ', chain_length, ':', db, ';2O'),
      lipid_subclass == '-SN1' ~ paste0(lipid_class, ' ', chain_length, ':', db),
      lipid_subclass == '-SN2' ~ paste0(lipid_class, ' ', chain_length, ':', db),
      chain_length == 0 & db == 0 ~ paste0(lipid_class),
      is.na(lipid_subclass) ~ paste0(lipid_class, ' ', chain_length, ':', db)
    )) %>% 
    select(abbr, lipid_class, chain_length, db, lipid_subclass)
  
  return(lipid_table)
}

# extract the abbreviation of lipid
lipid_c18_pos %>% pull(feature_name) %>% unique() %>% length()
lipid_c18_neg %>% pull(feature_name) %>% unique() %>% length()
lipid_hilic_pos %>% pull(feature_name) %>% unique() %>% length()
lipid_hilic_neg %>% pull(feature_name) %>% unique() %>% length()

abbr_c18_pos <- extract_lipid_abbr(lipid_c18_pos$name)
abbr_c18_neg <- extract_lipid_abbr(lipid_c18_neg$name)
abbr_hilic_pos <- extract_lipid_abbr(lipid_hilic_pos$name)
abbr_hilic_neg <- extract_lipid_abbr(lipid_hilic_neg$name)


# merge lipid annotation
lipid_c18_pos_table <- lipid_c18_pos %>% 
  bind_cols(abbr_c18_pos)
lipid_c18_neg_table <- lipid_c18_neg %>% 
  bind_cols(abbr_c18_neg)
lipid_hilic_pos_table <- lipid_hilic_pos %>% 
  bind_cols(abbr_hilic_pos)
lipid_hilic_neg_table <- lipid_hilic_neg %>% 
  bind_cols(abbr_hilic_neg)

dir.create('~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation', showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(lipid_c18_pos_table, 
                    path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_c18_pos_table.xlsx', 
                    format_headers = FALSE)
writexl::write_xlsx(lipid_c18_neg_table, 
                    path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_c18_neg_table.xlsx', 
                    format_headers = FALSE)
writexl::write_xlsx(lipid_hilic_pos_table, 
                    path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_hilic_pos_table.xlsx', 
                    format_headers = FALSE)
writexl::write_xlsx(lipid_hilic_neg_table, 
                    path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_hilic_neg_table.xlsx', 
                    format_headers = FALSE)




################################################################################
# lipid dereplication ----------------------------------------------------------
# extract_lipid_abbr ----------------------------------------------------------
extract_lipid_abbr <- function(x) {
  # browser()
  temp_label <- x %>% 
    str_split(' ') %>% 
    do.call(rbind, .)
  
  # extract lipid class
  temp_class <- temp_label[,1]
  
  # sum chain length
  temp_chain <- temp_label[,2] %>% 
    stringr::str_extract_all('\\d+\\:\\d+')
  
  temp_chain_table <- lapply(temp_chain, function(x){
    temp <- length(x)
    if (temp == 0){
      x <- c("0:0", "0:0", "0:0")
    } else if (temp == 1){
      x <- c(x, "0:0", "0:0")
    } else if (temp == 2){
      x <- c(x, "0:0")
    } else {
      x <- x
    }
    return(x)
  }) %>% 
    do.call(rbind, .)
  
  temp_chain_table <- temp_chain_table %>% 
    as.data.frame() %>% 
    separate(V1, into = c("chain_length_1", "db_1"), sep = ":") %>% 
    separate(V2, into = c("chain_length_2", "db_2"), sep = ":") %>%
    separate(V3, into = c("chain_length_3", "db_3"), sep = ":") %>%
    mutate_all(as.numeric) %>% 
    as_tibble() %>% 
    mutate(chain_length = chain_length_1 + chain_length_2 + chain_length_3) %>%
    mutate(db = db_1 + db_2 + db_3) %>% 
    select(chain_length, db)
  
  # extract O-, ;2O, -SN1, -SN2
  temp_subclass <- temp_label[,2] %>% 
    str_extract(pattern = 'O-|;2O|-SN1|-SN2')
  
  
  lipid_table <- temp_chain_table %>% 
    mutate(lipid_class = temp_class,
           lipid_subclass = temp_subclass) %>% 
    select(lipid_class, chain_length, db, lipid_subclass) %>% 
    mutate(abbr = case_when(
      lipid_subclass == 'O-' ~ paste0(lipid_class, ' O-', chain_length, ':', db),
      lipid_subclass == ';2O' ~ paste0(lipid_class, ' ', chain_length, ':', db, ';2O'),
      lipid_subclass == '-SN1' ~ paste0(lipid_class, ' ', chain_length, ':', db),
      lipid_subclass == '-SN2' ~ paste0(lipid_class, ' ', chain_length, ':', db),
      chain_length == 0 & db == 0 ~ paste0(lipid_class),
      is.na(lipid_subclass) ~ paste0(lipid_class, ' ', chain_length, ':', db)
    )) %>% 
    select(abbr, lipid_class, chain_length, db, lipid_subclass)
  
  return(lipid_table)
}

# retrieve the experimental lib match result ---------------------------------
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')

lipid_exp_result <- object_final %>% 
  extract_annotation_table() %>%
  filter(metabolon_class == 'Lipids') %>% 
  as_tibble()

writexl::write_xlsx(lipid_exp_result, 
                    path = '~/Project/00_IBD_project/Data/20250112_CD_complication_analysis/lipid_exp_result_250113.xlsx', 
                    format_headers = FALSE)


temp_lipid <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_exp_result_manual_241125.xlsx', sheet = 2)
temp_lipid <- temp_lipid %>% select(variable_id, abbr)

lipid_exp_result <- lipid_exp_result %>% 
  left_join(temp_lipid, by = 'variable_id')

temp_fa <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241104_metabolite_class_distribution_CD_nonIBD/annot_fatty_acids_manual_2_241104.xlsx')
temp_fa <- temp_fa %>% 
  select(variable_id, C_number:COOH)

lipid_exp_result <- lipid_exp_result %>%
  left_join(temp_fa, by = 'variable_id')

writexl::write_xlsx(lipid_exp_result, 
                    path = '~/Project/00_IBD_project/Data/20250112_CD_complication_analysis/lipid_exp_result_manual_check_250114.xlsx', 
                    format_headers = FALSE)



# read the lipid annotation table

lipid_exp_result2 <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20250112_CD_complication_analysis/lipid_exp_result_manual_check_250114_checked_1_240114.xlsx', sheet = 2)

lipid_exp_result2 <- lipid_exp_result2 %>% 
  select(variable_id, mz, rt, id, abbr, Compound.name, 
         formula, smiles, inchikey, adduct, mz_error, 
         rt_error, msms_score_forward, msms_score_reverse, msms_matched_frag, 
         confidence_level, class_adduct) %>% 
  rename(abbr_name = abbr,
         compound_name = Compound.name)

temp_abbr <- lipid_exp_result2$abbr_name %>% 
  extract_lipid_abbr()


lipid_exp_result2 <- lipid_exp_result2 %>%
  mutate(lipid_class = stringr::str_split(lipid_exp_result2$class_adduct, pattern = '_', simplify = TRUE)[,1],
         chain_length = temp_abbr$chain_length,
         db = temp_abbr$db, 
         source = 'exp_lib') %>% 
  select(variable_id:confidence_level, lipid_class:db, source, class_adduct)


lipid_c18_pos <- readxl::read_xlsx(path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_c18_pos_table_manual_241122.xlsx', sheet = 2) %>% 
  mutate(condition = 'c18_pos')
lipid_c18_neg <- readxl::read_xlsx(path = '~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_c18_neg_table_manual_241124.xlsx', sheet = 2) %>% 
  mutate(condition = 'c18_neg')
lipid_hilic_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_hilic_pos_table_manual_241124.xlsx', sheet = 2) %>% 
  mutate(condition = 'hilic_pos')
lipid_hilic_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20241120_CD_complicates/lipid_annotation/lipid_hilic_neg_table_manual_241124.xlsx', sheet = 2) %>% 
  mutate(condition = 'hilic_neg')


lipid_annot_table_msdial <- lipid_c18_pos %>% 
  bind_rows(lipid_c18_neg) %>% 
  bind_rows(lipid_hilic_pos) %>% 
  bind_rows(lipid_hilic_neg) %>% 
  mutate(comp_score = msms_matched_frag*msms_score_forward) %>% 
  arrange(name, desc(comp_score)) %>% 
  select(feature_name, mz:id, abbr, name:msms_matched_frag, comp_score, everything()) %>% 
  select(feature_name, mz, rt, id:adduct, mz_error, rt_error, msms_score_forward:condition) %>% 
  mutate(lipid_class = case_when(lipid_subclass == 'O-' ~ paste0('Ether', lipid_class),
                                 lipid_class == 'Tocopherol' ~ 'Vitamin E',
                                 lipid_class == 'SPB' ~ 'Sph',
                                 lipid_class == '25-hydroxycholecalciferol' ~ 'Vitamin D',
                                 TRUE ~ lipid_class)) %>% 
  select(-c('lipid_subclass', 'lipid_subclass', 'selected', 'condition')) %>% 
  rename(variable_id = feature_name,
         abbr_name = abbr,
         compound_name = name) %>% 
  mutate(confidence_level = 'Level2.3') %>% 
  select(variable_id, everything()) %>% 
  select(variable_id:msms_matched_frag, confidence_level, lipid_class:db) %>% 
  mutate(source = 'theo_lib') %>% 
  mutate(class_adduct = paste0(lipid_class, '_', adduct))





# merge lipid annotation:
# 1. define adduct & polarity
# 2. ordered according to the confidence level and the score (dp*n_frag)
# 3. remove duplicated lipid annotation, same feature


lipid_annot_table <- lipid_exp_result2 %>% 
  bind_rows(lipid_annot_table_msdial) %>% 
  filter(class_adduct %in% c(
    'Car_[M+H]+',
    'Cer_[M+H]+',
    'DG_[M+NH4]+',
    'EtherLPC_[M+H]+',
    'EtherPC_[M+H]+',
    'EtherPE_[M-H]-',
    'FAHFA_[M-H]-',
    'FA_[M-H]-',
    'HexCer_[M+CH3COO]-',
    'LNAPE_[M-H]-',
    'LPC_[M+H]+',
    'LPE_[M-H]-',
    'LPG_[M-H]-',
    'LPS_[M-H]-',
    'MG_[M+H]+',
    'OxFA_[M-H]-',
    'OxPE_[M-H]-',
    'PC_[M+H]+',
    'PE_[M-H]-',
    'PG_[M-H]-',
    'PI_[M-H]-',
    'PS_[M+H]+',
    'PhytoSph_[M+H]+',
    'SM_[M+H]+',
    'Sph_[M+H]+',
    'SphP_[M+H]+',
    'TG_[M+NH4]+',
    'Vitamin D_[M+H]+',
    'Vitamin E_[M+H]+'
  )) %>% 
  arrange(confidence_level, desc(msms_score_forward*msms_matched_frag), abbr_name) %>% 
  distinct(variable_id, .keep_all = TRUE) %>%
  distinct(abbr_name, .keep_all = TRUE) %>% 
  arrange(abbr_name)


writexl::write_xlsx(lipid_annot_table, 
                    path = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_dereplication_table_update_250401.xlsx', 
                    format_headers = FALSE)

save(lipid_annot_table, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_dereplication_table_250401.RData')


lipid_annot_table %>% count(lipid_class) %>% as.data.frame()
lipid_annot_table %>% count(confidence_level) %>% as.data.frame()
lipid_annot_table %>% count(class_adduct) %>% as.data.frame()

################################################################################
# Generate object_lipid --------------------------------------------------------

load('~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/lipid_dereplication_table_250401.RData')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_input_data_cleaning/04_object_c18_pos_outlier_removal.RData')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_input_data_cleaning/04_object_c18_neg_outlier_removal.RData')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/01_input_data_cleaning/04_object_hilic_pos_outlier_removal.RData')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/01_input_data_cleaning/04_object_hilic_neg_outlier_removal.RData')



temp_c18_pos <- object_c18_pos %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% lipid_annot_table$variable_id)

temp_c18_neg <- object_c18_neg %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% lipid_annot_table$variable_id)

temp_hilic_pos <- object_hilic_pos %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% lipid_annot_table$variable_id)

temp_hilic_neg <- object_hilic_neg %>% 
  activate_mass_dataset(what = 'variable_info') %>% 
  filter(variable_id %in% lipid_annot_table$variable_id)



object_lipid <- merge_mass_dataset(x = temp_c18_pos,
                                   y = temp_c18_neg,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))

object_lipid <- merge_mass_dataset(x = object_lipid,
                                   y = temp_hilic_pos,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))

object_lipid <- merge_mass_dataset(x = object_lipid,
                                   y = temp_hilic_neg,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))


colnames(object_lipid@sample_info)
object_lipid@sample_info <- object_lipid@sample_info %>% dplyr::select(sample_id:note)
object_lipid@variable_info <- object_lipid@variable_info %>% 
  dplyr::select(variable_id:rt)

object_lipid@sample_info_note <- object_lipid@sample_info_note %>% slice(1:21)
object_lipid@variable_info_note <- object_lipid@variable_info_note %>% slice(1:3)


lipid_annot_table <- lipid_annot_table %>% 
  arrange(match(variable_id, object_lipid@variable_info$variable_id))

object_lipid@annotation_table <- lipid_annot_table

save(object_lipid, 
     file = '~/Project/00_IBD_project/Data/20250401_CD_complication_analysis/object_lipid_250401.RData')


rm(list = ls());gc()
