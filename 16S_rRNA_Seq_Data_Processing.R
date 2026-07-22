setwd('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/')

################################################################################
# Functions --------------------------------------------------------------------

# Function: aggregate_taxa
# aggregates OTU counts to a specified taxonomic level

aggregate_taxa <- function(otu_table, taxa_table, level, relative = FALSE, remove_na = TRUE) {
  # Check if the specified level exists in the taxa data frame
  if (!level %in% colnames(taxa)) {
    stop("Specified taxonomic level does not exist in the taxa data frame.")
  }
  
  # Merge OTU table with taxa table to get taxonomic information
  merged_data <- otu_table %>%
    inner_join(taxa_table %>% select(sequence, all_of(level)), by = "sequence")
  
  # Optionally remove rows with NA or empty values in the specified taxonomic level
  if (remove_na) {
    merged_data <- merged_data %>%
      filter(!is.na(!!sym(level)) & !!sym(level) != "")
  }
  
  # Aggregate OTU counts by the specified taxonomic level
  aggregated_data <- merged_data %>%
    group_by(!!sym(level)) %>%
    summarise(across(-sequence, sum, na.rm = TRUE)) %>%
    ungroup()
  
  # If relative abundance is requested, convert counts to relative abundances
  if (relative) {
    aggregated_data <- aggregated_data %>%
      mutate(across(-!!sym(level), ~ .x / sum(.x)))
  }
  
  return(aggregated_data)
}




################################################################################
# Retrieve meta-data from 16S rRNA data - RISK cohort --------------------------

library(tidyverse)
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_250327.RData')

# read the 16S rRNA data metadata for RISK cohort
ibd_cohort_16s <- readxl::read_xlsx('~/Project/00_IBD_project/Data/00_meta_data/00_raw/risk_clinical_info/2023-09-18_Full Cohort_RISK/Full Cohort RISK Summary_lower_241002.xlsx')
ibd_cohort_16s <- ibd_cohort_16s %>% 
  select(deidentified_master_patient_id, visit_encounter_id:visit_month, age_at_encounter:gender, first_behavior:disease_journey, final_diagnosis, `16s`, antibiotics, wpcdai) %>% 
  filter(!is.na(`16s`)) %>% 
  separate_rows(`16s`, sep = ';') %>%
  mutate(`16s` = stringr::str_replace(`16s`, '\\.BAM', '')) %>%
  mutate(`16s` = stringr::str_replace(`16s`, ' ', '')) %>% 
  rename(data_16s = '16s') %>% 
  # count(final_diagnosis)
  mutate(final_diagnosis = case_when(
    final_diagnosis == "Crohn's Disease" ~ 'CD',
    final_diagnosis == "Ulcerative Colitis" ~ 'UC',
    final_diagnosis == 'Not IBD' ~ 'non_IBD',
    final_diagnosis == 'IBD Unclassified' ~ 'CD'
  )) %>% 
  mutate(phenotype_type2 = case_when(
    disease_journey == 'B1 -> B1' ~ 'B1',
    disease_journey == 'B1 -> B2' ~ 'B2',
    disease_journey == 'B1 -> B2 -> B2+B3' ~ 'B2+B3',
    disease_journey == 'B1 -> B2+B3' ~ 'B2+B3',
    disease_journey == 'B1 -> B3 -> B2+B3' ~ 'B2+B3',
    disease_journey == 'B1 -> B3' ~ 'B3',
    disease_journey == 'B3 -> B3' ~ 'B3'
  )) %>% 
  rename('phenotype_type' = 'final_diagnosis')

ibd_cohort_16s %>% pull('data_16s') %>% length()


# read the omics patient mapping file 
omics_data <- readr::read_csv('~/Project/00_IBD_project/Data/00_meta_data/00_raw/risk_clinical_info/2023-09-18_Full Cohort_RISK/data/RISK/22113_Omics_Patient_mapping_20230906_175053.txt')

omics_data <- omics_data %>% 
  dplyr::filter(REQUESTED.FILE.PATH == '/risk-16s')

colnames(omics_data) <- tolower(colnames(omics_data))

omics_data %>% 
  count(characteristics_bio_material)

omics_data <- omics_data %>% 
  select(deidentified_master_patient_id, raw.data.file.name, characteristics_bio_material) %>% 
  mutate(data_16s = stringr::str_replace(raw.data.file.name, '\\.BAM', ''))

ibd_cohort_16s_mapped <- ibd_cohort_16s %>% 
  inner_join(omics_data, by = c('data_16s'))

ibd_cohort_16s_mapped %>% 
  count(phenotype_type)

ibd_cohort_16s_mapped %>% 
  count(phenotype_type, characteristics_bio_material)

ibd_cohort_16s_mapped %>% 
  count(antibiotics)

# save(ibd_cohort_16s_mapped, file = '~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/ibd_cohort_16s_mapped_251022.RData')
save(ibd_cohort_16s_mapped, file = '~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/ibd_cohort_16s_mapped_251027.RData')

################################################################################
# 16S rRNA raw data processing -------------------------------------------------
library(dada2); packageVersion("dada2")

path <- "~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/" # CHANGE ME to the directory containing the fastq files after unzipping.

setwd("~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/")

list.files(path)


# RISK 16S rRNA data
# The overlapping paired-end reads were stitched together
fnFs <- sort(list.files(path, pattern="\\.fastq", full.names = TRUE, recursive = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "\\.fastq"), `[`, 1)


# Inspect read quality profiles
plotQualityProfile(fnFs[1:10])

# Filter and trim
filteredFs <- file.path(path, "filtered", paste0(sample.names, "_filtered.fastq.gz"))
names(filteredFs) <- sample.names


out <- filterAndTrim(fnFs, filteredFs, truncLen=c(120),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
dir.create('./00_intermidiate_data/', showWarnings = FALSE, recursive = TRUE)
save(out, file = './00_intermidiate_data/out_251022.RData')

# Learn the Error Rates
err <- learnErrors(filteredFs, multithread=TRUE)
save(err, file = './00_intermidiate_data/err_251022.RData')

# Visualize the estimated error rates
plotErrors(err, nominalQ=TRUE)

# Sample Inference
dada_result <- dada(filteredFs, err=err, multithread=TRUE)
save(dada_result, file = './00_intermidiate_data/dada_result_251022.RData')

# # Merge paired reads
# # Inspect the merger data.frame from the first sample

# Construct sequence table
seq_tab <- makeSequenceTable(dada_result)
dim(seq_tab)
save(seq_tab, file = './00_intermidiate_data/seqtab_251023.RData')


# Inspect distribution of sequence lengths
table(nchar(getSequences(seq_tab)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seq_tab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seq_tab)

class(seqtab_nochim)
save(seqtab_nochim, file = './00_intermidiate_data/seqtab_nochim_251023.RData')
write.csv(seqtab_nochim, file = './00_intermidiate_data/risk_16s_nochim_251023.csv')

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_result, getN), rowSums(seqtab_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write_csv(as.data.frame(track), file = './00_intermidiate_data/sequence_track_251023.csv')


# Assign taxonomy 
taxa <- assignTaxonomy(seqtab_nochim, "~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/silva_species_assignment_v138.1.fa.gz")


taxa_print <- taxa # Removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

save(taxa, file = './00_intermidiate_data/taxa_251023.RData')
write.csv(taxa_print, file = './00_intermidiate_data/taxa.nochim_251023.csv')

# # count the taxonomy levels
# taxa_print %>% as.data.frame() %>% count(Kingdom)
# taxa_print %>% as.data.frame() %>% count(Phylum)
# taxa_print %>% as.data.frame() %>% count(Class)
# taxa_print %>% as.data.frame() %>% count(Family)
# taxa_print %>% as.data.frame() %>% count(Genus)
# taxa_print %>% as.data.frame() %>% count(Species)
# 
# taxa %>% as.data.frame() %>% count(Order)


################################################################################
# PhyloSeq ---------------------------------------------------------------------

dir.create('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis', showWarnings = FALSE)
setwd('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis')

library(tidyverse)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

library(DECIPHER)
library(phangorn)
library(plotly)


load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/ibd_cohort_16s_mapped_251027.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/taxa_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/seqtab_nochim_251023.RData')


sample_table <- ibd_cohort_16s_mapped %>% 
  select(-deidentified_master_patient_id.y) %>% 
  dplyr::rename('deidentified_master_patient_id' = 'deidentified_master_patient_id.x') %>% 
  select(data_16s, everything())

temp_sample_name <- rownames(seqtab_nochim)

sample_table <- sample_table %>% 
  arrange(match(data_16s, temp_sample_name)) %>% 
  column_to_rownames(var = 'data_16s') %>%
  as.data.frame()

# create phyloseq-class object
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), 
               sample_data(sample_table), 
               tax_table(taxa))

save(ps, file = './00_intermidiate_data/phyloseq_object_251023.RData')


load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/phyloseq_object_251023.RData')

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- names(dna)
ps

# build the phylogenetic tree --------------------------------------------------

library(phyloseq)
library(DECIPHER)
library(phangorn)

# alighn the sequences
seqs <- Biostrings::DNAStringSet(dna)
alignment <- AlignSeqs(seqs, anchor = NA)

# Build the phynogenetic tree (Neighbor-Joining Tree)
phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm <- phangorn::dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm)

# Rooting the tree, required by UniFrac
treeNJ_rooted <- phangorn::midpoint(treeNJ)

# Adding the tree to phyloseq object
ps_prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ps_prop_with_tree <- merge_phyloseq(ps_prop, treeNJ_rooted)

save(ps_prop, file = './00_intermidiate_data/ps_prop_251023.RData')
save(ps_prop_with_tree, file = './00_intermidiate_data/ps_prop_with_tree_251023.RData')

ps_prop_with_tree@sam_data <- sample_data(sample_table)

load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/ps_prop_with_tree_251023.RData')

# PCoA ordination --------------------------------------------------------------
# calculate distance matrix and ordination 
ord_PCoA_unifrac <- ordinate(ps_prop_with_tree, 
                             method="PCoA", 
                             distance = "unifrac", 
                             weighted=FALSE)

save(ord_PCoA_unifrac, file = './00_intermidiate_data/ord_PCoA_unifrac_251023.RData')

load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/ord_PCoA_unifrac_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/alpha_index_251023.RData')


plot_ordination(ps_prop_with_tree, ord_PCoA_unifrac, color="phenotype_type", title="PCoA_Unifac")
plot_ordination(ps_prop_with_tree, ord_PCoA_unifrac, color="characteristics_bio_material", title="PCoA_Unifac")

pcoa_data_3d <- plot_ordination(ps_prop_with_tree, ord_PCoA_unifrac, axes = 1:3, justDF = TRUE)

# add alpha_index
alpha_index <- estimate_richness(ps, measures=c("Shannon", "Simpson","Chao1"))
# save(alpha_index, file = './00_intermidiate_data/alpha_index_251023.RData')

pcoa_data_3d <- pcoa_data_3d %>% 
  rownames_to_column('data_16s') %>% 
  left_join(alpha_index %>% rownames_to_column('data_16s'), by = 'data_16s') %>% 
  column_to_rownames('data_16s')




# characteristics_bio_material colors ------------------------------------------
library(scatterplot3d)

temp_colors_material <- pcoa_data_3d %>% 
  mutate(color = case_when(characteristics_bio_material == 'Ileum' ~ '#fdb562',
                           characteristics_bio_material == 'Rectum' ~ '#8cd4c8',
                           characteristics_bio_material == 'Stool' ~ '#ffed6f')) %>% 
  pull(color)

pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_material.pdf', width = 6, height = 6)
scatterplot3d(pcoa_data_3d$Axis.1, pcoa_data_3d$Axis.2, pcoa_data_3d$Axis.3, 
              color = temp_colors_material, 
              angle = 30,          
              pch = 19,
              box = FALSE,         
              grid = TRUE,         
              lty.hplot = 2,       
              scale.y = 0.7,       
              xlab = "PCoA1 (11.3%)", ylab = "PCoA2 (3.5%)", zlab = 'PCoA3 (2.8%)')

legend("topright", 
       legend = c('Ileum', 'Rectum', 'Stool'),
       col = c('#fdb562', '#8cd4c8', '#ffed6f'), 
       pch = 19, 
       bty = "n")
dev.off()


# disease status colors --------------------------------------------------------
temp_colors_phenotype <- pcoa_data_3d %>% 
  mutate(color = case_when(phenotype_type == 'CD' ~ 'tomato',
                           phenotype_type == 'non_IBD' ~ 'dodgerblue')) %>% 
  pull(color)


pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_disease_status.pdf', width = 6, height = 6)
scatterplot3d(pcoa_data_3d$Axis.1, pcoa_data_3d$Axis.2, pcoa_data_3d$Axis.3, 
              color = temp_colors_phenotype, 
              angle = 30,         
              pch = 19,
              box = FALSE,         
              grid = TRUE,         
              lty.hplot = 2,       
              scale.y = 0.7,       
              xlab = "PCoA1 (11.3%)", ylab = "PCoA2 (3.5%)", zlab = 'PCoA3 (2.8%)')

legend("topright", 
       legend = c('non-IBD', 'CD'),
       col = c('dodgerblue', 'tomato'), 
       pch = 19, 
       bty = "n")
dev.off()


# rainbow colorscale for alpha index -------------------------------------------

colors_vec <- c("blue", "cyan", "green", "yellow", "orange", "red", "red", "red", "darkred")


my_palette_func <- colorRampPalette(colors_vec)
n_colors <- 100

gradient_colors <- my_palette_func(n_colors)
color_var <- alpha_index$Chao1
color_indices <- as.numeric(cut(color_var, breaks = n_colors))
val_min <- min(color_var)
val_max <- max(color_var)
point_colors <- gradient_colors[color_indices]

pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_alpha_index.pdf', width = 6, height = 6)
scatterplot3d(pcoa_data_3d$Axis.1, pcoa_data_3d$Axis.2, pcoa_data_3d$Axis.3, 
              color = point_colors, 
              angle = 30,          
              pch = 19,
              box = FALSE,        
              grid = TRUE,         
              lty.hplot = 2,       
              scale.y = 0.7,       
              xlab = "PCoA1 (11.3%)", ylab = "PCoA2 (3.5%)", zlab = 'PCoA3 (2.8%)')

dev.off()


pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_alpha_index_legend.pdf', width = 6, height = 6)
# generate the legend for alpha index color scale
plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = 'Legend')
legend_image <- as.raster(matrix(rev(gradient_colors), ncol = 1))

rasterImage(legend_image, 0, 0, 0.4, 1)

axis(4, at = seq(0, 1, length.out = 5), 
     labels = round(seq(val_min, val_max, length.out = 5), 1), 
     las = 1) 
dev.off()

rm(list = ls());gc()


################################################################################
# stacked bar plot -------------------------------------------------------------

load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/ibd_cohort_16s_mapped_251027.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/taxa_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/seqtab_nochim_251023.RData')

taxa_table <- taxa %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'sequence') %>% 
  as_tibble()


# otu table 
otu_table <- seqtab_nochim %>% 
  as.data.frame()

rownames(otu_table) <- rownames(seqtab_nochim)

otu_table <- otu_table %>% 
  rotate_df() %>% 
  rownames_to_column(var = 'sequence') %>% 
  as_tibble()

# merged data
temp_data <- otu_table %>% 
  pivot_longer(-sequence, names_to = 'sample', values_to = 'count') %>% 
  left_join(taxa_table, by = c('sequence'))


temp_data <- temp_data %>% 
  left_join(ibd_cohort_16s_mapped, by = c('sample' = 'data_16s'))

temp_data_ileum <- temp_data %>% 
  dplyr::filter(characteristics_bio_material == 'Ileum')


temp_phylum <- temp_data_ileum %>% 
  filter(!is.na(Phylum)) %>% 
  count(Phylum) %>% 
  arrange(desc(n)) %>% 
  slice(1:6) %>% 
  pull(Phylum) %>% 
  unique()

temp_data2 <- temp_data_ileum %>%
  filter(!is.na(Phylum)) %>% 
  mutate(temp_phylum = case_when(
    Phylum %in% temp_phylum ~ Phylum,
    TRUE ~ 'Other'
  )) %>% 
  select(sample, phenotype_type, Phylum, temp_phylum, count) %>% 
  group_by(sample, temp_phylum) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()

# calculate the percentage for every sample and phylum
temp_data2 <- temp_data2 %>% 
  group_by(sample) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = count / total * 100) %>% 
  select(-total)


# plot stacked bar plot
temp_data3 <- temp_data2 %>% 
  left_join(ibd_cohort_16s_mapped, by = c('sample' = 'data_16s')) %>%
  select(sample, phenotype_type, temp_phylum, count, percentage) %>% 
  rename(phylum = temp_phylum)


# CD
temp_data3_cd <- temp_data3 %>% 
  dplyr::filter(phenotype_type == 'CD') %>% 
  arrange(match(phylum, temp_phylum), desc(percentage))

sample_index <- temp_data3_cd %>% pull(sample) %>% unique()

temp_data3_cd <- temp_data3_cd %>% 
  mutate(index = match(sample, sample_index))

temp_data3_cd$phylum <- factor(temp_data3_cd$phylum, levels = c(temp_phylum, 'Other'))

temp_plot_cd <- ggplot(temp_data3_cd) +
  geom_bar(aes(x = index, y = percentage, fill = phylum), stat = 'identity') +
  scale_fill_manual(values = c('Firmicutes' = '#8cd4c8',
                               'Proteobacteria' = '#fdb562',
                               'Bacteroidota' = '#ffed6f',
                               'Actinobacteriota' = '#b5de68',
                               'Fusobacteriota' = '#beb9da',
                               'Cyanobacteria' = '#fbcee5',
                               'Other' = '#fc8070')) + 
  ylab('Percentage (%)') +
  xlab('Sample index') +
  ZZWtool::ZZWTheme() +
  theme(margin = margin(0, 0, 0, 0))


# Non-IBD
temp_data3_non_IBD <- temp_data3 %>% 
  dplyr::filter(phenotype_type == 'non_IBD') %>% 
  arrange(match(phylum, temp_phylum), desc(percentage))

sample_index <- temp_data3_non_IBD %>% pull(sample) %>% unique()

temp_data3_non_IBD <- temp_data3_non_IBD %>% 
  mutate(index = match(sample, sample_index))

temp_data3_non_IBD$phylum <- factor(temp_data3_non_IBD$phylum, levels = c(temp_phylum, 'Other'))

temp_plot_non_IBD <- ggplot(temp_data3_non_IBD) +
  geom_bar(aes(x = index, y = percentage, fill = phylum), stat = 'identity') +
  scale_fill_manual(values = c('Firmicutes' = '#8cd4c8',
                               'Proteobacteria' = '#fdb562',
                               'Bacteroidota' = '#ffed6f',
                               'Actinobacteriota' = '#b5de68',
                               'Fusobacteriota' = '#beb9da',
                               'Cyanobacteria' = '#fbcee5',
                               'Other' = '#fc8070')) + 
  ylab('Percentage (%)') +
  xlab('Sample index') +
  ZZWtool::ZZWTheme() +
  theme(margin = margin(0, 0, 0, 0))


combined_data <- dplyr::bind_rows(temp_data3_cd, temp_data3_non_IBD)

combined_plot_facet <- ggplot(combined_data) +
  geom_bar(aes(x = index, y = percentage, fill = phylum), stat = 'identity') +
  facet_grid(~ phenotype_type, scales = 'free_x', space = 'free') +
  scale_fill_manual(values = c('Firmicutes' = '#8cd4c8',
                               'Proteobacteria' = '#fdb562',
                               'Bacteroidota' = '#ffed6f',
                               'Actinobacteriota' = '#b5de68',
                               'Fusobacteriota' = '#beb9da',
                               'Cyanobacteria' = '#fbcee5',
                               'Other' = '#fc8070')) + 
  ylab('Percentage (%)') +
  xlab('Sample index') +
  ZZWtool::ZZWTheme() +
  theme(margin = margin(0, 0, 0, 0),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = "bold")) 


ggsave(combined_plot_facet, 
       filename = '~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/Combined_stacked_phylum_ileum_260122.pdf', 
       width = 10, height = 6)

rm(list = ls());gc()





################################################################################
# differential taxonomies ------------------------------------------------------

load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/confident_otu_table_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/confident_taxa_table_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/confident_sample_table_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20251030_16S_metabolomics_correlation/00_intermidiate_data/object_ileum_251103.RData')
load('~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/result_16s_ileum_260121.RData')


# 16S rRNA data differential analysis with Maaslin2 ----------------------------
meta_table_ileum <- confident_sample_table_ileum %>% 
  select(data_16s:visit_encounter_id, phenotype_type, age, gender, race, use_antibiotics, wpcdai) %>%
  rename('ID' = 'data_16s',
         'subject_id' = 'deidentified_master_patient_id') %>% 
  mutate(age = as.numeric(age))

dir.create('~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv', showWarnings = FALSE, recursive = TRUE)
write_tsv(meta_table_ileum, 
          file = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/meta_table_ileum_260121.tsv')


otu_taxa_table_ileum <- confident_otu_table_ileum %>%
  mutate(asv = match(sequence, confident_taxa_table_ileum$sequence) %>% confident_taxa_table_ileum$asv[.]) %>%
  select(asv, everything(), -sequence) %>% 
  column_to_rownames('asv') %>%
  rotate_df() %>% 
  rownames_to_column('ID')

write_tsv(otu_taxa_table_ileum, 
          file = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/otu_taxa_table_ileum_260121.tsv')

library(Maaslin2)
fit_data <- Maaslin2(
  input_data = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/otu_taxa_table_ileum_260121.tsv', 
  input_metadata = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/meta_table_ileum_260121.tsv', 
  output = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/maaslin2_results',
  fixed_effects = c('phenotype_type', 'gender', 'age', 'use_antibiotics', 'race'),
  transform = 'LOG',
  # random_effects = c('subject_id'),
  reference = "phenotype_type,non_IBD;gender,Male;race,Caucasian;use_antibiotics,No",
  standardize = TRUE,
  save_models = TRUE)

# read results
result_16s_ileum <- readr::read_tsv('~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/maaslin2_results/all_results.tsv')
result_16s_ileum <- result_16s_ileum %>% 
  filter(metadata == 'phenotype_type' & value == 'CD') %>% 
  filter(qval < 0.2) %>% 
  left_join(confident_taxa_table_ileum, by = c('feature' = 'asv'))

save(result_16s_ileum, 
     file = '~/Project/00_IBD_project/Data/20260120_microbial_met_16S_rRNA_correlation/01_differential_analysis_16S_ileum_asv/result_16s_ileum_260121.RData')

rm(list = ls());gc()


# visualize the differential analysis results ----------------------------------
meta_table_ileum <- confident_sample_table_ileum %>% 
  select(data_16s:visit_encounter_id, phenotype_type, age, gender, race, use_antibiotics, wpcdai) %>%
  rename('ID' = 'data_16s',
         'subject_id' = 'deidentified_master_patient_id') %>% 
  mutate(age = as.numeric(age))

temp_data <- confident_otu_table_ileum %>% 
  dplyr::filter(sequence %in% result_16s_ileum$sequence) %>% 
  left_join(confident_taxa_table_ileum, by = 'sequence') %>%
  select(sequence, asv, Kingdom:Species, everything()) %>% 
  select(-c('sequence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% 
  column_to_rownames('asv') %>% 
  rotate_df() %>%  
  rownames_to_column('ID') %>%
  left_join(meta_table_ileum, by = 'ID')


summarize_availability_table <- lapply(result_16s_ileum$feature, function(x){
  temp_data_aim <- temp_data %>% 
    select(ID, subject_id:use_antibiotics, x)
  
  # count un-zero number, total number and percentage in different phenotype_type 
  temp_result <- temp_data_aim %>% 
    group_by(phenotype_type) %>%
    summarise(n = sum(get(x) > 0),
              n_total = n()) %>% 
    mutate(percentage = n / n_total * 100,
           feature = x) %>% 
    select(feature, everything())
  
  return(temp_result)
})

temp_data_availability_summary <- summarize_availability_table %>% 
  bind_rows() %>% 
  filter(!is.na(phenotype_type)) %>% 
  tidyr::pivot_wider(names_from = phenotype_type, values_from = c(n, n_total, percentage))

result_16s_ileum2 <- result_16s_ileum %>% 
  left_join(temp_data_availability_summary, by = 'feature') %>% 
  mutate(label = ifelse(coef > 0, 'increase', 'decrease')) %>% 
  arrange(coef) %>%
  mutate(index = seq(n())) %>% 
  mutate(y = 0.1)



# plot

library(ggplot2)
library(cowplot)

y_limits <- c(0.5, 25.5) 
y_expand <- c(0, 0)      

temp_plot1 <- ggplot(result_16s_ileum2, aes(x = index, y = coef)) +
  geom_segment(aes(x = index, xend = index, y = 0, yend = coef, color = label), size = 1) +
  geom_point(aes(color = label), size = 3) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c('increase' = "#fb7f72", 'decrease' = "#7fb1d3")) +
  scale_x_continuous(breaks = seq(1, 25, 1), 
                     labels = paste0(result_16s_ileum2$Genus, ' (', result_16s_ileum2$feature, ')'),
                     limits = y_limits, 
                     expand = y_expand) +
  xlab('Genus (ASV)') +
  ylab('Coefficient') +
  coord_flip(xlim = c(0.5, 25.5)) + 
  ZZWtool::ZZWTheme() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
        legend.position = c(0.8, 0.2))


temp_plot2 <- ggplot(result_16s_ileum2, aes(x = 1, y = index, fill = Phylum)) + 
  geom_tile(color = "black", lwd = 0.5, linetype = 1, height = 1) + 
  scale_fill_manual(values = c('Firmicutes' = '#8cd4c8',
                               'Proteobacteria' = '#fdb562',
                               'Bacteroidota' = '#ffed6f',
                               'Actinobacteriota' = '#b5de68',
                               'Fusobacteriota' = '#beb9da',
                               'Cyanobacteria' = '#fbcee5',
                               'Other' = '#fc8070')) +
  scale_y_continuous(limits = y_limits, expand = y_expand) +
  scale_x_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(legend.position = 'none') 



temp_plot3 <- ggplot(result_16s_ileum2) +
  geom_bar(aes(x = index, y = percentage_CD),
           fill = "#fb7f72", stat = 'identity', position = 'dodge', width = 0.8) +
  scale_x_continuous(breaks = seq(1, 25, 1), 
                     limits = y_limits, 
                     expand = y_expand,
                     labels = NULL) + 
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     labels = paste0(seq(0, 100, by = 20), '%'),
                     expand = c(0.01, 0.01)) +
  xlab(NULL) +
  ylab('Prevalence in CD (%)') +
  coord_flip(ylim = c(0, 100)) + 
  ZZWtool::ZZWTheme() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')


temp_plot4 <- ggplot(result_16s_ileum2) +
  geom_bar(aes(x = index, y = percentage_non_IBD),
           fill = "#7fb1d3", stat = 'identity', position = 'dodge', width = 0.8) +
  scale_x_continuous(breaks = seq(1, 25, 1), 
                     limits = y_limits, 
                     expand = y_expand,
                     labels = NULL) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     labels = paste0(seq(0, 100, by = 20), '%'),
                     expand = c(0.01, 0.01)) +
  xlab(NULL) +
  ylab('Prevalence in Non-IBD (%)') +
  coord_flip(ylim = c(0, 100)) +
  ZZWtool::ZZWTheme() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')


pp <- list(temp_plot1, temp_plot2, temp_plot3, temp_plot4)
temp_plot_merge <- plot_grid(plotlist = pp, 
                             ncol = 4, 
                             align = 'h', 
                             axis = 'bt', 
                             rel_widths = c(0.8, 0.05, 0.3, 0.3)) # 色块通常可以窄一点

print(temp_plot_merge)

ggsave(temp_plot_merge, 
       filename = '~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/differential_taxonmy_ileum_260122.pdf', 
       width = 12, height = 5)







################################################################################
# microbial dysbiosis index ----------------------------------------------------
load('~/Project/00_IBD_project/Data/20251104_correlation_metabolite_with_microbial_dysbiosis_inex/microbial_dysbiosis_index_251104.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/ord_PCoA_unifrac_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/alpha_index_251023.RData')
load('~/Project/00_IBD_project/Data/20251022_16S_rRNA_analysis/00_intermidiate_data/ps_prop_with_tree_251023.RData')


pcoa_data_3d <- plot_ordination(ps_prop_with_tree, ord_PCoA_unifrac, axes = 1:3, justDF = TRUE)

# add alpha_index
alpha_index <- estimate_richness(ps, measures=c("Shannon", "Simpson","Chao1"))
# save(alpha_index, file = './00_intermidiate_data/alpha_index_251023.RData')

pcoa_data_3d <- pcoa_data_3d %>%
  rownames_to_column('data_16s') %>%
  left_join(alpha_index %>% rownames_to_column('data_16s'), by = 'data_16s') %>%
  left_join(microbial_dysbiosis_index, by = c('data_16s' = 'sample_id')) %>% 
  filter(!is.na(microbial_dysbiosis_index) & !is.infinite(microbial_dysbiosis_index))

library(scatterplot3d)

colors_vec <- c("blue", "cyan", "green", "yellow", "orange", "red", "red", "red", "darkred")
my_palette_func <- colorRampPalette(colors_vec)
n_colors <- 100

gradient_colors <- my_palette_func(n_colors)
color_var <- pcoa_data_3d$microbial_dysbiosis_index
color_indices <- as.numeric(cut(color_var, breaks = n_colors))
val_min <- min(color_var)
val_max <- max(color_var)
point_colors <- gradient_colors[color_indices]

pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_microbial_dysbiosis_index_260123.pdf', width = 6, height = 6)
scatterplot3d(pcoa_data_3d$Axis.1, pcoa_data_3d$Axis.2, pcoa_data_3d$Axis.3,
              color = point_colors,
              angle = 30,        
              pch = 19,
              box = FALSE,         
              grid = TRUE,        
              lty.hplot = 2,      
              scale.y = 0.7,      
              xlab = "PCoA1 (11.3%)", ylab = "PCoA2 (3.5%)", zlab = 'PCoA3 (2.8%)')

dev.off()


pdf('~/Project/00_IBD_project/Data/20260122_16S_rRNA_data_analysis/PCoA_3D_scatterplot3d_dysbiosis_index_legend_260123.pdf', width = 6, height = 6)
# generate the legend for alpha index color scale ----------------------------
plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = 'Legend')
legend_image <- as.raster(matrix(rev(gradient_colors), ncol = 1))

rasterImage(legend_image, 0, 0, 0.4, 1)

axis(4, at = seq(0, 1, length.out = 5), 
     labels = round(seq(val_min, val_max, length.out = 5), 1), 
     las = 1) 
dev.off()

rm(list = ls());gc()