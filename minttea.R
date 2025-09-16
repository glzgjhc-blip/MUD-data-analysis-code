######function######
prep_add_col_prefix <- function(df, prefix) {
  sample_id_col <- df %>% select(sample_id__)
  df <- df %>% select(-sample_id__)
  colnames(df) <- paste0(prefix, "__", colnames(df))
  df <- bind_cols(sample_id_col, df)
  return(df)
}

# Performs TSS normalization
prep_convert_to_relative_abundances <- function(df) {
  sample_id_col <- df %>%
    select(sample_id__) %>%
    mutate(sample_id__ = as.character(sample_id__))
  
  df <- df %>%
    select(-sample_id__) %>%
    mutate_all(as.numeric) %>%
    mutate(totalSum = rowSums(across(everything()))) %>%
    mutate(across(.cols = -c(totalSum), ~ .x / totalSum)) %>%
    select(-totalSum) %>%
    bind_cols(sample_id_col) %>%
    relocate(sample_id__)
  return(df)
}

# Performs a CLR transformation on a given data table
prep_clr_transpose <- function(df, pseudo_count) {
  require(mixOmics)
  includes_sample_id <- ifelse('sample_id__' %in% colnames(df), TRUE, FALSE)
  tmp_df <- df
  if(includes_sample_id) tmp_df <- df %>% select(-sample_id__)
  
  # CLR transform values
  clr_df <- logratio.transfo(tmp_df, logratio = 'CLR', offset = pseudo_count)
  
  # Reformat as table and add sample ids if needed
  class(clr_df) <- "matrix"
  if (!includes_sample_id) return(data.frame(clr_df))
  return(bind_cols(data.frame(clr_df), df %>% select(sample_id__)))
}

# prep_load_metadata <- function(path) {
#   valid_labels <- c('H', config$disease_labels)
# 
#   metadata_df <- read.csv(file = path, sep = '\t', header = TRUE) %>%
#       as_tibble() %>%
#       rename(sample_id__ = 1) %>%
#       mutate(sample_id__ = as.character(sample_id__)) %>%
#       select(sample_id__, DiseaseState) %>%
#       filter(DiseaseState %in% valid_labels) %>%
#       mutate(DiseaseState = replace(DiseaseState, DiseaseState == 'H', 'healthy')) %>%
#       mutate(DiseaseState = replace(DiseaseState, DiseaseState != 'healthy', 'disease')) %>%
#       mutate_at(vars(DiseaseState), factor)
#   return(metadata_df)
# }


prep_load_pathways <- function(path, 
                               remove_non_bacterial = TRUE, 
                               remove_super_pwys = FALSE, 
                               convert_to_relab = TRUE) {
  # Load pathways and filter non-bacterial pathways and redundant (super-pathways) pathways
  pathways_df <- read.table(file = path, sep = '\t',
                            header = TRUE,
                            quote = "",
                            check.names = FALSE) 
  
  if (remove_non_bacterial) pathways_df <- filter_bacterial_pathways(pathways_df)
  if (remove_super_pwys) pathways_df <- filter_super_pathways(pathways_df)
  
  # Transpose pathways DF and correct column names to pathways names
  pathways_df <- pathways_df %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    as_tibble()
  
  names(pathways_df) <- unlist(pathways_df[1, ])
  
  pathways_df <- pathways_df[-(1:2), ]
  pathways_df <- rename(pathways_df, sample_id__ = 1)
  pathways_df <- prep_add_col_prefix(pathways_df, 'P')
  
  # Fix names
  names(pathways_df) <- make.names(names(pathways_df))
  
  if (convert_to_relab) pathways_df <- prep_convert_to_relative_abundances(pathways_df)
  return(pathways_df)
}

# prep_load_metagenome <- function(path, convert_to_relab = TRUE) {
#   # Read table, remove "description" column if exists (outputted by picrust)
#   df <- read.csv(file = path, header = TRUE, check.names = FALSE, sep = '\t') %>%
#     select(-any_of(c('description'))) 
#   
#   if ("function" %in% names(df)) 
#     df <- df %>% rename(`KO` = `function`)
#   
#   df <- df %>%
#     tibble::column_to_rownames(var = "KO") %>%
#     t() %>%
#     data.frame() %>%
#     tibble::rownames_to_column(var = "sample_id__") %>%
#     prep_add_col_prefix('G')
#   
#   if (convert_to_relab) df <- prep_convert_to_relative_abundances(df)
#   return(df)
# }

prep_load_taxonomy <- function(path, convert_to_relab = TRUE) {
  df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    select(-any_of(c('(UG)'))) %>%
    rename(sample_id__ = 1) %>%
    prep_add_col_prefix('T')
  
  # Fix names
  names(df) <- make.names(names(df))
  
  if (convert_to_relab) df <- prep_convert_to_relative_abundances(df)
  return(df)
}

prep_load_metabolites <- function(path, remove_unknown_metabolites = TRUE) {
  metabolites_df <- read.csv(file = path, sep = '\t', header = TRUE, check.names = FALSE) %>%
    rename(sample_id__ = 1) %>%
    mutate(sample_id__ = as.character(sample_id__)) %>%  
    prep_add_col_prefix('M')
  
  # We only keep identified metabolites to facilitate interpretation
  d <- basename(dirname(path))
  if (remove_unknown_metabolites & d == "esrd_wang_2020") {
    keep <- grep("_unknown$", colnames(metabolites_df), invert = TRUE, value = TRUE)
    log_debug(sprintf('Removed %d/%d metabolite features with no annotation',
                      ncol(metabolites_df)-length(keep),
                      ncol(metabolites_df)-1))
    metabolites_df <- metabolites_df[,keep]
  } else if (remove_unknown_metabolites & d %in% c("uc_franzosa_2019","cd_franzosa_2019")) {
    keep <- grep("\\: NA$", colnames(metabolites_df), invert = TRUE, value = TRUE) 
    log_debug(sprintf('Removed %d/%d metabolite features with no annotation',
                      ncol(metabolites_df)-length(keep),
                      ncol(metabolites_df)-1))
    metabolites_df <- metabolites_df[,keep]
  }
  
  # Fix names
  names(metabolites_df) <- make.names(names(metabolites_df))
  
  return(metabolites_df)
}



 prep_join_metadata <- function(dataset_df, metadata_df) {
   all_data <- dataset_df %>%
     inner_join(metadata_df, by = 'sample_id__')
   return(all_data)
 }
 

prep_join_features <- function(df1, df2) {
  joined_df <- df1 %>% 
    inner_join(df2, by = 'sample_id__') %>%
    relocate(sample_id__)
  return(joined_df)
}


prep_sanitize_dataset <- function(df, feature_set, rare_feature_cutoff = 0.15, mean_abundance_cutoff = NULL) {
  
  # Remove constant features
  nfeatures <- ncol(df)
  df <- df %>% select(where(~ n_distinct(.) > 1))
  log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) constant features for dataset %s',
                    (nfeatures - ncol(df)),
                    nfeatures,
                    as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                    feature_set))
  
  # Remove rare features (have <15% non-zero values)
  nfeatures <- ncol(df)
  non_zero_percentage <- colSums(df != 0) / nrow(df)
  rare_features <- names(non_zero_percentage[non_zero_percentage <= rare_feature_cutoff])
  if (length(rare_features)>0) df <- df %>% select(-all_of(rare_features))
  log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) low-prevalance features for dataset %s',
                    (nfeatures - ncol(df)),
                    nfeatures,
                    as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                    feature_set))
  
  # Remove rare features by mean abundance (only for TSS-normalized feature types)
  if (!is.null(mean_abundance_cutoff)) {
    # Verify the data looks like TSS-normalized data, i.e. each row sums to ~1.
    if (sum(round(rowSums(df %>% select(-sample_id__))) == 1) / nrow(df) < 0.95)
      log_error('Seems like the feature table is not TSS-normalized. Abundance filtering should not be performed')
    nfeatures <- ncol(df)
    mean_abundances <- apply(df %>% select(-sample_id__), 2, mean)
    rare_features <- names(mean_abundances[mean_abundances <= mean_abundance_cutoff])
    if (length(rare_features)>0) df <- df %>% select(-all_of(rare_features))
    log_debug(sprintf('Sanitizer: removed %d/%d (%d%%) low-abundance features for dataset %s',
                      (nfeatures - ncol(df)),
                      nfeatures,
                      as.integer(100 * (nfeatures - ncol(df)) / nfeatures),
                      feature_set))
  }
  
  return(df)
}
config <- config::get(file = "src/ml_pipeline/config.yml")

######start######
library(MintTea)
library(dplyr)
library(readxl)
source('src/analyses/plotting_functions.R')
source('src/ml_pipeline/clustering.R')
source('src/ml_pipeline/preprocessing.R')
#load data
df_metabolic <- as.data.frame(read_excel("代谢物一二级鉴定定量总表.xlsx", sheet = 1))
df_func <- as.data.frame(read.table("out_path.function.txt",header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F))
df_species <-  as.data.frame(read.table("species.txt",header=T, row.names=1))
metadata = read_xlsx("metadata.xlsx")
#pre
metadata <- metadata %>% rename(sample_id__ = 1)
#pre
rownames(df_func) <- sub(":.*", "", rownames(df_func))
df_func_t <- as.data.frame(t(df_func))
df_func_t$sample_id__ <- rownames(df_func_t)
df_func_t <- df_func_t %>% dplyr::relocate(sample_id__)
pathway <- prep_add_col_prefix(df_func_t, 'P')  # 给所有代谢物加前缀
#pre
df_metabolic <- df_metabolic %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "")
metabolic<- as.data.frame(df_metabolic[, c( 18:104)])
rownames(metabolic) <- df_metabolic[, 2]
metabolic <- as.data.frame(t(metabolic))
metabolic$sample_id__ <- rownames(metabolic)
metabolic <- metabolic %>% relocate(sample_id__)
metabolic <- prep_add_col_prefix(metabolic, 'S')  # 给所有代谢物加前缀
#pre
df_species <- df_species %>% rename(sample_id__ = 2)
df_species <- df_species[,-1]
species <- prep_add_col_prefix(df_species, 'T')  # 给所有代谢物加前缀
#common
common_ids <- Reduce(intersect, list(
  species$sample_id__,
  metabolic$sample_id__,
  pathway$sample_id__
))
species   <- species %>% filter(sample_id__ %in% common_ids)
metabolic <- metabolic %>% filter(sample_id__ %in% common_ids)
pathway   <- pathway %>% filter(sample_id__ %in% common_ids)

#metadata pre
metadata_filtered   <- metadata %>% filter(sample_id__ %in% common_ids)
metadata_T0_HC <- metadata_filtered %>%
  filter(Group %in% c("T0","HC")) %>%
  mutate(DiseaseState = ifelse(Group == "T0", 'disease', 'healthy'))
metadata_T1_HC <- metadata_filtered %>%
  filter(Group %in% c("T1","HC")) %>%
  mutate(DiseaseState = ifelse(Group == "T1", 'disease', 'healthy'))

#species pre
## Convert to relative abundances
species <- prep_convert_to_relative_abundances(species)
## Split 
species_T0_HC <- species %>% filter(sample_id__ %in% metadata_T0_HC$sample_id__)
species_T1_HC <- species %>% filter(sample_id__ %in% metadata_T1_HC$sample_id__)
## Remove rare features + constant features
species_T0_HC <- prep_sanitize_dataset(species_T0_HC, 'T', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)
species_T1_HC <- prep_sanitize_dataset(species_T1_HC, 'T', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)


#pathway pre
## Split 
pwy_T0_HC <- pathway %>% filter(sample_id__ %in% metadata_T0_HC$sample_id__)
pwy_T1_HC <- pathway %>% filter(sample_id__ %in% metadata_T1_HC$sample_id__)
## Remove rare features + constant features
pwy_T0_HC <- prep_sanitize_dataset(pwy_T0_HC, 'P', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)
pwy_T1_HC <- prep_sanitize_dataset(pwy_T1_HC, 'P', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)

#metabolic pre
# 手动去掉首列名中的隐藏字符
colnames(metabolic)[1] <- sub("^[^A-Za-z0-9]+", "", colnames(metabolic)[1])
colnames(metabolic) <- make.names(colnames(metabolic))
#取对数
mtb <- metabolic %>% mutate_if(is.numeric, log)
## Split 
mtb_T0_HC <- mtb %>% filter(sample_id__ %in% metadata_T0_HC$sample_id__)
mtb_T1_HC <- mtb %>% filter(sample_id__ %in% metadata_T1_HC$sample_id__)
## Remove rare features + constant features
# 定义函数：移除 CV < 0.005 的代谢物
remove_low_variation_features <- function(df, sample_id_col = "sample_id__", threshold = 0.1) {
  # 拆出 sample_id__
  id_col <- df[[sample_id_col]]
  
  # 计算 CV：按列（每个代谢物）
  df_features <- df %>% select(-all_of(sample_id_col))
  means <- colMeans(df_features, na.rm = TRUE)
  sds   <- apply(df_features, 2, sd, na.rm = TRUE)
  cvs   <- sds / means
  
  # 保留 CV >= threshold 的特征
  keep_features <- names(cvs)[cvs >= threshold]
  df_filtered <- df %>% select(all_of(sample_id_col), all_of(keep_features))
  
  return(df_filtered)
}

# 使用示例
mtb_T0_HC <- remove_low_variation_features(mtb_T0_HC, threshold = 0.1)
mtb_T1_HC <- remove_low_variation_features(mtb_T1_HC, threshold = 0.1)


metadata_T0_HC=metadata_T0_HC[,c(1,23)]
metadata_T1_HC=metadata_T1_HC[,c(1,23)]


# Cluster highly correlated features-----------------------------------------------
# 
# Paths for output files
clusters_T0_HC_file <- 'data/ml_input/metacardis_T0_HC/clusters.tsv'
clusters_T1_HC_file <- 'data/ml_input/metacardis_T1_HC/clusters.tsv'
final_T0_HC_file <- 'data/ml_input/metacardis_T0_HC/all_data.tsv'
final_T1_HC_file <- 'data/ml_input/metacardis_T1_HC/all_data.tsv'
source('src/ml_pipeline/clustering.R')
source('src/ml_pipeline/preprocessing.R')
cluster_init_clusters_table(clusters_T0_HC_file)
cluster_init_clusters_table(clusters_T1_HC_file)
# -----------------------------------------------

all_data_T0_HC <- metadata_T0_HC %>%
  select(sample_id__, DiseaseState) %>%
  inner_join(species_T0_HC, by = 'sample_id__') %>%
  inner_join(pwy_T0_HC, by = 'sample_id__') %>%
  inner_join(mtb_T0_HC, by = 'sample_id__') 

all_data_T0_HC_clustered <- cluster_cluster_features(
  dataset = all_data_T0_HC %>% select(-sample_id__),
  feature_set_type = 'T+P+S',
  cluster_type = 'clustering95',
  clusters_output = clusters_T0_HC_file
)
all_data_T0_HC <- bind_cols(all_data_T0_HC %>% select(sample_id__), all_data_T0_HC_clustered)

# ---
all_data_T1_HC <- metadata_T1_HC %>%
  select(sample_id__, DiseaseState) %>%
  inner_join(species_T1_HC, by = 'sample_id__') %>%
  inner_join(pwy_T1_HC, by = 'sample_id__') %>%
  inner_join(mtb_T1_HC, by = 'sample_id__') 

all_data_T1_HC_clustered <- cluster_cluster_features(
  dataset = all_data_T1_HC %>% select(-sample_id__),
  feature_set_type = 'T+P+S',
  cluster_type = 'clustering95',
  clusters_output = clusters_T1_HC_file
)
all_data_T1_HC <- bind_cols(all_data_T1_HC %>% select(sample_id__), all_data_T1_HC_clustered)

# ---------------------------------------------------
# Save to files
# ---------------------------------------------------

# Save for intermediate integration (one large table)
write_tsv(all_data_T0_HC, final_T0_HC_file)
write_tsv(all_data_T1_HC, final_T1_HC_file)

# Print summaries (numbers of cases/controls)
table(all_data_T0_HC$DiseaseState)
table(all_data_T1_HC$DiseaseState)

# 假设你的整合数据是 all_data_T0_HC，包含了 T__*, P__*, S__* 特征
# 还有 sample_id__ 列 和 disease_state 列
library(MintTea)
colnames(all_data_T0_HC) <- make.names(colnames(all_data_T0_HC), unique = TRUE)
colnames(all_data_T1_HC) <- make.names(colnames(all_data_T1_HC), unique = TRUE)
T0_HC_results <- MintTea(
 all_data_T0_HC,
  view_prefixes = c("T", "P", "S"),
  sample_id_column = "sample_id__",
  study_group_column = "DiseaseState",
  param_diablo_keepX = 10,
  param_sgcca_design = 0.7,
  param_n_repeats = 10,
  param_n_folds = 5,
  param_edge_thresholds = 0.8,
  n_evaluation_repeats = 3,
  n_evaluation_folds = 5,
  seed = 3407
)
T1_HC_results <- MintTea(
  all_data_T1_HC,
  view_prefixes = c("T", "P", "S"),
  sample_id_column = "sample_id__",
  study_group_column = "DiseaseState",
  param_diablo_keepX = 10,
  param_sgcca_design = 0.7,
  param_n_repeats = 10,
  param_n_folds = 5,
  param_edge_thresholds = 0.8,
  n_evaluation_repeats = 3,
  n_evaluation_folds = 5,
  seed = 3407
)
feature_type_color_map = 
  c("M" = "midnightblue", 
    "T" = "firebrick4", 
    "P" = "darkorange3", 
    "S" = "chartreuse4", 
    "All" = "grey30")
plot_module_with_igraph(T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]], prefix_colors = feature_type_color_map)
plot_module_with_igraph(T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module2"]], prefix_colors = feature_type_color_map)
plot_module_with_igraph(T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]], prefix_colors = feature_type_color_map)
plot_module_with_igraph(T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module2"]], prefix_colors = feature_type_color_map)
######
plot_modules_overview <- function(minttea_output,
                                  feature_type_color_map) {
  
  # ---- Strip 1: N features per dataset per module ---->
  tmp1 <- sens_analysis_modules %>%
    mutate(feature_type = substr(feature,1,1)) %>%
    mutate(module = as.numeric(gsub('module', '', module))) %>%
    group_by(dataset, module, feature_type) %>%
    summarise(N = n(), .groups = "drop")
  
  # Only plot modules with at least 2 types of features
  modules_to_plot <- modules_overview %>%
    filter(multi_view) %>%
    select(dataset, module, is_interesting) %>%
    mutate(module = as.numeric(gsub('module','',module)))
  
  tmp1 <- inner_join(tmp1, modules_to_plot, by = c('dataset','module'))
  
  tmp1$module2 <- factor(paste('Module', tmp1$module),
                         levels = paste('Module', max(tmp1$module):1))
  
  p1 <- ggplot(tmp1 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)),
               aes(x = module2, y = N, fill = feature_type)) +
    geom_rect(data = tmp1 %>%
                filter(dataset %in% dataset_order) %>%
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>%
                distinct(),
              aes(alpha = dataset %>% as.numeric() %% 2 == 0),
              xmin = -Inf,
              xmax = Inf,
              ymin = -Inf,
              ymax = Inf,
              fill = 'gray80',
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    new_scale("alpha") +
    geom_bar(aes(alpha = is_interesting, color = is_interesting),
             width = 0.8,
             position="stack",
             stat="identity") +
    scale_alpha_manual(values = c('FALSE' = 0.3, 'TRUE' = 0.9), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray70', 'TRUE' = 'black'), guide = "none") +
    geom_hline(yintercept = 0) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    scale_x_discrete(expand = c(0, 0.5)) +
    coord_flip() +
    theme_classic() +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    xlab(NULL) +
    ylab("No. of features\nof each type") +
    scale_fill_manual(name = "View", values = feature_type_color_map) +
    theme(panel.grid.major.x =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.minor.x = element_blank()) +
    theme(panel.grid.major.y =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.text.y = element_text(size = 9)) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text(size = 11)) +
    theme(strip.background = element_blank()) +
    theme(strip.placement = "outside") +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1)) +
    theme(panel.spacing.y = unit(6, "points"))
  
  if (hide_y_axis_text) p1 <- p1 + theme(axis.text.y = element_blank())
  
  # ---- Strip 3: AUC ---->
  tmp2 <- modules_overview %>%
    select(-is_interesting) %>%
    mutate(module = as.numeric(gsub('module','',module))) %>%
    inner_join(modules_to_plot, by = c('dataset','module'))
  
  tmp2$module2 <- factor(paste('Module', tmp2$module),
                         levels = levels(tmp1$module2))
  points_size <- ifelse(hide_y_axis_text, 3, 3.7)
  
  p2 <- ggplot(tmp2 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)),
               aes(x = module2)) +
    geom_rect(data = tmp1 %>%
                filter(dataset %in% dataset_order) %>%
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>%
                distinct(),
              aes(alpha = dataset %>% as.numeric() %% 2 == 0),
              xmin = -Inf,
              xmax = Inf,
              ymin = -Inf,
              ymax = Inf,
              fill = 'gray80',
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    geom_linerange(aes(ymax = mean_module_auc_shuffled + sd_module_auc_shuffled,
                       ymin = mean_module_auc_shuffled - sd_module_auc_shuffled),
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    geom_point(aes(y = mean_module_auc_shuffled),
               shape = 16, size = points_size - 0.5,
               color = 'grey60', alpha = 0.8) +
    geom_point(aes(y = mean_module_auc, fill = is_interesting, color = is_interesting),
               shape = 23,
               size = points_size, alpha = 0.9) +
    scale_fill_manual(values = c('FALSE' = '#D7E7ED', 'TRUE' = 'skyblue4'), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray60', 'TRUE' = 'black'), guide = "none") +
    scale_y_continuous(breaks = seq(0.5,1,0.1), expand = expansion(mult = c(0.05,0.1))) +
    scale_x_discrete(expand = c(0, 0.5)) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Module AUC") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(panel.grid.major.x =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.major.y =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(axis.text.y = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  if (show_rf) p2 <- p2 + geom_hline(aes(yintercept = mean_auc_rf), color = "goldenrod2", linewidth = 2, alpha = 0.7)
  
  # ---- Strip 2: Cross-view correlations ---->
  p3 <- ggplot(tmp2 %>%
                 filter(dataset %in% dataset_order) %>%
                 mutate(dataset = factor(dataset, levels = dataset_order)),
               aes(x = module2)) +
    geom_rect(data = tmp1 %>%
                filter(dataset %in% dataset_order) %>%
                mutate(dataset = factor(dataset, levels = dataset_order)) %>%
                select(dataset) %>%
                distinct(),
              aes(alpha = dataset %>% as.numeric() %% 2 == 0),
              xmin = -Inf, xmax = Inf, ymin = -Inf,
              ymax = Inf, fill = 'gray80', inherit.aes = FALSE) +
    scale_alpha_manual(values = c('FALSE' = 0, 'TRUE' = 0.3), guide = "none") +
    geom_linerange(aes(ymax = avg_spear_corr_shuffled + sd_spear_corr_shuffled,
                       ymin = avg_spear_corr_shuffled - sd_spear_corr_shuffled),
                   alpha = 0.4, linewidth = 2, color = "grey70") +
    geom_point(aes(y = avg_spear_corr_shuffled),
               shape = 16, size = points_size - 0.5, color = 'grey60', alpha = 0.8) +
    geom_point(aes(y = avg_spear_corr, fill = is_interesting, color = is_interesting),
               shape = 23,
               size = points_size,
               alpha = 0.9) +
    scale_fill_manual(values = c('FALSE' = '#D9B9AB', 'TRUE' = 'sienna4'), guide = "none") +
    scale_color_manual(values = c('FALSE' = 'gray60', 'TRUE' = 'black'), guide = "none") +
    scale_x_discrete(expand = c(0, 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.05,0.1))) +
    coord_flip() +
    theme_classic() +
    xlab(NULL) +
    ylab("Cross-view avg.\ncorrelation") +
    facet_grid(rows = vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    theme(panel.grid.major.x =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(panel.grid.major.y =
            element_line(linewidth = 0.5, color = "grey93")) +
    theme(axis.title.x = element_text(size = 11)) +
    theme(axis.text.y = element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(panel.spacing.y = unit(6, "points"))
  
  # Combine plots
  tmp_rel_widths <- c(8.1, 3, 3)
  if (hide_y_axis_text) tmp_rel_widths <- c(7, 3, 3)
  print(plot_grid(p1, p3, p2,
                  nrow = 1,
                  rel_widths = tmp_rel_widths,
                  align = 'h', axis = 'tb'))
}
plot_modules_overview(T1_HC_results,feature_type_color_map)

module_features <- list(
  T0_module1 = T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]]$module1$features,
  T0_module2 = T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]]$module2$features,
  T1_module1 = T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]]$module1$features,
  T1_module2 = T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]]$module2$features
)
library(tidyverse)

# 创建所有模块组合（不含重复、非自身对比）
module_pairs <- combn(names(module_features), 2, simplify = FALSE)

# 统一全集
global_universe <- unique(unlist(module_features))

# 修正后的 overlap 计算函数
get_overlap_stats <- function(pair) {
  m1 <- pair[1]
  m2 <- pair[2]
  f1 <- module_features[[m1]]
  f2 <- module_features[[m2]]
  
  # 用统一全集做背景
  shared_universe <- global_universe
  
  a <- length(base::intersect(f1, f2))
  b <- length(base::setdiff(f2, f1))
  c <- length(base::setdiff(f1, f2))
  d <- length(base::setdiff(shared_universe, base::union(f1, f2)))
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  
  tibble(
    module_left = m1,
    module_right = m2,
    n_intersection = a,
    n_left_only = c,
    n_right_only = b,
    neither = d,
    fisher_p = ft$p.value,
    odds_ratio = unname(ft$estimate)
  )
}
# 计算所有模块对的重叠并输出
overlap_results <- map_dfr(module_pairs, get_overlap_stats) %>%
  mutate(fisher_fdr = p.adjust(fisher_p, method = "fdr")) %>%
  filter(n_intersection >= 2) %>%
  arrange(fisher_fdr)

# 查看或保存结果
print(overlap_results)
write_csv(overlap_results, "Table_S7_module_overlap.csv")


####基于mint tea尝试做症状关联网络######
metadata_T1 <- metadata_filtered %>%
  filter(Group %in% c("T1")) 



all_data_T0=all_data_T0_HC %>%
  filter(DiseaseState %in% c("disease")) 

####相关####

library(WGCNA)
library(tidyverse)
library(readxl)

# ---- 0. 读取元数据 ----
metadata = read_xlsx("metadata_v2.xlsx")
metadata <- metadata %>% rename(sample_id__ = 1)
metadata_filtered <- metadata %>% filter(sample_id__ %in% common_ids)
metadata_T1 <- metadata_filtered %>% filter(Group == "T1")
metadata_T1 <- metadata_T1[, c(1,6,13:16,23:25)]  # 提取指定列

# ---- 1. 提取模块特征表达矩阵 ----
module_feats <- T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]][["features"]]
expr_mat <- all_data_T1 %>%
  select(sample_id__, all_of(module_feats)) %>%
  column_to_rownames("sample_id__")
expr_s <- scale(expr_mat)

# ---- 2. 计算 eigengenes ----
get_eigengene <- function(expr_mat, pattern) {
  feats <- grep(pattern, colnames(expr_mat), value = TRUE)
  if (length(feats) < 2) return(rep(NA, nrow(expr_mat)))
  expr_sub <- scale(expr_mat[, feats])
  prcomp(expr_sub)$x[, 1]
}

eig_all <- prcomp(expr_s)$x[, 1]
eig_S   <- get_eigengene(expr_s, "^S__")
eig_T   <- get_eigengene(expr_s, "^T__")
eig_P   <- get_eigengene(expr_s, "^P__")

df_eig <- data.frame(sample_id__ = rownames(expr_s),
                     eig_all = eig_all,
                     eig_S = eig_S,
                     eig_T = eig_T,
                     eig_P = eig_P)

# ---- 3. 合并 metadata 和 module 特征表达值 ----
df_expr <- expr_mat %>% rownames_to_column("sample_id__")
df <- df_eig %>%
  left_join(df_expr, by = "sample_id__") %>%
  left_join(metadata_T1, by = "sample_id__")

# ---- 4. 选择要做相关分析的症状变量列 ----
cor_vars <- colnames(metadata_T1)[-1]  # 除 sample_id__

# ---- 5. 定义要做相关分析的X变量（eigengene 和 module features）----
x_vars <- c("eig_all", "eig_S", "eig_T", module_feats)

# ---- 6. Spearman相关分析并FDR校正 ----
library(dplyr)

cor_results <- expand.grid(X = x_vars, Y = cor_vars, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    res = list(
      tryCatch(
        cor.test(df[[X]], df[[Y]], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
    ),
    R = if (!is.null(res)) as.numeric(res$estimate) else NA_real_,
    P = if (!is.null(res)) res$p.value else NA_real_
  ) %>%
  select(-res) %>%
  ungroup() %>%
  group_by(Y) %>%  # 👉 按照每个量表做FDR校正
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  ungroup()
cor_results
# ---- 7. 保存结果 ----
write.csv(cor_results, "T1_eigengene_module_spearman_results.csv", row.names = FALSE)
cor_results_T1=cor_results


# ---- 0. 读取元数据 ----
metadata = read_xlsx("metadata_v2.xlsx")
metadata <- metadata %>% rename(sample_id__ = 1)
metadata_filtered <- metadata %>% filter(sample_id__ %in% common_ids)
metadata_T0 <- metadata_filtered %>% filter(Group == "T0")
metadata_T0 <- metadata_T0[, c(1,6,13:16,23:25)]  # 提取指定列

# ---- 1. 提取模块特征表达矩阵 ----
module_feats <- T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]][["features"]]
expr_mat <- all_data_T0 %>%
  select(sample_id__, all_of(module_feats)) %>%
  column_to_rownames("sample_id__")
expr_s <- scale(expr_mat)

# ---- 2. 计算 eigengenes ----
get_eigengene <- function(expr_mat, pattern) {
  feats <- grep(pattern, colnames(expr_mat), value = TRUE)
  if (length(feats) < 2) return(rep(NA, nrow(expr_mat)))
  expr_sub <- scale(expr_mat[, feats])
  prcomp(expr_sub)$x[, 1]
}
expr_s=expr_s[,-20]
eig_all <- prcomp(expr_s)$x[, 1]
eig_S   <- get_eigengene(expr_s, "^S__")
eig_T   <- get_eigengene(expr_s, "^T__")
eig_P   <- get_eigengene(expr_s, "^P__")

df_eig <- data.frame(sample_id__ = rownames(expr_s),
                     eig_all = eig_all,
                     eig_S = eig_S,
                     eig_T = eig_T,
                     eig_P = eig_P)

# ---- 3. 合并 metadata 和 module 特征表达值 ----
df_expr <- expr_mat %>% rownames_to_column("sample_id__")
df <- df_eig %>%
  left_join(df_expr, by = "sample_id__") %>%
  left_join(metadata_T0, by = "sample_id__")

# ---- 4. 选择要做相关分析的症状变量列 ----
cor_vars <- colnames(metadata_T0)[-1]  # 除 sample_id__

# ---- 5. 定义要做相关分析的X变量（eigengene 和 module features）----
x_vars <- c("eig_all", "eig_S", "eig_T","eig_P", module_feats)

# ---- 6. Spearman相关分析并FDR校正 ----
library(dplyr)

cor_results <- expand.grid(X = x_vars, Y = cor_vars, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    res = list(
      tryCatch(
        cor.test(df[[X]], df[[Y]], method = "pearson", exact = FALSE),
        error = function(e) NULL
      )
    ),
    R = if (!is.null(res)) as.numeric(res$estimate) else NA_real_,
    P = if (!is.null(res)) res$p.value else NA_real_
  ) %>%
  select(-res) %>%
  ungroup() %>%
  group_by(Y) %>%  # 👉 按照每个量表做FDR校正
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  ungroup()

# ---- 7. 保存结果 ----
write.csv(cor_results, "T0_eigengene_module_spearman_results.csv", row.names = FALSE)
cor_results_T0=cor_results

library(ggplot2)
library(patchwork)

library(ggplot2)
library(dplyr)
library(forcats)

# 合并两个数据框，增加 Timepoint 列
plot_df <- bind_rows(
  cor_results_T0_plot %>% mutate(Timepoint = "T0"),
  cor_results_T1_plot %>% mutate(Timepoint = "T1")
)

library(ggplot2)
library(dplyr)
library(forcats)
# 假设你已经有 cor_results_T0 和 cor_results_T1

# 格式化函数
prepare_plot_df <- function(cor_res, timepoint_label) {
  cor_res %>%
    mutate(
      label = ifelse(FDR < 0.01, paste0("**\n", round(R, 2)),
                     ifelse(FDR < 0.05, paste0("*\n", round(R, 2)),
                            round(R, 2))),
      Timepoint = timepoint_label
    ) %>%
    select(X, Y, R, label, Timepoint)
}

# 应用到 T0 和 T1
cor_results_T0_plot <- prepare_plot_df(cor_results_T0, "T0")
cor_results_T1_plot <- prepare_plot_df(cor_results_T1, "T1")
library(ggplot2)
library(dplyr)

# ---- 准备绘图数据 ----
prepare_plot_df <- function(cor_res, timepoint_label) {
  cor_res %>%
    mutate(
      label = ifelse(FDR < 0.01, paste0("**\n", round(R, 2)),
                     ifelse(FDR < 0.05, paste0("*\n", round(R, 2)),
                            round(R, 2))),
      Timepoint = timepoint_label
    ) %>%
    select(X, Y, R, label, Timepoint)
}

# 合并 T0 和 T1
cor_results_T0_plot <- prepare_plot_df(cor_results_T0, "T0")
cor_results_T1_plot <- prepare_plot_df(cor_results_T1, "T1")
plot_df <- bind_rows(cor_results_T0_plot, cor_results_T1_plot)

# 创建行名，作为X轴（之前的X作为行）
plot_df <- plot_df %>%
  mutate(RowLabel = paste(Timepoint, X, sep = "_"))
# 假设 cor_results_combined 是合并好的 T0 和 T1 相关分析结果
# 且包含列：X（module or eigengene），Y（symptom），R，P，FDR，Timepoint

library(dplyr)
# 假设 cor_results_T0 和 cor_results_T1 都已经存在并格式一致
cor_results_combined <- bind_rows(
  cor_results_T0 %>% mutate(Timepoint = "T0"),
  cor_results_T1 %>% mutate(Timepoint = "T1")
)
# 1. 保留 FDR 显著的结果（比如 FDR < 0.1）
sig_threshold <- 0.05
sig_vars <- cor_results_combined %>%
  group_by(X) %>%
  filter(any(FDR < sig_threshold)) %>%
  pull(X) %>%
  unique()

# 2. 过滤掉 X 全部不显著的变量
filtered_cor_results <- cor_results_combined %>%
  filter(X %in% sig_vars)

# 3. 添加显著性标记（* 或 **）
plot_df <- filtered_cor_results %>%
  mutate(
    sig_star = case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(round(R, 2), sig_star),
    RowLabel = paste0(Timepoint, "_", X)  # 可选：用于绘图翻转轴
  )
sig_vars <- plot_df %>%
  group_by(RowLabel) %>%
  filter(any(FDR < sig_threshold)) %>%
  pull(RowLabel) %>%
  unique()
plot_df_1 <- plot_df %>%
  filter(RowLabel %in% sig_vars)
# 控制顺序
plot_df$RowLabel <- factor(plot_df$RowLabel, levels = rev(unique(plot_df$RowLabel)))
plot_df$Y <- factor(plot_df$Y, levels = unique(plot_df$Y))  # Y轴为症状变量

# ---- 绘图 ----
# ---- 2. 绘图 ----
p <- ggplot(plot_df_1, aes(x = RowLabel, y = Y, fill = R)) +
  geom_tile(color = "grey70", width = 1, height = 1) +
  geom_text(aes(label = label), size = 3, family = "Arial") +
  scale_fill_gradient2(
    low = "lightblue", mid = "white", high = "pink",
    midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  coord_fixed() +  # 保持正方形格子
  labs(title = "Correlation of Eigengenes and Features with Symptoms (T0 & T1)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 14, base_family = "arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )
p

# ---- 输出为 PDF ----
ggsave(filename = "Combined_Heatmap_T0_T1_with_TimeStripe.pdf",
       plot = p,
       device = cairo_pdf,
       width = 10, height = 6, units = "in")




pdf("Combined_Heatmap_T0_T1_with_TimeStripe.pdf", width = 12, height = 6)
p
dev.off()
#####中介分析####
# --- 载入库 ---
library(readxl)
library(mediation)
library(dplyr)
library(purrr)
library(broom)

# --- 读取数据和变量组合 ---
#dfT1 是前面热图绘制之前相关的df数据
dfT1=df
med_list <- read_excel("T1_mediate_analysis_list.xlsx")
dfT1backup=dfT1

#dfT1 <- dfT1 %>%
 #      dplyr::filter(!sample_id__ %in% c("SYHG162", "SYHG226", "SYHG173","SX034","SYHG311"))
# --- 初始化结果列表 ---
dfT1 <- dfT1 %>%
      dplyr::filter(!sample_id__ %in% c("SX045"))
results <- list()

for (i in 1:nrow(med_list)) {
  x <- med_list$X[i]
  m <- med_list$M[i]
  y <- med_list$Y[i]
  
  # 检查变量是否都在 dfT1 中
  if (!all(c(x, m, y) %in% colnames(dfT1))) {
    warning(paste("变量缺失:", x, m, y))
    next
  }
  
  dat <- dfT1[, c(x, m, y)] %>% drop_na()
  
  # 构建模型
  f_m <- as.formula(paste0("`", m, "` ~ `", x, "`"))
  f_y <- as.formula(paste0("`", y, "` ~ `", x, "` + `", m, "`"))
  model.m <- lm(f_m, data = dat)
  model.y <- lm(f_y, data = dat)
  
  # 提取路径 a 和 b 的系数和 P 值
  summary_m <- summary(model.m)
  summary_y <- summary(model.y)
  
  # a: X → M
  a_est <- summary_m$coefficients[2, 1]
  a_p <- summary_m$coefficients[2, 4]
  
  # b: M → Y（注意是第三个系数）
  b_est <- summary_y$coefficients[3, 1]
  b_p <- summary_y$coefficients[3, 4]
  
  # 中介分析
  x_col <- names(model.m$model)[2]
  m_col <- names(model.m$model)[1]
  
  med.out <- mediate(model.m, model.y,
                     treat = x_col,
                     mediator = m_col,
                     boot = TRUE, sims = 1000)
  
  res_df <- data.frame(
    X = x,
    M = m,
    Y = y,
    a = a_est,
    a_p = a_p,
    b = b_est,
    b_p = b_p,
    ACME = med.out$d0,
    ACME_p = med.out$d0.p,
    ADE = med.out$z0,
    ADE_p = med.out$z0.p,
    Total = med.out$tau.coef,
    Total_p = med.out$tau.p,
    Prop_mediated = med.out$n0,
    Prop_p = med.out$n0.p
  )
  
  results[[i]] <- res_df
}

results_df <- do.call(rbind, results)

write.csv(results_df, "All_Mediation_ResultsT1.csv", row.names = FALSE)
write.csv(results_df, "Significant_Mediation_ResultsT1.csv", row.names = FALSE)


#####T0中介######
#dfT1 是前面热图绘制之前相关的df数据

# --- 初始化结果列表 ---


library(mediation)
library(tidyverse)
library(readxl)

# 读取变量组合列表
med_list <- read_excel("T0_mediate_analysis_list.xlsx")
results <- list()
dfT0 <- dfT0 %>%
  +     dplyr::filter(!sample_id__ %in% c("SYHG296", "SYHG295", "SYHG163"))
for (i in 1:nrow(med_list)) {
  x <- med_list$X[i]
  m <- med_list$M[i]
  y <- med_list$Y[i]
  
  # 检查变量是否都在 dfT1 中
  if (!all(c(x, m, y) %in% colnames(dfT0))) {
    warning(paste("变量缺失:", x, m, y))
    next
  }
  
  dat <- dfT0[, c(x, m, y)] %>% drop_na()
  
  # 构建模型
  f_m <- as.formula(paste0("`", m, "` ~ `", x, "`"))
  f_y <- as.formula(paste0("`", y, "` ~ `", x, "` + `", m, "`"))
  model.m <- lm(f_m, data = dat)
  model.y <- lm(f_y, data = dat)
  
  # 提取路径 a 和 b 的系数和 P 值
  summary_m <- summary(model.m)
  summary_y <- summary(model.y)
  
  # a: X → M
  a_est <- summary_m$coefficients[2, 1]
  a_p <- summary_m$coefficients[2, 4]
  
  # b: M → Y（注意是第三个系数）
  b_est <- summary_y$coefficients[3, 1]
  b_p <- summary_y$coefficients[3, 4]
  
  # 中介分析
  x_col <- names(model.m$model)[2]
  m_col <- names(model.m$model)[1]
  
  med.out <- mediate(model.m, model.y,
                     treat = x_col,
                     mediator = m_col,
                     boot = TRUE, sims = 1000)
  
  res_df <- data.frame(
    X = x,
    M = m,
    Y = y,
    a = a_est,
    a_p = a_p,
    b = b_est,
    b_p = b_p,
    ACME = med.out$d0,
    ACME_p = med.out$d0.p,
    ADE = med.out$z0,
    ADE_p = med.out$z0.p,
    Total = med.out$tau.coef,
    Total_p = med.out$tau.p,
    Prop_mediated = med.out$n0,
    Prop_p = med.out$n0.p
  )
  
  results[[i]] <- res_df
}

# 合并结果
results_df <- bind_rows(results)


# --- 输出显著结果 ---
results_sig <- results_df %>%
  filter(ACME_p < 0.05 |  Total_p < 0.05)

# --- 保存结果 ---
write.csv(results_df, "All_Mediation_ResultsT0.csv", row.names = FALSE)
write.csv(results_sig, "Significant_Mediation_ResultsT0.csv", row.names = FALSE)



#####网络nodes上升下降标注，上接网络计算和变量提取########
# 第一步：提取变量名
target_vars <- module_features[["T1_module1"]]

# 第二步：确保这些变量名都在数据中
target_vars <- intersect(target_vars, colnames(all_data_T1_HC))

# 第三步：提取子数据集
data_sub <- all_data_T1_HC[, c("DiseaseState", target_vars)]

# 第四步：差异比较（t 检验）
library(dplyr)
library(tidyr)
library(purrr)

# 转长格式
long_data <- data_sub %>%
  pivot_longer(-DiseaseState, names_to = "Variable", values_to = "Value")

# 计算均值和标准差
summary_stats <- long_data %>%
  group_by(Variable, DiseaseState) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = DiseaseState,
    values_from = c(mean, sd),
    names_glue = "{.value}_{DiseaseState}"
  )

# 进行 Wilcoxon 检验（非参数检验）
wilcox_results <- long_data %>%
  group_by(Variable) %>%
  summarise(
    p_value = wilcox.test(Value ~ DiseaseState)$p.value,
    .groups = "drop"
  )

# 合并结果并整理输出
results <- summary_stats %>%
  left_join(wilcox_results, by = "Variable") %>%
  mutate(
    mean_sd_disease = paste0(round(mean_disease, 4), " ± ", round(sd_disease, 4)),
    mean_sd_healthy = paste0(round(mean_healthy, 4), " ± ", round(sd_healthy, 4)),
    p_value = signif(p_value, 3)
  ) %>%
  select(Variable, mean_sd_disease, mean_sd_healthy, p_value)

# 输出结果
print(results)
write.csv(results,"T1_nodes_tips.csv")


# 第一步：提取变量名
target_vars <- module_features[["T0_module1"]]

# 第二步：确保这些变量名都在数据中
target_vars <- intersect(target_vars, colnames(all_data_T0_HC))

# 第三步：提取子数据集
data_sub <- all_data_T0_HC[, c("DiseaseState", target_vars)]

# 第四步：差异比较（t 检验）
library(dplyr)
library(tidyr)
library(purrr)

# 转长格式
long_data <- data_sub %>%
  pivot_longer(-DiseaseState, names_to = "Variable", values_to = "Value")

# 计算均值和标准差
summary_stats <- long_data %>%
  group_by(Variable, DiseaseState) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = DiseaseState,
    values_from = c(mean, sd),
    names_glue = "{.value}_{DiseaseState}"
  )

# 进行 Wilcoxon 检验（非参数检验）
wilcox_results <- long_data %>%
  group_by(Variable) %>%
  summarise(
    p_value = wilcox.test(Value ~ DiseaseState)$p.value,
    .groups = "drop"
  )

# 合并结果并整理输出
results <- summary_stats %>%
  left_join(wilcox_results, by = "Variable") %>%
  mutate(
    mean_sd_disease = paste0(round(mean_disease, 4), " ± ", round(sd_disease, 4)),
    mean_sd_healthy = paste0(round(mean_healthy, 4), " ± ", round(sd_healthy, 4)),
    p_value = signif(p_value, 3)
  ) %>%
  select(Variable, mean_sd_disease, mean_sd_healthy, p_value)

# 输出结果
print(results)
write.csv(results,"T0_nodes_tips.csv")

#####网络数值可视化####
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 提取模块信息的函数
get_module_data <- function(results, label) {
  list(
    label = label,
    real_corr = mean(results[["inter_view_corr"]]),
    shuffled_corr = results[["shuffled_inter_view_corr"]],
    real_auroc = results[["auroc"]],
    shuffled_auroc = results[["shuffled_auroc"]]
  )
}

# 获取模块数据
m0 <- get_module_data(
  T0_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]],
  "T0 vs HC"
)
m1 <- get_module_data(
  T1_HC_results[["keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"]][["module1"]],
  "T1 vs HC"
)

# 准备图B数据：相关性
df_corr <- bind_rows(
  data.frame(label = m0$label, type = "shuffled", value = m0$shuffled_corr),
  data.frame(label = m0$label, type = "real", value = m0$real_corr),
  data.frame(label = m1$label, type = "shuffled", value = m1$shuffled_corr),
  data.frame(label = m1$label, type = "real", value = m1$real_corr)
)

# 准备图C数据：AUROC
df_auc <- bind_rows(
  data.frame(label = m0$label, type = "shuffled", value = m0$shuffled_auroc),
  data.frame(label = m0$label, type = "real", value = m0$real_auroc),
  data.frame(label = m1$label, type = "shuffled", value = m1$shuffled_auroc),
  data.frame(label = m1$label, type = "real", value = m1$real_auroc)
)

# baseline数据（平均±SD）
baseline_corr <- df_corr %>%
  filter(type == "shuffled") %>%
  group_by(label) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    ymin = mean - sd,
    ymax = mean + sd
  )

baseline_auc <- df_auc %>%
  filter(type == "shuffled") %>%
  group_by(label) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    ymin = mean - sd,
    ymax = mean + sd
  )

# 图B：模块相关性
p_corr <- ggplot(df_corr, aes(x = label, y = value, color = type, shape = type)) +
  geom_jitter(data = subset(df_corr, type == "shuffled"), width = 0.1, alpha = 0.4, color = "grey") +
  geom_point(data = subset(df_corr, type == "real"), size = 4, color = "red") +
  geom_errorbar(data = baseline_corr, inherit.aes = FALSE,
                aes(x = label, ymin = ymin, ymax = ymax),
                width = 0.2, color = "darkgray", linewidth = 1) +
  geom_point(data = baseline_corr, inherit.aes = FALSE,
             aes(x = label, y = mean),
             color = "darkgray", size = 3) +
  labs(y = "Cross-omics correlation", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text.y = element_text(size = 12, family = "Arial", face = "bold"),
    axis.title = element_text(size = 14, family = "Arial", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )

# 图C：模块AUC
p_auc <- ggplot(df_auc, aes(x = label, y = value, color = type, shape = type)) +
  geom_jitter(data = subset(df_auc, type == "shuffled"), width = 0.1, alpha = 0.4, color = "grey") +
  geom_point(data = subset(df_auc, type == "real"), size = 4, color = "blue") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkred") +
  geom_errorbar(data = baseline_auc, inherit.aes = FALSE,
                aes(x = label, ymin = ymin, ymax = ymax),
                width = 0.2, color = "darkgray", linewidth = 1) +
  geom_point(data = baseline_auc, inherit.aes = FALSE,
             aes(x = label, y = mean),
             color = "darkgray", size = 3) +
  labs(y = "Module AUROC", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text.y = element_text(size = 12, family = "Arial", face = "bold"),
    axis.title = element_text(size = 14, family = "Arial", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )

# 合并图
combined_plot <- p_corr / p_auc + plot_layout(heights = c(1, 1))

# 输出为PDF
ggsave(
  filename = "combined_plot.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  width = 6, height = 6, units = "in"
)
###
#####相关分析#####

