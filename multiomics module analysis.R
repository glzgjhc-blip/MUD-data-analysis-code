###############################################################################
# Multi-omics Module Discovery and Association with Clinical Symptoms
# Strictly follows Methods section
# Note: MetOrigin 2.0 and MetaboAnalyst KEGG enrichment are platform-based
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(stringr)
  library(MintTea)
  library(mediation)
  library(purrr)
})

set.seed(3407)

# =============================================================================
# 0. HELPER FUNCTIONS
# =============================================================================

prep_add_col_prefix <- function(df, prefix) {
  id_col <- df %>% select(sample_id__)
  df <- df %>% select(-sample_id__)
  colnames(df) <- paste0(prefix, "__", colnames(df))
  bind_cols(id_col, df)
}

prep_convert_to_relative_abundances <- function(df) {
  id_col <- df %>% select(sample_id__) %>% mutate(sample_id__ = as.character(sample_id__))
  df %>% select(-sample_id__) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(total = rowSums(across(everything()))) %>%
    mutate(across(-total, ~ .x / total)) %>%
    select(-total) %>%
    bind_cols(id_col) %>% relocate(sample_id__)
}

prep_sanitize_dataset <- function(df, rare_cutoff = 0.15, mean_cutoff = NULL) {
  df <- df %>% select(where(~ n_distinct(.) > 1))
  non_zero_pct <- colSums(df != 0) / nrow(df)
  rare <- names(non_zero_pct[non_zero_pct <= rare_cutoff])
  rare <- setdiff(rare, "sample_id__")
  if (length(rare) > 0) df <- df %>% select(-all_of(rare))
  if (!is.null(mean_cutoff)) {
    means <- colMeans(df %>% select(-sample_id__), na.rm = TRUE)
    low <- names(means[means <= mean_cutoff])
    if (length(low) > 0) df <- df %>% select(-all_of(low))
  }
  df
}

remove_low_cv_features <- function(df, threshold = 0.05) {
  id_col <- df$sample_id__
  feats <- df %>% select(-sample_id__)
  cvs <- apply(feats, 2, sd, na.rm = TRUE) / abs(apply(feats, 2, mean, na.rm = TRUE))
  keep <- names(cvs)[!is.na(cvs) & cvs >= threshold]
  df %>% select(sample_id__, all_of(keep))
}


# =============================================================================
# 1. MINTTEA MODULE DISCOVERY
#    Separate analyses: HC vs MUD-T0 and HC vs MUD-T1
#    sGCCA: 10 repeats x 5 folds, keepX = 10, design = 0.7, edge = 0.8
#    Evaluation: 3 repeats x 5 folds (null: 200 shuffled modules)
# =============================================================================


# Assemble integrated data tables
meta_T0 <- metadata_T0_HC %>% select(sample_id__, DiseaseState)
meta_T1 <- metadata_T1_HC %>% select(sample_id__, DiseaseState)

all_data_T0_HC <- meta_T0 %>%
  inner_join(species_T0_HC, by = "sample_id__") %>%
  inner_join(pwy_T0_HC, by = "sample_id__") %>%
  inner_join(mtb_T0_HC, by = "sample_id__")
colnames(all_data_T0_HC) <- make.names(colnames(all_data_T0_HC), unique = TRUE)

all_data_T1_HC <- meta_T1 %>%
  inner_join(species_T1_HC, by = "sample_id__") %>%
  inner_join(pwy_T1_HC, by = "sample_id__") %>%
  inner_join(mtb_T1_HC, by = "sample_id__")
colnames(all_data_T1_HC) <- make.names(colnames(all_data_T1_HC), unique = TRUE)

cat(sprintf("  T0 vs HC: %d samples x %d features\n",
            nrow(all_data_T0_HC), ncol(all_data_T0_HC) - 2))
cat(sprintf("  T1 vs HC: %d samples x %d features\n",
            nrow(all_data_T1_HC), ncol(all_data_T1_HC) - 2))

# Run MintTea
T0_HC_results <- MintTea(
  all_data_T0_HC, view_prefixes = c("T", "P", "S"),
  sample_id_column = "sample_id__", study_group_column = "DiseaseState",
  param_diablo_keepX = 10, param_sgcca_design = 0.7,
  param_n_repeats = 10, param_n_folds = 5,
  param_edge_thresholds = 0.8,
  n_evaluation_repeats = 3, n_evaluation_folds = 5, seed = 3407
)

T1_HC_results <- MintTea(
  all_data_T1_HC, view_prefixes = c("T", "P", "S"),
  sample_id_column = "sample_id__", study_group_column = "DiseaseState",
  param_diablo_keepX = 10, param_sgcca_design = 0.7,
  param_n_repeats = 10, param_n_folds = 5,
  param_edge_thresholds = 0.8,
  n_evaluation_repeats = 3, n_evaluation_folds = 5, seed = 3407
)

# Extract module features
param_key <- "keep_10//des_0.7//nrep_10//nfol_5//ncom_5//edge_0.8"
module_features <- list(
  T0_module1 = T0_HC_results[[param_key]]$module1$features,
  T1_module1 = T1_HC_results[[param_key]]$module1$features
)

# =============================================================================
# 2. MODULE FEATURE GROUP DIFFERENCES
#    Wilcoxon rank-sum on smoking-adjusted residuals + BH FDR
#    Sensitivity: linear model adjusting for smoking
# =============================================================================


test_module_features <- function(all_data, module_feats, label) {
  feats <- intersect(module_feats, colnames(all_data))
  data_sub <- all_data[, c("DiseaseState", feats)]
  results <- lapply(feats, function(v) {
    x <- data_sub[[v]]
    g <- data_sub$DiseaseState
    wt <- wilcox.test(x ~ g, exact = FALSE)
    data.frame(Feature = v, W = as.numeric(wt$statistic),
               P = wt$p.value, stringsAsFactors = FALSE)
  })
  out <- bind_rows(results)
  out$Q <- p.adjust(out$P, method = "BH")
  out$Module <- label
  out
}

diff_T0 <- test_module_features(all_data_T0_HC, module_features$T0_module1, "HC_vs_MUD-T0")
diff_T1 <- test_module_features(all_data_T1_HC, module_features$T1_module1, "HC_vs_MUD-T1")
diff_all <- bind_rows(diff_T0, diff_T1)

cat(sprintf("  T0 module: %d/%d significant (Q<0.05)\n",
            sum(diff_T0$Q < 0.05), nrow(diff_T0)))
cat(sprintf("  T1 module: %d/%d significant (Q<0.05)\n",
            sum(diff_T1$Q < 0.05), nrow(diff_T1)))

# =============================================================================
# 3. PC1 SCORES & SPEARMAN CORRELATIONS WITH CLINICAL MEASURES
#    SVD on standardized abundance matrix → PC1
#    Spearman correlations at T0 and T1 within MUD, FDR-corrected
# =============================================================================


# Read clinical metadata

compute_eigengene <- function(expr_mat, pattern) {
  feats <- grep(pattern, colnames(expr_mat), value = TRUE)
  if (length(feats) < 2) return(rep(NA, nrow(expr_mat)))
  prcomp(scale(expr_mat[, feats]))$x[, 1]
}

run_correlation_analysis <- function(all_data, module_feats, meta_clinical, timepoint) {
  # Filter MUD participants at this timepoint
  meta_tp <- meta_clinical %>% filter(Group == timepoint)
  mud_data <- all_data %>% filter(DiseaseState == "disease")

  # Module expression matrix
  feats <- intersect(module_feats, colnames(mud_data))
  expr_mat <- mud_data %>% select(sample_id__, all_of(feats)) %>%
    column_to_rownames("sample_id__")
  expr_s <- scale(expr_mat)

  # Compute PC1 scores (eigengenes)
  eig_all <- prcomp(expr_s, center = FALSE)$x[, 1]
  eig_T <- compute_eigengene(expr_s, "^T__")
  eig_P <- compute_eigengene(expr_s, "^P__")
  eig_S <- compute_eigengene(expr_s, "^S__")

  df_eig <- data.frame(sample_id__ = rownames(expr_s),
                        eig_all = eig_all, eig_T = eig_T,
                        eig_P = eig_P, eig_S = eig_S)

  # Merge with clinical scores
  df_expr <- expr_mat %>% rownames_to_column("sample_id__")
  symptom_cols <- intersect(c("PSQI", "BDI", "BAI", "ACSR score",
                               "ACSR craving factor score",
                               "ACSR affective factor score",
                               "ACSR somatic factor score"),
                             colnames(meta_tp))
  df <- df_eig %>%
    left_join(df_expr, by = "sample_id__") %>%
    left_join(meta_tp %>% select(sample_id__, all_of(symptom_cols)), by = "sample_id__")

  # Spearman correlations
  x_vars <- c("eig_all", "eig_T", "eig_P", "eig_S", feats)
  cor_results <- expand.grid(X = x_vars, Y = symptom_cols, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      res = list(tryCatch(cor.test(df[[X]], df[[Y]], method = "spearman", exact = FALSE),
                           error = function(e) NULL)),
      R = if (!is.null(res)) as.numeric(res$estimate) else NA_real_,
      P = if (!is.null(res)) res$p.value else NA_real_
    ) %>%
    select(-res) %>% ungroup() %>%
    group_by(Y) %>%
    mutate(FDR = p.adjust(P, method = "fdr")) %>%
    ungroup() %>%
    mutate(Timepoint = timepoint)

  list(cor = cor_results, df = df, eig = df_eig)
}

cor_T0 <- run_correlation_analysis(all_data_T0_HC, module_features$T0_module1, metadata, "T0")
cor_T1 <- run_correlation_analysis(all_data_T1_HC, module_features$T1_module1, metadata, "T1")

cor_combined <- bind_rows(cor_T0$cor, cor_T1$cor)

# =============================================================================
# 4. EXPLORATORY MEDIATION ANALYSIS
#    X = microbial species/pathway, M = metabolite, Y = clinical symptom
#    + Reverse mediation; Bootstrap 1000; Covariate-adjusted sensitivity
# =============================================================================


run_mediation <- function(df, med_list_file, covariates = NULL, label = "") {
  med_list <- read_excel(med_list_file)
  results <- list()
  for (i in 1:nrow(med_list)) {
    x <- med_list$X[i]; m <- med_list$M[i]; y <- med_list$Y[i]
    if (!all(c(x, m, y) %in% colnames(df))) next
    dat <- df[, c(x, m, y, covariates)] %>% drop_na()
    if (nrow(dat) < 10) next

    # Unadjusted models
    f_m <- as.formula(paste0("`", m, "` ~ `", x, "`",
                              if (!is.null(covariates)) paste0(" + ", paste0("`", covariates, "`", collapse = " + ")) else ""))
    f_y <- as.formula(paste0("`", y, "` ~ `", x, "` + `", m, "`",
                              if (!is.null(covariates)) paste0(" + ", paste0("`", covariates, "`", collapse = " + ")) else ""))
    model.m <- lm(f_m, data = dat)
    model.y <- lm(f_y, data = dat)

    x_col <- names(model.m$model)[2]
    m_col <- names(model.m$model)[1]
    med.out <- tryCatch(
      mediate(model.m, model.y, treat = x_col, mediator = m_col, boot = TRUE, sims = 1000),
      error = function(e) NULL
    )
    if (is.null(med.out)) next

    results[[length(results) + 1]] <- data.frame(
      X = x, M = m, Y = y, Adjusted = label,
      ACME = med.out$d0, ACME_p = med.out$d0.p,
      ADE = med.out$z0, ADE_p = med.out$z0.p,
      Total = med.out$tau.coef, Total_p = med.out$tau.p,
      Prop_mediated = med.out$n0, Prop_p = med.out$n0.p,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(results)
}

# Covariates for sensitivity analysis
adj_covs <- c("Age", "Sex", "BMI", "smoking", "alcohol", "meategg", "snacks",
              "comorbidity", "medication")

# T0 mediation (unadjusted + adjusted)
if (file.exists("T0_mediate_analysis_list.xlsx")) {
  med_T0_unadj <- run_mediation(cor_T0$df, "T0_mediate_analysis_list.xlsx", label = "Unadjusted")
  adj_covs_T0 <- intersect(adj_covs, colnames(cor_T0$df))
  med_T0_adj   <- run_mediation(cor_T0$df, "T0_mediate_analysis_list.xlsx",
                                 covariates = adj_covs_T0, label = "Adjusted")
}

# T1 mediation (unadjusted + adjusted)
if (file.exists("T1_mediate_analysis_list.xlsx")) {
  med_T1_unadj <- run_mediation(cor_T1$df, "T1_mediate_analysis_list.xlsx", label = "Unadjusted")
  adj_covs_T1 <- intersect(adj_covs, colnames(cor_T1$df))
  med_T1_adj   <- run_mediation(cor_T1$df, "T1_mediate_analysis_list.xlsx",
                                 covariates = adj_covs_T1, label = "Adjusted")
}

# =============================================================================
# 5. PAIRED WILCOXON SIGNED-RANK TESTS (T0 vs T1, 29 MUD)
#    Applied to individual module features + PC1 scores, BH FDR
# =============================================================================


# Build full (unfiltered) omic tables for paired analysis
species_full <- species
pathway_full <- pathway
mtb_full     <- mtb

build_paired <- function(omic_df, meta_df, omic_name) {
  merged <- meta_df %>%
    filter(Group %in% c("T0", "T1")) %>%
    inner_join(omic_df, by = "sample_id__")
  paired_ids <- merged %>%
    group_by(participantID) %>%
    filter(n() == 2, all(c("T0", "T1") %in% Group)) %>%
    pull(participantID) %>% unique()
  merged <- merged %>% filter(participantID %in% paired_ids) %>%
    arrange(participantID, Group)
  dt0 <- merged %>% filter(Group == "T0")
  dt1 <- merged %>% filter(Group == "T1")
  stopifnot(all(dt0$participantID == dt1$participantID))
  cat(sprintf("  %s: %d paired individuals\n", omic_name, length(paired_ids)))
  list(T0 = dt0, T1 = dt1, ids = paired_ids, n = length(paired_ids))
}

species_pair <- build_paired(species_full, metadata, "Species")
pathway_pair <- build_paired(pathway_full, metadata, "Pathway")
mtb_pair     <- build_paired(mtb_full, metadata, "Metabolite")

# Union of T0 and T1 module features
all_module_feats <- union(module_features$T0_module1, module_features$T1_module1)

# Paired Wilcoxon for each feature
paired_wilcox <- function(pd, target_feats, prefix) {
  feat_cols <- intersect(target_feats, colnames(pd$T0))
  if (length(feat_cols) == 0) return(NULL)
  results <- lapply(feat_cols, function(fc) {
    v0 <- pd$T0 %>% arrange(participantID) %>% pull(all_of(fc))
    v1 <- pd$T1 %>% arrange(participantID) %>% pull(all_of(fc))
    keep <- !is.na(v0) & !is.na(v1)
    if (sum(keep) < 3) return(data.frame(Feature = fc, W = NA, P = NA, n = sum(keep)))
    wt <- wilcox.test(v1[keep], v0[keep], paired = TRUE, exact = FALSE)
    data.frame(Feature = fc, W = as.numeric(wt$statistic),
               P = wt$p.value, n = sum(keep), stringsAsFactors = FALSE)
  })
  bind_rows(results)
}

mod_T <- all_module_feats[str_starts(all_module_feats, "T__")]
mod_P <- all_module_feats[str_starts(all_module_feats, "P__")]
mod_S <- all_module_feats[str_starts(all_module_feats, "S__")]

pw_species <- paired_wilcox(species_pair, mod_T, "Species")
pw_pathway <- paired_wilcox(pathway_pair, mod_P, "Pathway")
pw_mtb     <- paired_wilcox(mtb_pair, mod_S, "Metabolite")

pw_all <- bind_rows(pw_species, pw_pathway, pw_mtb)
pw_all$Q <- p.adjust(pw_all$P, method = "BH")

# =============================================================================
# 6. DELTA CORRELATIONS (ΔT1-T0 features vs Δ clinical measures)
#    Spearman + BH FDR
# =============================================================================

symptom_candidates <- c("PSQI", "BDI", "BAI", "ACSR score",
                         "ACSR craving factor score", "ACSR affective factor score",
                         "ACSR somatic factor score")
symptom_cols <- intersect(symptom_candidates, colnames(metadata))

# Compute deltas for module features
compute_deltas <- function(pd, target_feats) {
  feat_cols <- intersect(target_feats, colnames(pd$T0))
  if (length(feat_cols) == 0) return(NULL)
  mat_t0 <- pd$T0 %>% arrange(participantID) %>% select(all_of(feat_cols)) %>% as.matrix()
  mat_t1 <- pd$T1 %>% arrange(participantID) %>% select(all_of(feat_cols)) %>% as.matrix()
  delta <- as.data.frame(mat_t1 - mat_t0)
  colnames(delta) <- paste0("delta_", feat_cols)
  delta$participantID <- pd$T0 %>% arrange(participantID) %>% pull(participantID)
  delta %>% relocate(participantID)
}

delta_species <- compute_deltas(species_pair, mod_T)
delta_pathway <- compute_deltas(pathway_pair, mod_P)
delta_mtb     <- compute_deltas(mtb_pair, mod_S)

# Symptom deltas
meta_pair <- build_paired(metadata %>% select(sample_id__), metadata, "Metadata")
symptom_cols_avail <- intersect(symptom_cols, colnames(meta_pair$T0))
if (length(symptom_cols_avail) > 0) {
  sym_t0 <- meta_pair$T0 %>% arrange(participantID) %>%
    select(all_of(symptom_cols_avail)) %>% mutate(across(everything(), as.numeric))
  sym_t1 <- meta_pair$T1 %>% arrange(participantID) %>%
    select(all_of(symptom_cols_avail)) %>% mutate(across(everything(), as.numeric))
  delta_sym <- as.data.frame(as.matrix(sym_t1) - as.matrix(sym_t0))
  colnames(delta_sym) <- paste0("delta_", make.names(symptom_cols_avail))
  delta_sym$participantID <- meta_pair$T0 %>% arrange(participantID) %>% pull(participantID)
  delta_sym <- delta_sym %>% relocate(participantID)
}

# Merge all deltas
delta_all <- list(delta_species, delta_pathway, delta_mtb, delta_sym) %>%
  compact() %>% reduce(full_join, by = "participantID")

# Spearman correlations: Δ features vs Δ symptoms
delta_feat_cols <- grep("^delta_[TPS]__", colnames(delta_all), value = TRUE)
delta_sym_cols  <- grep("^delta_(?!.*__)", colnames(delta_all), value = TRUE, perl = TRUE)

delta_cor <- expand.grid(X = delta_feat_cols, Y = delta_sym_cols, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    res = list(tryCatch(cor.test(delta_all[[X]], delta_all[[Y]], method = "spearman", exact = FALSE),
                         error = function(e) NULL)),
    R = if (!is.null(res)) as.numeric(res$estimate) else NA_real_,
    P = if (!is.null(res)) res$p.value else NA_real_
  ) %>%
  select(-res) %>% ungroup() %>%
  mutate(FDR = p.adjust(P, method = "fdr"))
