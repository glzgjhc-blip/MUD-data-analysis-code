###############################################################################
# Identification and Validation of a Gut Microbiota Marker Panel for MUD
# Strictly follows Methods section
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(vegan)
  library(randomForest)
  library(e1071)
  library(pROC)
  library(caret)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

set.seed(3407)

# =============================================================================
# 0. DATA LOADING
# =============================================================================

# Abundance tables (rows = taxa, cols = samples) for each taxonomic level
read_abundance <- function(file, level_prefix) {
  df <- read.delim(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1] <- "Taxon"
  df <- df[, !colnames(df) %in% "All"]
  df$Taxon <- paste0(level_prefix, "_", df$Taxon)
  df
}

view(phylum)
view(class)
view(order)
view(family)
view(genus)
abundance_all <- bind_rows(phylum, class_, order_, family_, genus_)

# Metadata
view(metadata)
train_meta <- subset(metadata, site == "sy")  # Discovery cohort

# Align samples
train_samples <- intersect(train_meta$SampleID, colnames(abundance_all))
train_meta <- train_meta[train_meta$SampleID %in% train_samples, ]

# Build abundance matrix (rows = samples, cols = taxa)
abun_mat <- abundance_all[, c("Taxon", train_samples)] %>%
  column_to_rownames("Taxon") %>% as.matrix() %>% t()
abun_mat <- abun_mat[train_meta$SampleID, , drop = FALSE]

# Low-abundance filtering
prevalence <- colSums(abun_mat > 0.01) / nrow(abun_mat)
abun_mat <- abun_mat[, prevalence >= 0.10, drop = FALSE]

group_vec <- factor(train_meta$Group, levels = c("HC", "MUD"))


# =============================================================================
# 1. FEATURE SELECTION: 4 Methods + Consensus (>=3/4)
#    LEfSe, ALDEx2, DESeq2, Random Forest across all taxonomic levels
# =============================================================================

# --- Method thresholds ---
LEFSE_LDA_CUTOFF    <- 2
LEFSE_P_CUTOFF      <- 0.05
ALDEX_BH_CUTOFF     <- 0.1
ALDEX_EFFECT_CUTOFF <- 0.2
DESEQ_P_CUTOFF      <- 0.05
DESEQ_LOG2FC_CUTOFF <- 1
RF_GINI_CUTOFF      <- 0.2  # Mean Decrease Gini importance > 0.2

# --- LEfSe: Kruskal-Wallis + LDA effect size ---
method_lefse <- function(X, y) {
  selected <- character(0)
  overall_mean <- colMeans(X)
  for (j in seq_len(ncol(X))) {
    kw_p <- tryCatch(kruskal.test(X[, j] ~ y)$p.value, error = function(e) NA_real_)
    if (is.na(kw_p) || kw_p >= LEFSE_P_CUTOFF) next
    class_means <- tapply(X[, j], y, mean, na.rm = TRUE)
    lda_score <- log10(max(abs(class_means - overall_mean[j]), na.rm = TRUE) * 1e6 + 1)
    if (lda_score > LEFSE_LDA_CUTOFF) selected <- c(selected, colnames(X)[j])
  }
  selected
}

# --- ALDEx2: CLR + Wilcoxon + effect size (BH < 0.1, effect > 0.2) ---
method_aldex2 <- function(X, y) {
  X_clr <- t(apply(X + 0.5, 1, function(x) log(x) - mean(log(x))))
  pvals <- effects <- numeric(ncol(X_clr))
  for (j in seq_len(ncol(X_clr))) {
    pvals[j] <- tryCatch(wilcox.test(X_clr[, j] ~ y, exact = FALSE)$p.value,
                         error = function(e) NA_real_)
    g1 <- X_clr[y == levels(y)[1], j]; g2 <- X_clr[y == levels(y)[2], j]
    sd_pool <- sqrt((var(g1) + var(g2)) / 2)
    effects[j] <- if (!is.na(sd_pool) && sd_pool > 0) abs(mean(g2) - mean(g1)) / sd_pool else 0
  }
  p_bh <- p.adjust(pvals, method = "BH")
  colnames(X_clr)[which(!is.na(p_bh) & p_bh < ALDEX_BH_CUTOFF & effects > ALDEX_EFFECT_CUTOFF)]
}

# --- DESeq2: Wilcoxon + |log2FC| > 1, P < 0.05 ---
method_deseq2 <- function(X, y) {
  pvals <- log2fc <- numeric(ncol(X))
  for (j in seq_len(ncol(X))) {
    pvals[j] <- tryCatch(wilcox.test(X[, j] ~ y, exact = FALSE)$p.value,
                         error = function(e) NA_real_)
    g1 <- mean(X[y == levels(y)[1], j], na.rm = TRUE)
    g2 <- mean(X[y == levels(y)[2], j], na.rm = TRUE)
    log2fc[j] <- log2((g2 + 0.01) / (g1 + 0.01))
  }
  colnames(X)[which(!is.na(pvals) & pvals < DESEQ_P_CUTOFF & abs(log2fc) > DESEQ_LOG2FC_CUTOFF)]
}

# --- Random Forest: Mean Decrease Gini importance > 0.2 ---
method_rf <- function(X, y) {
  rf <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)
  gini <- importance(rf, type = 2)[, 1]
  names(gini)[gini > RF_GINI_CUTOFF]
}

# --- Run all 4 methods on full discovery cohort ---
cat("\n=== Feature Selection (4 methods) ===\n")
sel_lefse  <- method_lefse(abun_mat, group_vec)
sel_aldex2 <- method_aldex2(abun_mat, group_vec)
sel_deseq2 <- method_deseq2(abun_mat, group_vec)
sel_rf     <- method_rf(abun_mat, group_vec)


# Consensus: taxa selected by >= 3 of 4 methods
all_taxa <- colnames(abun_mat)
method_count <- (all_taxa %in% sel_lefse) + (all_taxa %in% sel_aldex2) +
  (all_taxa %in% sel_deseq2) + (all_taxa %in% sel_rf)
marker_panel <- all_taxa[method_count >= 3]


# =============================================================================
# 2. COMMUNITY-LEVEL VALIDATION: Marginal PERMANOVA on Marker Panel
#    adonis2(by = "margin") with MUD status + smoking + dietary covariates
# =============================================================================


# Marker taxa abundance matrix
marker_mat <- abun_mat[, marker_panel, drop = FALSE]
bc_marker <- vegdist(marker_mat, method = "bray")

# Align covariates
meta_aligned <- train_meta[match(rownames(marker_mat), train_meta$SampleID), ]
stopifnot(all(rownames(marker_mat) == meta_aligned$SampleID))

perm_marginal <- adonis2(
  bc_marker ~ Group + smoke + meategg + snakes,
  data = meta_aligned, permutations = 999, by = "margin"
)
print(perm_marginal)

# =============================================================================
# 3. INDIVIDUAL-TAXON SENSITIVITY ANALYSES FOR MARKER PANEL
#    (1) Residualized Wilcoxon + BH FDR
#    (2) Covariate-adjusted linear model + BH FDR
# =============================================================================


# Need tax_abun.txt for individual-taxon analysis with covariates
tax_abun <- read.delim("tax_abun.txt", header = TRUE, stringsAsFactors = FALSE)
meta_cov <- read.delim("metadata_0414.txt", header = TRUE, stringsAsFactors = FALSE)

df_tax <- merge(tax_abun, meta_cov[, c("SampleID", "Group", "smoking", "meat_egg", "salted_snack")],
                by = "SampleID")
df_tax$Group <- factor(df_tax$Group, levels = c("HC", "MUD"))

# Match marker names to tax_abun column names
features <- marker_panel
features_in_data <- features[features %in% colnames(df_tax)]
# Try fixing hyphens
if (length(features_in_data) < length(features)) {
  features_fixed <- gsub("-", ".", features)
  features_in_data <- features_fixed[features_fixed %in% colnames(df_tax)]
  if (length(features_in_data) > length(features[features %in% colnames(df_tax)])) {
    features <- features_fixed
  }
}
features <- features[features %in% colnames(df_tax)]

sensitivity_results <- lapply(features, function(taxon) {
  vals <- df_tax[[taxon]]
  mud <- vals[df_tax$Group == "MUD"]; hc <- vals[df_tax$Group == "HC"]
  
  # (1) Residualized Wilcoxon
  fit_null <- lm(df_tax[[taxon]] ~ df_tax$smoking + df_tax$meat_egg + df_tax$salted_snack,
                 na.action = na.exclude)
  resid_vals <- residuals(fit_null)
  wt <- wilcox.test(resid_vals ~ df_tax$Group, exact = FALSE)
  
  # (2) Covariate-adjusted linear model
  tmp <- data.frame(y = vals, group = df_tax$Group,
                    smoking = df_tax$smoking, meat_egg = df_tax$meat_egg,
                    salted_snack = df_tax$salted_snack)
  tmp <- tmp[complete.cases(tmp), ]
  fit_lm <- lm(y ~ group + smoking + meat_egg + salted_snack, data = tmp)
  lm_sum <- summary(fit_lm)
  gr <- grep("group", rownames(lm_sum$coefficients))
  lm_est <- if (length(gr) > 0) lm_sum$coefficients[gr[1], "Estimate"] else NA
  lm_p   <- if (length(gr) > 0) lm_sum$coefficients[gr[1], "Pr(>|t|)"] else NA
  
  data.frame(Taxon = taxon,
             Direction = ifelse(mean(mud, na.rm = TRUE) > mean(hc, na.rm = TRUE),
                                "MUD_enriched", "HC_enriched"),
             Wilcox_resid_P = signif(wt$p.value, 4),
             LM_adj_Estimate = round(lm_est, 4),
             LM_adj_P = signif(lm_p, 4),
             stringsAsFactors = FALSE)
})
sensitivity_results <- bind_rows(sensitivity_results)
sensitivity_results$Wilcox_resid_Q <- signif(p.adjust(sensitivity_results$Wilcox_resid_P, method = "BH"), 4)
sensitivity_results$LM_adj_Q       <- signif(p.adjust(sensitivity_results$LM_adj_P, method = "BH"), 4)


# =============================================================================
# 4. BOOTSTRAP STABILITY ANALYSIS (1000 iterations)
#    Stratified resample → 4 methods → consensus (>=3/4) → selection frequency
# =============================================================================


N_BOOT <- 1000
N_CORES <- max(1, parallel::detectCores() - 2)

run_one_bootstrap <- function(X, y, seed) {
  set.seed(seed)
  idx_hc  <- which(y == "HC")
  idx_mud <- which(y == "MUD")
  boot_idx <- c(sample(idx_hc, length(idx_hc), replace = TRUE),
                sample(idx_mud, length(idx_mud), replace = TRUE))
  X_boot <- X[boot_idx, , drop = FALSE]
  y_boot <- y[boot_idx]
  
  m1 <- tryCatch(method_lefse(X_boot, y_boot),  error = function(e) character(0))
  m2 <- tryCatch(method_aldex2(X_boot, y_boot), error = function(e) character(0))
  m3 <- tryCatch(method_deseq2(X_boot, y_boot), error = function(e) character(0))
  m4 <- tryCatch(method_rf(X_boot, y_boot),     error = function(e) character(0))
  
  taxa <- colnames(X)
  n_methods <- (taxa %in% m1) + (taxa %in% m2) + (taxa %in% m3) + (taxa %in% m4)
  list(consensus = taxa[n_methods >= 3],
       LEfSe = m1, ALDEx2 = m2, DESeq2 = m3, RF = m4)
}

cat(sprintf("  Running %d iterations on %d cores...\n", N_BOOT, N_CORES))
t_start <- Sys.time()

cl <- makeCluster(N_CORES)
registerDoParallel(cl)
boot_results <- foreach(
  i = 1:N_BOOT,
  .packages = c("randomForest"),
  .export = c("method_lefse", "method_aldex2", "method_deseq2", "method_rf",
              "run_one_bootstrap",
              "LEFSE_LDA_CUTOFF", "LEFSE_P_CUTOFF",
              "ALDEX_BH_CUTOFF", "ALDEX_EFFECT_CUTOFF",
              "DESEQ_P_CUTOFF", "DESEQ_LOG2FC_CUTOFF", "RF_GINI_CUTOFF")
) %dopar% {
  run_one_bootstrap(abun_mat, group_vec, seed = i)
}
stopCluster(cl)


# Summarize selection frequency
tabulate_freq <- function(results, field, taxa, n) {
  counts <- table(unlist(lapply(results, function(r) r[[field]])))
  sapply(taxa, function(t) if (t %in% names(counts)) as.numeric(counts[t]) / n * 100 else 0)
}

stability_table <- tibble(
  Taxon = all_taxa,
  In_Panel = all_taxa %in% marker_panel,
  LEfSe_pct     = round(tabulate_freq(boot_results, "LEfSe",     all_taxa, N_BOOT), 1),
  ALDEx2_pct    = round(tabulate_freq(boot_results, "ALDEx2",    all_taxa, N_BOOT), 1),
  DESeq2_pct    = round(tabulate_freq(boot_results, "DESeq2",    all_taxa, N_BOOT), 1),
  RF_pct        = round(tabulate_freq(boot_results, "RF",        all_taxa, N_BOOT), 1),
  Consensus_pct = round(tabulate_freq(boot_results, "consensus", all_taxa, N_BOOT), 1)
) %>% arrange(desc(Consensus_pct))

panel_stability <- stability_table %>% filter(In_Panel)

# =============================================================================
# 5. SVM CLASSIFIER: Polynomial Kernel
#    Internal: 10-fold CV on discovery cohort
#    External: independent dataset (PRJNA970410)
#    Single-level comparisons
# =============================================================================


# --- Load full abundance data for classification ---
rawdata <- read.delim("tax_abun_all.txt", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(t(rawdata))
raw$SampleID <- rownames(raw)

# Discovery cohort
sy_meta <- subset(metadata, site == "sy")
sy <- left_join(sy_meta[, 1:2], raw, by = "SampleID")
rownames(sy) <- sy$SampleID; sy$SampleID <- NULL

# External validation cohort (PRJNA970410)
wh_meta <- subset(metadata, site == "wh")
wh <- left_join(wh_meta[, 1:2], raw, by = "SampleID")
rownames(wh) <- wh$SampleID; wh$SampleID <- NULL

# Read marker panel feature list
panel_features <- readLines("f_123_124_134_234.txt")
panel_features <- trimws(panel_features[panel_features != ""])

# Prepare training and test data with panel features
prep_data <- function(df, features) {
  df_sub <- df[, c("Group", intersect(features, colnames(df)))]
  df_sub <- df_sub[, colSums(df_sub != 0, na.rm = TRUE) > 0 | colnames(df_sub) == "Group"]
  colnames(df_sub) <- gsub("-|\\[|\\]", "_", colnames(df_sub))
  df_sub$Group <- as.factor(df_sub$Group)
  df_sub
}

sy_data <- prep_data(sy, panel_features)
wh_data <- prep_data(wh, panel_features)

# Ensure same columns
shared_cols <- intersect(colnames(sy_data), colnames(wh_data))
sy_data <- sy_data[, shared_cols]
wh_data <- wh_data[, shared_cols]

# Normalize
feat_cols <- setdiff(colnames(sy_data), "Group")
sy_scaled <- sy_data
sy_scaled[, feat_cols] <- scale(sy_data[, feat_cols])
sy_scaled[is.na(sy_scaled)] <- 0

wh_scaled <- wh_data
wh_scaled[, feat_cols] <- scale(wh_data[, feat_cols])
wh_scaled[is.na(wh_scaled)] <- 0

# SVM hyperparameters (polynomial kernel, final parameters from grid search)
SVM_KERNEL <- "polynomial"
SVM_COST   <- 128
SVM_GAMMA  <- 10
SVM_DEGREE <- 1

# --- 5a. Internal 10-fold cross-validation ---

set.seed(3407)
folds <- createFolds(sy_scaled$Group, k = 10, list = TRUE, returnTrain = FALSE)

cv_preds <- data.frame(actual = character(), prob_MUD = numeric(), stringsAsFactors = FALSE)
for (fold_idx in folds) {
  train_fold <- sy_scaled[-fold_idx, ]
  test_fold  <- sy_scaled[fold_idx, ]
  m_cv <- svm(Group ~ ., data = train_fold, type = "C-classification",
              kernel = SVM_KERNEL, cost = SVM_COST, gamma = SVM_GAMMA,
              degree = SVM_DEGREE, probability = TRUE)
  pred_cv <- predict(m_cv, test_fold, probability = TRUE)
  prob_cv <- attr(pred_cv, "probabilities")
  cv_preds <- rbind(cv_preds, data.frame(
    actual = as.character(test_fold$Group),
    prob_MUD = as.numeric(prob_cv[, "MUD"]),
    predicted = as.character(pred_cv),
    stringsAsFactors = FALSE
  ))
}

roc_cv <- pROC::roc(cv_preds$actual, cv_preds$prob_MUD, levels = c("HC", "MUD"))
cv_tab <- table(Predicted = cv_preds$predicted, Actual = cv_preds$actual)
cv_accuracy <- sum(diag(cv_tab)) / sum(cv_tab)


# --- 5b. External validation ---

m_full <- svm(Group ~ ., data = sy_scaled, type = "C-classification",
              kernel = SVM_KERNEL, cost = SVM_COST, gamma = SVM_GAMMA,
              degree = SVM_DEGREE, probability = TRUE)

pred_ext <- predict(m_full, wh_scaled, probability = TRUE)
prob_ext <- attr(pred_ext, "probabilities")
roc_ext <- pROC::roc(wh_scaled$Group, as.numeric(prob_ext[, "MUD"]),
                     levels = c("HC", "MUD"), ci = TRUE)
ext_tab <- table(Predicted = pred_ext, Actual = wh_scaled$Group)
ext_accuracy <- sum(diag(ext_tab)) / sum(ext_tab)



# --- 5c. Single-level comparisons ---
cat("\n  Single-level SVM comparisons...\n")

# Split features by taxonomic prefix
level_prefixes <- c(Phylum = "p_", Class = "c_", Order = "o_", Family = "f_", Genus = "g_")

level_results <- data.frame()
roc_list_int <- list(Full_panel = roc_cv)
roc_list_ext <- list(Full_panel = roc_ext)

for (lev_name in names(level_prefixes)) {
  prefix <- level_prefixes[lev_name]
  lev_cols <- feat_cols[startsWith(feat_cols, prefix)]
  if (length(lev_cols) < 2) next
  
  sy_lev <- sy_scaled[, c("Group", lev_cols)]
  wh_lev <- wh_scaled[, c("Group", lev_cols)]
  
  # Internal CV
  cv_preds_lev <- data.frame(actual = character(), prob_MUD = numeric(),
                             predicted = character(), stringsAsFactors = FALSE)
  for (fold_idx in folds) {
    tryCatch({
      m_lev <- svm(Group ~ ., data = sy_lev[-fold_idx, ], type = "C-classification",
                   kernel = SVM_KERNEL, cost = SVM_COST, gamma = SVM_GAMMA,
                   degree = SVM_DEGREE, probability = TRUE)
      pred_lev <- predict(m_lev, sy_lev[fold_idx, ], probability = TRUE)
      prob_lev <- attr(pred_lev, "probabilities")
      cv_preds_lev <- rbind(cv_preds_lev, data.frame(
        actual = as.character(sy_lev$Group[fold_idx]),
        prob_MUD = as.numeric(prob_lev[, "MUD"]),
        predicted = as.character(pred_lev),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
  
  if (nrow(cv_preds_lev) == 0) next
  roc_int_lev <- pROC::roc(cv_preds_lev$actual, cv_preds_lev$prob_MUD, levels = c("HC", "MUD"))
  tab_int <- table(Predicted = cv_preds_lev$predicted, Actual = cv_preds_lev$actual)
  acc_int <- sum(diag(tab_int)) / sum(tab_int)
  
  # External
  m_lev_full <- svm(Group ~ ., data = sy_lev, type = "C-classification",
                    kernel = SVM_KERNEL, cost = SVM_COST, gamma = SVM_GAMMA,
                    degree = SVM_DEGREE, probability = TRUE)
  pred_ext_lev <- predict(m_lev_full, wh_lev, probability = TRUE)
  prob_ext_lev <- attr(pred_ext_lev, "probabilities")
  roc_ext_lev <- pROC::roc(wh_lev$Group, as.numeric(prob_ext_lev[, "MUD"]),
                           levels = c("HC", "MUD"), ci = TRUE)
  tab_ext <- table(Predicted = pred_ext_lev, Actual = wh_lev$Group)
  acc_ext <- sum(diag(tab_ext)) / sum(tab_ext)
  
  roc_list_int[[lev_name]] <- roc_int_lev
  roc_list_ext[[lev_name]] <- roc_ext_lev
  
  level_results <- rbind(level_results, data.frame(
    Panel = lev_name, N_features = length(lev_cols),
    Internal_AUC = round(as.numeric(roc_int_lev$auc), 3),
    Internal_Accuracy = round(acc_int, 3),
    External_AUC = round(as.numeric(roc_ext_lev$auc), 3),
    External_CI_lo = round(roc_ext_lev$ci[1], 3),
    External_CI_hi = round(roc_ext_lev$ci[3], 3),
    External_Accuracy = round(acc_ext, 3),
    stringsAsFactors = FALSE
  ))
}

# Add full panel row
level_results <- rbind(
  data.frame(Panel = "Full_panel", N_features = length(feat_cols),
             Internal_AUC = round(as.numeric(roc_cv$auc), 3),
             Internal_Accuracy = round(cv_accuracy, 3),
             External_AUC = round(as.numeric(roc_ext$auc), 3),
             External_CI_lo = round(roc_ext$ci[1], 3),
             External_CI_hi = round(roc_ext$ci[3], 3),
             External_Accuracy = round(ext_accuracy, 3),
             stringsAsFactors = FALSE),
  level_results
)

