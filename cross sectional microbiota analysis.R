###############################################################################
# Cross-sectional Analysis of Gut Microbiota Diversity and Dominant Taxa
# Strictly follows Methods: "Cross-sectional analysis of gut microbiota
# diversity and dominant taxa composition"
###############################################################################

library(vegan)
library(dplyr)
library(ggplot2)
library(ggsci)
library(patchwork)
library(gridExtra)

# =============================================================================
# 0. DATA LOADING
# =============================================================================

view(meta) 
meta <- subset(meta, site == "sy")

# Alpha diversity (precomputed by vegan)
view(alpha_raw)
view(alpha)

# OTU table (rarefied, for beta diversity)
view(otu.tab)
view(otu)

# Phylum and genus abundance tables (rows = taxa, cols = samples)
view(pt)
view(g)
pt <- pt[, colnames(pt) %in% meta$SampleID]
g  <- g[, colnames(g)  %in% meta$SampleID]

# Factor setup
alpha$Group <- factor(alpha$Group, levels = c("HC", "MUD"))
stopifnot(all(c("Group", "smoking", "meat_egg", "salted_snack") %in% colnames(alpha)))

# =============================================================================
# 1. ALPHA DIVERSITY
#    Primary: residualized Wilcoxon (regress out covariates, then rank-sum test)
#    Sensitivity: linear model with Group + covariates
# =============================================================================

# --- 1a. Primary analysis: residualized Wilcoxon ---
null_chao1   <- lm(chao1   ~ smoking + meat_egg + salted_snack, data = alpha)
null_shannon <- lm(shannon ~ smoking + meat_egg + salted_snack, data = alpha)

alpha$chao1_resid   <- residuals(null_chao1)
alpha$shannon_resid <- residuals(null_shannon)

w_chao1   <- wilcox.test(chao1_resid   ~ Group, data = alpha, exact = FALSE)
w_shannon <- wilcox.test(shannon_resid ~ Group, data = alpha, exact = FALSE)

# --- 1b. Sensitivity analysis: linear model ---
m_chao1   <- lm(chao1   ~ Group + smoking + meat_egg + salted_snack, data = alpha)
m_shannon <- lm(shannon ~ Group + smoking + meat_egg + salted_snack, data = alpha)

coef_chao1   <- coef(summary(m_chao1))["GroupMUD", ]
coef_shannon <- coef(summary(m_shannon))["GroupMUD", ]
ci_chao1     <- confint(m_chao1)["GroupMUD", ]
ci_shannon   <- confint(m_shannon)["GroupMUD", ]


# =============================================================================
# 2. BETA DIVERSITY
#    PCoA on Bray-Curtis dissimilarity
#    PERMANOVA (adonis2) with Group + smoking + meat_egg + salted_snack
# =============================================================================

otu_t <- t(otu)
bc_dist <- vegdist(otu_t, method = "bray")

# --- 2a. PERMANOVA with covariates ---
perm_res <- adonis2(bc_dist ~ Group + smoking + meat_egg + salted_snack,
                    data = meta, permutations = 999)
cat("\n=== Beta Diversity: PERMANOVA (Bray-Curtis, 999 permutations) ===\n")
print(perm_res)

# --- 2b. PCoA visualization ---
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(PCoA1 = pcoa$points[, 1],
                      PCoA2 = pcoa$points[, 2])
pcoa_df <- data.frame(meta, pcoa_df[match(meta$SampleID, rownames(pcoa_df)), ])

eig <- pcoa$eig
pct1 <- format(100 * eig[1] / sum(eig), digits = 4)
pct2 <- format(100 * eig[2] / sum(eig), digits = 4)

perm_label <- sprintf("PERMANOVA: R² = %.3f, P = %.3f",
                      perm_res["Group", "R2"], perm_res["Group", "Pr(>F)"])

# =============================================================================
# 3. DOMINANT TAXA COMPOSITION (Top 10 Phylum + F/B, Top 10 Genus)
#    Primary: residualized Wilcoxon + BH FDR
#    Sensitivity: linear model + BH FDR
# =============================================================================

# --- Helper: build sample-by-taxa data.frame with covariates ---
build_taxa_df <- function(taxa_mat, meta_df) {
  df <- as.data.frame(t(taxa_mat))
  df$SampleID <- rownames(df)
  df <- merge(df, meta_df[, c("SampleID", "Group", "smoking", "meat_egg", "salted_snack")],
              by = "SampleID", all.x = TRUE)
  taxa_cols <- setdiff(colnames(df), c("SampleID", "Group", "smoking", "meat_egg", "salted_snack"))
  df[, taxa_cols] <- lapply(df[, taxa_cols], as.numeric)
  df$Group <- factor(df$Group, levels = c("HC", "MUD"))
  df
}

# --- Helper: get top 10 taxa by mean relative abundance ---
get_top10 <- function(taxa_mat) {
  means <- sort(rowMeans(taxa_mat, na.rm = TRUE), decreasing = TRUE)
  names(means)[1:min(10, length(means))]
}

# --- Helper: run primary + sensitivity analyses on a set of taxa ---
compare_taxa <- function(df, taxa_names) {
  results <- lapply(taxa_names, function(tax) {
    vals <- df[[tax]]
    if (sum(!is.na(vals)) < 10 || var(vals, na.rm = TRUE) == 0) {
      return(data.frame(Taxon = tax, MUD_mean = NA, MUD_sd = NA,
                        HC_mean = NA, HC_sd = NA, Direction = NA,
                        Wilcox_resid_P = NA,
                        LM_adj_Estimate = NA, LM_adj_P = NA,
                        stringsAsFactors = FALSE))
    }
    
    mud_vals <- vals[df$Group == "MUD"]
    hc_vals  <- vals[df$Group == "HC"]
    
    # Primary: residualized Wilcoxon
    fit_null <- lm(vals ~ df$smoking + df$meat_egg + df$salted_snack, na.action = na.exclude)
    resid_vals <- residuals(fit_null)
    wt <- wilcox.test(resid_vals ~ df$Group, exact = FALSE)
    
    # Sensitivity: linear model with covariates
    tmp <- data.frame(y = vals, group = df$Group,
                      smoking = df$smoking, meat_egg = df$meat_egg,
                      salted_snack = df$salted_snack)
    tmp <- tmp[complete.cases(tmp), ]
    fit_lm <- lm(y ~ group + smoking + meat_egg + salted_snack, data = tmp)
    lm_sum <- summary(fit_lm)
    group_row <- grep("group", rownames(lm_sum$coefficients))
    lm_est <- ifelse(length(group_row) > 0, lm_sum$coefficients[group_row[1], "Estimate"], NA)
    lm_p   <- ifelse(length(group_row) > 0, lm_sum$coefficients[group_row[1], "Pr(>|t|)"], NA)
    
    data.frame(
      Taxon           = tax,
      MUD_mean        = round(mean(mud_vals, na.rm = TRUE), 4),
      MUD_sd          = round(sd(mud_vals, na.rm = TRUE), 4),
      HC_mean         = round(mean(hc_vals, na.rm = TRUE), 4),
      HC_sd           = round(sd(hc_vals, na.rm = TRUE), 4),
      Direction       = ifelse(mean(mud_vals, na.rm = TRUE) > mean(hc_vals, na.rm = TRUE),
                               "MUD_enriched", "HC_enriched"),
      Wilcox_resid_P  = signif(wt$p.value, 4),
      LM_adj_Estimate = round(lm_est, 4),
      LM_adj_P        = signif(lm_p, 4),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(results)
  # BH FDR correction
  out$Wilcox_resid_Q <- signif(p.adjust(out$Wilcox_resid_P, method = "BH"), 4)
  out$LM_adj_Q       <- signif(p.adjust(out$LM_adj_P, method = "BH"), 4)
  out$Consistent <- ifelse(
    (out$Wilcox_resid_Q < 0.05 & out$LM_adj_Q < 0.05) |
      (out$Wilcox_resid_Q >= 0.05 & out$LM_adj_Q >= 0.05),
    "Yes", "No")
  out
}

# --- 3a. Phylum level ---
ptax_df <- build_taxa_df(pt, meta)
top10_phylum <- get_top10(pt)

# Add F/B ratio
ptax_df$FB_ratio <- ptax_df$Firmicutes / (ptax_df$Bacteroidetes + 1e-6)
phylum_taxa <- c(top10_phylum, "FB_ratio")

phylum_results <- compare_taxa(ptax_df, phylum_taxa)

cat("\n=== Top 10 Phylum + F/B Ratio (Adjusted) ===\n")
print(phylum_results[, c("Taxon", "Direction", "Wilcox_resid_P", "Wilcox_resid_Q",
                         "LM_adj_Estimate", "LM_adj_P", "LM_adj_Q", "Consistent")])

# --- 3b. Genus level ---
gtax_df <- build_taxa_df(g, meta)
top10_genus <- get_top10(g)

genus_results <- compare_taxa(gtax_df, top10_genus)

cat("\n=== Top 10 Genus (Adjusted) ===\n")
print(genus_results[, c("Taxon", "Direction", "Wilcox_resid_P", "Wilcox_resid_Q",
                        "LM_adj_Estimate", "LM_adj_P", "LM_adj_Q", "Consistent")])

# --- 3c. Save supplementary tables ---
phylum_results$Level <- "Phylum"
genus_results$Level  <- "Genus"
all_taxa_results <- rbind(phylum_results, genus_results)
