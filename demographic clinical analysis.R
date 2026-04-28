###############################################################################
# Demographic and Clinical Data Analysis
# Strictly follows Methods: "Demographic and clinical data analysis"
###############################################################################

library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(readxl)

# =============================================================================
# 0. DATA LOADING
# =============================================================================
# --- Full cohort (103 MUD + 80 HC) cross-sectional data ---
# Columns expected: participantID, Group (MUD/HC),
#   Age, Sex (M/F), BMI,
#   meategg, milk, vegan, corn, snacks, drink, smoke (categorical: 0/1 or ordinal),
#   comorbidity_hypertension, comorbidity_dyslipidemia, medication_use,
#   GI_symptoms (0/1),
#   PSQI, BDI, BAI
view(df_cross) 

# --- Longitudinal MUD data (T0-T3, N=52 MUD participants) ---
# Columns expected: participantID, TimePoint (T0/T1/T2/T3),
#   Age, Sex, BMI,
#   meategg, milk, vegan, corn, snacks, drink, smoke,
#   PSQI, BDI, BAI,
#   ACSR_total, ACSR_affective, ACSR_craving, ACSR_somatic,
#   GI_symptoms
view(df_long)

# --- Multi-omics subset (29 MUD at T0 & T1, 29 HC at enrolment) ---
# Columns expected: participantID, Group (MUD-T0 / MUD-T1 / HC),
#   Age, Sex, BMI,
#   meategg, milk, vegan, corn, snacks, drink, smoke,
#   comorbidity_hypertension, comorbidity_dyslipidemia, medication_use,
#   PSQI, BDI, BAI, GI_symptoms,
#   ACSR_total, ACSR_affective, ACSR_craving, ACSR_somatic
view(df_omics)

# =============================================================================
# 1. CROSS-SECTIONAL COMPARISON: MUD vs HC (Full Cohort)
#    - t-test (normal) / Mann-Whitney U (non-normal) for continuous
#    - Chi-square for categorical
# =============================================================================

continuous_vars <- c("Age", "BMI", "PSQI", "BDI", "BAI")
categorical_vars <- c("Sex", "meategg", "milk", "vegan", "corn",
                      "snacks", "drink", "smoke",
                      "comorbidity_hypertension", "comorbidity_dyslipidemia",
                      "medication_use", "GI_symptoms")

# --- 1a. Continuous variables: normality test -> t-test or Mann-Whitney ---
cross_continuous_results <- lapply(continuous_vars, function(var) {
  mud_vals <- df_cross[[var]][df_cross$Group == "MUD"]
  hc_vals  <- df_cross[[var]][df_cross$Group == "HC"]
  
  # Shapiro-Wilk normality test (per group)
  shap_mud <- shapiro.test(mud_vals)$p.value
  shap_hc  <- shapiro.test(hc_vals)$p.value
  is_normal <- (shap_mud > 0.05) & (shap_hc > 0.05)
  
  if (is_normal) {
    test_res <- t.test(mud_vals, hc_vals)
    method   <- "t-test"
  } else {
    test_res <- wilcox.test(mud_vals, hc_vals)
    method   <- "Mann-Whitney U"
  }
  
  data.frame(
    Variable    = var,
    MUD_mean    = mean(mud_vals, na.rm = TRUE),
    MUD_sd      = sd(mud_vals, na.rm = TRUE),
    MUD_median  = median(mud_vals, na.rm = TRUE),
    MUD_IQR     = IQR(mud_vals, na.rm = TRUE),
    HC_mean     = mean(hc_vals, na.rm = TRUE),
    HC_sd       = sd(hc_vals, na.rm = TRUE),
    HC_median   = median(hc_vals, na.rm = TRUE),
    HC_IQR      = IQR(hc_vals, na.rm = TRUE),
    Method      = method,
    Statistic   = test_res$statistic,
    P_value     = test_res$p.value,
    stringsAsFactors = FALSE
  )
})
cross_continuous_results <- bind_rows(cross_continuous_results)

# --- 1b. Categorical variables: chi-square test ---
cross_categorical_results <- lapply(categorical_vars, function(var) {
  tbl <- table(df_cross$Group, df_cross[[var]])
  test_res <- chisq.test(tbl)
  data.frame(
    Variable  = var,
    Chi_sq    = test_res$statistic,
    Df        = test_res$parameter,
    P_value   = test_res$p.value,
    stringsAsFactors = FALSE
  )
})
cross_categorical_results <- bind_rows(cross_categorical_results)


# =============================================================================
# 2. LONGITUDINAL ANALYSIS WITHIN MUD (T0-T3, N=52)
#    - LMMs for continuous variables; likelihood ratio test (full vs reduced)
#    - GLMMs for categorical variables; likelihood ratio test
# =============================================================================

df_long$TimePoint <- factor(df_long$TimePoint, levels = c("T0", "T1", "T2", "T3"))

# --- 2a. Continuous variables: LMM + likelihood ratio test ---
long_continuous_vars <- c("Age", "BMI", "PSQI", "BDI", "BAI",
                          "ACSR_total", "ACSR_affective",
                          "ACSR_craving", "ACSR_somatic")

long_cont_results <- lapply(long_continuous_vars, function(var) {
  fml_full <- as.formula(paste0(var, " ~ TimePoint + (1 | participantID)"))
  fml_null <- as.formula(paste0(var, " ~ 1 + (1 | participantID)"))
  
  model_full <- lmer(fml_full, data = df_long, REML = FALSE)
  model_null <- lmer(fml_null, data = df_long, REML = FALSE)
  
  # Likelihood ratio test
  lrt <- anova(model_null, model_full)
  
  data.frame(
    Variable = var,
    Chisq    = lrt$Chisq[2],
    Df       = lrt$Df[2] - lrt$Df[1],
    P_value  = lrt[["Pr(>Chisq)"]][2],
    stringsAsFactors = FALSE
  )
})
long_cont_results <- bind_rows(long_cont_results)


# --- 2b. Categorical variables: GLMM + likelihood ratio test ---
long_categorical_vars <- c("meategg", "milk", "vegan", "corn",
                           "snacks", "drink", "smoke", "GI_symptoms")

long_cat_results <- lapply(long_categorical_vars, function(var) {
  fml_full <- as.formula(paste0(var, " ~ TimePoint + (1 | participantID)"))
  fml_null <- as.formula(paste0(var, " ~ 1 + (1 | participantID)"))
  
  model_full <- glmer(fml_full, data = df_long, family = binomial)
  model_null <- glmer(fml_null, data = df_long, family = binomial)
  
  lrt <- anova(model_null, model_full)
  
  data.frame(
    Variable = var,
    Chisq    = lrt$Chisq[2],
    Df       = lrt$Df[2] - lrt$Df[1],
    P_value  = lrt[["Pr(>Chisq)"]][2],
    stringsAsFactors = FALSE
  )
})
long_cat_results <- bind_rows(long_cat_results)


# =============================================================================
# 3. MULTI-OMICS SUBSET (29 MUD + 29 HC)
# =============================================================================

# --- 3a. Time-invariant variables: MUD vs HC (t-test / chi-square) ---
#     Use only one row per MUD participant (e.g., T0) since values are constant.
df_omics_static <- df_omics %>%
  filter(Group %in% c("MUD-T0", "HC")) %>%
  mutate(GroupBinary = ifelse(Group == "HC", "HC", "MUD"))

static_continuous <- c("Age", "BMI")
static_categorical <- c("Sex", "meategg", "milk", "vegan", "corn",
                        "snacks", "drink", "smoke",
                        "comorbidity_hypertension", "comorbidity_dyslipidemia",
                        "medication_use")

# Continuous: t-test
omics_static_cont <- lapply(static_continuous, function(var) {
  mud_vals <- df_omics_static[[var]][df_omics_static$GroupBinary == "MUD"]
  hc_vals  <- df_omics_static[[var]][df_omics_static$GroupBinary == "HC"]
  test_res <- t.test(mud_vals, hc_vals)
  data.frame(
    Variable  = var,
    MUD_mean  = mean(mud_vals, na.rm = TRUE),
    MUD_sd    = sd(mud_vals, na.rm = TRUE),
    HC_mean   = mean(hc_vals, na.rm = TRUE),
    HC_sd     = sd(hc_vals, na.rm = TRUE),
    T_stat    = test_res$statistic,
    P_value   = test_res$p.value,
    stringsAsFactors = FALSE
  )
})
omics_static_cont <- bind_rows(omics_static_cont)


# Categorical: chi-square
omics_static_cat <- lapply(static_categorical, function(var) {
  tbl <- table(df_omics_static$GroupBinary, df_omics_static[[var]])
  test_res <- chisq.test(tbl)
  data.frame(
    Variable = var,
    Chi_sq   = test_res$statistic,
    Df       = test_res$parameter,
    P_value  = test_res$p.value,
    stringsAsFactors = FALSE
  )
})
omics_static_cat <- bind_rows(omics_static_cat)




# --- 3b. Time-varying variables: GLMM across MUD-T0, MUD-T1, HC ---
#     Fixed effect: three-level Group; Random effect: participantID
df_omics$Group <- factor(df_omics$Group, levels = c("HC", "MUD-T0", "MUD-T1"))

time_varying_vars <- c("PSQI", "BDI", "BAI", "GI_symptoms")

omics_glmm_results <- lapply(time_varying_vars, function(var) {
  fml <- as.formula(paste0(var, " ~ Group + (1 | participantID)"))
  model <- lmer(fml, data = df_omics, REML = FALSE)
  
  # Type III Anova for overall Group effect
  anova_res <- car::Anova(model, type = 3)
  group_row <- anova_res["Group", ]
  
  # Pairwise fixed-effect estimates
  coef_tidy <- broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)")
  
  list(
    anova = data.frame(
      Variable = var,
      Chisq    = group_row$Chisq,
      Df       = group_row$Df,
      P_value  = group_row[["Pr(>Chisq)"]],
      stringsAsFactors = FALSE
    ),
    coefs = coef_tidy %>% mutate(Variable = var)
  )
})

omics_glmm_anova <- bind_rows(lapply(omics_glmm_results, `[[`, "anova"))
omics_glmm_coefs <- bind_rows(lapply(omics_glmm_results, `[[`, "coefs"))



# --- 3c. ACSR scores: GLMM, MUD participants only (T0 vs T1) ---
#     Fixed effect: TimePoint; Random effect: participantID
df_omics_mud <- df_omics %>%
  filter(Group %in% c("MUD-T0", "MUD-T1")) %>%
  mutate(TimePoint = ifelse(Group == "MUD-T0", "T0", "T1"),
         TimePoint = factor(TimePoint, levels = c("T0", "T1")))

acsr_vars <- c("ACSR_total", "ACSR_affective", "ACSR_craving", "ACSR_somatic")

acsr_results <- lapply(acsr_vars, function(var) {
  fml <- as.formula(paste0(var, " ~ TimePoint + (1 | participantID)"))
  model <- lmer(fml, data = df_omics_mud, REML = FALSE)
  
  # Fixed effect estimates (T1 vs T0)
  coef_tidy <- broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    mutate(Variable = var)
  
  # Descriptive stats
  desc <- df_omics_mud %>%
    group_by(TimePoint) %>%
    summarise(Mean = mean(.data[[var]], na.rm = TRUE),
              SD   = sd(.data[[var]], na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = TimePoint, values_from = c(Mean, SD))
  
  list(coefs = coef_tidy, desc = desc %>% mutate(Variable = var))
})

acsr_coefs <- bind_rows(lapply(acsr_results, `[[`, "coefs"))
acsr_desc  <- bind_rows(lapply(acsr_results, `[[`, "desc"))
