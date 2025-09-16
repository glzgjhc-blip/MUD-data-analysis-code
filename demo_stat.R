# 加载必要的包
library(lme4)
library(broom.mixed)
library(dplyr)
library(readr)

# 读取数据
data <- read_csv("lmm_data.csv")

# 设置因子
data$Group <- factor(data$Group, levels = c("HC", "T0", "T1"))

# 定义变量
vars <- c("PSQI", "BDI", "BAI")

# 循环建模
results <- lapply(vars, function(var) {
  model <- lmer(as.formula(paste0(var, " ~ Group + (1 | participantID)")), data = data)
  tidy_res <- broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    mutate(Variable = var)
  return(tidy_res)
})

# 合并所有结果
final_result <- bind_rows(results)
print(final_result)
library(car)

# 创建一个函数，返回整体 F 检验结果（Group 的主效应）
get_anova_p <- function(model) {
  anova_res <- car::Anova(model, type = 3)  # type=3 更严谨，适合非平衡数据
  return(anova_res["Group", ])
}

# 分别建立三个模型
model_psqi <- lmer(PSQI ~ Group + (1 | participantID), data = data)
model_bdi  <- lmer(BDI ~ Group + (1 | participantID), data = data)
model_bai  <- lmer(BAI ~ Group + (1 | participantID), data = data)

# 获取Group的主效应P值
anova_psqi <- get_anova_p(model_psqi)
anova_bdi  <- get_anova_p(model_bdi)
anova_bai  <- get_anova_p(model_bai)

# 合并显示
anova_results <- rbind(
  data.frame(Variable = "PSQI", anova_psqi),
  data.frame(Variable = "BDI", anova_bdi),
  data.frame(Variable = "BAI", anova_bai)
)
print(anova_results)