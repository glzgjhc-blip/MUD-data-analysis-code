# 安装依赖包
install.packages(c("readxl", "dplyr", "tidyr", "lme4", "broom.mixed"))
library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(broom.mixed)

# 读取数据
df <- read_excel("p2lifestylesxlsx.xlsx")

# 转换 GroupID 为数值型因子
df <- df %>%
  mutate(GroupNum = as.integer(factor(GroupID)))
#library(dplyr)
library(tidyr)

# 要填补的变量
vars <- c("meategg", "milk", "vegan", "corn", "snakes", "drink", "smoke")

# 定义函数：对单个变量按组内比例填补
fill_na_by_group_prop <- function(data, group_var, fill_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      x <- .x[[fill_var]]
      if (all(is.na(x))) {
        .x[[fill_var]] <- NA
      } else {
        prop <- prop.table(table(x, useNA = "no"))
        fill_values <- sample(names(prop), size = sum(is.na(x)), replace = TRUE, prob = prop)
        x[is.na(x)] <- as.numeric(fill_values)
        .x[[fill_var]] <- x
      }
      .x
    }) %>%
    ungroup()
}

# 循环填补所有变量
for (v in vars) {
  df <- fill_na_by_group_prop(df, group_var = "GroupID", fill_var = v)
}


# 创建统计表
freq_table <- df %>%
  pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value") %>%
  group_by(GroupID, Variable, Value) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(GroupID, Variable) %>%
  mutate(Total = sum(Count),
         Percent = round(Count / Total * 100, 1)) %>%
  ungroup()

# 查看统计结果
print(freq_table)



# 定义变量列表
vars <- c("meategg", "milk", "vegan", "corn", "snakes", "drink", "smoke")

# 用于存放结果
p_results <- data.frame(
  Variable = character(),
  Chisq = numeric(),
  Df = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# 循环每个变量，做模型对比
for (v in vars) {
  message("Processing variable: ", v)
  
  # 构造完整模型和空模型
  fml_full <- as.formula(paste0(v, " ~ GroupID + (1 | ID)"))
  fml_null <- as.formula(paste0(v, " ~ 1 + (1 | ID)"))
  
  # 拟合模型，处理出错情况
  try({
    model_full <- glmer(fml_full, data = df, family = gaussian)
    model_null <- glmer(fml_null, data = df, family = gaussian)
    
    # 似然比检验
    anova_res <- anova(model_null, model_full, test = "Chisq")
    
    # 取结果并保存
    p_results <- rbind(
      p_results,
      data.frame(
        Variable = v,
        Chisq = round(anova_res$Chisq[2], 3),
        Df = anova_res$Df[2] - anova_res$Df[1],
        P_value = signif(anova_res$`Pr(>Chisq)`[2], 4)
      )
    )
  }, silent = TRUE)
}

# 查看整体 P 值结果
print(p_results)


# 构造列联表（每行是类别，每列是组别）
tbl <- matrix(c(
  8, 6, 6, 6,     # 类别1
  4, 3, 2, 2,     # 类别2
  6, 6, 6, 4,     # 类别3
  34, 31, 31, 27  # 类别4
), nrow = 4, byrow = TRUE)

# 添加行列名称（可选）
rownames(tbl) <- c("Category1", "Category2", "Category3", "Category4")
colnames(tbl) <- c("Group1", "Group2", "Group3", "Group4")

# 查看表格
print(tbl)

# 执行卡方检验
result <- chisq.test(tbl)

# 输出结果
print(result)