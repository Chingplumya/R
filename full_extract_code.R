## ----setup, include=FALSE--------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----载入包-------------------------------------------------------------------------------------------------
library(survival)
library(tidyverse)
library(survminer)
library(broom)


## ----原始数据载入----------------------------------------------------------------------------------------------
df_raw <- lung


## ----inst列-----------------------------------------------------------------------------------------------
sum(is.na(df_raw$inst))
df <- df_raw %>%
  drop_na(inst)
table(df$inst) #查看分布


## ----time列-----------------------------------------------------------------------------------------------
sum(is.na(df$time)) #查看缺失值
sum(df$time <= 0) #查看小于等于0的异常值
summary(df$time) #查看描述性统计，排除极大极小值
df <- df %>%
  rename(survival_days = time) %>%
  drop_na(survival_days) %>%
  filter(survival_days > 0)


## ----status列---------------------------------------------------------------------------------------------
df <- df %>%
  drop_na(status) %>%
  mutate(status = status - 1)
unique(df$status) #检验是否仅为0/1，不存在其他异常值


## ----age列------------------------------------------------------------------------------------------------
sum(is.na(df$age))
summary(df$age)


## ----sex列------------------------------------------------------------------------------------------------
sum(is.na(df$sex))
df <- df %>%
  mutate(sex = sex - 1)
unique(df$sex)


## ----ph.ecog列--------------------------------------------------------------------------------------------
sum(is.na(df$ph.ecog))
unique(df$ph.ecog)
df <- df %>%
  mutate(
    ecog_mode = as.numeric(names(sort(table(ph.ecog), decreasing = TRUE, na.rm = TRUE))[1]),
    ph.ecog = ifelse(is.na(ph.ecog), ecog_mode, ph.ecog)
  )


## ----ph.karno列-------------------------------------------------------------------------------------------
sum(is.na(df$ph.karno))
summary(df$ph.karno)
unique(df$ph.karno)
df <- df %>%
  mutate(
    karno_mode = as.numeric(names(sort(table(ph.karno), decreasing = TRUE, na.rm = TRUE))[1]),
    ph.karno = ifelse(is.na(ph.karno), karno_mode, ph.karno)
  )


## ----karno列----------------------------------------------------------------------------------------------
sum(is.na(df$pat.karno))
summary(df$pat.karno)
unique(df$pat.karno)
df <- df %>%
  mutate(
    pat.karno = ifelse(
      is.na(pat.karno),
      median(pat.karno, na.rm = T),
      pat.karno
    )
  )


## ----meal.cal列-------------------------------------------------------------------------------------------
sum(is.na(df$meal.cal))
summary(df$meal.cal)


## ----wt.loss列--------------------------------------------------------------------------------------------
sum(is.na(df$wt.loss))
summary(df$wt.loss)


## ----相关性矩阵-----------------------------------------------------------------------------------------------
df_clean <- df # 备份
cor_vars <- df %>%
  select(age, sex, ph.ecog, ph.karno, pat.karno, meal.cal, wt.loss)

cor_mat <- cor(cor_vars, use = "pairwise.complete.obs")
print(cor_mat)


## ----wt.loss回归插补-----------------------------------------------------------------------------------------
fit_wt <- lm(
  wt.loss ~ age + sex + ph.ecog + pat.karno,
  data = df_clean
)
summary(fit_wt)
idx_wt <- which(is.na(df$wt.loss))
pred_wt <- predict(fit_wt, newdata = df[idx_wt,])
set.seed(123)
noise_wt <- sample(residuals(fit_wt), size = length(idx_wt), replace = T)
df$wt.loss[idx_wt] <- pred_wt + noise_wt
summary(df$wt.loss)


## ----meal.cal回归插补----------------------------------------------------------------------------------------
fit_cal <- lm(
  meal.cal ~ age + sex + ph.ecog + pat.karno + wt.loss,
  data = df
)
summary(fit_cal)
idx_cal <- which(is.na(df$meal.cal))
pred_cal <- predict(fit_cal, newdata = df[idx_cal,])
set.seed(123)
noise_cal <- sample(residuals(fit_cal), size = length(idx_cal), replace = T)
imputed_values <- pred_cal + noise_cal
# 物理约束：热量不能为负
df$meal.cal[idx_cal] <- pmax(0, imputed_values)
summary(df$meal.cal)


## ----数据整合------------------------------------------------------------------------------------------------
df_final <- df %>%
  select(inst, survival_days, status, age, sex, ph.ecog, pat.karno, meal.cal, wt.loss) %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), labels = c("Male", "Female")),
    ph.ecog = factor(ph.ecog, levels = 0:3, labels = c("0", "1", "2", "3"))
  )


## ----性别KM曲线----------------------------------------------------------------------------------------------
surv_obj <- Surv(df_final$survival_days, df_final$status) #拟合生存对象
fit_sex <- survfit(
  surv_obj ~ sex,
  data = df_final
)
ggsurvplot(fit_sex,
           data = df_final,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.3,
           title = "Survival Curve by Sex",
           xlab = "Time (Days)",
           palette = c("#2C7BB6", "#D7191C"))


## ----ECOG评分KM曲线------------------------------------------------------------------------------------------
df_final <- df_final %>%
  mutate(
    ph.ecog = as.numeric(as.character(ph.ecog)),
    # 创建一个新的分组变量
    ecog_group = case_when(
      ph.ecog == 0 ~ "0 (Asymptomatic)",
      ph.ecog == 1 ~ "1 (Symptomatic but ambulatory)",
      ph.ecog >= 2 ~ "2-3 (Bedridden/Severe)" # 合并 2 和 3
    ),
    # 转化为因子，并指定顺序
    ecog_group = factor(ecog_group, 
                        levels = c("0 (Asymptomatic)", "1 (Symptomatic but ambulatory)", "2-3 (Bedridden/Severe)"))
  )
fit_ecog <- survfit(
  surv_obj ~ ecog_group,
  data = df_final
)
e_KM <- ggsurvplot(
  fit_ecog,
  data = df_final,
  pval = TRUE,              
  conf.int = TRUE,          
  risk.table = TRUE,        
  risk.table.height = 0.3, 
  palette = c("#2C7BB6", "#FDAE61", "#D7191C"), 
  title = "Survival Curve by ECOG Performance Status (Grouped)",
  xlab = "Time (Days)",
  ylab = "Probability of Survival",
  legend.title = "ECOG Status",
  legend.labs = c(
    "0 (Asymptomatic)",
    "1 (Ambulatory)",
    "2-3 (Bedridden)"
  ), # 统一规范图例文本，避免自动换行截断
  pval.coord = c(60, 0.42), # 手动控制P值位置，不遮挡曲线
  # 曲线美化参数
  conf.int.alpha = 0.25,
  surv.line.size = 1.1,
  risk.table.y.text.col = TRUE,
  # 关键：基础主题放ggtheme，自定义样式放入独立theme参数，不会被覆盖
  ggtheme = theme_minimal()
)
# 修改上方生存曲线画布样式，增加边距、下移图例
e_KM$plot <- e_KM$plot +
  theme(
    # 顶部、右侧大幅增加留白，防止文字裁切
    plot.margin = margin(t = 30, r = 40, b = 20, l = 10),
    # 图例移到图表底部横向排列，彻底避开顶部裁切区域
    legend.position = "bottom",
    legend.box = "horizontal",
    # 文字美化
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11)
  )

# 同步修改下方风险表字体大小
e_KM$table <- e_KM$table +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    text = element_text(size = 10.5)
  )
print(e_KM)


## ----单因素Cox回归--------------------------------------------------------------------------------------------
candidates <- c("sex", "ecog_group", "age", "pat.karno", "meal.cal", "wt.loss")
univariate_results <- lapply(candidates, function(var) {
  formula <- as.formula(paste("surv_obj ~", var))
  cox_model <- coxph(formula, data = df_final)
  summary_cox <- summary(cox_model)
  
  p_val_global <- summary_cox$waldtest["pvalue"]
  
  coef_tab <- summary_cox$coefficients
  conf_tab <- summary_cox$conf.int
  n_rows <- nrow(coef_tab)
  
  if (n_rows == 1) {
    # 情况 A: 二分类 或 连续变量 (只有 1 行系数)
    hr_val <- exp(coef_tab[1, 1])
    lower_val <- conf_tab[1, 3]
    upper_val <- conf_tab[1, 4]
    display_hr <- round(hr_val, 3)
    display_ci <- paste0(round(lower_val, 3), "-", round(upper_val, 3))
  } else {
    # 情况 B: 多分类变量 (有多行系数)
    last_row <- n_rows
    
    hr_val <- exp(coef_tab[last_row, 1])
    lower_val <- conf_tab[last_row, 3]
    upper_val <- conf_tab[last_row, 4]
    
    display_hr <- round(hr_val, 3)
    display_ci <- paste0(round(lower_val, 3), "-", round(upper_val, 4))
  }
  
  # 提取关键指标
  data.frame(
    Variable = var,
    HR = display_hr,
    CI_95 = display_ci,
    P_Value = round(p_val_global, 4),
    Significant = ifelse(p_val_global < 0.05, "Yes", "No")
  )
})
# 合并结果并排序
results_df <- bind_rows(univariate_results) %>%
  arrange(P_Value)

print(results_df)


## ----多因素Cox回归--------------------------------------------------------------------------------------------
cox_multi <- coxph(surv_obj ~ age + sex + pat.karno + ecog_group, data = df_final)
summary_multi <- summary(cox_multi)
print(summary_multi)


## ----双因素Cox回归--------------------------------------------------------------------------------------------
fit <- coxph(surv_obj ~ sex + ecog_group, data = df_final)
print(summary(fit))


## ----森林图-------------------------------------------------------------------------------------------------
model_data <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>% # 去掉截距项
  mutate(
    label = case_when(
      term == "sexFemale" ~ "Female",
      term == "ecog_group1 (Symptomatic but ambulatory)" ~ "ECOG 1\n(Symptomatic)",
      term == "ecog_group2-3 (Bedridden/Severe)" ~ "ECOG 2-3\n(Severe)",
      TRUE ~ term
    ),
    # 格式化 P 值
    p_val_fmt = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)),
    # 格式化 HR 和 CI 显示文本
    hr_text = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high)
  )


p_final <- ggplot(model_data, aes(x = estimate, y = reorder(label, estimate))) +
  # 绘制置信区间线段
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#5C5C5C", size = 0.8) +
  # 绘制 HR 点 (方块)
  geom_point(shape = 22, fill = "#2C7BB6", color = "black", size = 5, stroke = 1) +
  # 添加垂直参考线 (HR=1)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  # 在右侧添加 P 值文本
  geom_text(aes(x = max(estimate) * 1.1, label = paste0("p=", p_val_fmt)), 
            hjust = 0, size = 4.5, fontface = "italic", color = "#D7191C") +
  # 在中间添加 HR 数值文本 (可选，如果不想太挤可以去掉这行)
  geom_text(aes(x = estimate, y = label, label = hr_text), vjust = -1.5, size = 3.5) +
  
  # 坐标轴和主题
  scale_x_continuous(limits = c(0, 4.5), breaks = c(0.5, 1, 2, 4), labels = c("0.5", "1.0", "2.0", "4.0")) +
  labs(
    title = "Multivariate Cox Regression: Final Model",
    x = "Hazard Ratio (95% CI)",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 20)),
    axis.text.y = element_text(face = "bold", color = "black"), # Y 轴标签加粗
    axis.text.x = element_text(color = "black"),
    panel.grid.major.y = element_blank(), # 去掉横向网格
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(), # 去掉 Y 轴线
    axis.ticks.y = element_blank() # 去掉 Y 轴刻度
  )

# 打印最终图表
print(p_final)

