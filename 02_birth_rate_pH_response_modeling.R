rm(list = ls())

library(tidyverse)
library(ggplot2)
library(mgcv)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(chngpt)

source("data/analyze_threshold_models.R", encoding = "UTF-8")         # 核心建模函数
source("data/threshold_analysis_helpers.R", encoding = "UTF-8")        # 辅助函数库

# 设置参数（全部可调）
# ------------------- 数据筛选条件 -------------------
min_unique_pH <- 1      # 每个分类单元至少需要多少个不同的pH值 
min_samples <- 1        # 每个分类单元至少需要的样本总数
k_gam <- 5               # GAM平滑项的基维度，控制曲线复杂度  # 一般取5-10，数据量大可适当增大
tax_level <- "tax_id"    # 分类水平，可选：Phylum, Class, Order, Family, Genus

# ------------------- 阈值模型估计参数 -------------------
ci_bootstrap_size <- 1000#1000  # 阈值置信区间bootstrap次数
gam_ci_boot <- 999         # 999 GAM置信带bootstrap次数（0表示不计算）# GAM整体趋势bootstrap次数 增大可提高置信带精度，但会增加计算时间
chngpt_est_method <- "smoothapprox"  # 阈值估计方法："grid"（网格搜索，稳定通用）或 "smoothapprox"（平滑近似，仅部分模型）
n_cores <- 6            # 并行计算核心数 # 若电脑支持可设为 detectCores()-1，但注意内存消耗

# ------------------- 显著性检验参数 -------------------
alpha_sig <- 0.05        # 显著性阈值，p值小于此值认为显著 # 常用0.05或0.01，可调
mc_n <- 1e5            # 阈值模型蒙特卡洛模拟次数 # 增大可提高p值精度（如1e5），但耗时增加；减小可加快速度（如1e4）
min_seg_points <- 5      # 分段线性回归每侧最少样本数 # 确保每侧有足够数据拟合直线，一般保持5-10
min_seg_unique <- 5      # 分段线性回归每侧最少自变量唯一值个数 # 确保每侧有足够变化，一般保持5-10
selection_method <- "p_then_AICc"  # 模型筛选方式 # "p_then_AICc"：先筛选显著模型，再按AICc选最佳  # "AICc"：直接在所有模型中按AICc选最佳（忽略显著性）

# 输出目录
output_dir <- paste0("./Fig/Threshold_Results_", tax_level)
fig_dir <- paste0(output_dir, "/Figures")
stats_dir <- paste0(output_dir, "/Statistics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 步骤1：加载数据并预处理
# =============================================================================
design <- read.table('data/design.txt', sep='\t', header=TRUE, row.names=1, comment.char='', check.names=FALSE)
Brith_rate_act <- read.csv("data/Brith_rate_indv_long_active.csv", header=TRUE, row.names=1, check.names=FALSE)

#要保留全部数据，非活的变成0
# 1. 获取所有唯一的 treatment 和 tax_id
treatments <- unique(Brith_rate_act$treatment)
tax_ids <- unique(Brith_rate_act$tax_id)
# 2. 生成所有组合（笛卡尔积）
complete_grid <- expand.grid(treatment = treatments, tax_id = tax_ids, stringsAsFactors = FALSE)
# 3. 左连接原始数据，填充缺失值为0
df_long_filled <- complete_grid %>%
  left_join(Brith_rate_act %>% select(treatment, tax_id, Brith_rate), by = c("treatment", "tax_id")) %>%
  mutate(Brith_rate = ifelse(is.na(Brith_rate), 0, Brith_rate))
# 4. 提取每个 tax_id 的唯一分类信息（去重）
tax_info <- Brith_rate_act %>%
  select(tax_id, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  distinct()
# 5. 合并分类信息，得到最终长表
df_final <- df_long_filled %>%
  left_join(tax_info, by = "tax_id")

#活性非活都保留
df <- merge(design, df_final, by = "treatment", all.x = TRUE)

# 过滤Unassigned
df <- df %>%
  filter(
    !is.na(.data[[tax_level]]),
    !tolower(trimws(.data[[tax_level]])) %in% c("unassigned", "unclassified", "unknown", "")
  )

# 去除缺失值
df_clean <- df %>%
  filter(!is.na(pH), !is.na(Brith_rate), !is.na(.data[[tax_level]]))

# 查看各分类单元统计
tax_stats <- df_clean %>%
  group_by(across(all_of(tax_level))) %>%
  summarise(
    n_total = n(),
    n_unique_pH = length(unique(pH)),
    pH_min = min(pH, na.rm = TRUE),
    pH_max = max(pH, na.rm = TRUE),
    pH_range = pH_max - pH_min,
    .groups = 'drop'
  ) %>%
  arrange(desc(n_unique_pH))

print(tax_stats)
write.table(tax_stats, file.path(stats_dir, "all_taxa.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)  

# 筛选符合条件的分类单元
#当用全部数据（包括非活性）时，这里的筛选不起作用
selected_taxa <- tax_stats %>%
  filter(n_unique_pH >= min_unique_pH,
         n_total >= min_samples,
         pH_range >= 0) %>%
  pull(.data[[tax_level]])

cat("\n符合筛选条件的", tax_level, "数量:", length(selected_taxa), "\n")


# 保留筛选后的数据
df_clean <- df_clean %>% filter(.data[[tax_level]] %in% selected_taxa)

# 读取分类学信息（门映射）——可选
tax_file <- './data/OTU/GTDB_taxonomy_adj.txt'
if (file.exists(tax_file)) {
  tax <- read.table(tax_file, sep='\t', header=TRUE, comment.char='', row.names=1)
  if (!(tax_level %in% colnames(tax))) stop(paste("分类学文件中没有", tax_level, "列"))
  tax_unique <- tax[!duplicated(tax[[tax_level]]), ]
  phylum_map <- setNames(tax_unique$Phylum, tax_unique[[tax_level]])
} else {
  phylum_map <- c()
}

# =============================================================================
# 步骤2：对每个分类单元运行 analyze_threshold_models，保存结果到列表
# =============================================================================
res_list <- list()

for (i in seq_along(selected_taxa)) {
  tax_name <- selected_taxa[i]
  cat("正在分析", i, "/", length(selected_taxa), ":", tax_name, "\n")
  
  tax_data <- df_clean %>% filter(.data[[tax_level]] == tax_name)
  
  # 基本检查
  if (nrow(tax_data) < min_samples || length(unique(tax_data$pH)) < min_unique_pH) {
    cat("  样本量不足，跳过\n")
    next
  }
  
  res <- tryCatch(
    analyze_threshold_models(
      data = tax_data,
      x_col = "pH",
      y_col = "Brith_rate",
      family = "gaussian",
      fit_linear = TRUE,
      fit_quadratic = TRUE,
      fit_gam = TRUE,
      threshold_types = c("segmented", "stegmented", "step"),
      overall_trend = "best", #整体趋势线来源，可选 "best", "gam", "none"，一般用 "best" 使整体线等于最佳模型
      gam_k = k_gam,
      gam_bs = "tp", #GAM 平滑基类型，默认 "tp"
      gam_ci_boot = gam_ci_boot, #GAM 置信带 bootstrap 次数
      seed = 123, #随机种子，保证可重复性	设为 123 + i 使每个分类种子不同
      n_cores = n_cores,
      parallel_mode = "auto",
      chngpt_var_type = "bootstrap", #阈值模型方差估计方法，常用 "bootstrap"
      ci_bootstrap_size = ci_bootstrap_size, #阈值置信区间 bootstrap 次数
      chngpt_est_method = chngpt_est_method, #阈值估计方法，"grid" 或 "smoothapprox"
      chngpt_search_lower = min(tax_data$pH), #阈值搜索范围的下限和上限	设为当前数据的最小/最大 pH，无需手动干预
      chngpt_search_upper = max(tax_data$pH),
      chngpt_grid_length = 100, #网格搜索时的网格点数	默认100，精度足够
      alpha = 0.05, #置信区间显著性水平，默认0.05
      report_segment_fit = TRUE, #是否返回阈值前后段描述性R²
      grid_n = 500, #预测网格点数，用于绘图数据
      do_chngpt_test = FALSE #是否执行阈值效应检验	设为 FALSE，因为我们在步骤3单独计算p值
    ),
    error = function(e) { cat("  出错:", e$message, "\n"); NULL }
  )
  
  if (!is.null(res)) {
    res_list[[tax_name]] <- res
  }
}

# 可选：保存结果列表到文件，以便后续加载（避免重复运行）
saveRDS(res_list, file = file.path(stats_dir, "res_list.rds"))
cat("\n成功运行", length(res_list), "个分类单元，结果已保存至:", file.path(stats_dir, "res_list.rds"), "\n")

# =============================================================================
# 步骤3：为每个模型计算p值，并添加显著性标记
# =============================================================================
# 如果是从文件加载，可以跳过步骤2直接加载
# res_list <- readRDS(file.path(stats_dir, "res_list.rds"))

all_model_tables <- list()

for (tax_name in names(res_list)) {
  cat("计算p值:", tax_name, "\n")
  res <- res_list[[tax_name]]
  
  # 从原始数据集中重新提取该分类的数据（包含原始列名 pH 和 Brith_rate）
  tax_data <- df_clean %>% filter(.data[[tax_level]] == tax_name)
  if (nrow(tax_data) == 0) {
    warning(tax_name, " 在 df_clean 中找不到数据，跳过")
    next
  }
  
  model_table <- compute_model_pvalues(res, tax_data, tax_name, mc_n, alpha_sig)
  all_model_tables[[tax_name]] <- model_table
}

# 合并所有模型表并保存
if (length(all_model_tables) > 0) {
  all_models_df <- bind_rows(all_model_tables)
  all_models_df <- all_models_df %>% select(Tax_Name, Model, Class, Family, EstMethod, N, LogLik, K, AIC, AICc, R2_gaussian, p_value, Significant, everything())
  write.table(all_models_df, file.path(stats_dir, "All_Model_Comparisons.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n模型比较总表已保存至:", file.path(stats_dir, "All_Model_Comparisons.tsv"), "\n")
}

#all_models_df <- read.table(file.path(stats_dir, "All_Model_Comparisons.tsv"),sep = "\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

#all_model_tables <- split(all_models_df, all_models_df$Tax_Name)# 关键一步：按 Tax_Name 拆回 list

# =============================================================================
# 步骤4：对每个分类单元选择最佳模型（基于p值和AICc）
# =============================================================================
best_models_info <- list()

for (tax_name in names(res_list)) {
  model_table <- all_model_tables[[tax_name]]
  best_info <- select_best_model_by_criteria(model_table, selection_method)
  best_models_info[[tax_name]] <- best_info
}

# =============================================================================
# 步骤5：为每个最佳模型生成图形并保存到对应文件夹
# =============================================================================
# 初始化图形收集列表
threshold_plots <- list()
threshold_other_plots <- list()
quadgam_plots <- list()
linear_plots <- list()
nonsig_plots <- list()

for (tax_name in names(res_list)) {
  cat("处理图形:", tax_name, "\n")
  
  res <- res_list[[tax_name]]
  best_info <- best_models_info[[tax_name]]
  model_table <- all_model_tables[[tax_name]]
  
  # 从原始数据中提取当前分类的数据（包含原始列名 pH 和 Brith_rate）
  tax_data <- df_clean %>% filter(.data[[tax_level]] == tax_name)
  if (nrow(tax_data) == 0) {
    warning("在 df_clean 中找不到 ", tax_name, " 的数据，跳过")
    next
  }
  
  phylum_name <- ifelse(tax_name %in% names(phylum_map), phylum_map[tax_name], "Unknown")
  
  # 处理无显著模型的情况
  if (!best_info$any_significant && selection_method == "p_then_AICc") {
    p <- ggplot(tax_data, aes(x = pH, y = Brith_rate)) +
      geom_point(color = "gray70", alpha = 0.8, size = 2) +
      labs(title = paste0(tax_name, " (", phylum_name, ")"),
           subtitle = "No significant model",
           x = "pH", y = "Brith_rate") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"))
    
    nonsig_plots[[tax_name]] <- p
    sub_dir <- "NonSignificant"
    best_model_name <- "None"
    
    # 保存图形
    plot_path <- paste0(fig_dir, "/", sub_dir, "/", tax_name, "_", best_model_name, ".pdf")
    dir.create(paste0(fig_dir, "/", sub_dir), showWarnings = FALSE, recursive = TRUE)
    ggsave(plot_path, p, width = 6, height = 6, device = "pdf")
    
  } else if (best_info$any_significant) {
    # 有显著模型，创建图形
    plot_info <- create_model_plot(
      best_model_name = best_info$best_model_name,
      best_row = best_info$best_row,
      res = res,
      tax_data = tax_data,      # 使用正确的数据框
      tax_name = tax_name,
      phylum_name = phylum_name,
      min_seg_points = min_seg_points,
      min_seg_unique = min_seg_unique
    )
    
    p <- plot_info$p
    sub_dir <- plot_info$sub_dir
    model_type <- plot_info$model_type
    best_model_name <- best_info$best_model_name
    
    # 将图形存入对应列表
    if (model_type == "Threshold_SinglePeak") {
      threshold_plots[[tax_name]] <- p
    } else if (model_type == "Threshold_Other") {
      threshold_other_plots[[tax_name]] <- p
    } else if (model_type == "Quadratic/GAM") {
      quadgam_plots[[tax_name]] <- p
    } else if (model_type == "Linear") {
      linear_plots[[tax_name]] <- p
    }
    
    # 保存单个图形
    plot_path <- paste0(fig_dir, "/", sub_dir, "/", tax_name, "_", best_model_name, ".pdf")
    dir.create(paste0(fig_dir, "/", sub_dir), showWarnings = FALSE, recursive = TRUE)
    ggsave(plot_path, p, width = 6, height = 6, device = "pdf")
  }
}

# 保存组合图
graphics.off()
save_combined_plots(threshold_plots, fig_dir, "Threshold_SinglePeak")
save_combined_plots(threshold_other_plots, fig_dir, "Threshold_Other")
save_combined_plots(quadgam_plots, fig_dir, "QuadraticGAM")
save_combined_plots(linear_plots, fig_dir, "Linear")
save_combined_plots(nonsig_plots, fig_dir, "NonSignificant")

# =============================================================================
# 步骤6：汇总每个分类的统计信息
# =============================================================================
stats_summary <- data.frame()

for (tax_name in names(res_list)) {
  res <- res_list[[tax_name]]
  best_info <- best_models_info[[tax_name]]
  model_table <- all_model_tables[[tax_name]]
  
  tax_data <- res$raw_data
  phylum_name <- ifelse(tax_name %in% names(phylum_map), phylum_map[tax_name], "Unknown")
  
  # 如果是无显著模型，则相关信息为NA
  if (!best_info$any_significant && selection_method == "p_then_AICc") {
    best_model_name <- "None"
    model_type <- "NonSignificant"
    optimum_x <- NA
    left_slope <- NA; right_slope <- NA; left_r2 <- NA; right_r2 <- NA
  } else {
    # 从之前的信息中获取（需要从 plot_info 中提取，但此处没有 plot_info，需要从 res 中提取）
    # 实际上，我们可以在步骤5中同时收集统计信息，但为了分离，我们重新计算一次
    # 简单起见，我们从模型表中提取阈值信息（如果是阈值模型）
    best_model_name <- best_info$best_model_name
    if (grepl("^chngpt_", best_model_name)) {
      best_row <- model_table[model_table$Model == best_model_name, , drop = FALSE]
      optimum_x <- best_row$Threshold
      left_slope <- best_row$Slope_Pre
      right_slope <- best_row$Slope_Post
      left_r2 <- if (!is.null(res$segment_fit)) res$segment_fit$R2_pre else NA
      right_r2 <- if (!is.null(res$segment_fit)) res$segment_fit$R2_post else NA
      
      if (!is.na(left_slope) && !is.na(right_slope) && left_slope > 0 && right_slope < 0) {
        model_type <- "Threshold_SinglePeak"
      } else {
        model_type <- "Threshold_Other"
      }
    } else if (best_model_name %in% c("Quadratic", "GAM")) {
      # 对于二次/GAM，需要重新计算最优值和分段回归，但我们已经有了 plot_info 中的属性
      # 为简化，这里从模型表中无法直接得到，我们暂时设为 NA
      optimum_x <- NA
      left_slope <- NA; right_slope <- NA; left_r2 <- NA; right_r2 <- NA
      model_type <- "Quadratic/GAM"
    } else if (best_model_name == "Linear") {
      optimum_x <- NA
      left_slope <- NA; right_slope <- NA; left_r2 <- NA; right_r2 <- NA
      model_type <- "Linear"
    } else {
      optimum_x <- NA
      left_slope <- NA; right_slope <- NA; left_r2 <- NA; right_r2 <- NA
      model_type <- "Unknown"
    }
  }
  
  stats_row <- data.frame(
    Tax_Level = tax_level,
    Tax_Name = tax_name,
    Phylum = phylum_name,
    n_total = nrow(tax_data),
    n_unique_pH = length(unique(tax_data$pH)),
    pH_min = min(tax_data$pH),
    pH_max = max(tax_data$pH),
    pH_range = max(tax_data$pH) - min(tax_data$pH),
    best_model = best_model_name,
    model_type = model_type,
    any_significant = best_info$any_significant,
    optimum_pH = optimum_x,
    left_slope = left_slope,
    right_slope = right_slope,
    left_r2 = left_r2,
    right_r2 = right_r2,
    stringsAsFactors = FALSE
  )
  stats_summary <- rbind(stats_summary, stats_row)
}

# 保存详细统计表
write.table(stats_summary, file.path(stats_dir, "Detailed_Statistics.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 生成汇总统计
summary_stats <- data.frame(
  Parameter = c(
    "分析分类水平", "最小pH唯一值", "最小样本数",
    "模型筛选方式",
    paste0("原始", tax_level, "数量"), paste0("初步筛选", tax_level, "数量"),
    "成功分析总数", "单峰阈值模型数", "其他阈值模型数", "二次/GAM模型数", "线性模型数", "无显著模型数",
    "单峰阈值模型占比", "其他阈值模型占比", "二次/GAM占比", "线性占比", "无显著模型占比",
    "平均最适pH（仅阈值+二次/GAM）"
  ),
  Value = c(
    tax_level, min_unique_pH, min_samples,
    selection_method,
    nrow(tax_stats), length(selected_taxa),
    nrow(stats_summary),
    sum(stats_summary$model_type == "Threshold_SinglePeak", na.rm = TRUE),
    sum(stats_summary$model_type == "Threshold_Other", na.rm = TRUE),
    sum(stats_summary$model_type == "Quadratic/GAM", na.rm = TRUE),
    sum(stats_summary$model_type == "Linear", na.rm = TRUE),
    sum(stats_summary$model_type == "NonSignificant", na.rm = TRUE),
    paste0(round(mean(stats_summary$model_type == "Threshold_SinglePeak", na.rm = TRUE)*100, 1), "%"),
    paste0(round(mean(stats_summary$model_type == "Threshold_Other", na.rm = TRUE)*100, 1), "%"),
    paste0(round(mean(stats_summary$model_type == "Quadratic/GAM", na.rm = TRUE)*100, 1), "%"),
    paste0(round(mean(stats_summary$model_type == "Linear", na.rm = TRUE)*100, 1), "%"),
    paste0(round(mean(stats_summary$model_type == "NonSignificant", na.rm = TRUE)*100, 1), "%"),
    round(mean(stats_summary$optimum_pH, na.rm = TRUE), 3)
  )
)

write.table(summary_stats, file.path(stats_dir, "Summary_Statistics.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
print(summary_stats)

# =============================================================================
# 步骤7：生成模型类型分布图和最适pH分布图
# =============================================================================
model_type_dist <- ggplot(stats_summary, aes(x = model_type, fill = model_type)) +
  geom_bar(alpha = 0.7) +
  scale_fill_manual(values = c(
    "Threshold_SinglePeak" = "#4CAF50",
    "Threshold_Other" = "#FF9800",
    "Quadratic/GAM" = "#2196F3",
    "Linear" = "#F44336",
    "NonSignificant" = "#9E9E9E"
  )) +
  labs(title = paste("Model Type Distribution by", tax_level),
       x = "Model Type", y = paste("Number of", tax_level, "Units")) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(fig_dir, "/Model_Type_Distribution.pdf"), model_type_dist, width = 8, height = 5, device = "pdf")

if (any(!is.na(stats_summary$optimum_pH))) {
  optimum_dist <- ggplot(filter(stats_summary, !is.na(optimum_pH)), aes(x = optimum_pH)) +
    geom_histogram(bins = 15, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean(stats_summary$optimum_pH, na.rm = TRUE),
               color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste("Distribution of Optimum pH (Threshold + Quadratic/GAM) by", tax_level),
         subtitle = paste("Mean =", round(mean(stats_summary$optimum_pH, na.rm = TRUE), 2),
                          " | Median =", round(median(stats_summary$optimum_pH, na.rm = TRUE), 2),
                          " | n =", sum(!is.na(stats_summary$optimum_pH))),
         x = "Optimum pH", y = paste("Number of", tax_level, "Units")) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(paste0(fig_dir, "/Optimum_pH_Distribution.pdf"), optimum_dist, width = 8, height = 6, device = "pdf")
}

cat("\n所有分析完成！结果已保存在:", output_dir, "\n")

#
