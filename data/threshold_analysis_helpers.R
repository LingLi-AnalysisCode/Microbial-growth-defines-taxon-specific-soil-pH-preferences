# =============================================================================
# 阈值分析辅助函数库（可调节参数版，统一绘图风格：标题居中，副标题含模型和最优pH，左右段标注）
# 包含：两段线性回归、模型p值计算、组合图保存、阈值模型绘图、二次/GAM模型绘图、线性模型绘图
# 依赖包：tidyverse, ggplot2, mgcv, chngpt, gridExtra
# LingLi 2026-03-18
# =============================================================================

# 确保必要包已安装
required_packages <- c("tidyverse", "ggplot2", "mgcv", "chngpt", "gridExtra")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("请先安装包:", pkg))
  }
}

# =============================================================================
# 1. 两段线性回归函数（参数可调节）
# 功能：在给定的分割点处，对数据左右两侧分别拟合线性回归，返回斜率、截距、R²、p值等
# 参数：
#   data        : 数据框，包含自变量和因变量
#   x_var       : 自变量列名（字符串）
#   y_var       : 因变量列名（字符串）
#   split_point : 分割点数值
#   min_points  : 每侧最少样本数（默认5）
#   min_unique_x: 每侧最少自变量唯一值个数（默认5）
# 返回值：
#   列表，包含 left 和 right 两个子列表，每个子列表有 slope, intercept, r2, p_value, n
#   若某侧数据不足，则该侧为 NULL
# =============================================================================
fit_two_segment_regression <- function(data, x_var, y_var, split_point, 
                                       min_points = 5, min_unique_x = 5) {
  data$segment <- ifelse(data[[x_var]] <= split_point, "left", "right")
  results <- list(left = NULL, right = NULL)
  
  # 左侧回归
  left_data <- data %>% filter(segment == "left")
  if (nrow(left_data) >= min_points && length(unique(left_data[[x_var]])) >= min_unique_x) {
    left_lm <- lm(as.formula(paste(y_var, "~", x_var)), data = left_data)
    left_summary <- summary(left_lm)
    left_p <- tryCatch(coef(left_summary)[2, 4], error = function(e) NA)
    results$left <- list(
      slope = coef(left_lm)[2],
      intercept = coef(left_lm)[1],
      r2 = left_summary$r.squared,
      p_value = left_p,
      n = nrow(left_data)
    )
  }
  
  # 右侧回归
  right_data <- data %>% filter(segment == "right")
  if (nrow(right_data) >= min_points && length(unique(right_data[[x_var]])) >= min_unique_x) {
    right_lm <- lm(as.formula(paste(y_var, "~", x_var)), data = right_data)
    right_summary <- summary(right_lm)
    right_p <- tryCatch(coef(right_summary)[2, 4], error = function(e) NA)
    results$right <- list(
      slope = coef(right_lm)[2],
      intercept = coef(right_lm)[1],
      r2 = right_summary$r.squared,
      p_value = right_p,
      n = nrow(right_data)
    )
  }
  
  return(results)
}

# =============================================================================
# 2. 获取模型显著性p值函数（可调节模拟次数）
# 功能：根据模型类型提取或计算p值
# 参数：
#   model      : 模型对象（lm, gam, chngptm 等）
#   data       : 原始数据框（用于阈值模型检验）
#   model_name : 模型名称（字符串）
#   x_col      : 自变量列名（默认 "pH"）
#   y_col      : 因变量列名（默认 "Brith_rate"）
#   mc_n       : 蒙特卡洛模拟次数（仅用于阈值模型，默认 50000）
# 返回值：
#   数值型p值，若无法计算则返回 NA
# =============================================================================
get_model_pvalue <- function(model, data, model_name, x_col = "pH", y_col = "Brith_rate", mc_n = 50000) {
  tryCatch({
    if (grepl("^Linear", model_name)) {
      s <- summary(model)
      return(s$coefficients[2, 4])   # 斜率 p 值
    } else if (grepl("^Quadratic", model_name)) {
      s <- summary(model)
      return(s$coefficients[3, 4])   # 二次项 p 值
    } else if (grepl("^GAM", model_name)) {
      s <- summary(model)
      return(s$s.pv[1])               # 平滑项 p 值
    } else if (grepl("^chngpt_", model_name)) {
      thr_type <- sub("^chngpt_", "", model_name)
      test_res <- chngpt::chngpt.test(
        formula.null = as.formula(paste(y_col, "~ 1")),
        formula.chngpt = as.formula(paste("~", x_col)),
        data = data,
        family = "gaussian",
        type = thr_type,
        test.statistic = "lr",
        p.val.method = "MC",
        mc.n = mc_n
      )
      if (!is.null(test_res) && !is.null(test_res$p.value)) {
        return(test_res$p.value)
      } else {
        warning("阈值模型 p 值返回为空")
        return(NA)
      }
    } else {
      return(NA)
    }
  }, error = function(e) {
    warning("get_model_pvalue 错误: ", e$message)
    NA
  })
}

# =============================================================================
# 3. 组合图保存函数
# 功能：将多个 ggplot 对象按指定行列分页组合，保存为 PDF
# 参数：
#   plot_list     : 包含 ggplot 对象的列表
#   output_dir    : 输出目录（将在此目录下创建 category_name 子文件夹）
#   category_name : 类别名称，用于文件夹命名
#   rows_per_page : 每页行数（默认 2）
#   cols_per_page : 每页列数（默认 3）
# 返回值：无
# =============================================================================
save_combined_plots <- function(plot_list, output_dir, category_name, 
                                rows_per_page = 2, cols_per_page = 3) {
  if (length(plot_list) == 0) {
    message("类别 ", category_name, " 无图形，跳过组合图")
    return()
  }
  
  n_per_page <- rows_per_page * cols_per_page
  n_plots <- length(plot_list)
  n_pages <- ceiling(n_plots / n_per_page)
  
  message("\n创建 ", category_name, " 分页图形，共 ", n_pages, " 页...")
  
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * n_per_page + 1
    end_idx <- min(page * n_per_page, n_plots)
    page_plots <- plot_list[start_idx:end_idx]
    
    # 补空白图
    if (length(page_plots) < n_per_page) {
      n_blank <- n_per_page - length(page_plots)
      for (i in 1:n_blank) {
        page_plots[[length(page_plots) + 1]] <- ggplot2::ggplot() + ggplot2::theme_void()
      }
    }
    
    # 动态调整图形尺寸
    if (n_per_page <= 6) {
      fig_width <- 10; fig_height <- 8; title_size <- 14
    } else if (n_per_page <= 12) {
      fig_width <- 12; fig_height <- 10; title_size <- 12
    } else {
      fig_width <- 15; fig_height <- 12; title_size <- 10
    }
    
    # 组合图形
    combined_plot <- gridExtra::grid.arrange(
      grobs = page_plots,
      ncol = cols_per_page,
      nrow = rows_per_page,
      top = grid::textGrob(
        paste("Threshold Analysis - ", category_name, " - Page ", page, " of ", n_pages),
        gp = grid::gpar(fontsize = title_size, fontface = "bold")
      )
    )
    
    # 创建输出目录
    combined_dir <- file.path(output_dir, category_name, "Combined_Pages")
    dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
    file_path <- file.path(combined_dir, paste0("All_", category_name, "_Page_", page, ".pdf"))
    
    # 手动打开 PDF 设备，绘制，关闭
    grDevices::pdf(file_path, width = fig_width, height = fig_height)
    grid::grid.draw(combined_plot)
    grDevices::dev.off()
    
    message("    已保存第 ", page, " 页: ", file_path)
  }
}

# =============================================================================
# 4. 阈值模型绘图函数（增加左右段线性回归标注，与 GAM 风格一致）
# =============================================================================
# 修正版：阈值模型绘图函数（使用选定的最佳模型，而非全局 AICc 最佳）
# 参数：
#   tax_data      : 当前分类单元的原始数据（包含 pH 和 Brith_rate）
#   tax_name      : 分类单元名称
#   phylum_name   : 门名称（用于标题）
#   best_model_name: 最佳模型名称（如 "chngpt_segmented"）
#   best_row      : 最佳模型所在行（必须包含 Threshold, Thr_LB, Thr_UB, R2_gaussian, p_value）
#   thr_fit       : 对应的阈值模型对象（res$models[[best_model_name]]）
#   grid_n        : 预测网格点数（默认 500）
# 返回值：ggplot 对象，并附加 optimum_x, left_slope, right_slope, left_r2, right_r2 属性
# =============================================================================
# =============================================================================
# 修正版：阈值模型绘图函数（绘制左右分段线性回归线，而非模型预测线）
# 参数：
#   tax_data      : 当前分类单元的原始数据（包含 pH 和 Brith_rate）
#   tax_name      : 分类单元名称
#   phylum_name   : 门名称（用于标题）
#   best_model_name: 最佳模型名称（如 "chngpt_segmented"）
#   best_row      : 最佳模型所在行（必须包含 Threshold, Thr_LB, Thr_UB, R2_gaussian, p_value）
#   thr_fit       : 对应的阈值模型对象（res$models[[best_model_name]]，仅用于提取阈值，不再用于绘图）
#   grid_n        : 预测网格点数（默认 500，实际未使用，保留兼容性）
# 返回值：ggplot 对象，并附加 optimum_x, left_slope, right_slope, left_r2, right_r2 属性
# =============================================================================
plot_threshold_model <- function(tax_data, tax_name, phylum_name, best_model_name, 
                                 best_row, thr_fit, grid_n = 500) {
  
  # 提取阈值和置信区间
  thr_value <- best_row$Threshold
  thr_lb    <- best_row$Thr_LB
  thr_ub    <- best_row$Thr_UB
  
  # 计算左右分段线性回归（基于阈值）
  segment_fits <- fit_two_segment_regression(
    data = tax_data,
    x_var = "pH",
    y_var = "Brith_rate",
    split_point = thr_value,
    min_points = 3,
    min_unique_x = 3
  )
  
  # 准备左右回归线的数据点（用于 geom_line）
  left_line <- NULL
  right_line <- NULL
  
  if (!is.null(segment_fits$left)) {
    x_left <- seq(min(tax_data$pH), thr_value, length.out = 100)
    y_left <- segment_fits$left$intercept + segment_fits$left$slope * x_left
    left_line <- data.frame(x = x_left, y = y_left)
  }
  
  if (!is.null(segment_fits$right)) {
    x_right <- seq(thr_value, max(tax_data$pH), length.out = 100)
    y_right <- segment_fits$right$intercept + segment_fits$right$slope * x_right
    right_line <- data.frame(x = x_right, y = y_right)
  }
  
  # 绘图
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = tax_data, 
                        ggplot2::aes(x = pH, y = Brith_rate),
                        color = "gray70", alpha = 0.8, size = 2) +
    ggplot2::geom_vline(xintercept = thr_value, 
                        color = "#4DBBD5FF", linetype = "dashed", linewidth = 0.8) +
    ggplot2::annotate("rect", xmin = thr_lb, xmax = thr_ub, ymin = -Inf, ymax = Inf,
                      fill = "#4DBBD5FF", alpha = 0.08)
  
  # 添加左侧回归线（如果存在）
  if (!is.null(left_line)) {
    p <- p + ggplot2::geom_line(data = left_line,
                                ggplot2::aes(x = x, y = y),
                                color = "#2E86AB", linewidth = 1.2, linetype = "solid")
  }
  
  # 添加右侧回归线（如果存在）
  if (!is.null(right_line)) {
    p <- p + ggplot2::geom_line(data = right_line,
                                ggplot2::aes(x = x, y = y),
                                color = "#A23B72", linewidth = 1.2, linetype = "solid")
  }
  
  # 添加标题和标签
  p <- p + ggplot2::labs(
    title = paste0(tax_name, " (", phylum_name, ")"),
    subtitle = paste("Best model:", best_model_name,
                     " | Optimum pH =", round(thr_value, 2)),
    x = "pH", y = "Brith_rate"
  ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"))
  
  # 添加整体 R² 和 p 值（左上角）
  overall_r2 <- best_row$R2_gaussian
  overall_p  <- best_row$p_value
  p <- p + ggplot2::annotate(
    "text", 
    x = min(tax_data$pH) + 0.1,
    y = max(tax_data$Brith_rate) * 0.95,
    label = paste0("Overall R² = ", round(overall_r2, 3),
                   "\np = ", format.pval(overall_p, digits = 3)),
    hjust = 0, vjust = 1, size = 3, color = "black"
  )
  
  # 添加左右分段标注（左下角、右下角）
  y_min <- min(tax_data$Brith_rate)
  y_range <- diff(range(tax_data$Brith_rate))
  
  if (!is.null(segment_fits$left)) {
    left_label <- paste0("Left: slope=", round(segment_fits$left$slope, 3),
                         ", R²=", round(segment_fits$left$r2, 3),
                         ", p=", format.pval(segment_fits$left$p_value, digits = 2))
    p <- p + ggplot2::annotate(
      "text", 
      x = min(tax_data$pH) + 0.1, 
      y = y_min + 0.1 * y_range,
      label = left_label,
      hjust = 0, vjust = 0, size = 3, color = "#2E86AB"
    )
  }
  
  if (!is.null(segment_fits$right)) {
    right_label <- paste0("Right: slope=", round(segment_fits$right$slope, 3),
                          ", R²=", round(segment_fits$right$r2, 3),
                          ", p=", format.pval(segment_fits$right$p_value, digits = 2))
    p <- p + ggplot2::annotate(
      "text", 
      x = max(tax_data$pH) - 0.1, 
      y = y_min + 0.1 * y_range,
      label = right_label,
      hjust = 1, vjust = 0, size = 3, color = "#A23B72"
    )
  }
  
  # 附加属性
  attr(p, "optimum_x")   <- thr_value
  attr(p, "left_slope")  <- if (!is.null(segment_fits$left)) segment_fits$left$slope else NA
  attr(p, "right_slope") <- if (!is.null(segment_fits$right)) segment_fits$right$slope else NA
  attr(p, "left_r2")     <- if (!is.null(segment_fits$left)) segment_fits$left$r2 else NA
  attr(p, "right_r2")    <- if (!is.null(segment_fits$right)) segment_fits$right$r2 else NA
  
  return(p)
}

# =============================================================================
# 5. 二次/GAM模型绘图函数（统一标注位置）
# 功能：绘制二次或GAM模型图，包含最优值竖线、左右线性段和标注，并添加居中副标题
# 参数：
#   tax_data      : 当前分类单元的原始数据
#   tax_name      : 分类单元名称
#   phylum_name   : 门名称（用于标题）
#   best_model_name: 最佳模型名称
#   best_row      : 最佳模型所在行（含 R² 和 p 值）
#   model_curve   : 模型曲线数据框（含 x 和 y 列）
#   optimum_x     : 最优值（顶点）
#   segment_fits  : 分段线性回归结果（来自 fit_two_segment_regression）
# 返回值：
#   ggplot 对象，并附加 optimum_x, left_slope, right_slope, left_r2, right_r2 作为属性
# =============================================================================
plot_quadgam_model <- function(tax_data, tax_name, phylum_name, best_model_name, 
                               best_row, model_curve, optimum_x, segment_fits) {
  
  overall_r2 <- best_row$R2_gaussian
  overall_p <- best_row$p_value
  left_slope <- ifelse(!is.null(segment_fits$left), segment_fits$left$slope, NA)
  right_slope <- ifelse(!is.null(segment_fits$right), segment_fits$right$slope, NA)
  left_r2 <- ifelse(!is.null(segment_fits$left), segment_fits$left$r2, NA)
  right_r2 <- ifelse(!is.null(segment_fits$right), segment_fits$right$r2, NA)
  
  p <- ggplot() +
    geom_point(data = tax_data, aes(x = pH, y = Brith_rate),
               color = "gray70", alpha = 0.8, size = 2) +
    geom_line(data = model_curve, aes(x = x, y = y),
              color = "black", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = optimum_x, color = "black",
               linetype = "dashed", linewidth = 1) +
    labs(title = paste0(tax_name, " (", phylum_name, ")"),
         subtitle = paste("Best model:", best_model_name,
                          " | Optimum pH =", round(optimum_x, 2)),
         x = "pH", y = "Brith_rate") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"))
  
  # 添加左右分段线
  if (!is.null(segment_fits$left)) {
    left_x <- seq(min(tax_data$pH), optimum_x, length.out = 100)
    left_y <- segment_fits$left$intercept + segment_fits$left$slope * left_x
    p <- p + geom_line(data = data.frame(x = left_x, y = left_y),
                       aes(x = x, y = y), color = "#2E86AB", linewidth = 1.2)
  }
  if (!is.null(segment_fits$right)) {
    right_x <- seq(optimum_x, max(tax_data$pH), length.out = 100)
    right_y <- segment_fits$right$intercept + segment_fits$right$slope * right_x
    p <- p + geom_line(data = data.frame(x = right_x, y = right_y),
                       aes(x = x, y = y), color = "#A23B72", linewidth = 1.2)
  }
  
  # 添加整体 R² 和 p 值（左上角）
  p <- p + annotate("text", x = min(tax_data$pH) + 0.1,
                    y = max(tax_data$Brith_rate) * 0.95,
                    label = paste0("Overall R² = ", round(overall_r2, 3),
                                   ", p = ", format.pval(overall_p, digits = 3)),
                    hjust = 0, vjust = 1, size = 3, color = "black")
  
  # 添加左右侧标注（左下角、右下角）
  y_min <- min(tax_data$Brith_rate)
  y_range <- diff(range(tax_data$Brith_rate))
  
  # 左侧标注
  if (!is.null(segment_fits$left)) {
    left_label <- paste0("Left: slope=", round(left_slope, 3),
                         ", R²=", round(left_r2, 3),
                         ", p=", format.pval(segment_fits$left$p_value, digits = 2))
    p <- p + annotate("text", 
                      x = min(tax_data$pH) + 0.1, 
                      y = y_min + 0.1 * y_range,
                      label = left_label,
                      hjust = 0, vjust = 0, size = 3, color = "#2E86AB")
  }
  
  # 右侧标注
  if (!is.null(segment_fits$right)) {
    right_label <- paste0("Right: slope=", round(right_slope, 3),
                          ", R²=", round(right_r2, 3),
                          ", p=", format.pval(segment_fits$right$p_value, digits = 2))
    p <- p + annotate("text", 
                      x = max(tax_data$pH) - 0.1, 
                      y = y_min + 0.1 * y_range,
                      label = right_label,
                      hjust = 1, vjust = 0, size = 3, color = "#A23B72")
  }
  
  # 将信息附加为属性
  attr(p, "optimum_x") <- optimum_x
  attr(p, "left_slope") <- left_slope
  attr(p, "right_slope") <- right_slope
  attr(p, "left_r2") <- left_r2
  attr(p, "right_r2") <- right_r2
  
  return(p)
}

# =============================================================================
# 6. 线性模型绘图函数
# 功能：绘制线性回归图，包含回归线和标注，并添加居中副标题（无最优pH）
# 参数：
#   tax_data      : 当前分类单元的原始数据
#   tax_name      : 分类单元名称
#   phylum_name   : 门名称（用于标题）
#   best_row      : 最佳模型所在行（含 R² 和 p 值）
#   slope         : 回归斜率
# 返回值：ggplot 对象
# =============================================================================
plot_linear_model <- function(tax_data, tax_name, phylum_name, best_row, slope) {
  
  overall_r2 <- best_row$R2_gaussian
  overall_p <- best_row$p_value
  
  p <- ggplot(tax_data, aes(x = pH, y = Brith_rate)) +
    geom_point(color = "gray70", alpha = 0.8, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    labs(title = paste0(tax_name, " (", phylum_name, ")"),
         subtitle = "Best model: Linear",
         x = "pH", y = "Brith_rate") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"))
  
  p <- p + annotate("text",
                    x = min(tax_data$pH) + 0.1,
                    y = max(tax_data$Brith_rate) * 0.95,
                    label = paste0("Slope = ", round(slope, 3),
                                   "\nR² = ", round(overall_r2, 3),
                                   ", p = ", format.pval(overall_p, digits = 3)),
                    hjust = 0, vjust = 1, size = 3, color = "black")
  
  return(p)
}

# =============================================================================
# 7. 为所有模型计算p值并更新模型表（高级函数）
# 参数：
#   res           : analyze_threshold_models 返回的对象
#   tax_data      : 当前分类单元的原始数据
#   tax_name      : 分类单元名称
#   mc_n          : 蒙特卡洛模拟次数（默认50000）
#   alpha_sig     : 显著性阈值（默认0.05）
# 返回值：
#   更新后的模型表（已添加 p_value 和 Significant 列），同时将原表存入全局列表 all_model_tables
#   注意：all_model_tables 需在主环境中定义（通过 <<- 赋值）
# =============================================================================
compute_model_pvalues <- function(res, tax_data, tax_name, mc_n = 50000, alpha_sig = 0.05) {
  model_pvals <- sapply(res$model_table$Model, function(mn) {
    mod <- res$models[[mn]]
    if (is.null(mod)) return(NA)
    get_model_pvalue(mod, tax_data, mn, x_col = "pH", y_col = "Brith_rate", mc_n = mc_n)
  })
  res$model_table$p_value <- model_pvals
  res$model_table$Significant <- ifelse(is.na(model_pvals), FALSE, model_pvals < alpha_sig)
  res$model_table$Tax_Name <- tax_name
  all_model_tables[[tax_name]] <<- res$model_table   # 存入全局收集列表
  return(res$model_table)
}

# =============================================================================
# 8. 根据筛选方式选择最佳模型（高级函数）
# 参数：
#   model_table     : 包含 p_value 和 Significant 列的模型表
#   selection_method: 筛选方式，可选 "p_then_AICc" 或 "AICc"
# 返回值：
#   列表，包含：
#     best_model_name : 最佳模型名称（若无则 "None"）
#     best_row        : 最佳模型所在行（数据框）
#     any_significant : 是否存在候选模型（布尔值）
#     candidate_models: 候选模型集（可能为空）
# =============================================================================
select_best_model_by_criteria <- function(model_table, selection_method = "p_then_AICc") {
  if (selection_method == "p_then_AICc") {
    candidate_models <- model_table[model_table$Significant == TRUE, , drop = FALSE]
  } else {
    candidate_models <- model_table
  }
  any_significant <- nrow(candidate_models) > 0
  
  if (any_significant) {
    # 按AICc排序（若AICc全NA则用AIC）
    if (all(is.na(candidate_models$AICc))) {
      candidate_models <- candidate_models[order(candidate_models$AIC), ]
    } else {
      candidate_models <- candidate_models[order(candidate_models$AICc), ]
    }
    best_model_name <- candidate_models$Model[1]
    best_row <- candidate_models[1, ]
  } else {
    best_model_name <- "None"
    best_row <- NULL
  }
  
  return(list(
    best_model_name = best_model_name,
    best_row = best_row,
    any_significant = any_significant,
    candidate_models = candidate_models
  ))
}

# =============================================================================
# 9. 根据最佳模型创建图形并提取统计信息（高级函数）
# 参数：
#   best_model_name : 最佳模型名称
#   best_row        : 最佳模型所在行
#   res             : analyze_threshold_models 返回的对象
#   tax_data        : 当前分类单元的原始数据
#   tax_name        : 分类单元名称
#   phylum_name     : 门名称
#   min_seg_points  : 分段线性每侧最少样本数（默认5）
#   min_seg_unique  : 分段线性每侧最少自变量唯一值个数（默认5）
# 返回值：
#   列表，包含：
#     p           : ggplot 对象
#     optimum_x   : 最优pH（若存在）
#     left_slope  : 左侧斜率
#     right_slope : 右侧斜率
#     left_r2     : 左侧R²
#     right_r2    : 右侧R²
#     sub_dir     : 子文件夹名称（用于保存）
#     model_type  : 模型类型标签
# =============================================================================
create_model_plot <- function(best_model_name, best_row, res, tax_data, tax_name, 
                              phylum_name, min_seg_points = 5, min_seg_unique = 5) {
  
  # 初始化返回值
  out <- list(
    p = NULL,
    optimum_x = NA,
    left_slope = NA,
    right_slope = NA,
    left_r2 = NA,
    right_r2 = NA,
    sub_dir = NA,
    model_type = NA
  )
  
  if (grepl("^chngpt_", best_model_name)) {
    # ---------- 阈值模型（使用新版 plot_threshold_model）----------
    # 获取对应的阈值模型对象
    thr_fit <- res$models[[best_model_name]]
    
    # 调用新版函数（不再需要 thr_rows 和 segment_fit）
    p <- plot_threshold_model(
      tax_data = tax_data,
      tax_name = tax_name,
      phylum_name = phylum_name,
      best_model_name = best_model_name,
      best_row = best_row,
      thr_fit = thr_fit,
      grid_n = 500
    )
    
    # 从图形属性中提取统计信息（新版函数已附加）
    out$optimum_x   <- attr(p, "optimum_x")
    out$left_slope  <- attr(p, "left_slope")
    out$right_slope <- attr(p, "right_slope")
    out$left_r2     <- attr(p, "left_r2")
    out$right_r2    <- attr(p, "right_r2")
    
    # 判断是否为单峰（左升右降）
    if (!is.na(out$left_slope) && !is.na(out$right_slope) && 
        out$left_slope > 0 && out$right_slope < 0) {
      out$sub_dir <- "Threshold_SinglePeak"
      out$model_type <- "Threshold_SinglePeak"
    } else {
      out$sub_dir <- "Threshold_Other"
      out$model_type <- "Threshold_Other"
    }
    out$p <- p
    
  } else if (best_model_name %in% c("Quadratic", "GAM")) {
    # ---------- 二次或GAM模型 ----------
    if (best_model_name == "Quadratic") {
      quad_model <- res$models$Quadratic
      coefs <- coef(quad_model)
      a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
      optimum_x <- -b / (2 * c)
      optimum_x <- max(min(tax_data$pH), min(max(tax_data$pH), optimum_x))
      x_seq <- seq(min(tax_data$pH), max(tax_data$pH), length.out = 200)
      y_seq <- a + b*x_seq + c*x_seq^2
      model_curve <- data.frame(x = x_seq, y = y_seq)
    } else { # GAM
      gam_model <- res$models$GAM
      x_seq <- seq(min(tax_data$pH), max(tax_data$pH), length.out = 500)
      y_pred <- predict(gam_model, newdata = data.frame(x = x_seq))
      max_idx <- which.max(y_pred)
      optimum_x <- x_seq[max_idx]
      model_curve <- data.frame(x = x_seq, y = y_pred)
    }
    
    segment_fits <- fit_two_segment_regression(
      data = tax_data,
      x_var = "pH",
      y_var = "Brith_rate",
      split_point = optimum_x,
      min_points = min_seg_points,
      min_unique_x = min_seg_unique
    )
    
    p <- plot_quadgam_model(tax_data, tax_name, phylum_name, best_model_name, best_row,
                            model_curve, optimum_x, segment_fits)
    out$optimum_x <- attr(p, "optimum_x")
    out$left_slope <- attr(p, "left_slope")
    out$right_slope <- attr(p, "right_slope")
    out$left_r2 <- attr(p, "left_r2")
    out$right_r2 <- attr(p, "right_r2")
    out$sub_dir <- "QuadraticGAM"
    out$model_type <- "Quadratic/GAM"
    out$p <- p
    
  } else if (best_model_name == "Linear") {
    # ---------- 线性模型 ----------
    lin_model <- res$models$Linear
    slope <- coef(lin_model)[2]
    p <- plot_linear_model(tax_data, tax_name, phylum_name, best_row, slope)
    out$p <- p
    out$sub_dir <- "Linear"
    out$model_type <- "Linear"
  }
  
  return(out)
}