rm(list=ls())

# =============================================================================
# 生理类群划分 + 按类群复制图形（不重新绘图）
# =============================================================================

library(tidyverse)
library(ggplot2)

# 1. 设置路径 ----------------------------------------------------------------
output_dir <- "./Fig"   # 请修改为实际路径
stats_dir <- file.path(output_dir, "Statistics")
fig_dir <- file.path(output_dir, "Figures")                # 原始图形所在目录
new_fig_dir <- file.path(output_dir, "Figures_by_Physiology1")  # 新文件夹，按类群存放
dir.create(new_fig_dir, recursive = TRUE, showWarnings = FALSE) #建文件夹

# 文件路径
res_list_file <- file.path(stats_dir, "res_list.rds")
stats_file <- file.path(stats_dir, "Detailed_Statistics.tsv")
models_file <- file.path(stats_dir, "All_Model_Comparisons.tsv")

# 2. 读取数据 ----------------------------------------------------------------
if (!file.exists(res_list_file)) stop("res_list.rds not found")
if (!file.exists(stats_file)) stop("Detailed_Statistics.tsv not found")
if (!file.exists(models_file)) stop("All_Model_Comparisons.tsv not found")

res_list <- readRDS(res_list_file)
stats <- read.table(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
all_models <- read.table(models_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 3. 重新计算活性pH范围（仅 Brith_rate > 0）---------------------------------
active_ranges <- tibble(Tax_Name = character(), pH_min = double(), pH_max = double())
for (tax_name in names(res_list)) {
  res <- res_list[[tax_name]]
  if (is.null(res)) next
  tax_data <- res$raw_data
  names(tax_data) <- c("pH", "Brith_rate")
  active <- tax_data %>% filter(Brith_rate > 0)
  if (nrow(active) > 0) {
    pH_min <- min(active$pH)
    pH_max <- max(active$pH)
  } else {
    pH_min <- NA; pH_max <- NA
  }
  active_ranges <- rbind(active_ranges, data.frame(Tax_Name = tax_name, pH_min = pH_min, pH_max = pH_max))
}

# 4. 合并信息，构建最终统计表 ------------------------------------------------
# 准备 all_models 中的关键列，重命名避免冲突
all_models_selected <- all_models %>%
  select(Tax_Name, Model, R2_gaussian, p_value, Slope_Pre, Slope_Post, Threshold) %>%
  rename(
    overall_r2 = R2_gaussian,
    overall_p = p_value,
    left_slope_new = Slope_Pre,
    right_slope_new = Slope_Post,
    thr = Threshold
  )

stats_joined <- stats %>%
  left_join(all_models_selected, by = c("Tax_Name" = "Tax_Name", "best_model" = "Model")) %>%
  mutate(
    overall_r2 = overall_r2,
    overall_p = overall_p,
    left_slope = ifelse(is.na(left_slope_new), left_slope, left_slope_new),
    right_slope = ifelse(is.na(right_slope_new), right_slope, right_slope_new),
    optimum_pH = ifelse(is.na(optimum_pH) & !is.na(thr), thr, optimum_pH)
  ) %>%
  select(-left_slope_new, -right_slope_new, -thr)

# 合并活性范围，并明确使用活性范围覆盖 pH_min 和 pH_max
stats_final <- stats_joined %>%
  left_join(active_ranges, by = "Tax_Name") %>%
  mutate(
    pH_min = ifelse(is.finite(pH_min.y), pH_min.y, NA_real_),
    pH_max = ifelse(is.finite(pH_max.y), pH_max.y, NA_real_),
    pH_range = pH_max - pH_min
  ) %>%
  select(-pH_min.x, -pH_max.x, -pH_min.y, -pH_max.y)

# 5. 为二次/GAM模型补充最适pH（如果缺失）--------------------------------------
get_optimum_quadratic_safe <- function(model, data) {
  if (is.null(model)) return(NA)
  if (!inherits(model, "lm")) return(NA)
  coefs <- tryCatch(coef(model), error = function(e) NULL)
  if (is.null(coefs) || length(coefs) < 3) return(NA)
  a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
  if (c >= 0) return(NA)  # 开口向上，无峰值
  vertex <- -b / (2 * c)
  if (vertex < min(data$pH) || vertex > max(data$pH)) return(NA)
  return(vertex)
}

get_optimum_gam_safe <- function(model, data, n = 500) {
  if (is.null(model)) return(NA)
  if (!inherits(model, "gam")) return(NA)
  x_seq <- tryCatch(seq(min(data$pH), max(data$pH), length.out = n), error = function(e) NULL)
  if (is.null(x_seq)) return(NA)
  y_pred <- tryCatch(predict(model, newdata = data.frame(x = x_seq)), error = function(e) NULL)
  if (is.null(y_pred)) return(NA)
  max_idx <- which.max(y_pred)
  return(x_seq[max_idx])
}

for (i in 1:nrow(stats_final)) {
  if (stats_final$model_type[i] == "Quadratic/GAM" && is.na(stats_final$optimum_pH[i])) {
    tax_name <- stats_final$Tax_Name[i]
    res <- res_list[[tax_name]]
    if (is.null(res)) next
    tax_data <- res$raw_data
    names(tax_data) <- c("pH", "Brith_rate")
    best_model <- stats_final$best_model[i]
    if (best_model == "Quadratic") {
      quad <- res$models$Quadratic
      if (!is.null(quad)) {
        opt <- get_optimum_quadratic_safe(quad, tax_data)
        if (!is.na(opt)) stats_final$optimum_pH[i] <- opt
      }
    } else if (best_model == "GAM") {
      gam <- res$models$GAM
      if (!is.null(gam)) {
        opt <- get_optimum_gam_safe(gam, tax_data)
        if (!is.na(opt)) stats_final$optimum_pH[i] <- opt
      }
    }
  }
}

# 5.1 为线性模型补充斜率（因为合并时可能丢失）-----------------------------------
for (i in 1:nrow(stats_final)) {
  if (stats_final$best_model[i] == "Linear" && is.na(stats_final$left_slope[i])) {
    tax_name <- stats_final$Tax_Name[i]
    res <- res_list[[tax_name]]
    if (!is.null(res) && !is.null(res$models$Linear)) {
      lm_model <- res$models$Linear
      slope <- coef(lm_model)[2]
      stats_final$left_slope[i] <- slope
      stats_final$right_slope[i] <- slope
      cat("已补充线性斜率:", tax_name, "斜率 =", round(slope, 4), "\n")
    }
  }
}

source("./data/threshold_analysis_helpers.R") 

# 6. 生理类群划分（先范围后模型特征，不考虑p值）-----------------------------------------------
acid_low <- 5.5
alk_high <- 8.0
min_slope_abs <- 1e-6

stats_final$physiological_group <- NA_character_
stats_final$classification_notes <- NA_character_

for (i in 1:nrow(stats_final)) {
  row <- stats_final[i, ]
  
  # 无效 pH 范围 → pH-insensitive
  if (!is.finite(row$pH_min) || !is.finite(row$pH_max)) {
    stats_final$physiological_group[i] <- "pH-insensitive"
    stats_final$classification_notes[i] <- "No active data (Brith_rate>0)"
    next
  }
  
  # 无显著模型 → pH-insensitive
  if (row$model_type == "NonSignificant") {
    stats_final$physiological_group[i] <- "pH-insensitive"
    stats_final$classification_notes[i] <- "No significant model"
    next
  }
  
  pH_min <- row$pH_min
  pH_max <- row$pH_max
  
  # ===== 第一步：根据活性范围判断专性 =====
  if (pH_max <= acid_low) {
    stats_final$physiological_group[i] <- "Obligate acidophile"
    stats_final$classification_notes[i] <- paste0("活性范围完全≤5.5 (", pH_min, "-", pH_max, ")")
    next
  }
  if (pH_min >= alk_high) {
    stats_final$physiological_group[i] <- "Obligate alkaliphile"
    stats_final$classification_notes[i] <- paste0("活性范围完全≥8.0 (", pH_min, "-", pH_max, ")")
    next
  }
  if (pH_min >= acid_low && pH_max <= alk_high) {
    stats_final$physiological_group[i] <- "Obligate neutrophile"
    stats_final$classification_notes[i] <- paste0("活性范围完全在5.5-8.0 (", pH_min, "-", pH_max, ")")
    next
  }
  
  # ===== 第二步：范围跨越边界，根据模型特征判断兼性 =====
  
  # ----- 2.1 单峰阈值或二次/GAM有最适pH -----
  if ((row$model_type == "Threshold_SinglePeak" || row$model_type == "Quadratic/GAM") && !is.na(row$optimum_pH)) {
    opt <- row$optimum_pH
    if (opt < acid_low) {
      grp <- "Facultative acidophile"
      note <- paste0("最适pH=", round(opt,2), ", 活性范围跨越5.5")
    } else if (opt > alk_high) {
      grp <- "Facultative alkaliphile"
      note <- paste0("最适pH=", round(opt,2), ", 活性范围跨越8.0")
    } else {
      grp <- "Facultative neutrophile"
      note <- paste0("最适pH=", round(opt,2), ", 活性范围超出5.5-8.0")
    }
    stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
    next
  }
  
  # ----- 2.2 线性模型（仅根据斜率正负）-----
  if (row$model_type == "Linear" && !is.na(row$left_slope)) {
    slope <- row$left_slope
    if (abs(slope) < min_slope_abs) {
      grp <- "pH-insensitive"
      note <- paste0("线性斜率≈0 (", round(slope,6), "), 活性范围跨区间")
    } else if (slope > 0) {
      grp <- "Facultative alkaliphile"
      note <- paste0("线性正斜率 (", round(slope,3), "), 活性范围跨区间")
    } else {
      grp <- "Facultative acidophile"
      note <- paste0("线性负斜率 (", round(slope,3), "), 活性范围跨区间")
    }
    stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
    next
  }
  
  # ----- 2.3 阈值Other模型（根据斜率组合、跳变及显著侧）-----
  if (row$model_type == "Threshold_Other" && !is.na(row$left_slope) && !is.na(row$right_slope)) {
    left <- row$left_slope
    right <- row$right_slope
    # 极小斜率归零
    if (abs(left) < min_slope_abs) left <- 0
    if (abs(right) < min_slope_abs) right <- 0
    
    thr <- row$optimum_pH  # 阈值位置（仅用于左正右负或左负右正时作为极值点）
    
    # 获取左右截距（用于计算阈值处的预测值）
    tax_name <- row$Tax_Name
    res_tmp <- res_list[[tax_name]]
    left_pred <- NA
    right_pred <- NA
    if (!is.null(res_tmp)) {
      tax_data <- res_tmp$raw_data
      names(tax_data) <- c("pH", "Brith_rate")
      seg_fits <- fit_two_segment_regression(tax_data, "pH", "Brith_rate", thr, min_points = 3, min_unique_x = 3)
      if (!is.null(seg_fits$left)) {
        left_pred <- seg_fits$left$intercept + left * thr
      }
      if (!is.null(seg_fits$right)) {
        right_pred <- seg_fits$right$intercept + right * thr
      }
    }
    
    # 情况1：左正右负（先增后减，有峰值）→ 视为单峰，使用阈值作为最适pH
    if (left > 0 && right < 0) {
      if (!is.na(thr)) {
        if (thr < acid_low) {
          grp <- "Facultative acidophile"
          note <- paste0("阈值Other, 左正右负(峰值), 最适pH=", round(thr,2), " 在酸性区")
        } else if (thr > alk_high) {
          grp <- "Facultative alkaliphile"
          note <- paste0("阈值Other, 左正右负(峰值), 最适pH=", round(thr,2), " 在碱性区")
        } else {
          grp <- "Facultative neutrophile"
          note <- paste0("阈值Other, 左正右负(峰值), 最适pH=", round(thr,2), " 在中性区")
        }
      } else {
        grp <- "Others"
        note <- "阈值Other, 左正右负但无阈值位置"
      }
      stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
      next
    }
    
    # 情况2：左负右正（先减后增，有谷底）→ 直接归为 Others
    if (left < 0 && right > 0) {
      grp <- "Others"
      note <- "阈值Other, 左负右正(谷底), 无明确酸碱偏好"
      stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
      next
    }
    
    # 情况3：左正右正（双正斜率）或 左0右正 / 左正右0
    if ((left > 0 && right > 0) || (left == 0 && right > 0) || (left > 0 && right == 0)) {
      if (!is.na(left_pred) && !is.na(right_pred) && left_pred > right_pred) {
        grp <- "Others"
        note <- paste0("阈值Other, 双正斜率但有向下的跳变 (左侧预测值", round(left_pred,3), " > 右侧", round(right_pred,3), ")")
      } else {
        grp <- "Facultative alkaliphile"
        note <- "阈值Other, 整体递增趋势"
      }
      stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
      next
    }
    
    # 情况4：左负右负（双负斜率）或 左0右负 / 左负右0
    if ((left < 0 && right < 0) || (left == 0 && right < 0) || (left < 0 && right == 0)) {
      if (!is.na(left_pred) && !is.na(right_pred) && left_pred < right_pred) {
        grp <- "Others"
        note <- paste0("阈值Other, 双负斜率但有向上的跳变 (左侧预测值", round(left_pred,3), " < 右侧", round(right_pred,3), ")")
      } else {
        grp <- "Facultative acidophile"
        note <- "阈值Other, 整体递减趋势"
      }
      stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
      next
    }
    
    # 情况5：两侧均为0（显著但平缓）→ 归为 Others（pH-insensitive已用于非显著模型）
    if (left == 0 && right == 0) {
      grp <- "Others"
      note <- "阈值Other, 模型显著但两侧斜率均为0"
      stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
      next
    }
    
    # 情况6：其他组合 → Others
    grp <- "Others"
    note <- "阈值Other, 其他斜率组合"
    stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
    next
  }
  
  # ----- 2.4 二次/GAM无最适pH（根据整体斜率）-----
  if (row$model_type == "Quadratic/GAM" && is.na(row$optimum_pH)) {
    tax_name <- row$Tax_Name
    res_tmp <- res_list[[tax_name]]
    slope <- NA
    if (!is.null(res_tmp)) {
      tax_data <- res_tmp$raw_data
      names(tax_data) <- c("pH", "Brith_rate")
      if (nrow(tax_data) >= 3) {
        lm_overall <- lm(Brith_rate ~ pH, data = tax_data)
        slope <- coef(lm_overall)[2]
      }
    }
    if (is.na(slope) || abs(slope) < min_slope_abs) {
      grp <- "pH-insensitive"
      note <- "二次/GAM无峰值, 整体斜率≈0"
    } else if (slope > 0) {
      grp <- "Facultative alkaliphile"
      note <- "二次/GAM无峰值, 整体上升"
    } else {
      grp <- "Facultative acidophile"
      note <- "二次/GAM无峰值, 整体下降"
    }
    stats_final[i, c("physiological_group", "classification_notes")] <- c(grp, note)
    next
  }
  
  # ----- 2.5 其他未覆盖情况 -----
  stats_final[i, c("physiological_group", "classification_notes")] <- c("Others", "Significant but not classified")
}

# 7. 保存带生理类群的统计表 ----------------------------------------------------
out_stats <- file.path(stats_dir, "Detailed_Statistics_with_Physiology.tsv")
write.table(stats_final, out_stats, sep = "\t", row.names = FALSE, quote = FALSE)

# 8. 打印各类群数量 ----------------------------------------------------------
cat("\n生理类群统计:\n")
print(table(stats_final$physiological_group))

# =============================================================================
# Manual refinement of physiological pH classifications
#
# Following the initial model-based classification, growth-response curves
# were visually inspected for each taxon. Taxa with inaccurate initial
# assignments were manually reclassified based on their growth patterns.
# The manually refined classifications were saved in
# "Detailed_Statistics_with_Physiology_adjust.tsv".
#
# Organize growth-response plots according to the refined classifications.
# =============================================================================

# =============================================================================
# 按手动调整的生理类群复制图形
# =============================================================================
rm(list = ls())
library(tidyverse)

# ---------------------------- 1. 设置路径 ------------------------------------
# 源图形目录（使用您重新生成后的图形目录）
source_fig_dir <- "./Fig/Figures_OriginalSlope"
# 目标根目录（按生理类群存放）
target_base_dir <- "./Fig/Figures_by_Physiology_Adjusted"
# 统计表文件（包含手动调整后的生理类群列）
stats_file <- "./Fig/Statistics/Detailed_Statistics_with_Physiology_adjust.tsv"
# 注意：该统计表应已包含 Tax_Name, best_model, physiological_group_adjust 列

# 如果您的统计表是其他文件，请修改
if (!file.exists(stats_file)) {
  stop(paste("统计表不存在:", stats_file))
}

# ---------------------------- 2. 读取统计表 ----------------------------------
stats <- read.table(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8")

# 检查必要列
required_cols <- c("Tax_Name", "best_model", "physiological_group_adjust")
missing <- required_cols[!required_cols %in% names(stats)]
if (length(missing) > 0) {
  stop(paste("统计表中缺少列:", paste(missing, collapse = ", ")))
}

# 只保留有有效类群的行（非空、非NA）
stats_valid <- stats %>%
  filter(!is.na(physiological_group_adjust), 
         physiological_group_adjust != "", 
         best_model != "None")

cat("共有", nrow(stats_valid), "个分类单元有有效的生理类群\n")

# ---------------------------- 3. 列出源目录中所有 PDF 文件 -------------------
if (!dir.exists(source_fig_dir)) {
  stop(paste("源图形目录不存在:", source_fig_dir))
}
all_pdf <- list.files(source_fig_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
all_pdf_names <- basename(all_pdf)
cat("源目录中共有", length(all_pdf), "个 PDF 文件\n")

# ---------------------------- 4. 创建目标子文件夹 ----------------------------
groups <- unique(stats_valid$physiological_group_adjust)
for (g in groups) {
  dir.create(file.path(target_base_dir, g), recursive = TRUE, showWarnings = FALSE)
}
cat("创建了", length(groups), "个类群文件夹\n")

# ---------------------------- 5. 匹配并复制 ----------------------------------
success_count <- 0
fail_list <- c()

for (i in 1:nrow(stats_valid)) {
  tax_name <- stats_valid$Tax_Name[i]
  best_model <- stats_valid$best_model[i]
  group <- stats_valid$physiological_group_adjust[i]
  
  # 构建期望的文件名（根据您生成图形时的命名规则）
  # 注意：您的图形文件名格式是 "OTU1_chngpt_stegmented.pdf" 等（Tax_Name + "_" + best_model + ".pdf"）
  expected_name <- paste0(tax_name, "_", best_model, ".pdf")
  
  # 精确匹配
  matched <- all_pdf[all_pdf_names == expected_name]
  
  # 如果精确匹配失败，尝试模糊匹配（以防文件名有额外字符）
  if (length(matched) == 0) {
    pattern <- paste0("^", tax_name, ".*", best_model, "\\.pdf$")
    matched <- all_pdf[grepl(pattern, all_pdf_names, ignore.case = TRUE)]
    if (length(matched) == 0) {
      # 最后尝试只匹配分类名（有些图形可能单独命名）
      pattern_tax <- paste0("^", tax_name, ".*\\.pdf$")
      matched <- all_pdf[grepl(pattern_tax, all_pdf_names, ignore.case = TRUE)]
      if (length(matched) > 1) {
        # 多个匹配时，优先选择包含 best_model 的
        matched2 <- matched[grepl(best_model, basename(matched), ignore.case = TRUE)]
        if (length(matched2) > 0) matched <- matched2
        else matched <- matched[1]  # 取第一个
      }
    }
  }
  
  if (length(matched) == 0) {
    fail_list <- c(fail_list, tax_name)
    next
  }
  
  src_file <- matched[1]
  dest_file <- file.path(target_base_dir, group, expected_name)
  
  # 复制（覆盖已存在的文件）
  file.copy(src_file, dest_file, overwrite = TRUE)
  success_count <- success_count + 1
  
  if (i %% 50 == 0) cat("已处理", i, "个分类...\n")
}

# ---------------------------- 6. 报告结果 ------------------------------------
cat("\n复制完成！\n")
cat("成功复制:", success_count, "/", nrow(stats_valid), "\n")
if (length(fail_list) > 0) {
  cat("未找到对应图形文件的分类 (", length(fail_list), "个):\n")
  print(head(fail_list, 20))
  # 可选：将失败的列表保存到文件
  writeLines(fail_list, file.path(target_base_dir, "missing_files.txt"))
} else {
  cat("所有分类均已成功复制！\n")
}
