# =============================================================================
# 步骤0：准备工作（只需运行一次）
# =============================================================================
rm(list = ls())
library(dplyr)
library(tidyr)

# 设置参数（全部可调）
# ------------------- 数据筛选条件 -------------------
min_unique_pH <- 1      # 每个分类单元至少需要多少个不同的pH值 
min_samples <- 1        # 每个分类单元至少需要的样本总数
tax_level <- "tax_id"    # 分类水平，可选：Phylum, Class, Order, Family, Genus
number_otus <- 1              # 至少包含多少个不同的OTU，OTU水平

# 输出目录
output_dir <- "./Fig/Threshold_Results_tax_id"
stats_dir <- paste0(output_dir, "/Statistics")
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

# ========== 每个科至少包含几个不同的OTU ==========
otu_count_per_family <- df_clean %>%
  group_by(!!sym(tax_level)) %>%
  summarise(n_otus = n_distinct(tax_id), .groups = 'drop')

families_with_enough_otus <- otu_count_per_family %>%
  filter(n_otus >= number_otus) %>%   
  pull(!!sym(tax_level))

df_clean <- df_clean %>%
  filter(!!sym(tax_level) %in% families_with_enough_otus)

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

# 筛选符合条件的分类单元
#当用全部数据（包括非活性）时，这里的筛选不起作用
selected_taxa <- tax_stats %>%
  filter(n_unique_pH >= min_unique_pH,
         n_total >= min_samples,
         pH_range >= 0) %>%
  pull(.data[[tax_level]])

cat("\n符合筛选条件的", tax_level, "数量:", length(selected_taxa), "\n")

tax_stats <- tax_stats[which(tax_stats$tax_id %in% selected_taxa),]

write.table(tax_stats, file.path(stats_dir, "all_taxa_filtered.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)  

# 1. 文件路径
stats_file_Physiology <- file.path(stats_dir, "Detailed_Statistics_with_Physiology_adjust.tsv")
# 2. 读取数据 ----------------------------------------------------------------
if (!file.exists(stats_file_Physiology)) stop("Detailed_Statistics_with_Physiology_adjust.tsv not found")
stats_Physiology <- read.table(stats_file_Physiology, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")

# 保留筛选后的数据
stats_Physiology <- stats_Physiology %>% filter(Tax_Name %in% selected_taxa)

# 找出 selected_taxa 中不在 stats_Physiology 里的部分
missing_taxa <- setdiff(selected_taxa, stats_Physiology$Tax_Name) #应该是只有两个点，做不了回归等

print(stats_Physiology)
write.table(stats_Physiology, file.path(stats_dir, "Detailed_Statistics_with_Physiology_adjust_filtered.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)  

# ============================================================================
# 生理类群过滤：基于活性样本比例的验证
# ============================================================================
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

# ---------------------------- 1. 路径设置 ------------------------------------
#output_dir <- "./Fig/非活性调0/Threshold_Results_Family"
output_dir <- "./Fig/Threshold_Results_tax_id"
stats_dir <- file.path(output_dir, "Statistics")
fig_dir <- file.path(output_dir, "Figures")

stats_file <- file.path(stats_dir, "Detailed_Statistics_with_Physiology_adjust_filtered.tsv")
res_list_file <- file.path(stats_dir, "res_list.rds")

if (!file.exists(stats_file)) stop("找不到 Detailed_Statistics_with_Physiology_adjust_filtered.tsv")
if (!file.exists(res_list_file)) stop("找不到 res_list.rds")

# ---------------------------- 2. 读取数据 ------------------------------------
stats_final <- read.table(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,encoding = "UTF-8")  # 或 "GBK"

res_list <- readRDS(res_list_file)

# 定义 pH 区间阈值
acid_high <- 5.5
alk_low <- 8.0

# ========== 可调节参数（修改这里即可全局调整）==========
min_prop_obligate <- 0.5      # 专性类群的主区间最小活性比例
min_prop_facultative <- 0.5   # 兼性类群的主区间最小活性比例
min_prop_insensitive <- 0.5   # pH-insensitive 类群整体最小活性比例
zero_tolerance <- 1e-6        # 判断“无活性”的容忍值（严格为0时设1e-6）
# =====================================================

# # 定义过滤条件（使用上面的参数）
conditions <- list(
  "Obligate acidophile" = list(
    check = function(acid, neut, alk, overall) acid >= min_prop_obligate && (neut + alk) < zero_tolerance,
    note = paste0("酸性区活性≥", min_prop_obligate, "，中性+碱性区无活性")
  ),
  "Facultative acidophile" = list(
    check = function(acid, neut, alk, overall) acid >= min_prop_facultative && acid > neut,
    note = paste0("酸性区活性≥", min_prop_facultative, "且大于中性区")
  ),
  "Obligate neutrophile" = list(
    check = function(acid, neut, alk, overall) neut >= min_prop_obligate && (acid + alk) < zero_tolerance,
    note = paste0("中性区活性≥", min_prop_obligate, "，酸+碱区无活性")
  ),
  "Facultative neutrophile" = list(
    check = function(acid, neut, alk, overall) neut >= min_prop_facultative && (acid > zero_tolerance | alk > zero_tolerance),
    note = paste0("中性区活性≥", min_prop_facultative, "，且酸或碱区有活性")
  ),
  "Obligate alkaliphile" = list(
    check = function(acid, neut, alk, overall) alk >= min_prop_obligate && (acid + neut) < zero_tolerance,
    note = paste0("碱性区活性≥", min_prop_obligate, "，酸+中性区无活性")
  ),
  "Facultative alkaliphile" = list(
    check = function(acid, neut, alk, overall) alk >= min_prop_facultative && alk > neut,
    note = paste0("碱性区活性≥", min_prop_facultative, "且大于中性区")
  ),
  "pH-insensitive" = list(
    check = function(acid, neut, alk, overall) overall >= min_prop_insensitive,
    note = paste0("整体活性比例≥", min_prop_insensitive)
  ),
  "Others" = list(
    check = function(acid, neut, alk, overall) TRUE,
    note = "不限制"
  )
)

# ---------------------------- 过滤计算----------------------------
filter_results <- data.frame()
for (i in 1:nrow(stats_final)) {
  tax_name <- stats_final$Tax_Name[i]
  phys_group <- stats_final$physiological_group_adjust[i]
  if (is.na(phys_group) || phys_group == "") next
  
  res <- res_list[[tax_name]]
  if (is.null(res)) {
    warning("res_list 中缺失: ", tax_name)
    next
  }
  raw <- res$raw_data
  names(raw) <- c("pH", "Brith_rate")
  
  raw$active <- raw$Brith_rate > 0
  overall_active_prop <- mean(raw$active)
  raw$zone <- cut(raw$pH, breaks = c(-Inf, acid_high, alk_low, Inf),
                  labels = c("acid", "neutral", "alkaline"), right = FALSE)
  
  zone_summary <- raw %>%
    group_by(zone) %>%
    summarise(total = n(), active = sum(active), .groups = "drop") %>%
    mutate(prop = active / total)
  
  acid_prop <- ifelse("acid" %in% zone_summary$zone, zone_summary$prop[zone_summary$zone == "acid"], 0)
  neut_prop <- ifelse("neutral" %in% zone_summary$zone, zone_summary$prop[zone_summary$zone == "neutral"], 0)
  alk_prop  <- ifelse("alkaline" %in% zone_summary$zone, zone_summary$prop[zone_summary$zone == "alkaline"], 0)
  
  if (phys_group %in% names(conditions)) {
    cond <- conditions[[phys_group]]
    pass <- cond$check(acid_prop, neut_prop, alk_prop, overall_active_prop)
    note <- ifelse(pass, "通过", paste0("不通过: ", cond$note))
  } else {
    pass <- FALSE
    note <- "未知生理类群"
  }
  
  filter_results <- rbind(filter_results, data.frame(
    Tax_Name = tax_name,
    Physiological_Group = phys_group,
    Acid_Prop = acid_prop,
    Neutral_Prop = neut_prop,
    Alkaline_Prop = alk_prop,
    Overall_Prop = overall_active_prop,  
    Pass = pass,
    Filter_Note = note,
    stringsAsFactors = FALSE
  ))
}

# 查看通过情况
table(filter_results$Pass)
table(filter_results$Physiological_Group, filter_results$Pass)

# 保存结果
out_filter <- file.path(stats_dir, "Physiology_Filter_Results.tsv")
write.table(filter_results, out_filter, sep = "\t", row.names = FALSE, fileEncoding = "GBK")
cat("过滤结果已保存至:", out_filter, "\n")

# ---------------------------- 5. 可视化：各区间活性比例分面图 ----------------
# 准备长格式数据用于绘图
plot_data <- filter_results %>%
  select(Tax_Name, Physiological_Group, Acid_Prop, Neutral_Prop, Alkaline_Prop) %>%
  pivot_longer(cols = c(Acid_Prop, Neutral_Prop, Alkaline_Prop),
               names_to = "Zone", values_to = "Proportion") %>%
  mutate(Zone = factor(Zone, levels = c("Acid_Prop", "Neutral_Prop", "Alkaline_Prop"),
                       labels = c("Acid (<5.5)", "Neutral (5.5-8.0)", "Alkaline (≥8.0)")))

# 按生理类群分面，绘制箱线图或点图
p <- ggplot(plot_data, aes(x = Zone, y = Proportion, fill = Zone)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
  facet_wrap(~ Physiological_Group, scales = "free_y") +
  labs(title = "活性样本比例分布按生理类群",
       x = "pH 区间", y = "活性样本比例") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
p

out_plot <- file.path(fig_dir, "Active_Proportion_by_Physiology.pdf")
ggsave(out_plot, p, width = 12, height = 8, device = "pdf")
cat("可视化图已保存至:", out_plot, "\n")

# ---------------------------- 6. 列出不通过的分类（供手动复查） --------------
failed_list <- filter_results %>% filter(!Pass) %>% pull(Tax_Name)
if (length(failed_list) > 0) {
  cat("\n以下", length(failed_list), "个分类单元未通过过滤条件，建议复查：\n")
  print(failed_list)
} else {
  cat("\n所有分类均通过过滤条件！\n")
}

library(ggplot2)
library(dplyr)

# 仅保留通过筛选的OTU
dat_plot <- filter_results %>%
  filter(Pass %in% c(TRUE, "TRUE", "True"))

# 固定类别顺序，与Fig. 1的生理分类顺序保持一致
group_order <- c(
  "Obligate acidophile",
  "Facultative acidophile",
  "Obligate neutrophile",
  "Facultative neutrophile",
  "Obligate alkaliphile",
  "Facultative alkaliphile",
  "Others",
  "pH-insensitive"
)

# 图中显示名称
group_labels <- c(
  "Obligate acidophile"     = "Obligate acidophile",
  "Facultative acidophile"  = "Facultative acidophile",
  "Obligate neutrophile"    = "Obligate neutrophile",
  "Facultative neutrophile" = "Facultative neutrophile",
  "Obligate alkaliphile"    = "Obligate alkaliphile",
  "Facultative alkaliphile" = "Facultative alkaliphile",
  "Others"                  = "Other responders",
  "pH-insensitive"          = "pH-insensitive"
)

# 汇总数量和比例
plot_data <- dat_plot %>%
  count(Physiological_Group, name = "Number") %>%
  mutate(
    Percentage = Number / sum(Number) * 100,
    Physiological_Group = factor(
      Physiological_Group,
      levels = rev(group_order)
    ),
    Label = sprintf("%d  (%.1f%%)", Number, Percentage)
  )

# 查看汇总结果
print(plot_data)
stopifnot(sum(plot_data$Number) == 2589)

# 配色
# 同一pH偏好的专性和兼性类型使用同一色系
group_colours  <- c(
  "Obligate acidophile"      = "#D62728",   # 深红
  "Facultative acidophile"   = "#FF7F0E",   # 橙
  "Obligate neutrophile"     = "#3A8C6E",   # 中绿
  "Facultative neutrophile"  = "#8DCC6B",   # 黄绿
  "Obligate alkaliphile"     = "#1F77B4",   # 深蓝
  "Facultative alkaliphile"  = "#17BECF",   # 天蓝
  "pH-insensitive"           = "#7F7F7F",    # 中灰
  "Other responders"         = "#FFFFFF"
)

# 绘图
p <- ggplot(
  plot_data,
  aes(x = Number, y = Physiological_Group, fill = Physiological_Group)
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = Label),
    hjust = -0.08,
    size = 3.5
  ) +
  scale_fill_manual(
    values = group_colours,
    labels = group_labels
  ) +
  scale_y_discrete(labels = group_labels) +
  scale_x_continuous(
    limits = c(0, max(plot_data$Number) * 1.28),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Number of OTUs",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(
      colour = "black",
      size = 9
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 9
    ),
    axis.title.x = element_text(
      colour = "black",
      size = 10,
      margin = margin(t = 7)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(8, 25, 8, 8)
  )

p

# 保存为矢量PDF
ggsave(
  "Fig/Fig_S1_pH_strategy_distribution.pdf",
  plot = p,
  width = 170,
  height = 105,
  units = "mm",
  device = cairo_pdf
)

library(ggplot2)
library(dplyr)

# 筛选通过的OTU，并统一类别名称
dat_plot <- filter_results %>%
  filter(Pass %in% c(TRUE, "TRUE", "True")) %>%
  mutate(
    Physiological_Group = recode(
      Physiological_Group,
      "Other" = "Others",
      "Other responders" = "Others"
    )
  )

group_order <- c(
  "Obligate acidophile",
  "Facultative acidophile",
  "Obligate neutrophile",
  "Facultative neutrophile",
  "Obligate alkaliphile",
  "Facultative alkaliphile",
  "Others",
  "pH-insensitive"
)

group_labels <- c(
  "Obligate acidophile"     = "Obligate acidophile",
  "Facultative acidophile"  = "Facultative acidophile",
  "Obligate neutrophile"    = "Obligate neutrophile",
  "Facultative neutrophile" = "Facultative neutrophile",
  "Obligate alkaliphile"    = "Obligate alkaliphile",
  "Facultative alkaliphile" = "Facultative alkaliphile",
  "Others"                  = "Other responders",
  "pH-insensitive"          = "pH-insensitive"
)

group_colours <- c(
  "Obligate acidophile"     = "#D62728",
  "Facultative acidophile"  = "#FF7F0E",
  "Obligate neutrophile"    = "#3A8C6E",
  "Facultative neutrophile" = "#8DCC6B",
  "Obligate alkaliphile"    = "#1F77B4",
  "Facultative alkaliphile" = "#17BECF",
  "Others"                  = "#000000",
  "pH-insensitive"          = "#BDBDBD"
)

# 汇总数据
# 不再手动计算ypos
plot_data <- dat_plot %>%
  count(Physiological_Group, name = "Number") %>%
  mutate(
    Physiological_Group = factor(
      Physiological_Group,
      levels = group_order
    )
  ) %>%
  arrange(Physiological_Group) %>%
  mutate(
    Percentage = Number / sum(Number) * 100,
    Label = sprintf("%.1f%%", Percentage),
    
    # 黑色扇区使用白色文字
    Text_colour = ifelse(
      Physiological_Group == "pH-insensitive",
      "white",
      "black"
    )
  )

# 检查数据与顺序
print(
  plot_data %>%
    select(Physiological_Group, Number, Percentage, Label)
)

stopifnot(sum(plot_data$Number) == 2589)

# 绘制甜甜圈图
p_donut <- ggplot(
  plot_data,
  aes(
    x = 2,
    y = Number,
    fill = Physiological_Group,
    group = Physiological_Group
  )
) +
  geom_col(
    width = 1,
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = Label,
      colour = Text_colour,
      group = Physiological_Group
    ),
    position = position_stack(vjust = 0.5),
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  annotate(
    "text",
    x = 0.5,
    y = 0,
    label = "n = 2,589\nOTUs",
    size = 4
  ) +
  scale_fill_manual(
    values = group_colours,
    breaks = group_order,
    labels = group_labels[group_order],
    drop = FALSE
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

p_donut

ggsave(
  "Fig/Fig_S1_donut_chart.pdf",
  plot = p_donut,
  width = 180,
  height = 120,
  units = "mm",
  device = cairo_pdf
)


# ============================================================================
# 根据过滤结果重新组织图形：Pass/Fail -> Physiological_Group
# ============================================================================
rm(list = ls())

library(dplyr)
library(stringr)

# ---------------------------- 1. 路径设置 ------------------------------------
output_dir <- "./Fig"
#output_dir <- "./Fig/非活性调0/Threshold_Results_Family"
stats_dir <- file.path(output_dir, "Statistics")
fig_by_phys_dir <- file.path(output_dir, "Figures_by_Physiology")        # 原始按类群分好的图形
target_base_dir <- file.path(output_dir, "Figures_by_Physiology_Filtered") # 新根目录

# 输入文件：过滤结果表（您已经生成的 stats_final）
stats_file <- file.path(stats_dir, "Physiology_Filter_Results.tsv") # 请确认文件名
if (!file.exists(stats_file)) {
  # 如果文件不存在，尝试直接使用内存中的 stats_final（如果 R 会话未关闭）
  if (!exists("stats_final")) stop("请提供 stats_final 数据框或正确的文件路径")
} else {
  stats_final <- read.table(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")
}

# ---------------------------- 2. 创建目标文件夹结构 --------------------------
# 为每个 Pass 状态和每个生理类群创建子文件夹
pass_status <- c("Pass", "Fail")
phys_groups <- unique(stats_final$Physiological_Group)
phys_groups <- phys_groups[!is.na(phys_groups) & phys_groups != ""]

for (status in pass_status) {
  for (grp in phys_groups) {
    dir.create(file.path(target_base_dir, status, grp), recursive = TRUE, showWarnings = FALSE)
  }
}

# ---------------------------- 3. 复制文件 ------------------------------------
success <- 0
fail <- c()

for (i in 1:nrow(stats_final)) {
  tax <- stats_final$Tax_Name[i]
  grp <- stats_final$Physiological_Group[i]
  pass <- stats_final$Pass[i]
  status <- ifelse(pass, "Pass", "Fail")
  
  # 在原始目录中搜索（递归，不限定子目录）
  # 方法1：精确匹配文件名（可能包含下划线或点）
  pattern <- paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", tax), "(_|\\.)")
  all_files <- list.files(fig_by_phys_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  
  # 方法2：如果没找到，使用包含匹配（更宽松）
  if (length(all_files) == 0) {
    all_files <- list.files(fig_by_phys_dir, pattern = tax, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }
  
  if (length(all_files) == 0) {
    fail <- c(fail, tax)
    next
  }
  
  src_file <- all_files[1]  # 取第一个匹配
  dest_file <- file.path(target_base_dir, status, grp, basename(src_file))
  file.copy(src_file, dest_file, overwrite = TRUE)
  success <- success + 1
}

cat("成功复制:", success, "/", nrow(stats_final), "\n")
if (length(fail) > 0) {
  cat("仍缺失:", length(fail), "\n")
  print(fail)
}



