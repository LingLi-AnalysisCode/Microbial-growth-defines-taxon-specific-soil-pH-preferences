
#design表
rm(list = ls())
# ============================================================
# 0. 加载R包
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(agricolae)

# ============================================================
# 1. 文件路径
# ============================================================

design_file <- "./design.txt"
TPM_file <- "./data/MAG_TPM.tsv"
MAG_group_file <- "./Process/kegg_decoder.csv"

output_dir <- "./Process"
figure_dir <- "./Fig"

dir.create(output_dir,showWarnings = FALSE,recursive = TRUE)
dir.create(figure_dir,showWarnings = FALSE, recursive = TRUE)
# ============================================================
# 2. 读取design
# treatment1已调整为与TPM列名一致
# 以design为唯一的样品和pH依据
# ============================================================
design <- read.table(design_file,sep = "\t",header = TRUE,row.names = 1, comment.char = "",check.names = FALSE, quote = "", stringsAsFactors = FALSE)
# 保留原始行名
design$design_id <- rownames(design)

# 检查必要列
required_design_columns <- c("site","fertilizer1","replicate", "treatment", "treatment1", "pH")
missing_design_columns <- setdiff( required_design_columns,colnames(design))
if (length(missing_design_columns) > 0) {
  stop(
    paste0("design文件缺少以下列：",paste(missing_design_columns, collapse = ", ") )
  )
}

# treatment1就是TPM中的样品名
design <- design %>%
  mutate(
    TPM_sample = as.character(treatment1),
    pH = as.numeric(pH),
    replicate = as.integer(replicate)
  )

# 检查pH
if (anyNA(design$pH)) {
  stop("design中的pH列存在NA或非数值内容。")
}

# 检查treatment1
if (anyNA(design$TPM_sample) ||
    any(design$TPM_sample == "")) {
  stop("design中的treatment1存在空值。")
}

# 检查重复样品名
if (anyDuplicated(design$TPM_sample) > 0) {
  
  duplicated_samples <- design %>%
    count(TPM_sample) %>%
    filter(n > 1)
  
  print(duplicated_samples)
  
  stop("design中的treatment1存在重复样品名。")
}

cat(
  "删除CS后，design样品数量：",
  nrow(design),
  "\n"
)

print(
  table(design$site)
)

# 新design删除CS后应为48个样品
if (nrow(design) != 48) {
  warning(
    paste0(
      "删除CS后design中有",
      nrow(design),
      "个样品，不是预期的48个，请检查。"
    )
  )
}

# 保存本次分析实际使用的design
write.csv(design %>% select(design_id,TPM_sample,site,fertilizer1,replicate,treatment, pH ),
          file.path(output_dir,"design_used_for_MAG_weighted_pH_noCS.csv"), row.names = FALSE)
#TPM表
# ============================================================
# 3. 读取MAG TPM矩阵
# ============================================================
TPM <- read.table(TPM_file,sep = "\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "",stringsAsFactors = FALSE)

if (!"Genome" %in% colnames(TPM)) {
  stop("MAG_TPM.tsv中没有找到Genome列。")
}

cat("原始TPM矩阵：",nrow(TPM),"行 ×",ncol(TPM),"列\n")

# ============================================================
# 4. 检查design样品与TPM列名是否匹配
# ============================================================
missing_in_TPM <- setdiff(
  design$TPM_sample,
  colnames(TPM)
)

if (length(missing_in_TPM) > 0) {
  
  print(missing_in_TPM)
  
  stop(
    paste0(
      "以下design样品没有在TPM矩阵中找到：",
      paste(
        missing_in_TPM,
        collapse = ", "
      )
    )
  )
}

cat(
  "design中的",
  nrow(design),
  "个非CS样品均已在TPM矩阵中找到。\n"
)

# TPM中存在、但本次不使用的样品
unused_TPM_samples <- setdiff(
  setdiff(
    colnames(TPM),
    "Genome"
  ),
  design$TPM_sample
)

cat(
  "TPM中未进入本次分析的样品数量：",
  length(unused_TPM_samples),
  "\n"
)

print(unused_TPM_samples)

write.csv(
  data.frame(
    TPM_sample = unused_TPM_samples
  ),
  file.path(
    output_dir,
    "unused_TPM_samples_noCS.csv"
  ),
  row.names = FALSE
)
# ============================================================
# 5. 只保留design中的48个样品
# ============================================================

TPM_used <- TPM %>%
  
  # 删除unmapped
  filter(
    Genome != "unmapped",
    !is.na(Genome),
    Genome != ""
  ) %>%
  
  # 只使用design中存在的48个样品
  select(
    Genome,
    all_of(design$TPM_sample)
  )

cat(
  "用于分析的TPM矩阵：",
  nrow(TPM_used),
  "个MAG ×",
  ncol(TPM_used) - 1,
  "个样品\n"
)

if (ncol(TPM_used) - 1 != nrow(design)) {
  stop("TPM_used样品数量与design样品数量不一致。")
}

#pH生态型表
# ============================================================
# 6. 读取MAG的pH生态型
# ============================================================
library(readr)
MAG_group <- read_csv(MAG_group_file, show_col_types = FALSE)
# 删除可能存在的索引列
MAG_group <- MAG_group %>% select(-any_of( c("...1", "X","Unnamed: 0" )))

required_group_columns <- c(  "MAG",  "pH_group")

missing_group_columns <- setdiff(
  required_group_columns,
  colnames(MAG_group)
)

if (length(missing_group_columns) > 0) {
  stop(
    paste0(
      "kegg_decoder.csv缺少以下列：",
      paste(
        missing_group_columns,
        collapse = ", "
      )
    )
  )
}

MAG_group <- MAG_group %>%
  select(
    MAG,
    pH_group
  ) %>%
  filter(
    !is.na(MAG),
    MAG != "",
    !is.na(pH_group),
    pH_group != ""
  ) %>%
  distinct()

# 检查同一个MAG是否被分到多个组
duplicated_MAG_group <- MAG_group %>%
  count(MAG) %>%
  filter(n > 1)

if (nrow(duplicated_MAG_group) > 0) {
  
  print(duplicated_MAG_group)
  
  stop("部分MAG对应多个pH_group，请检查。")
}

# ============================================================
# 7. 设置pH生态型顺序
# ============================================================

pH_order <- c(
  "Acidophile-associated MAG",
  "Acid–neutrophile MAG",
  "Neutrophile-associated MAG",
  "Neutrophile–alkaliphile MAG",
  "Alkaliphile-associated MAG",
  "Broad pH-range MAG"
)

unknown_groups <- setdiff(
  unique(MAG_group$pH_group),
  pH_order
)

if (length(unknown_groups) > 0) {
  stop(
    paste0(
      "以下pH_group名称不在预设分类中：",
      paste(
        unknown_groups,
        collapse = ", "
      )
    )
  )
}

MAG_group$pH_group <- factor(
  MAG_group$pH_group,
  levels = pH_order
)

cat("\nMAG生态型数量：\n")
print(
  table(MAG_group$pH_group)
)
#MAG丰度与生态型合并
# ============================================================
# 8. 检查MAG名称匹配
# ============================================================
MAG_in_TPM <- unique(TPM_used$Genome)
MAG_in_group <- unique(MAG_group$MAG)

common_MAG <- intersect(
  MAG_in_TPM,
  MAG_in_group
)

classified_MAG_not_in_TPM <- setdiff(
  MAG_in_group,
  MAG_in_TPM
)

TPM_MAG_without_group <- setdiff(
  MAG_in_TPM,
  MAG_in_group
)

cat(
  "\n具有pH生态型且在TPM中找到的MAG数量：",
  length(common_MAG),
  "\n"
)

cat(
  "有pH生态型但未在TPM中找到的MAG数量：",
  length(classified_MAG_not_in_TPM),
  "\n"
)

cat(
  "TPM中存在但未划分pH生态型的MAG数量：",
  length(TPM_MAG_without_group),
  "\n"
)

if (length(common_MAG) == 0) {
  stop("TPM文件与MAG生态型文件之间没有匹配的MAG。")
}

write.csv(
  data.frame(
    MAG = classified_MAG_not_in_TPM
  ),
  file.path(
    output_dir,
    "classified_MAG_not_found_in_TPM.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    MAG = TPM_MAG_without_group
  ),
  file.path(
    output_dir,
    "TPM_MAG_without_pH_group.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. 将TPM转换成长表，并加入样品pH和MAG生态型
# ============================================================

TPM_long <- TPM_used %>%
  
  # 只保留已经分类的MAG
  filter(
    Genome %in% common_MAG
  ) %>%
  
  pivot_longer(
    cols = -Genome,
    names_to = "TPM_sample",
    values_to = "TPM"
  ) %>%
  
  mutate(
    TPM = as.numeric(TPM),
    TPM = replace_na(TPM, 0)
  ) %>%
  
  left_join(
    design %>%
      select(
        TPM_sample,
        design_id,
        site,
        fertilizer1,
        replicate,
        treatment,
        pH
      ),
    by = "TPM_sample"
  ) %>%
  
  rename(
    MAG = Genome
  ) %>%
  
  left_join(
    MAG_group,
    by = "MAG"
  )

# 检查pH是否全部匹配
missing_pH_samples <- TPM_long %>%
  filter(is.na(pH)) %>%
  distinct(TPM_sample)

if (nrow(missing_pH_samples) > 0) {
  
  print(missing_pH_samples)
  
  stop("部分TPM样品没有匹配到pH。")
}

# 检查每个MAG是否均有48条记录
MAG_sample_number <- TPM_long %>%
  count(MAG)

incorrect_MAG_records <- MAG_sample_number %>%
  filter(n != nrow(design))

if (nrow(incorrect_MAG_records) > 0) {
  
  print(incorrect_MAG_records)
  
  stop("部分MAG的样品记录数不是48。")
}
定义丰度加权函数
# ============================================================
# 10. 定义丰度加权函数
# ============================================================
weighted_pH_center_function <- function( pH, abundance) {
  valid <- is.finite(pH) &
    is.finite(abundance) &
    abundance >= 0
  
  pH <- pH[valid]
  abundance <- abundance[valid]
  
  if (
    length(pH) == 0 ||
    sum(abundance) <= 0
  ) {
    return(NA_real_)
  }
  
  weighted.mean(
    x = pH,
    w = abundance
  )
}


weighted_pH_SD_function <- function(
    pH,
    abundance
) {
  
  valid <- is.finite(pH) &
    is.finite(abundance) &
    abundance >= 0
  
  pH <- pH[valid]
  abundance <- abundance[valid]
  
  if (
    length(pH) == 0 ||
    sum(abundance) <= 0
  ) {
    return(NA_real_)
  }
  
  center <- weighted.mean(
    x = pH,
    w = abundance
  )
  
  sqrt(
    sum(
      abundance *
        (pH - center)^2
    ) /
      sum(abundance)
  )
}


effective_sample_number_function <- function(
    abundance
) {
  
  abundance <- abundance[
    is.finite(abundance) &
      abundance > 0
  ]
  
  if (
    length(abundance) == 0 ||
    sum(abundance) <= 0
  ) {
    return(NA_real_)
  }
  
  sum(abundance)^2 /
    sum(abundance^2)
}


peak_pH_function <- function(
    pH,
    abundance
) {
  
  valid <- is.finite(pH) &
    is.finite(abundance) &
    abundance >= 0
  
  pH <- pH[valid]
  abundance <- abundance[valid]
  
  if (
    length(pH) == 0 ||
    sum(abundance) <= 0
  ) {
    return(NA_real_)
  }
  
  pH[
    which.max(abundance)
  ]
}


detected_min_pH_function <- function(
    pH,
    abundance
) {
  
  keep <- is.finite(pH) &
    is.finite(abundance) &
    abundance > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  min(
    pH[keep]
  )
}


detected_max_pH_function <- function(
    pH,
    abundance
) {
  
  keep <- is.finite(pH) &
    is.finite(abundance) &
    abundance > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  max(
    pH[keep]
  )
}

#计算每个MAG的丰度加权pH中心
# ============================================================
# 11. 计算每个MAG的丰度加权pH中心
# ============================================================
MAG_weighted_pH <- TPM_long %>%
  group_by(
    MAG,
    pH_group
  ) %>%
  summarise(
    
    # 主要指标：丰度加权pH中心
    weighted_pH_center =
      weighted_pH_center_function(
        pH = pH,
        abundance = TPM
      ),
    
    # 丰度沿pH梯度的加权标准差
    weighted_pH_SD =
      weighted_pH_SD_function(
        pH = pH,
        abundance = TPM
      ),
    
    # 丰度实际分布的有效样品数
    effective_sample_n =
      effective_sample_number_function(
        TPM
      ),
    
    # 48个样品中的总TPM
    total_TPM = sum(
      TPM,
      na.rm = TRUE
    ),
    
    # TPM大于0的样品数
    n_detected = sum(
      TPM > 0,
      na.rm = TRUE
    ),
    
    # TPM不小于1的样品数
    n_detected_TPM1 = sum(
      TPM >= 1,
      na.rm = TRUE
    ),
    
    # 最大TPM
    max_TPM = max(
      TPM,
      na.rm = TRUE
    ),
    
    # 丰度最高样品对应的pH
    peak_abundance_pH =
      peak_pH_function(
        pH = pH,
        abundance = TPM
      ),
    
    # 检出范围的最低pH
    detected_pH_min =
      detected_min_pH_function(
        pH = pH,
        abundance = TPM
      ),
    
    # 检出范围的最高pH
    detected_pH_max =
      detected_max_pH_function(
        pH = pH,
        abundance = TPM
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    detected_pH_range =
      detected_pH_max -
      detected_pH_min
  )

# 检查是否有总TPM为0的MAG
zero_abundance_MAG <- MAG_weighted_pH %>%
  filter(
    is.na(weighted_pH_center)
  )

if (nrow(zero_abundance_MAG) > 0) {
  warning(
    paste0(
      nrow(zero_abundance_MAG),
      "个MAG在48个样品中总TPM为0，无法计算加权pH中心。"
    )
  )
}

# 保存每个MAG的结果
write.csv(
  MAG_weighted_pH,
  file.path(
    output_dir,
    "MAG_abundance_weighted_pH_center_noCS.csv"
  ),
  row.names = FALSE
)

#按生态型汇总每个MAG的加权pH中心
# ============================================================
# 12. 按生态型汇总每个MAG的加权pH中心
# 每个MAG作为一个观测值
# ============================================================
MAG_group_summary <- MAG_weighted_pH %>%
  group_by(
    pH_group
  ) %>%
  summarise(
    
    n_MAG = sum(
      !is.na(weighted_pH_center)
    ),
    
    mean_weighted_pH = mean(
      weighted_pH_center,
      na.rm = TRUE
    ),
    
    median_weighted_pH = median(
      weighted_pH_center,
      na.rm = TRUE
    ),
    
    sd_weighted_pH = sd(
      weighted_pH_center,
      na.rm = TRUE
    ),
    
    se_weighted_pH = ifelse(
      n_MAG > 1,
      sd_weighted_pH /
        sqrt(n_MAG),
      NA_real_
    ),
    
    min_weighted_pH = min(
      weighted_pH_center,
      na.rm = TRUE
    ),
    
    max_weighted_pH = max(
      weighted_pH_center,
      na.rm = TRUE
    ),
    
    mean_weighted_pH_SD = mean(
      weighted_pH_SD,
      na.rm = TRUE
    ),
    
    median_weighted_pH_SD = median(
      weighted_pH_SD,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(MAG_group_summary)

write.csv(MAG_group_summary,
          file.path(output_dir, "MAG_ecotype_weighted_pH_summary_noCS.csv" ), row.names = FALSE)

# ============================================================
# 13. 计算每个生态型合并丰度后的总体加权pH中心
#
# 先将同一生态型所有MAG在每个样品中的TPM相加，
# 再计算该生态型整体的加权pH中心。
# ============================================================

ecotype_abundance_by_sample <- TPM_long %>%
  group_by(
    pH_group,
    TPM_sample,
    pH
  ) %>%
  summarise(
    ecotype_TPM = sum(
      TPM,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ecotype_pooled_weighted_pH <- ecotype_abundance_by_sample %>%
  group_by(
    pH_group
  ) %>%
  summarise(
    
    pooled_weighted_pH_center =
      weighted_pH_center_function(
        pH = pH,
        abundance = ecotype_TPM
      ),
    
    pooled_weighted_pH_SD =
      weighted_pH_SD_function(
        pH = pH,
        abundance = ecotype_TPM
      ),
    
    pooled_effective_sample_n =
      effective_sample_number_function(
        ecotype_TPM
      ),
    
    total_ecotype_TPM = sum(
      ecotype_TPM,
      na.rm = TRUE
    ),
    
    n_detected_samples = sum(
      ecotype_TPM > 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    MAG_group %>%
      count(
        pH_group,
        name = "n_MAG"
      ),
    by = "pH_group"
  )

print(ecotype_pooled_weighted_pH)

write.csv(
  ecotype_pooled_weighted_pH,
  file.path(
    output_dir,
    "pH_ecotype_pooled_abundance_weighted_center_noCS.csv"
  ),
  row.names = FALSE
)

#检验
# ============================================================
# 14. Kruskal-Wallis整体检验
# ============================================================
stat_data <- MAG_weighted_pH %>%
  filter(
    !is.na(weighted_pH_center),
    !is.na(pH_group)
  ) %>%
  droplevels()

KW_result <- kruskal.test(
  weighted_pH_center ~ pH_group,
  data = stat_data
)

print(KW_result)

KW_table <- data.frame(
  method = KW_result$method,
  statistic = unname(
    KW_result$statistic
  ),
  df = unname(
    KW_result$parameter
  ),
  p_value = KW_result$p.value
)

write.csv(
  KW_table,
  file.path(
    output_dir,
    "MAG_weighted_pH_Kruskal_Wallis.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. agricolae事后比较
# 各生态型两两比较，组内使用FDR校正
# ============================================================

posthoc_result <- agricolae::kruskal(
  y = stat_data$weighted_pH_center,
  trt = stat_data$pH_group,
  alpha = 0.05,
  p.adj = "fdr",
  group = TRUE,
  console = TRUE
)

posthoc_raw <- posthoc_result$groups

posthoc_groups <- data.frame(
  pH_group = rownames(posthoc_raw),
  mean_value = posthoc_raw[[1]],
  significance_group = posthoc_raw$groups,
  stringsAsFactors = FALSE
)

posthoc_groups$pH_group <- factor(
  posthoc_groups$pH_group,
  levels = pH_order
)

posthoc_groups <- posthoc_groups %>%
  arrange(pH_group)

print(posthoc_groups)

write.csv(
  posthoc_groups,
  file.path(
    output_dir,
    "MAG_weighted_pH_posthoc_FDR_groups.csv"
  ),
  row.names = FALSE
)

#作图
# ============================================================
# 16. 设置绘图颜色
# ============================================================
pH_colors <- c(
  "Acidophile-associated MAG"  = "#D55E00",
  "Acid–neutrophile MAG"       = "#BDC609",
  "Neutrophile-associated MAG" = "#009E73",
  "Neutrophile–alkaliphile MAG"= "#56B4E9",
  "Alkaliphile-associated MAG" = "#0072B2",
  "Broad pH-range MAG"         = "#CC79A7")

# ============================================================
# 17. 构建横轴标签和显著性字母位置
# ============================================================

group_n <- stat_data %>%
  count(pH_group, .drop = FALSE ) %>%
  mutate( axis_label = paste0( pH_group, "\n(n = ", n, ")" ))

axis_labels <- setNames( group_n$axis_label, as.character( group_n$pH_group))

y_padding <- 0.04 * diff( range( design$pH, na.rm = TRUE ) )

letter_position <- stat_data %>%
  group_by( pH_group ) %>%
  summarise( y_position = max( weighted_pH_center,na.rm = TRUE  ) + y_padding, .groups = "drop"
  ) %>%
  left_join( posthoc_groups %>% select( pH_group, significance_group ),  by = "pH_group"  )

# ============================================================
# 18. 绘制每个MAG的丰度加权pH中心
#
# 彩色点：每个MAG
# 箱线图：同一生态型内MAG的分布
# 白色菱形：该生态型合并全部MAG后的总体加权中心
# 字母：FDR校正后的事后比较
# ============================================================

p <- ggplot(stat_data,aes(x = pH_group,y = weighted_pH_center,color = pH_group,fill = pH_group)) +
  # 生理分类边界，只作为参考
  geom_hline( yintercept = 5.5,linetype = 2, color = "grey65", linewidth = 0.45 ) +
  geom_hline( yintercept = 8.0,linetype = 2,color = "grey65",linewidth = 0.45 ) +
  # 箱线图
  geom_boxplot(width = 0.55,alpha = 0.18, outlier.shape = NA,linewidth = 0.6 ) +
  # 每个MAG
  geom_jitter( width = 0.14, height = 0, size = 2.8, alpha = 0.90 ) +
  # 生态型合并丰度后的总体加权中心
  geom_point(data = ecotype_pooled_weighted_pH,aes( x = pH_group,y = pooled_weighted_pH_center ),
             inherit.aes = FALSE,shape = 23,size = 3.5,stroke = 0.8,color = "black", fill = "white" ) +
  # 显著性字母
  geom_text( data = letter_position,aes(  x = pH_group, y = y_position, label = significance_group  ),
             inherit.aes = FALSE,color = "black", size = 4) +
  scale_color_manual(values = pH_colors, drop = FALSE ) +
  scale_fill_manual( values = pH_colors,drop = FALSE ) +
  scale_x_discrete(labels = axis_labels,drop = FALSE ) +
  scale_y_continuous(expand = expansion( mult = c( 0.04,0.13))) +
  labs(
    x = NULL,
    y = "Abundance-weighted pH center"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 35, hjust = 1, vjust = 1,  size = 9 ),
    axis.text.y = element_text( size = 9 ),
    axis.title.y = element_text(  size = 11 ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p)
ggsave(filename = file.path(figure_dir, "MAG_abundance_weighted_pH_center_noCS.pdf"), plot = p, width = 6, height = 6)

# ============================================================
# 19. 单独绘制每个生态型的合并丰度加权中心
# 这是描述性结果，每个生态型只有一个值
# ============================================================
p_pooled <- ggplot( ecotype_pooled_weighted_pH, aes(x = pH_group, y = pooled_weighted_pH_center,color = pH_group)) +
  geom_hline( yintercept = 5.5, linetype = 2,color = "grey65", linewidth = 0.45 ) +
  geom_hline(yintercept = 8.0, linetype = 2,  color = "grey65",  linewidth = 0.45  ) +
  geom_point(  size = 4  ) +
  geom_text(  aes( label = round( pooled_weighted_pH_center,  2 ) ),nudge_y = 0.12, color = "black", size = 3.5) +
  scale_color_manual( values = pH_colors, drop = FALSE) +
  scale_x_discrete( labels = axis_labels,drop = FALSE ) +
  scale_y_continuous( expand = expansion( mult = c( 0.05,0.12 ) ) ) +
  labs( x = NULL,y = "Pooled abundance-weighted pH center"  ) +
  theme_bw() +
  theme(   axis.text.x = element_text( angle = 35, hjust = 1, vjust = 1, size = 9 ),
           panel.grid.major.x = element_blank(),
           panel.grid.minor =   element_blank(),
           legend.position =   "none"  )

print(p_pooled)
ggsave( filename = file.path( figure_dir, "pH_ecotype_pooled_weighted_center_noCS.pdf" ),  plot = p_pooled,  width = 6,  height = 6)

# ============================================================
# 20. 分析完成后的控制台信息
# ============================================================

cat("\n================ 分析完成 ================\n")

cat( "实际使用样品数：",nrow(design), "\n")
cat( "实际参与计算的MAG数：",nrow(MAG_weighted_pH), "\n")
cat( "Kruskal-Wallis P值：", KW_result$p.value,"\n")

cat(
  "\n主要输出文件：\n",
  "1. Process/design_used_for_MAG_weighted_pH_noCS.csv\n",
  "2. Process/MAG_abundance_weighted_pH_center_noCS.csv\n",
  "3. Process/MAG_ecotype_weighted_pH_summary_noCS.csv\n",
  "4. Process/pH_ecotype_pooled_abundance_weighted_center_noCS.csv\n",
  "5. Process/MAG_weighted_pH_Kruskal_Wallis.csv\n",
  "6. Process/MAG_weighted_pH_posthoc_FDR_groups.csv\n",
  "7. Fig/MAG_abundance_weighted_pH_center_noCS.pdf\n",
  "8. Fig/pH_ecotype_pooled_weighted_center_noCS.pdf\n"
)

