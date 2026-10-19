rm(list=ls())

subphylum_colors <- c(
  "Actinomycetota"           = "#F685F9",
  "Planctomycetota"          = "#CD00FF",
  "Acidobacteriota"          = "#7570B3",
  "Chloroflexota"            = "#E7298A", 
  "Bacillota"                = "#66A61E",
  "Bacteroidota"             = "#E6AB02",  
  "Gemmatimonadota"          = "#A6761D",  
  "Armatimonadota"           = "#0000DB", 
  "Verrucomicrobiota"        = "#1F78B4", 
  "Cyanobacteriota"          = "#00C3FF", 
  "BRC1"                     = "#FB8072",
  "Nitrospirota"             = "#80B1D3",
  "Rhodothermota"            = "#FDB462",
  "candidate division WPS-2" = "#B3DE69",
  "Latescibacteria"          = "#FCCDE5",
  "Thermomicrobiota"         = "#BC80BD",
  "Deinococcota"             = "#8DD3C7",
  "Abditibacteriota"         = "#CCEBC5",
  "Hydrogenedentes"          = "#FFED6F",
  "Poribacteria"             = "#8C510A",
  "Alphaproteobacteria"      = "#D73027",
  "Betaproteobacteria"       = "#4575B4",
  "Gammaproteobacteria"      = "#F46D43", 
  "Deltaproteobacteria"      = "#542788", 
  "Oligoflexia"              = "#7FBC41",
  "Unassigned"               = "#A3A7A8"
)

library(tidyr)
library(dplyr)

stats <- read.csv("Process/pH_response_tax_classification.csv",  header = TRUE)

# 删除 Family 为 "Unassigned" 的行，然后计算每个科内各生理类型的比例
# 1. 获取长表（含 n 和 prop）
family_long <- stats %>%
  filter(Family != "Unassigned") %>%
  count(Phylum, Class,Family, Physiological_Group, name = "n") %>%
  group_by(Phylum, Class,Family) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 2. 计算每个科的总数（作为单独数据框）
family_total <- family_long %>%
  group_by(Phylum,Class, Family) %>%
  summarise(total_n = sum(n), .groups = "drop")

# 3. 将比例转为宽表
#长数据转宽数据
#同时保留比例和计数    
family_wide <- family_long %>%
  select(Phylum,Class, Family, Physiological_Group, prop, n) %>%
  pivot_wider(
    names_from = Physiological_Group,
    values_from = c(prop, n), # 填充值
    values_fill = 0  # 缺失值填0
  ) 

# 4. 合并总数
family_wide_with_total <- family_wide %>%
  left_join(family_total, by = c("Phylum","Class", "Family")) %>%
  relocate(total_n, .after = Family)   # 把总数放在科名后面，方便查看

write.csv(family_wide_with_total, "Process/pH_response_wide_Family.csv", row.names = FALSE)


rm(list=ls())
stats <- read.csv("Process/pH_response_wide_Family.csv",  header = TRUE)
# ============================================
# 1. 加载必要的 R 包
# ============================================
library(dplyr)
library(pheatmap)

# ============================================
# 2. 创建 Subphylum 列（如果尚未创建）
# ============================================
stats <- stats %>%
  mutate(Subphylum = ifelse(
    Phylum == "Pseudomonadota",
    Class,      # 变形菌门用 Class 作为亚门级别
    Phylum
  ))

# ============================================
# 3. 选择并重命名比例列（简化列名）
# ============================================
stats_plot <- stats %>%
  select(Subphylum, Family,
         Obligate.acidophile      = n_Obligate.acidophile,
         Facultative.acidophile   = n_Facultative.acidophile,
         Obligate.neutrophile     = n_Obligate.neutrophile,
         Facultative.neutrophile  = n_Facultative.neutrophile,
         Obligate.alkaliphile     = n_Obligate.alkaliphile,
         Facultative.alkaliphile  = n_Facultative.alkaliphile,
         pH.insensitive           = n_pH.insensitive#,
         # Others                   = n_Others
  )

# ============================================
# 4. 定义热图列的顺序（固定顺序，不聚类）
# ============================================
col_order <- c(
  "Obligate.acidophile", "Facultative.acidophile",
  "Obligate.neutrophile", "Facultative.neutrophile",
  "Obligate.alkaliphile", "Facultative.alkaliphile",
  "pH.insensitive"#, "Others"
)

# ============================================
# 5. 排序：先按 Subphylum 字母顺序，再按 Family 字母顺序
# ============================================
stats_plot <- stats_plot %>%
  mutate(Subphylum = as.character(Subphylum),
         Family = as.character(Family)) %>%
  arrange(Subphylum, Family)

# ============================================
# 6. 将 Subphylum 转换为因子，水平按字母排序（保证图注顺序正确）
# ============================================
stats_plot$Subphylum <- factor(stats_plot$Subphylum,
                               levels = sort(unique(stats_plot$Subphylum)))

# ============================================
# 7. 构建数值矩阵（行名为 Family）
# ============================================
heat_mat <- stats_plot %>%
  select(all_of(col_order)) %>%
  as.matrix()
rownames(heat_mat) <- stats_plot$Family

# ============================================
# 8. 构建行注释数据框（Subphylum）
# ============================================
annotation_row <- data.frame(Subphylum = stats_plot$Subphylum)
rownames(annotation_row) <- rownames(heat_mat)

# ============================================
# 9. 计算 Subphylum 之间的行分隔位置
# ============================================
subphylum_counts <- table(stats_plot$Subphylum)
gaps_row <- cumsum(subphylum_counts)
gaps_row <- gaps_row[-length(gaps_row)]   # 去掉最后一个（结尾不需要分隔线）

# ============================================
# 10. 自定义颜色：为每个 Subphylum 指定颜色
#     注意：颜色向量的名称必须与因子水平（字母顺序的亚门）完全一致
# ============================================
# 获取排序后的亚门名称（即因子水平）
sorted_subphyla <- levels(stats_plot$Subphylum)

# 在这里根据您的实际亚门名称和喜欢的颜色进行自定义
#读取自定义的颜色subphylum_colors <- c("Actinomycetota" = "#F685F9","Planctomycetota"= "#CD00FF",……)
# 按字母顺序重排颜色向量，确保与因子水平顺序一致
subphylum_colors <- subphylum_colors[sorted_subphyla]

# 检查是否所有亚门都已分配颜色，如果缺失则报错
if (any(is.na(subphylum_colors))) {
  missing <- sorted_subphyla[is.na(subphylum_colors)]
  stop("以下亚门缺少颜色定义：", paste(missing, collapse = ", "))
}

# 将颜色向量放入 ann_colors 列表中（必须命名为 Subphylum）
ann_colors <- list(Subphylum = subphylum_colors)

# ============================================
# 11. 绘制热图
# ============================================
pdf(file = "./Fig/Taxonomy_pH_response_Brith_rate_classification_Z_Score_Family_heatmap.pdf", width = 5, height = 16)
pheatmap(heat_mat,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         scale = "row",                # 对每行（每个科）进行标准化
         cluster_rows = FALSE,         # 不聚类行（按给定的顺序）
         cluster_cols = FALSE,         # 不聚类列（按 col_order 顺序）
         gaps_row = gaps_row,          # 亚门之间添加分隔线
         #color = colorRampPalette(c("#7DA7D5", "#FDFDFD", "#D83434"))(100),
         color = colorRampPalette(c("#8C5100", "#F5F5F5", "#35978F"))(100),
         border_color = FALSE,
         fontsize_row = 8,
         fontsize_col = 10,
         angle_col = 45)
dev.off() 