#查看不同MAG代谢通路完整性的差异（只关注筛选出来的MAG）
rm(list = ls())
library(dplyr)

Meta_data <- read.table("./data/kegg_decoder_results_fixed.txt", 
                        sep = "\t", header = TRUE, row.names = 1, 
                        stringsAsFactors = FALSE,
                        check.names = FALSE)   # 关键参数，保留原始列名

Ecotypes <- read.csv("./Process/Ecotypes_MAGtrait_group.csv",  row.names=1,header = TRUE) #

data=merge(Ecotypes[,c("MAG","pH_group")],Meta_data,by.x="MAG",by.y="row.names",all.x=T)
data <- subset(data, !(MAG %in% c("MAG_652","MAG_416","MAG_635"))) #这三个是古菌
write.csv(data,"./Process/kegg_decoder.csv")

#筛选不同组有差异的代谢通路
rm(list = ls())
## 读取数据
dat <- read_csv("./Process/kegg_decoder.csv", show_col_types = FALSE)
index <- dat %>% select(-any_of(c("...1", "X", "Unnamed: 0")))

## 建立分组列
index$treatment <- factor(index$pH_group)

## 加载检验函数
source("data/boxplerk.R")

## 代谢通路列
pathway_cols <- setdiff(
  names(index),
  c("MAG", "pH_group", "treatment")
)

## 非参数检验
l <- list()

for (nm in pathway_cols) {
  
  X <- as.numeric(index[[nm]])
  
  ## 全部相同，例如全部为0，则跳过
  if (length(unique(X[!is.na(X)])) < 2) {
    next
  }
  
  a <- boxplerk(
    X,
    index$treatment,
    ylab = nm,
    xlab = "pH group",
    bcol = "bisque",
    p.adj = "fdr",
    las = 1
  )
  ##用Kruskal-Wallis事后比较字母标记组间显著性
  stat <- a$comparison
  stat$p_value <- a$p.value
  stat$index <- nm
  stat$treatment <- rownames(stat)
  
  l[[nm]] <- stat
}

## 合并结果
Result <- do.call(rbind, l)

write.csv(Result,"./Process/MAG_metabolic_pathway_非参.csv")

#筛选具有差异的代谢通路
# 筛选显著（去除 NA）
sig_result <- Result[!is.na(Result$p_value) & Result$p_value < 0.05, ]

# 提取显著代谢物的列名（去重）
sig_metabolites <- unique(sig_result$index)

cat("显著差异的代谢物数量：", length(sig_metabolites), "\n")

# 提取子集数据：保留 MAG, pH_group 和显著代谢物列
sub_data <- dat[, c("MAG", "pH_group", sig_metabolites), drop = FALSE]
write.csv(sub_data,"./Process/MAG_metabolic_pathway_显著差异.csv")

#差异代谢通路气泡图
rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)
#--------------------------------------------------
# 1. 读取显著差异通路数据
#--------------------------------------------------
dat <- read_csv("./Process/MAG_metabolic_pathway_显著差异.csv", show_col_types = FALSE)
sub_data <- dat %>% select( -any_of( c("...1", "X", "Unnamed: 0") ) )

#--------------------------------------------------
# 2. 读取通路分类
#--------------------------------------------------
Class <- read.table("./data/kegg_decoder_headers_classification_table.txt",sep = "\t",
                    header = TRUE,row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)

# 如果Header_Name被读成了行名，重新转回普通列
if (!"Header_Name" %in% colnames(Class)) {
  Class <- Class %>%
    rownames_to_column("Header_Name")
}

#--------------------------------------------------
# 3. 确定通路列
#--------------------------------------------------
valid_pathways <- setdiff(colnames(sub_data), c("MAG", "pH_group"))

# 确保通路为数值型
sub_data[valid_pathways] <- lapply(sub_data[valid_pathways], function(x) as.numeric(as.character(x)))

#--------------------------------------------------
# 4. 转成长表
#--------------------------------------------------
long_data <- sub_data %>%
  pivot_longer(
    cols = all_of(valid_pathways),
    names_to = "Pathway",
    values_to = "Completeness"
  )

#--------------------------------------------------
# 5. 按pH组汇总
#--------------------------------------------------
plot_data <- long_data %>%
  group_by(pH_group, Pathway) %>%
  summarise(n_MAG = n(),
            # 包括0在内的平均通路完整性
            mean_completeness =  mean(Completeness, na.rm = TRUE),
            median_completeness =  median(Completeness, na.rm = TRUE),
            # 检测到该通路组成的MAG比例
            prevalence = mean(Completeness > 0, na.rm = TRUE),
            n_present = sum(Completeness > 0, na.rm = TRUE),
            .groups = "drop" )

#--------------------------------------------------
# 6. 添加Level1和Level2分类
#--------------------------------------------------
plot_data <- plot_data %>%
  left_join(
    Class %>%
      select(Header_Name,Level1_Name,Level2_Name),
    by = c("Pathway" = "Header_Name")) %>%
  mutate(Level1 = ifelse(is.na(Level1_Name),"Unknown",Level1_Name),
         Level2 = ifelse(is.na(Level2_Name),"Unknown",Level2_Name))

#--------------------------------------------------
# 7. 设置pH组顺序
#--------------------------------------------------
pH_order <- c(
  "Acidophile-associated MAG",
  "Acid–neutrophile MAG",
  "Neutrophile-associated MAG",
  "Neutrophile–alkaliphile MAG",
  "Alkaliphile-associated MAG",
  "Broad pH-range MAG")

plot_data$pH_group <- factor(plot_data$pH_group, levels = pH_order)
#--------------------------------------------------
# 8. 通路排序
# 按Level1、Level2以及最大组间差异排序
#--------------------------------------------------
pathway_order <- plot_data %>%
  group_by(Pathway,Level1,Level2) %>%
  summarise(group_range = max(mean_completeness) - min(mean_completeness),.groups = "drop" ) %>%
  arrange(Level1, Level2, desc(group_range), Pathway) %>%
  pull(Pathway)

# 反向，使第一条显示在图的上方
plot_data$Pathway <- factor(
  plot_data$Pathway,levels = rev(unique(pathway_order)))

#--------------------------------------------------
# 9. 给横轴添加各组MAG数量
#--------------------------------------------------
group_n <- sub_data %>%
  count(pH_group) %>%
  mutate(axis_label = paste0(pH_group, "\n(n = ", n,")"))
axis_labels <- setNames(group_n$axis_label, group_n$pH_group)

#--------------------------------------------------
# 10. 绘制气泡热图
#--------------------------------------------------
p <- ggplot(plot_data, aes(x = pH_group, y = Pathway)) +
  geom_point(aes(size = prevalence,color = mean_completeness)) +
  # 平均完整性固定为0–1
  scale_color_gradientn(
    #   colors = c("#8C5100","#A86E1F","#C89F65","#DCCFBF", "#8FC5B4", "#35978F"),
    colors = c("#8C5100", "#B8E2DC", "#35978F"), #"#8C5100", "#F5F5F5", "#35978F" # "#F7FBFF","#9ECAE1","#3182BD","#08519C"
    limits = c(0, 1),
    oob = scales::squish,
    name = "Mean pathway\ncompleteness") +
  # 气泡大小表示检出率
  scale_size_continuous(
    limits = c(0, 1),
    range = c(1, 8),
    breaks = c(0.25,0.50,0.75,1.00),
    labels = percent_format(accuracy = 1),
    name = "MAG prevalence") +
  scale_x_discrete(labels = axis_labels) +
  # 按Level1分块
  facet_grid(Level1+Level2 ~ .,scales = "free_y",space = "free_y",switch = "y") +
  labs(x = NULL,y = NULL,
       #  title = "Differential metabolic pathways among pH-associated MAG groups",
       # subtitle ="Color represents mean pathway completeness; point size represents the proportion of MAGs with completeness > 0"
  ) +
  theme_bw() +
  theme(
    panel.grid.major  = element_line(color = "grey90",linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(angle = 35,hjust = 1,vjust = 1,size = 9),
    axis.text.y       = element_text(size = 8),
    # strip.text.y.left = element_text(angle = 0,face = "bold",size = 9),
    strip.background  = element_rect(fill = "grey95",color = "grey75"),
    legend.position   = "right",
    plot.title        = element_text(face = "bold",size = 13),
    plot.subtitle     = element_text(size = 10))
p
#--------------------------------------------------
# 11. 保存
#--------------------------------------------------
figure_height <- max( 7, 3 + 0.28 * length(valid_pathways))
ggsave("./Fig/MAG_significant_pathways_bubble_heatmap.pdf", plot = p,width = 8,height = figure_height,limitsize = FALSE)

#整体代谢通路合并MAG热图
rm(list = ls())
setwd("G:/1 实验/实验：不同点qSIP个体CUE/")
# ============================================
# 0. 加载包
# ============================================
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# ============================================
# 1. 读取全部通路数据
# ============================================
data <- read_csv("./Process/kegg_decoder.csv",show_col_types = FALSE)
data <- data %>% select(-any_of(c("...1", "X", "Unnamed: 0")))

# ============================================
# 2. 读取显著通路数据
# 用于给全部通路标记Significant/Not significant
# ============================================
sig_data <- read_csv("./Process/MAG_metabolic_pathway_显著差异.csv",show_col_types = FALSE)
sig_data <- sig_data %>% select(-any_of(c("...1", "X", "Unnamed: 0")))

# 如果显著结果文件是MAG × pathway的宽表
if (all(c("MAG", "pH_group") %in% colnames(sig_data))) {
  significant_pathways <- setdiff(colnames(sig_data),c("MAG", "pH_group"))
  # 如果显著结果文件是boxplerk的长表，并包含index列
} else if ("index" %in% colnames(sig_data)) {
  significant_pathways <- unique(sig_data$index)
} else {
  stop("无法从显著结果文件识别通路名称。请确认文件中包含MAG/pH_group列或index列。")
}

# ============================================
# 3. 读取通路分类表
# ============================================
Class <- read.table("./宏基因组/代谢通路完整性/kegg_decoder_headers_classification_table.txt",
                    sep = "\t",header = TRUE,row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)

# 如果Header_Name被设置为行名，则恢复为普通列
if (!"Header_Name" %in% colnames(Class)) {
  Class <- Class %>% rownames_to_column("Header_Name") }

# ============================================
# 4. 确定全部代谢通路列
# ============================================
pathway_cols <- setdiff(colnames(data), c("MAG", "pH_group"))

# 转为数值
data[pathway_cols] <- lapply(data[pathway_cols],
                             function(x) { as.numeric(as.character(x)) } )

# 检查NA
if (anyNA(data[pathway_cols])) {
  warning("部分通路完整性存在NA；计算均值时将忽略NA。请确认NA是否具有实际含义。" )}

# ============================================
# 5. 识别全部为0的通路
# ============================================
all_zero_pathways <- pathway_cols[ colSums(data[pathway_cols], na.rm = TRUE ) == 0 ]

# 不显著通路要保留；
# 但全部为0的通路没有绘图信息，建议不放入热图
plot_pathways <- setdiff(pathway_cols, all_zero_pathways)

# 单独输出全零通路，附表中仍然提供
write_csv(data.frame( Pathway = all_zero_pathways, Note = "Completeness was zero in all MAGs" ),
          "./Process/All_zero_metabolic_pathways.csv")

# ============================================
# 6. 转成长表
# ============================================
long_data <- data %>%
  pivot_longer(
    cols = all_of(plot_pathways),
    names_to = "Pathway",
    values_to = "Completeness"
  )

# ============================================
# 7. 计算每个pH组的平均完整性和检出率
# ============================================
group_summary <- long_data %>%
  group_by(Pathway, pH_group) %>%
  summarise(n_MAG = n(),
            mean_completeness = mean(Completeness,na.rm = TRUE),
            median_completeness = median(Completeness,na.rm = TRUE),
            prevalence = mean(Completeness > 0,na.rm = TRUE),
            n_present = sum(Completeness > 0,na.rm = TRUE),
            .groups = "drop"  )

# ============================================
# 8. 设置pH组顺序
# ============================================
pH_order <- c(
  "Acidophile-associated MAG",
  "Acid–neutrophile MAG",
  "Neutrophile-associated MAG",
  "Neutrophile–alkaliphile MAG",
  "Alkaliphile-associated MAG",
  "Broad pH-range MAG")

group_summary$pH_group <- factor(group_summary$pH_group, levels = pH_order)

# ============================================
# 9. 转换成 Pathway × pH_group矩阵
# ============================================
mean_mat <- group_summary %>%
  select(Pathway,pH_group,mean_completeness) %>%
  pivot_wider(
    names_from = pH_group,
    values_from = mean_completeness,
    values_fill = 0
  ) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

# 按指定pH组顺序排列
mean_mat <- mean_mat[ ,pH_order,drop = FALSE]

# ============================================
# 10. 构建通路分类和显著性信息
# ============================================
pathway_info <- data.frame(
  Pathway = rownames(mean_mat),
  stringsAsFactors = FALSE) %>%
  left_join(
    Class %>%
      select(Header_Name,Level1_Name,Level2_Name), by = c("Pathway" = "Header_Name")) %>%
  mutate(Level1 = ifelse(is.na(Level1_Name),"Unknown",Level1_Name),
         Level2 = ifelse(is.na(Level2_Name),"Unknown",Level2_Name),
         Significance = ifelse(Pathway %in% significant_pathways,"Significant","Not significant"))

# ============================================
# 11. 计算每条通路的最大组间差异
# 用于通路排序
# ============================================
range_data <- data.frame(Pathway = rownames(mean_mat),
                         group_range = apply(mean_mat, 1, function(x) {max(x, na.rm = TRUE) - min(x, na.rm = TRUE) } ),
                         stringsAsFactors = FALSE)

pathway_info <- pathway_info %>%
  left_join(range_data, by = "Pathway")

# ============================================
# 12. 通路排序
#
# 先Level1
# 再Level2
# 再把显著通路放在该分类前面
# 再按组间差异大小排序
# ============================================
pathway_order <- pathway_info %>%
  arrange(Level1,Level2,desc(Significance == "Significant"),
          desc(group_range),Pathway) %>%
  pull(Pathway)

mean_mat <- mean_mat[pathway_order,  , drop = FALSE]
pathway_info <- pathway_info %>%
  slice(match(pathway_order,Pathway))
# ============================================
# 13. 构建热图行注释
#
# Level2用于排序，不再使用颜色注释，
# 避免出现非常庞大的Level2图例
# ============================================
# ============================================
# 行注释
# ============================================
row_anno <- data.frame(
  Level1 = pathway_info$Level1,
  Level2 = pathway_info$Level2,
  Significance = pathway_info$Significance,
  stringsAsFactors = FALSE
)

rownames(row_anno) <- pathway_info$Pathway

row_anno$Level1 <- factor(
  row_anno$Level1,
  levels = unique(row_anno$Level1)
)

row_anno$Level2 <- factor(
  row_anno$Level2,
  levels = unique(row_anno$Level2)
)

row_anno$Significance <- factor(
  row_anno$Significance,
  levels = c("Significant", "Not significant")
)

# ============================================
# Level1颜色
# ============================================

level1_levels <- levels(row_anno$Level1)
n1 <- length(level1_levels)

if (n1 <= 12) {
  
  level1_colors <- setNames(
    brewer.pal(max(3, n1), "Set3")[seq_len(n1)],
    level1_levels
  )
  
} else {
  
  level1_colors <- setNames(
    hcl.colors(n1, "Dynamic"),
    level1_levels
  )
}

# ============================================
# Level2颜色
# ============================================

level2_levels <- levels(row_anno$Level2)
n2 <- length(level2_levels)

if (n2 <= 12) {
  
  level2_colors <- setNames(
    brewer.pal(max(3, n2), "Pastel1")[seq_len(n2)],
    level2_levels
  )
  
} else {
  
  level2_colors <- setNames(
    hcl.colors(n2, "Dynamic"),
    level2_levels
  )
}

# ============================================
# 显著性颜色
# ============================================

significance_colors <- c(
  "Significant" = "#D55E00",
  "Not significant" = "#D9D9D9"
)

ann_colors <- list(
  Level1 = level1_colors,
  Level2 = level2_colors,
  Significance = significance_colors
)

# ============================================
# Level2分隔线
# ============================================

level12_group <- interaction(
  row_anno$Level1,
  row_anno$Level2,
  drop = TRUE,
  lex.order = TRUE
)

level12_run <- rle(
  as.character(level12_group)
)

gaps_row <- cumsum(level12_run$lengths)

if (length(gaps_row) > 1) {
  gaps_row <- gaps_row[-length(gaps_row)]
} else {
  gaps_row <- NULL
}

# ============================================
# 16. 给列名增加各组MAG数量
# ============================================
group_n <- data %>%
  count(pH_group)

column_labels <- sapply(
  pH_order,
  function(g) {
    current_n <- group_n$n[group_n$pH_group == g]
    paste0(g,"\n(n = ",current_n,")")
  }
)

# ============================================
# 17. 热图颜色
#
# 保持与你显著通路气泡图一致
# ============================================
heat_colors <- colorRampPalette(c("#8C5100", "#F5F5F5", "#35978F"))(100)

# 固定通路完整性范围为0–1
heat_breaks <- seq(0,1,length.out = 101)

# ============================================
# 18. 自动设置图片高度
# ============================================
figure_height <- max(12,0.18 * nrow(mean_mat))
# ============================================
# 19. 绘制全部通路平均完整性热图
# ============================================

pdf("./Fig/All_pathways_pH_group_mean_heatmap.pdf",width = 11,height = figure_height)

pheatmap(
  mean_mat,
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  color = heat_colors,
  breaks = heat_breaks,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 9,
  angle_col = 45,
  labels_col = column_labels,
  main = paste0(
    "Mean metabolic pathway completeness across pH groups\n",
    "Significant and non-significant pathways (n = ", nrow(mean_mat), ")"  )
)

dev.off()

# ============================================
# 20. 输出全部通路组统计附表
# ============================================

all_pathway_table <- group_summary %>%
  left_join(
    pathway_info %>%
      select(Pathway,Level1,Level2,Significance,group_range),
    by = "Pathway") %>%
  arrange(Level1,Level2,
          desc(Significance == "Significant"),
          desc(group_range),Pathway,pH_group)

write_csv(all_pathway_table,"./Process/All_pathways_pH_group_summary.csv")


# ============================================
#功能MAG流行率热图  显著差异的8条通路
# ============================================

rm(list = ls())
pkgs <- c("data.table","dplyr","tidyr","tibble","stringr","ggplot2","patchwork","scales")
for (pkg in pkgs) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(123)
metadata_file <- "./Process/Bacterial_MAG_metadata.rds"
ko_file <- "./data/kegg_function_MAG_01matrix.tsv"
mapping_file <- "./data/KEGG_KO_Tab_Table.txt"
analysis_dir <- "./Process"
figure_dir <- "./Fig"

dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(metadata_file)) stop("没有找到元数据：", metadata_file)
if (!file.exists(ko_file)) stop("没有找到KO矩阵：", ko_file)
if (!file.exists(mapping_file)) stop("没有找到KEGG映射表：", mapping_file)

pH_order <- c("Acidophile-associated MAG","Acid–neutrophile MAG","Neutrophile-associated MAG",
              "Neutrophile–alkaliphile MAG","Alkaliphile-associated MAG","Broad pH-range MAG")
pH_short <- c("Acidophile-associated MAG" = "Acid","Acid–neutrophile MAG" = "A-Neu","Neutrophile-associated MAG" = "Neu",
              "Neutrophile–alkaliphile MAG" = "N-Alk","Alkaliphile-associated MAG" = "Alk","Broad pH-range MAG" = "Broad")

pathway_colors <- c("alpha-amylase" = "#1B9E77","starch/glycogen synthesis" = "#66A61E","TCA Cycle" = "#377EB8",
                    "Mixed acid: PEP to Succinate via OAA, malate & fumarate" = "#4DAF4A","riboflavin biosynthesis" = "#E6AB02",
                    "sulfite dehydrogenase" = "#7570B3","ammonia oxidation (amo/pmmo)" = "#D95F02",
                    "Cytochrome b6/f complex" = "#E7298A")

metadata <- readRDS(metadata_file) %>% mutate(MAG = trimws(as.character(MAG)), 
                                              pH_group = factor(trimws(as.character(pH_group)), levels = pH_order)) %>% arrange(pH_group, MAG)

if (anyDuplicated(metadata$MAG)) stop("元数据中存在重复MAG。")
if (any(is.na(metadata$pH_group))) stop("存在无法匹配到六个pH生态型的MAG。")

mag_ids <- metadata$MAG
cat("参与分析的MAG数量：", length(mag_ids), "\n")
print(table(metadata$pH_group))

# ============================================================
# 1. 定义8条显著功能的局部相关KO集合
# ============================================================
P <- function(pathway, pathway_order, step_order, step, kos, scope) {
  tibble(Significant_pathway = pathway, Pathway_order = pathway_order, Step_order = step_order,
         KO_order = seq_along(kos), Functional_step = step, KO = kos, Scope = scope)
}

local_panel_raw <- bind_rows(
  P("alpha-amylase", 1, 1, "Alpha-amylase", c("K01176","K07405","K05343"), "Direct alpha-amylase orthologs"),
  
  P("starch/glycogen synthesis", 2, 1, "UDP-glucose activation", "K00963", "UDP-glucose precursor route"),
  P("starch/glycogen synthesis", 2, 2, "UDP-glucose-dependent glucan synthesis",  c("K00693","K00750","K16150","K16153","K13679","K20812"), "Glucan synthase alternatives"),
  P("starch/glycogen synthesis", 2, 3, "ADP-glucose activation", "K00975", "ADP-glucose precursor route"),
  P("starch/glycogen synthesis", 2, 4, "ADP-glucose-dependent starch synthesis",c("K00703","K13679","K20812"), "Starch synthase alternatives"),
  P("starch/glycogen synthesis", 2, 5, "Glucan branching", c("K00700","K16149"), "Branching enzyme alternatives"),
  
  P("TCA Cycle", 3, 1, "Oxaloacetate and acetyl-CoA to citrate", c("K01647","K05942","K01659"), "Citrate synthase alternatives"),
  P("TCA Cycle", 3, 2, "Citrate to isocitrate", c("K01681","K27802","K01682"), "Aconitate hydratase alternatives or subunits"),
  P("TCA Cycle", 3, 3, "Isocitrate to 2-oxoglutarate", c("K00031","K00030"), "Isocitrate dehydrogenase alternatives"),
  P("TCA Cycle", 3, 4, "2-Oxoglutarate to succinyl-CoA", c("K00164","K00658","K01616","K00382","K00174","K00175","K00176","K00177"), "Oxidative and ferredoxin-dependent alternatives"),
  P("TCA Cycle", 3, 5, "Succinyl-CoA to succinate", c("K01902","K01903","K01899","K01900","K18118"), "Succinyl-CoA synthetase alternatives"),
  P("TCA Cycle", 3, 6, "Succinate and fumarate interconversion", c("K00234","K00235","K00236","K00237","K25801","K00239","K00240","K00241","K00242","K18859","K18860","K00244","K00245","K00246","K00247"), "Succinate dehydrogenase and fumarate reductase alternatives"),
  P("TCA Cycle", 3, 7, "Fumarate to malate", c("K01676","K01679","K01677","K01678"), "Fumarate hydratase alternatives"),
  P("TCA Cycle", 3, 8, "Malate to oxaloacetate", c("K00026","K00025","K00024","K00116"), "Malate dehydrogenase alternatives"),
  
  P("Mixed acid: PEP to Succinate via OAA, malate & fumarate", 4, 1, "PEP to oxaloacetate", c("K01596","K20370","K01610"), "PEP carboxylation alternatives"),
  P("Mixed acid: PEP to Succinate via OAA, malate & fumarate", 4, 2, "Oxaloacetate to malate", c("K00025","K00026","K00024","K00051","K00116"), "Malate dehydrogenase alternatives"),
  P("Mixed acid: PEP to Succinate via OAA, malate & fumarate", 4, 3, "Malate to fumarate", c("K01676","K01677","K01678","K01679"), "Fumarate hydratase alternatives"),
  P("Mixed acid: PEP to Succinate via OAA, malate & fumarate", 4, 4, "Fumarate to succinate", c("K00234","K00235","K00236","K00237","K25801","K00239","K00240","K00241","K00242","K18859","K18860","K00244","K00245","K00246","K00247"), "Fumarate reductase and succinate dehydrogenase alternatives"),
  
  P("riboflavin biosynthesis", 5, 1, "GTP-derived precursor formation", c("K01497","K14652","K01498","K00082","K11752","K22912","K20860","K20861","K20862","K21063","K21064"), "Early riboflavin precursor reactions"),
  P("riboflavin biosynthesis", 5, 2, "Ribulose-phosphate-derived precursor", c("K02858","K14652"), "Second precursor route"),
  P("riboflavin biosynthesis", 5, 3, "Lumazine formation", "K00794", "Lumazine synthase"),
  P("riboflavin biosynthesis", 5, 4, "Riboflavin formation", "K00793", "Riboflavin synthase"),
  P("riboflavin biosynthesis", 5, 5, "Riboflavin activation to FMN or FAD", c("K20884","K22949","K11753"), "Downstream flavin activation"),
  
  P("sulfite dehydrogenase", 6, 1, "Cytochrome-dependent sulfite dehydrogenase SorAB", c("K05301","K00386"), "SorAB sulfite oxidation system"),
  P("sulfite dehydrogenase", 6, 2, "Quinone-dependent sulfite dehydrogenase SoeABC", c("K21307","K21308","K21309"), "SoeABC sulfite oxidation system"),
  
  P("ammonia oxidation (amo/pmmo)", 7, 1, "Ammonia or methane monooxygenase subunit A", c("K10944","K28504"), "Ambiguous pmoA/amoA and ammonia-specific amoA"),
  P("ammonia oxidation (amo/pmmo)", 7, 2, "Ammonia or methane monooxygenase subunit B", "K10945", "Shared amoB or pmoB ortholog"),
  P("ammonia oxidation (amo/pmmo)", 7, 3, "Ammonia or methane monooxygenase subunit C", "K10946", "Shared amoC or pmoC ortholog"),
  P("ammonia oxidation (amo/pmmo)", 7, 4, "Hydroxylamine to nitrite", "K10535", "Hydroxylamine oxidation step"),
  
  P("Cytochrome b6/f complex", 8, 1, "Core cytochrome b6/f complex", c("K02635","K02637","K02634","K02636","K02642","K02643","K03689","K02640"), "Cytochrome b6/f complex subunits")
)

# 同一个KO在同一条显著通路中可能出现在多个步骤中，将步骤合并，避免后续重复连接
local_panel <- local_panel_raw %>%
  group_by(Significant_pathway, Pathway_order, KO) %>%
  summarise(Step_order = min(Step_order), KO_order = min(KO_order), Functional_step = paste(unique(Functional_step), collapse = " / "), Scope = paste(unique(Scope), collapse = "; "), .groups = "drop") %>%
  arrange(Pathway_order, Step_order, KO_order, KO)

# ============================================================
# 2. 从KEGG映射表补充KO名称
# ============================================================
kegg_mapping_raw <- fread(mapping_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE, colClasses = "character")
if (!all(c("KO_ID","KO_Description") %in% names(kegg_mapping_raw))) stop("KEGG映射表缺少KO_ID或KO_Description列。")

ko_info <- kegg_mapping_raw %>% transmute(KO = trimws(as.character(KO_ID)), KO_Description = trimws(as.character(KO_Description))) %>% filter(!is.na(KO), KO != "") %>% distinct(KO, .keep_all = TRUE)

local_panel <- local_panel %>%
  left_join(ko_info, by = "KO") %>%
  mutate(KO_name = str_trim(str_extract(KO_Description, "^[^;]+")), KO_name = ifelse(is.na(KO_name) | KO_name == "", Functional_step, KO_name), KO_label = paste0(KO, ": ", str_trunc(KO_name, 34))) %>%
  arrange(Pathway_order, Step_order, KO_order, KO)

write.csv(local_panel, file.path(analysis_dir, "Significant_pathway_local_KO_definition.csv"), row.names = FALSE)

# ============================================================
# 3. 读取MAG × KO的0/1矩阵
# ============================================================

ko_dat <- fread(ko_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
colnames(ko_dat)[1] <- "KO"
ko_dat$KO <- trimws(as.character(ko_dat$KO))

if (any(is.na(ko_dat$KO) | ko_dat$KO == "")) stop("KO矩阵第一列存在空KO编号。")

missing_MAG <- setdiff(mag_ids, colnames(ko_dat))
if (length(missing_MAG) > 0) stop("KO矩阵缺少以下MAG：", paste(missing_MAG, collapse = ", "))

if (anyDuplicated(ko_dat$KO)) {
  ko_dat <- ko_dat %>% group_by(KO) %>% summarise(across(all_of(mag_ids), ~max(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
}

original_kos <- unique(ko_dat$KO)

ko_mat <- t(as.matrix(ko_dat[, mag_ids, drop = FALSE]))
storage.mode(ko_mat) <- "numeric"
rownames(ko_mat) <- mag_ids
colnames(ko_mat) <- ko_dat$KO

if (anyNA(ko_mat)) stop("KO矩阵中存在NA。")
if (any(!ko_mat %in% c(0,1))) stop("KO矩阵中存在0和1以外的数值。")

cat("原始KO矩阵MAG数量：", nrow(ko_mat), "\n")
cat("原始KO矩阵KO数量：", ncol(ko_mat), "\n")

# ============================================================
# 4. 补齐局部面板中不存在于原始KO矩阵的KO
# ============================================================

panel_kos <- unique(local_panel$KO)
missing_panel_kos <- setdiff(panel_kos, colnames(ko_mat))

if (length(missing_panel_kos) > 0) {
  zero_mat <- matrix(0, nrow = nrow(ko_mat), ncol = length(missing_panel_kos), dimnames = list(rownames(ko_mat), missing_panel_kos))
  ko_mat <- cbind(ko_mat, zero_mat)
}

sub_mat <- ko_mat[mag_ids, panel_kos, drop = FALSE]

cat("局部面板KO数量：", length(panel_kos), "\n")
cat("原矩阵中不存在、补为全0的KO数量：", length(missing_panel_kos), "\n")

if (length(missing_panel_kos) > 0) {
  write.csv(data.frame(KO = missing_panel_kos), file.path(analysis_dir, "Local_panel_KOs_absent_from_original_matrix.csv"), row.names = FALSE)
}

# ============================================================
# 5. 生成MAG层面的KO有无长表
# ============================================================

ko_presence_base <- as.data.frame(sub_mat) %>% rownames_to_column("MAG") %>% pivot_longer(-MAG, names_to = "KO", values_to = "Presence") %>% left_join(metadata %>% select(MAG, pH_group), by = "MAG")

# 一个KO可能属于多条显著通路，所以按KO与局部通路定义进行多对多连接
ko_presence <- ko_presence_base %>% left_join(local_panel, by = "KO", relationship = "many-to-many") %>% mutate(pH_group = factor(pH_group, levels = pH_order)) %>% arrange(Pathway_order, Step_order, KO_order, pH_group, MAG)

# ============================================================
# 6. 计算六组KO检出率
# ============================================================

ko_group_prevalence <- ko_presence %>%
  group_by(Significant_pathway, Pathway_order, Functional_step, Step_order, KO, KO_order, KO_label, Scope, pH_group) %>%
  summarise(n_MAG = n(), n_present = sum(Presence), Prevalence = mean(Presence), .groups = "drop") %>%
  mutate(pH_group = factor(pH_group, levels = pH_order)) %>%
  arrange(Pathway_order, Step_order, KO_order, pH_group)

# ============================================================
# 7. KO质控表
# 修正点：first(KO)确保每组只返回一个TRUE/FALSE
# ============================================================

ko_qc <- ko_presence %>%
  group_by(Significant_pathway, Pathway_order, Functional_step, Step_order, KO, KO_order, KO_label, Scope) %>%
  summarise(KO_in_original_matrix = first(KO) %in% original_kos, Detected_in_any_MAG = any(Presence > 0), n_MAG_present = sum(Presence), Overall_prevalence = mean(Presence), .groups = "drop") %>%
  arrange(Pathway_order, Step_order, KO_order)

write.csv(ko_presence, file.path(analysis_dir, "Significant_pathway_local_KO_MAG_presence.csv"), row.names = FALSE)
write.csv(ko_group_prevalence, file.path(analysis_dir, "Significant_pathway_local_KO_group_prevalence.csv"), row.names = FALSE)
write.csv(ko_qc, file.path(analysis_dir, "Significant_pathway_local_KO_QC.csv"), row.names = FALSE)

print(ko_qc %>% count(Significant_pathway, KO_in_original_matrix, Detected_in_any_MAG))

# ============================================================
# 8. 绘制纵向排列的正方形KO检出率热图
# ============================================================

draw_local_KO_heatmap <- function(group_data, include_all_zero = FALSE, output_name) {
  plot_data <- group_data
  
  if (!include_all_zero) {
    detected_kos <- ko_qc %>% filter(Detected_in_any_MAG) %>% select(Significant_pathway, KO)
    plot_data <- plot_data %>% semi_join(detected_kos, by = c("Significant_pathway","KO"))
  }
  
  pathways_used <- names(pathway_colors)[names(pathway_colors) %in% unique(plot_data$Significant_pathway)]
  if (length(pathways_used) == 0) stop("没有可用于绘图的通路或KO。")
  
  plot_list <- vector("list", length(pathways_used))
  panel_heights <- numeric(length(pathways_used))
  
  for (i in seq_along(pathways_used)) {
    pathway_name <- pathways_used[i]
    sub_dat <- plot_data %>% filter(Significant_pathway == pathway_name) %>% arrange(Step_order, KO_order, KO)
    
    row_order <- sub_dat %>%
      distinct(Step_order, KO_order, Functional_step, KO, KO_label) %>%
      arrange(Step_order, KO_order, KO) %>%
      mutate(Row_label = paste0(Step_order, ". ", KO_label)) %>%
      pull(Row_label)
    
    sub_dat <- sub_dat %>% mutate(Row_label = paste0(Step_order, ". ", KO_label), Row_label = factor(Row_label, levels = rev(row_order)))
    panel_heights[i] <- n_distinct(sub_dat$Row_label)
    base_col <- unname(pathway_colors[pathway_name])
    
    p_i <- ggplot(sub_dat, aes(x = pH_group, y = Row_label, alpha = Prevalence)) +
      geom_tile(fill = base_col, color = "grey20", linewidth = 0.35, width = 0.95, height = 0.95) +
      scale_x_discrete(labels = pH_short, limits = pH_order, drop = FALSE, expand = c(0,0)) +
      scale_y_discrete(drop = FALSE, expand = c(0,0)) +
      scale_alpha_continuous(limits = c(0,1), range = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = percent_format(accuracy = 1), name = "MAG prevalence") +
      coord_fixed(ratio = 1, clip = "off") +
      labs(title = pathway_name, x = NULL, y = NULL) +
      theme_bw() +
      theme(panel.grid = element_blank(), panel.border = element_rect(color = "grey20", linewidth = 0.6), plot.title = element_text(face = "bold", size = 10, color = base_col, margin = margin(b = 3)), axis.text.y = element_text(size = 6.8, color = "black"), axis.ticks = element_blank(), plot.margin = margin(2,6,2,6), legend.position = "right")
    
    if (i < length(pathways_used)) {
      p_i <- p_i + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
      p_i <- p_i + theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9, color = "black"))
    }
    
    plot_list[[i]] <- p_i
  }
  
  panel_heights <- pmax(panel_heights, 1)
  
  combined <- wrap_plots(plot_list, ncol = 1, heights = panel_heights, guides = "collect") +
    plot_annotation(title = "Local KO composition of significant pathways across pH ecotypes", subtitle = "Only KOs directly related to each significant functional pathway are included") &
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9))
  
  total_rows <- sum(panel_heights)
  figure_height <- max(10, total_rows * 0.35 + 3)
  
  ggsave(file.path(figure_dir, output_name), combined, width = 10.5, height = figure_height, limitsize = FALSE, useDingbats = FALSE)
  combined
}

# 正文候选：只显示43个MAG中至少检出一次的KO
p_detected <- draw_local_KO_heatmap(ko_group_prevalence, include_all_zero = FALSE, 
                                    output_name = "Significant_pathway_local_KO_heatmap_detected.pdf")

# 补充材料：显示完整局部KO面板，包括全0 KO
p_complete <- draw_local_KO_heatmap(ko_group_prevalence, include_all_zero = TRUE, 
                                    output_name = "Significant_pathway_local_KO_heatmap_complete_panel.pdf")

cat("\n================ 分析完成 ================\n")
cat("局部KO定义：", file.path(analysis_dir, "Significant_pathway_local_KO_definition.csv"), "\n")
cat("MAG层面KO有无：", file.path(analysis_dir, "Significant_pathway_local_KO_MAG_presence.csv"), "\n")
cat("六组KO检出率：", file.path(analysis_dir, "Significant_pathway_local_KO_group_prevalence.csv"), "\n")
cat("KO质控表：", file.path(analysis_dir, "Significant_pathway_local_KO_QC.csv"), "\n")
cat("正文候选热图：", file.path(figure_dir, "Significant_pathway_local_KO_heatmap_detected.pdf"), "\n")
cat("完整面板热图：", file.path(figure_dir, "Significant_pathway_local_KO_heatmap_complete_panel.pdf"), "\n")
