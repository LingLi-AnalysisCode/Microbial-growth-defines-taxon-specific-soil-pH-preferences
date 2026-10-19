rm(list = ls())
library(dplyr)

# ---------------------------- 1. 路径设置 ------------------------------------
#Thr_dir <- "./Fig/非活性调0/Threshold_Results_Family"
Thr_dir <- "./Fig"
stats_dir <- file.path(Thr_dir, "Statistics")
Meta_dir <- "./data"

# 输入文件：过滤结果表（您已经生成的 stats_final）
stats <- file.path(stats_dir, "Physiology_Filter_Results.tsv") # 请确认文件名
stats <- read.table(stats, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")

OTU_MAG <- file.path(Meta_dir, "otu_to_mag_mapping.tsv") # 请确认文件名
OTU_MAG <- read.table(OTU_MAG, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data=merge(stats,OTU_MAG,by.x="Tax_Name",by.y="OTU_ID",all=F)
# 提取从开头到第二个下划线之前的内容（即 MAG_数字）
data$Best_MAG <- sub("^(MAG_\\d+)_.*", "\\1", data$Best_MAG)
# 查看修改后的列
head(data$Best_MAG)
colnames(data)[colnames(data) == "Best_MAG"] <- "MAG"

data <- data %>%  filter(Pass == TRUE)

tapply(data$MAG, data$Physiological_Group, function(x) length(unique(x)))
write.csv(data,"./Process/Physiological_Group_MAG_mach_OTU.csv")

rm(list=ls())
library(tidyr)
library(dplyr)

stats <- read.csv("./Process/Physiological_Group_MAG_mach_OTU.csv",  row.names=1,header = TRUE) #54个独有MAG

# 统计每个 MAG 在每个生理类型下的频次
stats_wide <- stats %>%
  group_by(MAG, Physiological_Group) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Physiological_Group, 
              values_from = Count, 
              values_fill = 0)
# 查看结果
head(stats_wide)

MAG_tax <- read.table("./data/gtdbtk.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
MAG_data <- merge(stats_wide,MAG_tax,by="MAG",all.x=T)

MAG_trait <- read.table("./data/summarize.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
MAG_data <- merge(MAG_data,MAG_trait,by="MAG",all.x=T)

MAG_trait1 <- read.table("./data/gRodon_MAG.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
MAG_data1 <- merge(MAG_data,MAG_trait1,by.x="MAG",by.y="ID",all.x=T)

MAG_data1$Maximum_growth_rate=1/MAG_data1$d #最大生长速率是d的倒数
MAG_data1$Genome_Size=MAG_data1$Genome_Size/1000000 #换算成Mbp

write.csv(MAG_data1, "Process/MAG_tree_active.csv", row.names = FALSE)


rm(list=ls())
library(tidyr)
library(dplyr)

stats <- read.csv("Process/MAG_tree_active.csv", header = TRUE) #54个独有MAG

sub_data <- stats %>%
  mutate(
    acid = (Facultative.acidophile + Obligate.acidophile) > 0,
    neut = (Facultative.neutrophile + Obligate.neutrophile) > 0,
    alk  = (Facultative.alkaliphile + Obligate.alkaliphile) > 0
  ) %>%
  mutate(
    pH_group = case_when(
      # 单一偏好
      acid & !neut & !alk ~ "Acidophile-associated MAG",
      !acid & neut & !alk ~ "Neutrophile-associated MAG",
      !acid & !neut & alk  ~ "Alkaliphile-associated MAG",
      # 两种组合
      acid & neut & !alk ~ "Acid–neutrophile MAG",
      !acid & neut & alk  ~ "Neutrophile–alkaliphile MAG",
      # 包含酸+碱 → 宽pH范围
      acid & alk ~ "Broad pH-range MAG",
      # 什么都没有 → 标记为Unclassified
      TRUE ~ "Unclassified"
    )
  ) %>%
  # 删除Unclassified行
  filter(pH_group != "Unclassified")

# 查看分组计数
table(sub_data$pH_group)
write.csv(sub_data,"./Process/Ecotypes_MAGtrait_group.csv")

data <- subset(sub_data, !(MAG %in% c("MAG_652","MAG_416","MAG_635"))) #这三个是古菌
table(data$pH_group)
write.csv(data ,"./Process/Ecotypes_MAGtrait_group_Bac.csv")

rm(list=ls())
setwd('G:/1 实验/实验：不同点qSIP个体CUE/')
sub_data <- read.csv("./Process/Ecotypes_MAGtrait_group_Bac.csv", header = TRUE) #46个独有MAG

library(ggradar)
library(tidyverse)
#"Genome_Size","GC_Content","Completeness","Contamination","CDS","CUBHE","Maximum_growth_rate"
#按 Physiological_Group 聚合数据（这里用中位数）
df <- sub_data %>%
  group_by(pH_group) %>%
  summarise(across(c(Genome_Size, GC_Content, CDS,CUBHE, Maximum_growth_rate),
                   median, na.rm = TRUE)) %>%
  as.data.frame()
#归一化
df <- df %>%
  mutate(across(-pH_group, ~ (.x - min(.x)) / (max(.x) - min(.x))))
############## 雷达图 ------------
p2 <- ggradar(
  df,
  #最内圈（圆心）对应的数值。中间圈对应的数值（“平均值”）。最外圈对应的数值。
  grid.min = 0, grid.mid = 0.5, grid.max = 1, 
  values.radar = c("0%", "50%", "100%"),#这三条网格线旁边显示的标签文字，长度必须为 3。
  #是否分别显示三条网格线的标签（即 values.radar 中的文字），字的大小
  label.gridline.min = F, 
  label.gridline.mid = F,
  label.gridline.max = F,
  grid.label.size = 3,  
  
  #分别对应内圈、中圈、外圈的线型 和 颜色
  grid.line.width = 1,
  gridline.min.linetype = "longdash",
  gridline.mid.linetype = "longdash",
  gridline.max.linetype = "longdash",
  gridline.min.colour = "grey",
  gridline.mid.colour = "#FF0000",#"#FF0000" #"grey"
  gridline.max.colour = "grey", 
  
  axis.label.offset = 1.15,
  axis.label.size = 3,
  axis.line.colour = "grey",  
  
  # 雷达图线的粗细和颜色：
  group.line.width = 1,
  group.point.size = 5,
  group.colours = c(
    "Acidophile-associated MAG"      = "#D55E00",   # 橙红
    "Acid–neutrophile MAG"           = "#BDC609",   # 橙黄 ##E69F00 #BDC609
    "Neutrophile-associated MAG"     = "#009E73",   # 草绿
    "Neutrophile–alkaliphile MAG"    = "#56B4E9",   # 天蓝
    "Alkaliphile-associated MAG"     = "#0072B2",   # 深蓝
    "Broad pH-range MAG"             = "#CC79A7"   ), # 紫粉 #
  # 背景边框线颜色：
  background.circle.colour = "white", #"#D7D6D1"
  background.circle.transparency = 0.2,
  legend.position = "none" #"none", # "right"  
)+
  theme(plot.background = element_blank(),
        panel.background = element_blank())
p2
ggsave("./Fig/Ecotypes_MAGtrait_雷达图_legend.pdf", p2, width = 8, height = 7)

############## 圆环注释+文本注释 ------------
library(ggplot2)

# 注释数据：
tmp <- data.frame(x = rep(1, 5),
                  y = rep(1, 5),
                  group = colnames(df)[-1])
# 绘图：
p1 <- ggplot()+
  # 圆环：
  geom_bar(data = tmp, aes(x, y, fill = group), stat = "identity", position = "dodge")+
  # 文本注释： 
  # geom_text(data = tmp,aes(x, y, label = gsub("[.]", " "), group = group),
  #           color = "white",
  #          position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#e64b35", "#4dbbd5", "#00a087","#3c5488","#f39b7f"))+
  ylim(-5.5,2)+
  # 0.63是计算得来，5个色块，第一个色块的正中心要对准0的位置，
  # 所以2pi/10=0.628即为第一个色块左边界的位置
  coord_polar(start = -0.63)+
  theme_void()+
  theme(legend.position = "none")

p1

######## 拼图 ---------------------
library(patchwork)
p3=p1 + inset_element(p2, left = 0, bottom = 0, right = 0.99, top = 0.99)
p3
ggsave("./Fig/Ecotypes_MAGtrait_雷达图.pdf", p3, width = 8, height = 7)
