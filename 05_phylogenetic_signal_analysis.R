
library(reshape2)
library(ggplot2)
library(vegan)
library(dplyr)

Thr_dir <- "./Fig/Threshold_Results_tax_id"
stats_dir <- file.path(Thr_dir, "Statistics")
Seq_dir <- "./data"

# 输入文件：过滤结果表（您已经生成的 stats_final）
stats <- file.path(stats_dir, "Physiology_Filter_Results.tsv") # 请确认文件名
stats <- read.table(stats, sep = "\t", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")
stats <- stats %>%  filter(Pass == TRUE)
#stats <- subset(stats, Physiological_Group != "Others")

tax <- file.path(Seq_dir, "pick_taxonomy.txt") 
tax <- read.table(tax, sep='\t', header=T, comment.char='', row.names=1)
spe=merge(tax,stats,by.x="row.names",by.y="Tax_Name",all=F)
write.csv(spe, "Process/pH_response_tax_classification.csv", row.names = FALSE)


rm(list=ls())#清除工作历史
library("ctv")
library(geiger)
library(picante)
library(phytools)
set.seed(123)

data= read.csv("Process/pH_response_tax_classification.csv", header=T, row.names=1, comment.char="")
# 将 class 列转换为数字编码
data$class_numeric <- as.integer(factor(data$Physiological_Group))
data1=data[,c("Physiological_Group","class_numeric")]

#读取树文件
library(ape)
tree<- read.tree("data/0.Galaxy117-[FastTree.nwk].nwk")

anoleData=data1[which(rownames(data1) %in% tree[["tip.label"]]),]#匹配树
otu <- data.frame(t(anoleData))
##使用 picante 包计算 MNTD 和 NTI
library(picante)
#修剪系统发育树以仅包含群落数据集中存在的物种或具有非缺失性状数据的物种
anoleTree <- prune.sample(otu, tree)
name.check(anoleTree,anoleData)
anoleSize <- anoleData[, 2]
names(anoleSize)=rownames(anoleData)

phylosig(anoleTree,anoleSize,method = "K",test = T)
phylosig(anoleTree,anoleSize,method = "lambda",test=T)


#分组
library(phytools)

# 2. 确保数据顺序与树标签一致
anoleData <- anoleData[anoleTree$tip.label, ]

# 3. 获取所有分组
groups <- unique(anoleData$Physiological_Group)

# 4. 创建结果表格
results_KL <- data.frame(
  Group = groups,
  K = NA,
  K_pvalue = NA,
  Lambda = NA,
  Lambda_pvalue = NA
)

# 5. 循环计算每个分类的 K 和 Lambda
for (i in seq_along(groups)) {
  g <- groups[i]
  
  # 构造二值性状：属于该组 = 1，不属于 = 0
  trait <- ifelse(anoleData$Physiological_Group == g, 1, 0)
  names(trait) <- rownames(anoleData)   # 务必给向量命名
  
  # 计算 Blomberg's K（test=TRUE 做随机化检验）
  k_res <- phylosig(anoleTree, trait, method = "K", test = TRUE)
  
  # 计算 Pagel's λ（test=TRUE 做似然比检验）
  l_res <- phylosig(anoleTree, trait, method = "lambda", test = TRUE)
  
  # 存入结果
  results_KL$K[i] <- k_res$K
  results_KL$K_pvalue[i] <- k_res$P
  results_KL$Lambda[i] <- l_res$lambda
  results_KL$Lambda_pvalue[i] <- l_res$P
}

# 打印结果
print(results_KL)
