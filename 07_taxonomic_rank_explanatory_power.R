rm(list=ls())#清除工作历史

data= read.csv("Process/pH_response_tax_classification.csv", header=T, row.names=1, comment.char="")
# 将 class 列转换为数字编码
data$class_numeric <- as.integer(factor(data$Physiological_Group))

library(picante)
library(ape)
library(factoextra)
library(car)
library(nlme)
set.seed(123)

data_model <- lme(class_numeric~1, 
                  random = ~1|Phylum/Class/Order/Family/Genus,#
                  data = data,na.action= na.omit)
summary(data_model)
data_model_v=varcomp(data_model,T,T)
plot(data_model_v)

data_model_total=as.matrix(data_model_v)
data_model_total=as.data.frame(data_model_total)
data_model_total$varia =rownames(data_model_total)

data_model_total$varia  = factor(data_model_total$varia, levels=rev(rownames(data_model_total)))#rev：倒序


write.csv(data_model_total,"./process/lme_total.csv")

data_model_total$varia  = factor(data_model_total$varia, levels=rownames(data_model_total))#rev：倒序

library("ggplot2")
q1 <- ggplot(data_model_total, aes(x = varia, y = V1)) +
  geom_point(alpha = 0.8, size = 2, shape = 21) +
  geom_line(aes(group = 1), color = "#000000") +   # 添加 group = 1
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    strip.text = element_text(size=12,face="plain",color="black"),
    text=element_text(size=12,face="plain",color="black"),
    axis.title=element_text(size=12,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.background = element_rect(fill=alpha("white",0)))+
  xlab(expression(""))+
  ylab(expression("Proportion variation explained" ))

q1
ggsave("./Fig/Proportion variation explained.pdf", q1, width = 5, height = 4)