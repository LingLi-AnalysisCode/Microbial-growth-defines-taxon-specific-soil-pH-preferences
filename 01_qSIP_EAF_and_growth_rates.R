rm(list=ls())
library(phyloseq)
library(qsip)

sam1 <- read.table('./data/samdata.txt', sep='\t', header=T, row.names=1, comment.char='', check.names=F)
sam1$isotope=as.factor(sam1$isotope)

otu <- read.table('./data/pick_otutab.txt', sep='\t', header=T, comment.char='', check.names=F)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]
tax <- read.table('./data/GTDB_taxonomy_adj.txt', sep='\t', header=T, comment.char='', row.names=1)


sam1=sam1[which(sam1$fertilizer %in% c("NJXCK","NJXNPK","NJXOM","NJXNPKM")), ]

#删除第0天没有的OTU
library(vegan)

#江西
otu_CK=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXCK" & sam1$day == "0"), ]))],method = "pa")) >= 2),
           which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXCK"), ]))]#
otu_NPK=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXNPK" & sam1$day == "0"), ]))],method = "pa")) >= 2),
            which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXNPK"), ]))]#
otu_OM=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXOM" & sam1$day == "0"), ]))],method = "pa")) >= 2),
           which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXOM"), ]))]#
otu_NPKM=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXNPKM" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "NJXNPKM"), ]))]#
otu_CKNPK=merge(otu_CK,otu_NPK,by="row.names",all=T)
otu_CKNPKOM=merge(otu_CKNPK,otu_OM,by.x="Row.names",by.y="row.names",all=T)
otu=merge(otu_CKNPKOM,otu_NPKM,by.x="Row.names",by.y="row.names",all=T)

rownames(otu) <- otu$Row.names
otu <- otu[,-1]
otu[is.na(otu)] <- 0
write.csv(otu,"./process/filter_OTUtable_xinjiangxi.csv")
tax=tax[rownames(otu),]

otu <- as.matrix(otu)
tax <- as.matrix(tax)

dat <- phyloseq(otu_table(otu, taxa_are_rows=T),
                tax_table(tax),
                sample_data(sam1))

#save phyloseq data here and make a new file for qSIP analysis
save(dat, file='process/phyloseq_data_xinjiangxi.RData')

load('process/phyloseq_data_xinjiangxi.RData')
#specify qsip-related dat_cona from phyloseq object, create phylosip object
qsip_dat <- specify_qsip(dat,
                         density='density_g_ml',
                         abund='qpcr_abund',
                         rep_id='sample_id',
                         rep_num='replicate',
                         rep_group='fertilizer',
                         iso='18O',
                         iso_trt='isotope',
                         timepoint='day')                         

#specify filtering level:
#and in more than 5 of the fractions that make up a replicate
qsip_dat@qsip@filter_levels <- create_filters(replicate = 2, fraction = 5, soft=1)##几个重复，几层中出现#参考ISME方法

qsip_dat1_1 <- calc_excess(qsip_dat,
                           separate_label=T,
                           max_label = 0.98,#
                           filter=T,
                           correction=T)
#extract excess atom fraction (EAF)                           
library(reshape2)
qsip_enrich1 <- dcast(qsmelt(qsip_dat1_1, include=c('excess'), regex=T),
                      tax_id~sample_id+fertilizer,value.var="atom_excess")

#95%CI
qsip_dat1_2 <- calc_excess(qsip_dat,
                           #percent = FALSE,##是否将结果展示为百分比，即EAF值表示为**%，默认为F
                           ci_method = c("bootstrap"),
                           ci = 0.95, iters = 999,
                           #offset_taxa = 0.1,##用于计算偏移校正值的类群的百分比，在correction = T时有效
                           max_label = 0.98,##默认为1，同位素最高掺入量，受实验中提供的原子富集百分比的限制，原子富集量
                           recalc = TRUE,##是否重新计算WAD和分子量值。默认为T，使用bootstrapped将自动重新计算所有值。
                           separate_light = FALSE,##是否将三个重复的WAD取均值，FALSE为取均值
                           separate_label = TRUE,
                           filter = TRUE,##是不是要执行上一步的过滤
                           correction = T)

atom_excess_ci_l <- qsip_dat1_2@qsip[["atom_excess_ci_l"]]
atom_excess <- qsip_dat1_2@qsip[["atom_excess"]]
atom_excess_ci_u <- qsip_dat1_2@qsip[["atom_excess_ci_u"]]
colnames(atom_excess_ci_l)=c("NJXCK_ci_l","NJXNPK_ci_l","NJXNPKM_ci_l","NJXOM_ci_l")
colnames(atom_excess)=c("NJXCK_mean","NJXNPK_mean","NJXNPKM_mean","NJXOM_mean")
colnames(atom_excess_ci_u)=c("NJXCK_ci_u","NJXNPK_ci_u","NJXNPKM_ci_u","NJXOM_ci_u")

qsip_enrich2=as.data.frame(cbind(atom_excess_ci_l, atom_excess, atom_excess_ci_u))
atom_excess1=merge(qsip_enrich1,qsip_enrich2,by.x="tax_id",by.y="row.names",all=F)

write.csv(atom_excess1,'process/EAF_xinjiangxi.csv')

par(mfrow=c(1,1), mar=c(2,2,1,1)+.1)
hist(atom_excess)

# per-capita growth
qsip_dat2_1 <- calc_pop(qsip_dat ,
                        separate_label=T,
                        max_label = 0.98,
                        filter=T,
                        correction=T,
                        match_replicate=T)##F为利用在time0的丰度，在每组重复中取平均值；T为一一匹配
qsip_rate1 <- dcast(melt(qsmelt(qsip_dat2_1, include=c('rate'), regex=T),ID=c("tax_id")),
                    tax_id~variable+sample_id)

#calculate population growth and death
qsip_dat2_2 <- calc_pop(qsip_dat, 
                        ci_method = c("bootstrap"),
                        ci = 0.95,
                        separate_label=T, 
                        max_label = 0.98,
                        filter=T, 
                        correction=T, 
                        match_replicate=T) 

birth_rate_ci_l <- qsip_dat2_2@qsip[["birth_rate_ci_l"]]
birth_rate <- qsip_dat2_2@qsip[["birth_rate"]]
birth_rate_ci_u <- qsip_dat2_2@qsip[["birth_rate_ci_u"]]                      

death_rate_ci_l <- qsip_dat2_2@qsip[["death_rate_ci_l"]]
death_rate <- qsip_dat2_2@qsip[["death_rate"]]
death_rate_ci_u <- qsip_dat2_2@qsip[["death_rate_ci_u"]]    

growth_rate_ci_l <- qsip_dat2_2@qsip[["growth_rate_ci_l"]]
growth_rate <- qsip_dat2_2@qsip[["growth_rate"]]
growth_rate_ci_u <- qsip_dat2_2@qsip[["growth_rate_ci_u"]]   

colnames(birth_rate_ci_l)=c("NJXCK_birth_ci_l","NJXNPK_birth_ci_l","NJXNPKM_birth_ci_l","NJXOM_birth_ci_l")
colnames(birth_rate)=c("NJXCK_birth_mean","NJXNPK_birth_mean","NJXNPKM_birth_mean","NJXOM_birth_mean")
colnames(birth_rate_ci_u)=c("NJXCK_birth_ci_u","NJXNPK_birth_ci_u","NJXNPKM_birth_ci_u","NJXOM_birth_ci_u")

colnames(death_rate_ci_l)=c("NJXCK_death_ci_l","NJXNPK_death_ci_l","NJXNPKM_death_ci_l","NJXOM_death_ci_l")
colnames(death_rate)=c("NJXCK_death_mean","NJXNPK_death_mean","NJXNPKM_death_mean","NJXOM_death_mean")
colnames(death_rate_ci_u)=c("NJXCK_death_ci_u","NJXNPK_death_ci_u","NJXNPKM_death_ci_u","NJXOM_death_ci_u")

colnames(growth_rate_ci_l)=c("NJXCK_growth_ci_l","NJXNPK_growth_ci_l","NJXNPKM_growth_ci_l","NJXOM_growth_ci_l")
colnames(growth_rate)=c("NJXCK_growth_mean","NJXNPK_growth_mean","NJXNPKM_growth_mean","NJXOM_growth_mean")
colnames(growth_rate_ci_u)=c("NJXCK_growth_ci_u","NJXNPK_growth_ci_u","NJXNPKM_growth_ci_u","NJXOM_growth_ci_u")

qsip_rate2=as.data.frame(cbind(birth_rate_ci_l, birth_rate, birth_rate_ci_u,
                               death_rate_ci_l,death_rate,death_rate_ci_u,
                               growth_rate_ci_l,growth_rate,growth_rate_ci_u))

rates=merge(qsip_rate1,qsip_rate2,by="tax_id",by.x="tax_id",by.y="row.names",all=F)

write.csv(rates,'process/rates_xinjiangxi.csv')
par(mfrow=c(1,3), mar=c(2,2,1,1)+.1)
hist(birth_rate); hist(death_rate); hist(growth_rate)

#RG
sam1 <- read.table('./data/samdata.txt', sep='\t', header=T, row.names=1, comment.char='', check.names=F)
sam1$isotope=as.factor(sam1$isotope)

otu <- read.table('./data/pick_otutab.txt', sep='\t', header=T, comment.char='', check.names=F)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]
tax <- read.table('./data/GTDB_taxonomy_adj.txt', sep='\t', header=T, comment.char='', row.names=1)

sam1=sam1[which(sam1$fertilizer %in% c("RGLN","RGMN","RGHN")), ]

#单独跑，对OTU进行两步筛选再合并(第0天3个重复至少2个存在数值)
#删除第0天没有的OTU
library(vegan)
#如皋
otu_N0=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGLN" & sam1$day == "0"), ]))],method = "pa")) >= 2),
           which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGLN"), ]))]#15472
otu_N180=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGMN" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGMN"), ]))]#15719
otu_N360=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGHN" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "RGHN"), ]))]#16759
otu_N0N180=merge(otu_N0,otu_N180,by="row.names",all=T)
otu=merge(otu_N0N180,otu_N360,by.x="Row.names",by.y="row.names",all=T)

rownames(otu) <- otu$Row.names
otu <- otu[,-1]
otu[is.na(otu)] <- 0

write.csv(otu,"./process/filter_OTUtable_rugao.csv")

tax=tax[rownames(otu),]

otu <- as.matrix(otu)
tax <- as.matrix(tax)

dat <- phyloseq(otu_table(otu, taxa_are_rows=T),
                tax_table(tax),
                sample_data(sam1))

#save phyloseq data here and make a new file for qSIP analysis
save(dat, file='process/phyloseq_data_rugao.RData')

load('process/phyloseq_data_rugao.RData')
#specify qsip-related dat_cona from phyloseq object, create phylosip object
qsip_dat <- specify_qsip(dat,
                         density='density_g_ml',
                         abund='qpcr_abund',
                         rep_id='sample_id',
                         rep_num='replicate',
                         rep_group='fertilizer',
                         iso='18O',
                         iso_trt='isotope',
                         timepoint='day')                         

#specify filtering level:
#and in more than 5 of the fractions that make up a replicate
qsip_dat@qsip@filter_levels <- create_filters(replicate = 2, fraction = 5, soft=1)##几个重复，几层中出现#参考ISME方法

qsip_dat1_1 <- calc_excess(qsip_dat,
                           separate_label=T,
                           max_label = 0.98,#
                           filter=T,
                           correction=T)
#extract excess atom fraction (EAF)                           
library(reshape2)
qsip_enrich1 <- dcast(qsmelt(qsip_dat1_1, include=c('excess'), regex=T),
                      tax_id~sample_id+fertilizer,value.var="atom_excess")

#95%CI
qsip_dat1_2 <- calc_excess(qsip_dat,
                           #percent = FALSE,##是否将结果展示为百分比，即EAF值表示为**%，默认为F
                           ci_method = c("bootstrap"),
                           ci = 0.95, iters = 999,
                           #offset_taxa = 0.1,##用于计算偏移校正值的类群的百分比，在correction = T时有效
                           max_label = 0.98,##默认为1，同位素最高掺入量，受实验中提供的原子富集百分比的限制，原子富集量
                           recalc = TRUE,##是否重新计算WAD和分子量值。默认为T，使用bootstrapped将自动重新计算所有值。
                           separate_light = FALSE,##是否将三个重复的WAD取均值，FALSE为取均值
                           separate_label = TRUE,
                           filter = TRUE,##是不是要执行上一步的过滤
                           correction = T)

atom_excess_ci_l <- qsip_dat1_2@qsip[["atom_excess_ci_l"]]
atom_excess <- qsip_dat1_2@qsip[["atom_excess"]]
atom_excess_ci_u <- qsip_dat1_2@qsip[["atom_excess_ci_u"]]

colnames(atom_excess_ci_l)=c("RGHN_ci_l","RGLN_ci_l","RGMN_ci_l")
colnames(atom_excess)=c("RGHN_mean","RGLN_mean","RGMN_mean")
colnames(atom_excess_ci_u)=c("RGHN_ci_u","RGLN_ci_u","RGMN_ci_u")

qsip_enrich2=as.data.frame(cbind(atom_excess_ci_l, atom_excess, atom_excess_ci_u))

atom_excess1=merge(qsip_enrich1,qsip_enrich2,by.x="tax_id",by.y="row.names",all=F)

write.csv(atom_excess1,'process/EAF_rugao.csv')

par(mfrow=c(1,1), mar=c(2,2,1,1)+.1)
hist(atom_excess)

# per-capita growth
qsip_dat2_1 <- calc_pop(qsip_dat ,
                        separate_label=T,
                        max_label = 0.98,
                        filter=T,
                        correction=T,
                        match_replicate=T)##F为利用在time0的丰度，在每组重复中取平均值；T为一一匹配
qsip_rate1 <- dcast(melt(qsmelt(qsip_dat2_1, include=c('rate'), regex=T),ID=c("tax_id")),
                    tax_id~variable+sample_id)

#calculate population growth and death
qsip_dat2_2 <- calc_pop(qsip_dat, 
                        ci_method = c("bootstrap"),
                        ci = 0.95,
                        separate_label=T, 
                        max_label = 0.98,
                        filter=T, 
                        correction=T, 
                        match_replicate=T) 

birth_rate_ci_l <- qsip_dat2_2@qsip[["birth_rate_ci_l"]]
birth_rate <- qsip_dat2_2@qsip[["birth_rate"]]
birth_rate_ci_u <- qsip_dat2_2@qsip[["birth_rate_ci_u"]]                      

death_rate_ci_l <- qsip_dat2_2@qsip[["death_rate_ci_l"]]
death_rate <- qsip_dat2_2@qsip[["death_rate"]]
death_rate_ci_u <- qsip_dat2_2@qsip[["death_rate_ci_u"]]    

growth_rate_ci_l <- qsip_dat2_2@qsip[["growth_rate_ci_l"]]
growth_rate <- qsip_dat2_2@qsip[["growth_rate"]]
growth_rate_ci_u <- qsip_dat2_2@qsip[["growth_rate_ci_u"]]   

colnames(birth_rate_ci_l)=c("RGHN_birth_ci_l","RGLN_birth_ci_l","RGMN_birth_ci_l")
colnames(birth_rate)=c("RGHN_birth_mean","RGLN_birth_mean","RGMN_birth_mean")
colnames(birth_rate_ci_u)=c("RGHN_birth_ci_u","RGLN_birth_ci_u","RGMN_birth_ci_u")

colnames(death_rate_ci_l)=c("RGHN_death_ci_l","RGLN_death_ci_l","RGMN_death_ci_l")
colnames(death_rate)=c("RGHN_death_mean","RGLN_death_mean","RGMN_death_mean")
colnames(death_rate_ci_u)=c("RGHN_death_ci_u","RGLN_death_ci_u","RGMN_death_ci_u")

colnames(growth_rate_ci_l)=c("RGHN_growth_ci_l","RGLN_growth_ci_l","RGMN_growth_ci_l")
colnames(growth_rate)=c("RGHN_growth_mean","RGLN_growth_mean","RGMN_growth_mean")
colnames(growth_rate_ci_u)=c("RGHN_growth_ci_u","RGLN_growth_ci_u","RGMN_growth_ci_u")

qsip_rate2=as.data.frame(cbind(birth_rate_ci_l, birth_rate, birth_rate_ci_u,
                               death_rate_ci_l,death_rate,death_rate_ci_u,
                               growth_rate_ci_l,growth_rate,growth_rate_ci_u))

rates=merge(qsip_rate1,qsip_rate2,by="tax_id",by.x="tax_id",by.y="row.names",all=F)

write.csv(rates,'process/rates_rugao.csv')

par(mfrow=c(1,3), mar=c(2,2,1,1)+.1)
hist(birth_rate); hist(death_rate); hist(growth_rate)

#HB
sam1 <- read.table('./data/samdata.txt', sep='\t', header=T, row.names=1, comment.char='', check.names=F)
sam1$isotope=as.factor(sam1$isotope)

otu <- read.table('./data/pick_otutab.txt', sep='\t', header=T, comment.char='', check.names=F)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]
tax <- read.table('./data/GTDB_taxonomy_adj.txt', sep='\t', header=T, comment.char='', row.names=1)

sam1=sam1[which(sam1$fertilizer %in% c("HBCK","HBLN","HBMN","HBHN")), ]

#单独跑，对OTU进行两步筛选再合并(第0天3个重复至少2个存在数值)
#删除第0天没有的OTU
library(vegan)
#草原
otu_HBCK=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBCK" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBCK"), ]))]#15472
otu_HBLN=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBLN" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBLN"), ]))]#15719
otu_HBMN=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBMN" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBMN"), ]))]#16759
otu_HBHN=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBHN" & sam1$day == "0"), ]))],method = "pa")) >= 3),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "HBHN"), ]))]#6807

otu_CKLN=merge(otu_HBCK,otu_HBLN,by="row.names",all=T)
otu_CKLNMN=merge(otu_CKLN,otu_HBMN,by.x="Row.names",by.y="row.names",all=T)
otu=merge(otu_CKLNMN,otu_HBHN,by.x="Row.names",by.y="row.names",all=T)

rownames(otu) <- otu$Row.names
otu <- otu[,-1]
otu[is.na(otu)] <- 0  
write.csv(otu,"./process/filter_OTUtable_haibei.csv")
tax=tax[rownames(otu),]

otu <- as.matrix(otu)
tax <- as.matrix(tax)

dat <- phyloseq(otu_table(otu, taxa_are_rows=T),
                tax_table(tax),
                sample_data(sam1))

#save phyloseq data here and make a new file for qSIP analysis
save(dat, file='process/phyloseq_data_haibei.RData')

load('process/phyloseq_data_haibei.RData')
#specify qsip-related dat_cona from phyloseq object, create phylosip object
qsip_dat <- specify_qsip(dat,
                         density='density_g_ml',
                         abund='qpcr_abund',
                         rep_id='sample_id',
                         rep_num='replicate',
                         rep_group='fertilizer',
                         iso='18O',
                         iso_trt='isotope',
                         timepoint='day')                         

#specify filtering level:
#for each ecosystem (the grouping factor), taxa must occur in more than 2 of the 3 replicates 
#and in more than 5 of the fractions that make up a replicate
qsip_dat@qsip@filter_levels <- create_filters(replicate = 2, fraction = 5, soft=1)##几个重复，几层中出现#参考ISME方法
qsip_dat1_1 <- calc_excess(qsip_dat,
                           separate_label=T,
                           max_label = 0.45,
                           filter=T,
                           correction=T)
#extract excess atom fraction (EAF)                           
library(reshape2)
qsip_enrich1 <- dcast(qsmelt(qsip_dat1_1, include=c('excess'), regex=T),
                      tax_id~sample_id+fertilizer,value.var="atom_excess")

#95%CI
qsip_dat1_2 <- calc_excess(qsip_dat,
                           #percent = FALSE,##是否将结果展示为百分比，即EAF值表示为**%，默认为F
                           ci_method = c("bootstrap"),
                           ci = 0.95, iters = 999,
                           #offset_taxa = 0.1,##用于计算偏移校正值的类群的百分比，在correction = T时有效
                           max_label = 0.45,##默认为1，同位素最高掺入量，受实验中提供的原子富集百分比的限制，原子富集量
                           recalc = TRUE,##是否重新计算WAD和分子量值。默认为T，使用bootstrapped将自动重新计算所有值。
                           separate_light = FALSE,##是否将三个重复的WAD取均值，FALSE为取均值
                           separate_label = TRUE,
                           filter = TRUE,##是不是要执行上一步的过滤
                           correction = T)

atom_excess_ci_l <- qsip_dat1_2@qsip[["atom_excess_ci_l"]]
atom_excess <- qsip_dat1_2@qsip[["atom_excess"]]
atom_excess_ci_u <- qsip_dat1_2@qsip[["atom_excess_ci_u"]]

colnames(atom_excess_ci_l)=c("HBCK_ci_l","HBHN_ci_l","HBLN_ci_l","HBMN_ci_l")
colnames(atom_excess)=c("HBCK_mean","HBHN_mean","HBLN_mean","HBMN_mean")
colnames(atom_excess_ci_u)=c("HBCK_ci_u","HBHN_ci_u","HBLN_ci_u","HBMN_ci_u")

qsip_enrich2=as.data.frame(cbind(atom_excess_ci_l, atom_excess, atom_excess_ci_u))

atom_excess1=merge(qsip_enrich1,qsip_enrich2,by.x="tax_id",by.y="row.names",all=F)

write.csv(atom_excess1,'process/EAF_haibei.csv')

par(mfrow=c(1,1), mar=c(2,2,1,1)+.1)
hist(atom_excess)

# per-capita growth
qsip_dat2_1 <- calc_pop(qsip_dat ,
                        separate_label=T,
                        max_label = 0.45,
                        filter=T,
                        correction=T,
                        match_replicate=T)##F为利用在time0的丰度，在每组重复中取平均值；T为一一匹配
qsip_rate1 <- dcast(melt(qsmelt(qsip_dat2_1, include=c('rate'), regex=T),ID=c("tax_id")),
                    tax_id~variable+sample_id)

#calculate population growth and death
qsip_dat2_2 <- calc_pop(qsip_dat, 
                        ci_method = c("bootstrap"),
                        ci = 0.95,
                        separate_label=T, 
                        max_label = 0.45,
                        filter=T, 
                        correction=T, 
                        match_replicate=T) 

birth_rate_ci_l <- qsip_dat2_2@qsip[["birth_rate_ci_l"]]
birth_rate <- qsip_dat2_2@qsip[["birth_rate"]]
birth_rate_ci_u <- qsip_dat2_2@qsip[["birth_rate_ci_u"]]                      

death_rate_ci_l <- qsip_dat2_2@qsip[["death_rate_ci_l"]]
death_rate <- qsip_dat2_2@qsip[["death_rate"]]
death_rate_ci_u <- qsip_dat2_2@qsip[["death_rate_ci_u"]]    

growth_rate_ci_l <- qsip_dat2_2@qsip[["growth_rate_ci_l"]]
growth_rate <- qsip_dat2_2@qsip[["growth_rate"]]
growth_rate_ci_u <- qsip_dat2_2@qsip[["growth_rate_ci_u"]]   

colnames(birth_rate_ci_l)=c("HBCK_birth_ci_l","HBHN_birth_ci_l","HBLN_birth_ci_l","HBMN_birth_ci_l")
colnames(birth_rate)=c("HBCK_birth_mean","HBHN_birth_mean","HBLN_birth_mean","HBMN_birth_mean")
colnames(birth_rate_ci_u)=c("HBCK_birth_ci_u","HBHN_birth_ci_u","HBLN_birth_ci_u","HBMN_birth_ci_u")

colnames(death_rate_ci_l)=c("HBCK_death_ci_l","HBHN_death_ci_l","HBLN_death_ci_l","HBMN_death_ci_l")
colnames(death_rate)=c("HBCK_death_mean","HBHN_death_mean","HBLN_death_mean","HBMN_death_mean")
colnames(death_rate_ci_u)=c("HBCK_death_ci_u","HBHN_death_ci_u","HBLN_death_ci_u","HBMN_death_ci_u")

colnames(growth_rate_ci_l)=c("HBCK_growth_ci_l","HBHN_growth_ci_l","HBLN_growth_ci_l","HBMN_growth_ci_l")
colnames(growth_rate)=c("HBCK_growth_mean","HBHN_growth_mean","HBLN_growth_mean","HBMN_growth_mean")
colnames(growth_rate_ci_u)=c("HBCK_growth_ci_u","HBHN_growth_ci_u","HBLN_growth_ci_u","HBMN_growth_ci_u")

qsip_rate2=as.data.frame(cbind(birth_rate_ci_l, birth_rate, birth_rate_ci_u,
                               death_rate_ci_l,death_rate,death_rate_ci_u,
                               growth_rate_ci_l,growth_rate,growth_rate_ci_u))

rates=merge(qsip_rate1,qsip_rate2,by="tax_id",by.x="tax_id",by.y="row.names",all=F)

write.csv(rates,'process/rates_haibei.csv')
par(mfrow=c(1,3), mar=c(2,2,1,1)+.1)
hist(birth_rate); hist(death_rate); hist(growth_rate)

#XJ
sam1 <- read.table('./data/samdata.txt', sep='\t', header=T, row.names=1, comment.char='', check.names=F)
sam1$isotope=as.factor(sam1$isotope)

otu <- read.table('./data/pick_otutab.txt', sep='\t', header=T, comment.char='', check.names=F)
rownames(otu) <- otu$`#OTU ID`
otu <- otu[,-1]
tax <- read.table('./data/GTDB_taxonomy_adj.txt', sep='\t', header=T, comment.char='', row.names=1)

sam1=sam1[which(sam1$fertilizer %in% c("XJCK","XJNPK","XJNPKHS","XJNPKLS","XJNPKM")), ]

#单独跑，对OTU进行两步筛选再合并(第0天3个重复至少2个存在数值)
#删除第0天没有的OTU
library(vegan)
#新疆
otu_XJCK=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJCK" & sam1$day == "0"), ]))],method = "pa")) >= 2),
             which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJCK"), ]))]#15472
otu_XJNPK=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPK" & sam1$day == "0"), ]))],method = "pa")) >= 2),
              which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPK"), ]))]#15719
otu_XJNPKHS=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKHS" & sam1$day == "0"), ]))],method = "pa")) >= 2),
                which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKHS"), ]))]#16759
otu_XJNPKLS=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKLS" & sam1$day == "0"), ]))],method = "pa")) >= 3),
                which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKLS"), ]))]#6807
otu_XJNPKM=otu[which(rowSums(decostand(otu[,which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKM" & sam1$day == "0"), ]))],method = "pa")) >= 3),
               which(colnames(otu) %in% rownames(sam1[which(sam1$fertilizer == "XJNPKM"), ]))]#6807

otu_CKNPK=merge(otu_XJCK,otu_XJNPK,by="row.names",all=T)
otu_CKNPKHS=merge(otu_CKNPK,otu_XJNPKHS,by.x="Row.names",by.y="row.names",all=T)
otu_CKNPKHSLS=merge(otu_CKNPKHS,otu_XJNPKLS,by.x="Row.names",by.y="row.names",all=T)
otu=merge(otu_CKNPKHSLS,otu_XJNPKM,by.x="Row.names",by.y="row.names",all=T)

rownames(otu) <- otu$Row.names
otu <- otu[,-1]
otu[is.na(otu)] <- 0  
write.csv(otu,"./process/filter_OTUtable_xinjiang.csv")
tax=tax[rownames(otu),]
otu <- as.matrix(otu)
tax <- as.matrix(tax)

dat <- phyloseq(otu_table(otu, taxa_are_rows=T),
                tax_table(tax),
                sample_data(sam1))

#save phyloseq data here and make a new file for qSIP analysis
save(dat, file='process/phyloseq_data_xinjiang.RData')

load('process/phyloseq_data_xinjiang.RData')
#specify qsip-related dat_cona from phyloseq object, create phylosip object
qsip_dat <- specify_qsip(dat,
                         density='density_g_ml',
                         abund='qpcr_abund',
                         rep_id='sample_id',
                         rep_num='replicate',
                         rep_group='fertilizer',
                         iso='18O',
                         iso_trt='isotope',
                         timepoint='day')                         

#specify filtering level:
#for each ecosystem (the grouping factor), taxa must occur in more than 2 of the 3 replicates 
#and in more than 5 of the fractions that make up a replicate
qsip_dat@qsip@filter_levels <- create_filters(replicate = 2, fraction = 5, soft=1)##几个重复，几层中出现#参考ISME方法
qsip_dat1_1 <- calc_excess(qsip_dat,
                           separate_label=T,
                           max_label = 0.98,
                           filter=T,
                           correction=T)
#extract excess atom fraction (EAF)                           
library(reshape2)
qsip_enrich1 <- dcast(qsmelt(qsip_dat1_1, include=c('excess'), regex=T),
                      tax_id~sample_id+fertilizer,value.var="atom_excess")

#95%CI
qsip_dat1_2 <- calc_excess(qsip_dat,
                           #percent = FALSE,##是否将结果展示为百分比，即EAF值表示为**%，默认为F
                           ci_method = c("bootstrap"),
                           ci = 0.95, iters = 999,
                           #offset_taxa = 0.1,##用于计算偏移校正值的类群的百分比，在correction = T时有效
                           max_label = 0.98,##默认为1，同位素最高掺入量，受实验中提供的原子富集百分比的限制，原子富集量
                           recalc = TRUE,##是否重新计算WAD和分子量值。默认为T，使用bootstrapped将自动重新计算所有值。
                           separate_light = FALSE,##是否将三个重复的WAD取均值，FALSE为取均值
                           separate_label = TRUE,
                           filter = TRUE,##是不是要执行上一步的过滤
                           correction = T)

atom_excess_ci_l <- qsip_dat1_2@qsip[["atom_excess_ci_l"]]
atom_excess <- qsip_dat1_2@qsip[["atom_excess"]]
atom_excess_ci_u <- qsip_dat1_2@qsip[["atom_excess_ci_u"]]

colnames(atom_excess_ci_l)=c("XJCK_ci_l","XJNPK_ci_l","XJNPKHS_ci_l","XJNPKLS_ci_l","XJNPKM_ci_l")
colnames(atom_excess)=c("XJCK_mean","XJNPK_mean","XJNPKHS_mean","XJNPKLS_mean","XJNPKM_mean")
colnames(atom_excess_ci_u)=c("XJCK_ci_u","XJNPK_ci_u","XJNPKHS_ci_u","XJNPKLS_ci_u","XJNPKM_ci_u")

qsip_enrich2=as.data.frame(cbind(atom_excess_ci_l, atom_excess, atom_excess_ci_u))

atom_excess1=merge(qsip_enrich1,qsip_enrich2,by.x="tax_id",by.y="row.names",all=F)

write.csv(atom_excess1,'process/EAF_xinjiang.csv')

par(mfrow=c(1,1), mar=c(2,2,1,1)+.1)
hist(atom_excess)

# per-capita growth
qsip_dat2_1 <- calc_pop(qsip_dat ,
                        separate_label=T,
                        max_label = 0.98,
                        filter=T,
                        correction=T,
                        match_replicate=T)##F为利用在time0的丰度，在每组重复中取平均值；T为一一匹配
qsip_rate1 <- dcast(melt(qsmelt(qsip_dat2_1, include=c('rate'), regex=T),ID=c("tax_id")),
                    tax_id~variable+sample_id)

#calculate population growth and death
qsip_dat2_2 <- calc_pop(qsip_dat, 
                        ci_method = c("bootstrap"),
                        ci = 0.95,
                        separate_label=T, 
                        max_label = 0.98,
                        filter=T, 
                        correction=T, 
                        match_replicate=T) 

birth_rate_ci_l <- qsip_dat2_2@qsip[["birth_rate_ci_l"]]
birth_rate <- qsip_dat2_2@qsip[["birth_rate"]]
birth_rate_ci_u <- qsip_dat2_2@qsip[["birth_rate_ci_u"]]                      

death_rate_ci_l <- qsip_dat2_2@qsip[["death_rate_ci_l"]]
death_rate <- qsip_dat2_2@qsip[["death_rate"]]
death_rate_ci_u <- qsip_dat2_2@qsip[["death_rate_ci_u"]]    

growth_rate_ci_l <- qsip_dat2_2@qsip[["growth_rate_ci_l"]]
growth_rate <- qsip_dat2_2@qsip[["growth_rate"]]
growth_rate_ci_u <- qsip_dat2_2@qsip[["growth_rate_ci_u"]]   

colnames(birth_rate_ci_l)=c("XJCK_birth_ci_l","XJNPK_birth_ci_l","XJNPKHS_birth_ci_l","XJNPKLS_birth_ci_l","XJNPKM_birth_ci_l")
colnames(birth_rate)=c("XJCK_birth_mean","XJNPK_birth_mean","XJNPKHS_birth_mean","XJNPKLS_birth_mean","XJNPKM_birth_mean")
colnames(birth_rate_ci_u)=c("XJCK_birth_ci_u","XJNPK_birth_ci_u","XJNPKHS_birth_ci_u","XJNPKLS_birth_ci_u","XJNPKM_birth_ci_u")

colnames(death_rate_ci_l)=c("XJCK_death_ci_l","XJNPK_death_ci_l","XJNPKHS_death_ci_l","XJNPKLS_death_ci_l","XJNPKM_death_ci_l")
colnames(death_rate)=c("XJCK_death_mean","XJNPK_death_mean","XJNPKHS_death_mean","XJNPKLS_death_mean","XJNPKM_death_mean")
colnames(death_rate_ci_u)=c("XJCK_death_ci_u","XJNPK_death_ci_u","XJNPKHS_death_ci_u","XJNPKLS_death_ci_u","XJNPKM_death_ci_u")

colnames(growth_rate_ci_l)=c("XJCK_growth_ci_l","XJNPK_growth_ci_l","XJNPKHS_growth_ci_l","XJNPKLS_growth_ci_l","XJNPKM_growth_ci_l")
colnames(growth_rate)=c("XJCK_growth_mean","XJNPK_growth_mean","XJNPKHS_growth_mean","XJNPKLS_growth_mean","XJNPKM_growth_mean")
colnames(growth_rate_ci_u)=c("XJCK_growth_ci_u","XJNPK_growth_ci_u","XJNPKHS_growth_ci_u","XJNPKLS_growth_ci_u","XJNPKM_growth_ci_u")

qsip_rate2=as.data.frame(cbind(birth_rate_ci_l, birth_rate, birth_rate_ci_u,
                               death_rate_ci_l,death_rate,death_rate_ci_u,
                               growth_rate_ci_l,growth_rate,growth_rate_ci_u))

rates=merge(qsip_rate1,qsip_rate2,by="tax_id",by.x="tax_id",by.y="row.names",all=F)

write.csv(rates,'process/rates_xinjiang.csv')
par(mfrow=c(1,3), mar=c(2,2,1,1)+.1)
hist(birth_rate); hist(death_rate); hist(growth_rate)






