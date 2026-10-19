#pH适应机制分成几个功能模块

rm(list = ls())

library(data.table)
library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(ape)
library(patchwork)

set.seed(123)
input_dir <- "./Process"
analysis_dir <- "./Process"
figure_dir <- "./Fig"
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. 读取元数据并删除古菌
# ============================================================

pH_order <- c("Acidophile-associated MAG", "Acid–neutrophile MAG", "Neutrophile-associated MAG", "Neutrophile–alkaliphile MAG", "Alkaliphile-associated MAG", "Broad pH-range MAG")
color_dict <- c("Acidophile-associated MAG" = "#D55E00", "Acid–neutrophile MAG" = "#BDC609", "Neutrophile-associated MAG" = "#009E73", "Neutrophile–alkaliphile MAG" = "#56B4E9", "Alkaliphile-associated MAG" = "#0072B2", "Broad pH-range MAG" = "#CC79A7")

metadata <- fread(file.path(input_dir, "Ecotypes_MAGtrait_group_Bac.csv"), check.names = FALSE, data.table = FALSE)
metadata <- metadata %>% select(-any_of(c("Unnamed: 0", "...1", "X","V"))) %>% mutate(MAG = trimws(MAG), Domain = trimws(Domain), pH_group = trimws(pH_group))

if (anyDuplicated(metadata$MAG)) stop("元数据中存在重复MAG。")
if (any(is.na(metadata$pH_group))) stop("部分细菌MAG没有正确匹配到六种pH生态型。")
if (nrow(metadata) != 43) warning("删除古菌后，细菌MAG数量不是43，请检查元数据。")

write.csv(metadata, file.path(analysis_dir, "Bacterial_MAG_metadata.csv"), row.names = FALSE)
saveRDS(metadata, file.path(analysis_dir, "Bacterial_MAG_metadata.rds"))


#pH稳态机制
#把两类证据联合起来：
#关键基因是否存在；
#相关通路或模块完整性是否提高。
#重点包括：
#F-type、V-type和A-type ATPase；
#Na⁺/H⁺反向转运系统；
#K⁺稳态；
#氨基酸脱羧；
#精氨酸脱亚氨基途径；
#尿素利用；
#谷氨酸和谷氨酰胺代谢；
#氧化胁迫；
#膜与细胞壁合成；
#渗透保护物质。
#数据读取
rm(list = ls())

pkgs <- c("data.table","dplyr","tidyr","tibble","stringr","ggplot2","patchwork","scales","purrr")
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

pH_order <- c("Acidophile-associated MAG","Acid–neutrophile MAG","Neutrophile-associated MAG","Neutrophile–alkaliphile MAG",
              "Alkaliphile-associated MAG","Broad pH-range MAG")
pH_short <- c("Acidophile-associated MAG" = "Acid","Acid–neutrophile MAG" = "A-Neu","Neutrophile-associated MAG" = "Neu",
              "Neutrophile–alkaliphile MAG" = "N-Alk","Alkaliphile-associated MAG" = "Alk","Broad pH-range MAG" = "Broad")

category_order <- c("Proton-translocating ATPases","Na+/H+ antiport systems",
                    "Amino-acid decarboxylation systems","Alkali-generating systems")

category_colors <- c("Proton-translocating ATPases"="#377EB8",
                     "Na+/H+ antiport systems"="#1B9E77",
                     "Amino-acid decarboxylation systems"="#D95F02",
                     "Alkali-generating systems"="#984EA3")

metadata <- readRDS(metadata_file) %>% 
  mutate(MAG = trimws(as.character(MAG)), pH_group = factor(trimws(as.character(pH_group)), levels = pH_order)) %>% 
  arrange(pH_group, MAG)
if (anyDuplicated(metadata$MAG)) stop("元数据中存在重复MAG。")
if (any(is.na(metadata$pH_group))) stop("存在无法匹配六种pH生态型的MAG。")

mag_ids <- metadata$MAG
cat("参与分析的细菌MAG数量：", length(mag_ids), "\n")
print(table(metadata$pH_group))

# ============================================================
# 1. 定义pH稳态机制、替代路线、功能组件和KO
# 同一Component内的多个KO表示替代KO；同一Route内的不同Component均需存在
# ============================================================

R <- function(category,mechanism,mechanism_order,route,component_order,component,kos,interpretation,rule="AND") {
  tibble(Category=category,
         Mechanism=mechanism,
         Mechanism_order=mechanism_order,
         Route=route,
         Component_order=component_order,
         Component=component,
         KO_order=seq_along(kos),
         KO=kos,
         Interpretation=interpretation,
         Rule=rule)
}

panel <- bind_rows(
  
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",1,"Subunit a","K02108","Membrane proton-conducting subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",2,"Subunit b","K02109","Membrane stator subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",3,"Subunit c","K02110","Membrane proton-conducting rotor subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",4,"Subunit alpha","K02111","Catalytic-sector structural subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",5,"Subunit beta","K02112","Catalytic ATPase subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",6,"Subunit delta","K02113","Peripheral stalk subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",7,"Subunit epsilon","K02114","Rotary regulatory subunit","coverage"),
  R("Proton-translocating ATPases","F-type H+-transporting ATP synthase",1,"F-type H+-transporting ATP synthase",8,"Subunit gamma","K02115","Central stalk subunit","coverage"),
  
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",1,"Subunit A","K02117","Catalytic subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",2,"Subunit B","K02118","Non-catalytic nucleotide-binding subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",3,"Subunit C","K02119","Stator subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",4,"Subunit D","K02120","Central stalk subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",5,"Subunit E","K02121","Peripheral stalk subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",6,"Subunit F","K02122","Central stalk subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",7,"Subunit G/H","K02107","Peripheral stalk subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",8,"Subunit I","K02123","Membrane ion-translocating subunit","coverage"),
  R("Proton-translocating ATPases","V/A-type H+-translocating ATPase",2,"V/A-type H+-translocating ATPase",9,"Subunit K","K02124","Membrane rotor subunit","coverage"),
  
  R("Na+/H+ antiport systems","NhaA-type Na+/H+ antiporter",3,"NhaA-type Na+/H+ antiporter",1,"NhaA","K03313","Na+/H+ antiporter, NhaA family","single"),
  R("Na+/H+ antiport systems","NhaB-type Na+/H+ antiporter",4,"NhaB-type Na+/H+ antiporter",1,"NhaB","K03314","Na+/H+ antiporter, NhaB family","single"),
  R("Na+/H+ antiport systems","NhaC-type Na+/H+ antiporter",5,"NhaC-type Na+/H+ antiporter",1,"NhaC","K03315","Na+/H+ antiporter, NhaC family","single"),
  R("Na+/H+ antiport systems","NhaP-type Na+/H+ antiporter",6,"NhaP-type Na+/H+ antiporter",1,"NhaP","K03316","NhaP-type Na+/H+ antiporter","single"),
  R("Na+/H+ antiport systems","NhaS-type Na+/H+ antiporter",7,"NhaS-type Na+/H+ antiporter",1,"NhaS","K26730","NhaS-type monovalent cation/H+ antiporter","single"),
  
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",1,"MrpA","K05565","Multicomponent Na+/H+ antiporter subunit A","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",2,"MrpB","K05566","Multicomponent Na+/H+ antiporter subunit B","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",3,"MrpC","K05567","Multicomponent Na+/H+ antiporter subunit C","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",4,"MrpD","K05568","Multicomponent Na+/H+ antiporter subunit D","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",5,"MrpE","K05569","Multicomponent Na+/H+ antiporter subunit E","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",6,"MrpF","K05570","Multicomponent Na+/H+ antiporter subunit F","coverage"),
  R("Na+/H+ antiport systems","Mrp Na+/H+ antiporter complex",8,"Mrp Na+/H+ antiporter complex",7,"MrpG","K05571","Multicomponent Na+/H+ antiporter subunit G","coverage"),
  
  R("Amino-acid decarboxylation systems","Glutamate decarboxylase (GAD) system",9,"GadAB-GadC",1,"Glutamate decarboxylase","K01580","GadA/GadB glutamate decarboxylase","AND"),
  R("Amino-acid decarboxylation systems","Glutamate decarboxylase (GAD) system",9,"GadAB-GadC",2,"Glutamate/GABA antiporter","K20265","GadC glutamate/GABA antiporter","AND"),
  
  R("Amino-acid decarboxylation systems","Arginine decarboxylase (ADC) system",10,"AdiA-AdiC",1,"Arginine decarboxylase","K01583","Alternative arginine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Arginine decarboxylase (ADC) system",10,"AdiA-AdiC",1,"Arginine decarboxylase","K01584","Alternative arginine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Arginine decarboxylase (ADC) system",10,"AdiA-AdiC",1,"Arginine decarboxylase","K01585","Alternative arginine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Arginine decarboxylase (ADC) system",10,"AdiA-AdiC",1,"Arginine decarboxylase","K02626","Alternative arginine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Arginine decarboxylase (ADC) system",10,"AdiA-AdiC",2,"Arginine/agmatine antiporter","K03759","AdiC arginine/agmatine antiporter","AND"),
  
  R("Amino-acid decarboxylation systems","Lysine decarboxylase (LDC) system",11,"CadA-CadB",1,"Lysine decarboxylase","K01582","Alternative lysine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Lysine decarboxylase (LDC) system",11,"CadA-CadB",1,"Lysine decarboxylase","K28682","Alternative lysine decarboxylase","OR"),
  R("Amino-acid decarboxylation systems","Lysine decarboxylase (LDC) system",11,"CadA-CadB",2,"Lysine/cadaverine antiporter","K03757","CadB lysine/cadaverine antiporter","AND"),
  
  R("Amino-acid decarboxylation systems","Ornithine decarboxylase system",12,"SpeF-PotE",1,"Ornithine decarboxylase","K01581","SpeF ornithine decarboxylase","AND"),
  R("Amino-acid decarboxylation systems","Ornithine decarboxylase system",12,"SpeF-PotE",2,"Ornithine/putrescine antiporter","K03756","PotE ornithine/putrescine antiporter","AND"),
  
  R("Alkali-generating systems","Arginine deiminase system",13,"Arc pathway",1,"Arginine deiminase","K01478","Arginine deiminase","AND"),
  R("Alkali-generating systems","Arginine deiminase system",13,"Arc pathway",2,"Ornithine carbamoyltransferase","K00611","Ornithine carbamoyltransferase family","AND"),
  R("Alkali-generating systems","Arginine deiminase system",13,"Arc pathway",3,"Carbamate kinase","K00926","Carbamate kinase","AND"),
  R("Alkali-generating systems","Arginine deiminase system",13,"Arc pathway",4,"Arginine/ornithine antiporter","K03758","Arginine/ornithine antiporter","AND"),
  
  R("Alkali-generating systems","Urease system",14,"Three-subunit urease",1,"Urease alpha subunit","K01428","UreC","AND"),
  R("Alkali-generating systems","Urease system",14,"Three-subunit urease",2,"Urease beta subunit","K01429","UreB","AND"),
  R("Alkali-generating systems","Urease system",14,"Three-subunit urease",3,"Urease gamma subunit","K01430","UreA","AND"),
  R("Alkali-generating systems","Urease system",14,"Single-chain urease",1,"Single-chain urease","K01427","Single-polypeptide urease alternative","OR")
  
) %>% mutate(Category=factor(Category,levels=category_order)) %>% 
  arrange(Mechanism_order,Route,Component_order,KO_order)

#用KEGG描述表核对KO名称
# ============================================================
# 2. 用KEGG描述表核对KO名称
# ============================================================

mapping_raw <- fread(mapping_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE, colClasses = "character")
if (!all(c("KO_ID","KO_Description") %in% names(mapping_raw))) stop("KEGG_KO_Tab_Table.txt缺少KO_ID或KO_Description列。")

ko_info <- mapping_raw %>% transmute(KO = trimws(KO_ID), KO_Description = trimws(KO_Description)) %>% filter(!is.na(KO), KO != "") %>% distinct(KO, .keep_all = TRUE)

panel <- panel %>% left_join(ko_info, by = "KO") %>% mutate(KO_name = str_trim(str_extract(KO_Description, "^[^;]+")), KO_name = ifelse(is.na(KO_name) | KO_name == "", Component, KO_name), KO_label = paste0(KO, ": ", str_trunc(KO_name, 30)))

write.csv(panel, file.path(analysis_dir, "pH_homeostasis_KO_definition.csv"), row.names = FALSE)

# ============================================================
# 3. 读取MAG × KO 0/1矩阵
# ============================================================

ko_dat <- fread(ko_file, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
colnames(ko_dat)[1] <- "KO"
ko_dat$KO <- trimws(as.character(ko_dat$KO))

if (any(is.na(ko_dat$KO) | ko_dat$KO == "")) stop("KO矩阵中存在空KO编号。")
missing_MAG <- setdiff(mag_ids, colnames(ko_dat))
if (length(missing_MAG) > 0) stop("KO矩阵缺少MAG：", paste(missing_MAG, collapse = ", "))

if (anyDuplicated(ko_dat$KO)) ko_dat <- ko_dat %>% group_by(KO) %>% summarise(across(all_of(mag_ids), ~max(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

original_kos <- unique(ko_dat$KO)
ko_mat <- t(as.matrix(ko_dat[, mag_ids, drop = FALSE]))
storage.mode(ko_mat) <- "numeric"
rownames(ko_mat) <- mag_ids
colnames(ko_mat) <- ko_dat$KO

if (anyNA(ko_mat)) stop("KO矩阵存在NA。")
if (any(!ko_mat %in% c(0,1))) stop("KO矩阵存在0和1以外的数值。")

panel_kos <- unique(panel$KO)
missing_panel_kos <- setdiff(panel_kos, colnames(ko_mat))

if (length(missing_panel_kos) > 0) {
  zero_mat <- matrix(0, nrow = nrow(ko_mat), ncol = length(missing_panel_kos), dimnames = list(rownames(ko_mat), missing_panel_kos))
  ko_mat <- cbind(ko_mat, zero_mat)
}

write.csv(data.frame(KO = panel_kos, KO_in_original_matrix = panel_kos %in% original_kos), 
          file.path(analysis_dir, "pH_homeostasis_KO_matrix_check.csv"), row.names = FALSE)

#计算KO、组件、路线和机制层面的完整度
# ============================================================
# 4. 计算KO、组件、路线和机制层面的完整度
# ============================================================
sub_mat <- ko_mat[mag_ids, panel_kos, drop = FALSE]

ko_presence_base <- as.data.frame(sub_mat) %>% rownames_to_column("MAG") %>% pivot_longer(-MAG, names_to = "KO", values_to = "KO_presence") %>% left_join(metadata %>% select(MAG, pH_group), by = "MAG")

ko_presence <- ko_presence_base %>% left_join(panel, by = "KO", relationship = "many-to-many") %>% arrange(Mechanism_order, Route, Component_order, KO_order, pH_group, MAG)

component_mag <- ko_presence %>%
  group_by(MAG, pH_group, Category, Mechanism, Mechanism_order, Route, Component, Component_order) %>%
  summarise(Component_presence = max(KO_presence), KO_detected = paste(KO[KO_presence > 0], collapse = ";"), .groups = "drop")

route_mag <- component_mag %>% group_by(MAG,pH_group,Category,Mechanism,Mechanism_order,Route) %>%
  summarise(n_components=n(),n_components_present=sum(Component_presence),
            Route_coverage=mean(Component_presence),Route_complete=all(Component_presence==1),
            .groups="drop")
##route_mag = 所有路线层面 (route_mag 是所有可能路线的信息)
##system_mag = 每个系统选最佳路线层面 (每个系统选择一个“最佳路线”后的代表完整度)

system_mag <- route_mag %>% group_by(MAG,pH_group,Category,Mechanism,Mechanism_order) %>%
  arrange(desc(Route_coverage),Route) %>%
  slice(1) %>%
  ungroup() %>%
  rename(Best_route = Route,System_coverage = Route_coverage,System_complete = Route_complete) %>%
  mutate(High_coverage = System_coverage >= 0.75)

write.csv(ko_presence, file.path(analysis_dir, "pH_homeostasis_KO_MAG_presence.csv"), row.names = FALSE)
write.csv(component_mag, file.path(analysis_dir, "pH_homeostasis_component_MAG_presence.csv"), row.names = FALSE)
write.csv(route_mag, file.path(analysis_dir, "pH_homeostasis_route_MAG_coverage.csv"), row.names = FALSE)
write.csv(system_mag, file.path(analysis_dir, "pH_homeostasis_system_MAG_coverage.csv"), row.names = FALSE)

#六组描述性汇总
# ============================================================
# 5. 六组描述性汇总
# ============================================================
system_group_summary <- system_mag %>%
  group_by(Category, Mechanism, Mechanism_order, pH_group) %>%
  summarise(n_MAG = n(), Mean_coverage = mean(System_coverage), Median_coverage = median(System_coverage), Complete_n = sum(System_complete), Complete_prevalence = mean(System_complete), High_coverage_n = sum(High_coverage), High_coverage_prevalence = mean(High_coverage), .groups = "drop") %>%
  arrange(Mechanism_order, pH_group)

component_group_summary <- component_mag %>%
  group_by(Category, Mechanism, Mechanism_order, Route, Component, Component_order, pH_group) %>%
  summarise(n_MAG = n(), n_present = sum(Component_presence), Prevalence = mean(Component_presence), .groups = "drop") %>%
  arrange(Mechanism_order, Component_order, pH_group)

write.csv(system_group_summary, file.path(analysis_dir, "pH_homeostasis_system_group_summary.csv"), row.names = FALSE)
write.csv(component_group_summary, file.path(analysis_dir, "pH_homeostasis_component_group_summary.csv"), row.names = FALSE)

#根据正态性和方差齐性选择统计方法
# ============================================================
# pH-homeostasis system statistics
# ============================================================
library(dplyr)
source("./data/boxplerk.R")
source("./data/boxplert4.R")

pairwise_results <- list()
analyze_homeostasis <- function(dat){
  dat <- dat %>% filter(!is.na(System_coverage))
  Category <- unique(dat$Category)
  Mechanism <- unique(dat$Mechanism)
  Mechanism_order <- unique(dat$Mechanism_order)
  
  if(length(unique(dat$System_coverage)) <= 1){
    return(tibble(Category=Category,Mechanism=Mechanism,Mechanism_order=Mechanism_order,Method="No_variation",Shapiro_P=NA,Bartlett_P=NA,P_value=NA))
  }
  
  sub_data <- data.frame(System_coverage=dat$System_coverage,pH_group=dat$pH_group)
  shapiro_p <- tryCatch(shapiro.test(sub_data$System_coverage)$p.value,error=function(e) NA)
  bartlett_p <- tryCatch(bartlett.test(System_coverage~pH_group,data=sub_data)$p.value,error=function(e) NA)
  
  if(!is.na(shapiro_p) & !is.na(bartlett_p) & shapiro_p>0.05 & bartlett_p>0.05){
    #参数检验
    a <- tryCatch(boxplert(sub_data$System_coverage,sub_data$pH_group,
                           ylab=Mechanism,xlab="pH group",bcol="bisque",p.adj="fdr",las=1),error=function(e) NULL)
    if(is.null(a)){
      return(tibble(Category=Category,Mechanism=Mechanism,Mechanism_order=Mechanism_order,
                    Method="boxplert_failed",Shapiro_P=shapiro_p,Bartlett_P=bartlett_p,P_value=NA))
    }
    stat <- a$comparison
    stat$p_value <- a$p.value
    stat$shapiro_p <- a$shapiro_p
    stat$bartlett_p <- a$bartlett_p
    stat$treatment <- rownames(stat)
    rownames(stat)<-NULL
    stat$Method <- "ANOVA_LSD_FDR"
  }else{
    #非参检验
    a <- tryCatch(boxplerk(sub_data$System_coverage,sub_data$pH_group,
                           ylab=Mechanism,xlab="pH group",bcol="bisque",p.adj="fdr",las=1),error=function(e) NULL)
    if(is.null(a)){
      return(tibble(Category=Category,Mechanism=Mechanism,Mechanism_order=Mechanism_order,
                    Method="boxplerk_failed",Shapiro_P=shapiro_p,Bartlett_P=bartlett_p,P_value=NA))
    }
    stat <- a$comparison
    stat$p_value <- a$p.value
    stat$treatment <- rownames(stat)
    rownames(stat)<-NULL
    stat$Method <- "Kruskal_Wallis_FDR"
  }
  
  stat$Category <- Category
  stat$Mechanism <- Mechanism
  stat$Mechanism_order <- Mechanism_order
  pairwise_results[[paste(Category,Mechanism,sep="_")]] <<- stat
  return(tibble(Category=Category,Mechanism=Mechanism,Mechanism_order=Mechanism_order,Method=unique(stat$Method),
                Shapiro_P=shapiro_p,Bartlett_P=bartlett_p,P_value=unique(stat$p_value)[1]))
}

statistics_summary <- system_mag %>% group_by(Category,Mechanism,Mechanism_order) %>% 
  group_split() %>% lapply(analyze_homeostasis) %>% bind_rows()
pairwise_summary <- bind_rows(pairwise_results)

write.csv(statistics_summary,file.path(analysis_dir,"pH_homeostasis_global_statistics.csv"),row.names=FALSE)
write.csv(pairwise_summary,file.path(analysis_dir,"pH_homeostasis_pairwise_FDR_results.csv"),row.names=FALSE)

print(statistics_summary)
print(pairwise_summary)
#系统层面气泡图
# ============================================================
#  系统层面气泡图
# 颜色为平均组件覆盖度，大小为≥75%覆盖MAG比例
# ============================================================
mechanism_levels <- system_group_summary %>% distinct(Mechanism, Mechanism_order) %>% arrange(Mechanism_order) %>% pull(Mechanism)

bubble_dat <- system_group_summary %>% mutate(Mechanism = factor(Mechanism, levels = rev(mechanism_levels)), Category = factor(Category, levels = category_order))

p_bubble <- ggplot(bubble_dat, aes(x = pH_group, y = Mechanism)) +
  geom_point(aes(size = High_coverage_prevalence, fill = Mean_coverage), shape = 21, color = "grey25", stroke = 0.35) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(labels = pH_short, drop = FALSE) +
  scale_size_continuous(range = c(1.5, 8), limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = percent_format(accuracy = 1), name = "MAG prevalence\n(coverage ≥75%)") +
  scale_fill_gradient2(low = "#A6611A", mid = "#F5F5F5", high = "#018571", midpoint = 0.5, limits = c(0,1), name = "Mean system\ncoverage") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.placement = "outside", strip.background = element_rect(fill = "grey95", color = "grey50"), strip.text.y.left = element_text(angle = 0, face = "bold", size = 8), axis.text.x = element_text(angle = 35, hjust = 1, size = 8), axis.text.y = element_text(size = 8), legend.position = "right")
p_bubble

ggsave(file.path(figure_dir, "pH_homeostasis_system_bubble.pdf"), p_bubble, width = 180/25.4, height = 145/25.4, units = "in", useDingbats = FALSE)

#组件层面正方形热图
# ============================================================
# 10. 组件层面正方形热图
# ============================================================
#component_detected <- component_mag %>% group_by(Mechanism, Component) %>% summarise(Detected = any(Component_presence > 0), .groups = "drop")
#component_plot_dat <- component_group_summary %>% semi_join(component_detected %>% filter(Detected), by = c("Mechanism","Component"))
#mechanisms_used <- mechanism_levels[mechanism_levels %in% unique(component_plot_dat$Mechanism)]
component_plot_dat <- component_group_summary
mechanisms_used <- mechanism_levels

component_plot_list <- vector("list", length(mechanisms_used))
component_heights <- numeric(length(mechanisms_used))

for (i in seq_along(mechanisms_used)) {
  mechanism_name <- mechanisms_used[i]
  sub_dat <- component_plot_dat %>% filter(Mechanism == mechanism_name) %>% arrange(Component_order)
  component_order <- sub_dat %>% distinct(Component_order, Component) %>% arrange(Component_order) %>% pull(Component)
  sub_dat$Component <- factor(sub_dat$Component, levels = rev(component_order))
  category_name <- as.character(first(sub_dat$Category))
  base_col <- unname(category_colors[category_name])
  component_heights[i] <- n_distinct(sub_dat$Component)
  
  p_i <- ggplot(sub_dat, aes(x = pH_group, y = Component, alpha = Prevalence)) +
    geom_tile(fill = base_col, color = "grey25", linewidth = 0.3, width = 0.95, height = 0.95) +
    scale_x_discrete(labels = pH_short, limits = pH_order, drop = FALSE, expand = c(0,0)) +
    scale_y_discrete(drop = FALSE, expand = c(0,0)) +
    scale_alpha_continuous(limits = c(0,1), range = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = percent_format(accuracy = 1), name = "MAG prevalence") +
    coord_fixed(ratio = 1) +
    labs(title = mechanism_name, x = NULL, y = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "grey25"), plot.title = element_text(face = "bold", color = base_col, size = 9), axis.text.y = element_text(size = 7), axis.ticks = element_blank(), legend.position = "right", plot.margin = margin(2,5,2,5))
  
  if (i < length(mechanisms_used)) p_i <- p_i + theme(axis.text.x = element_blank()) else p_i <- p_i + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
  component_plot_list[[i]] <- p_i
}

component_heights <- pmax(component_heights, 1)

p_component <- wrap_plots(component_plot_list, ncol = 1, heights = component_heights, guides = "collect") +
  plot_annotation(title = "Components of pH-homeostasis systems across pH ecotypes", subtitle = "Color intensity indicates the proportion of MAGs carrying each functional component") &
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
p_component
component_height <- max(10, sum(component_heights) * 0.33 + 3)
ggsave(file.path(figure_dir, "pH_homeostasis_component_heatmap.pdf"), p_component, width = 180/25.4, height = component_height, units = "in", limitsize = FALSE, useDingbats = FALSE)

