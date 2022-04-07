setwd("/Users/zhaoy2/Desktop/Copied_20220309/Added_after_copy/Source_data_and_codes")
### Survival plot:

mean_final <- read.csv("index_new.csv", row.names=1)
metadata_x <- read.delim("group_info_all_3groups_pathology.txt", sep="\t", header=T)
head(metadata_x)
merged <- merge(mean_final, metadata_x, by="SampleID", all.x=F, all.y=F)

survival_index <- read.csv("survival_index.csv")
head(survival_index)
head(merged)
nrow(survival_index)
survival_index_new <- merge(survival_index, merged[,c(1:4)], by.x="RNA.LC", by.y="SampleID", all.x=T, all.y=F)
head(survival_index_new)
m <- median(survival_index_new$index)

for (i in 1:nrow(survival_index_new)){
  if (survival_index_new[i,"index"] >= m){
    survival_index_new[i,"index_group"] <- "high"
  } else {survival_index_new[i,"index_group"] <- "low"}
  if (survival_index_new[i,"index"] >= 0){
    survival_index_new[i,"index_group_0"] <- "high"
  } else {survival_index_new[i,"index_group_0"] <- "low"}
}


library(survminer)
library(gridExtra)
library(survival)
ggsurvplot(survfit(Surv(as.numeric(survival_index_new$RFS),as.logical(as.numeric(survival_index_new$RFS.E)))~as.factor(survival_index_new$index_group),data = survival_index_new),
           #title and legend of survival plot
           title = "",
           #font.title = c(14, "bold", "darkred"),
           subtitle = "",
           #font.subtitle = c(14, "plain", "darkred"),
           #legend
           legend.title = "Tumor Progressive Index",
           legend.labs = c("High", "Low"),
           #legend = c(0.84,0.92),
           font.legend = c(12,"plain","black"),
           #x and y
           xlab = "Time (months)",
           font.x = c(14,"bold","black"),
           ylab = "Recurrence-free survival",
           font.y = c(14,"bold","black"),
           conf.int.fill = "yellow", 
           axes.offset = T,
           surv.scale = "percent",
           break.time.by = c(12),
           font.tickslab = c(12, "bold", "black"),
           #p value and confidence interval
           pval = TRUE,
           pval.size = c(5),
           #conf.int.fill = c("#9F0B16", "#0C689F"),
           conf.int = F,
           conf.int.style = "step",
           conf.int.alpha = 0.7,
           #censor and survival curve
           censor = F,
           censor.shape = 124,
           censor.size = c(3),
           linetype = "solid", 
           surv.median.line = "none", 
           #risk plot
           risk.table = T,
           #risk.table.fontsize = c(3.7,"bold","black"),
           font.table = c(3.2,"bold","black"),
           tables.y.text = F,
           tables.y.text.col = T,
           tables.height = 0.27,
           risk.table.title = "",
           risk.table.pos = "in",
           risk.table.col = "darkblue",
           #censor plot
           ncensor.plot = F,
           ncensor.plot.title = "",
           ncensor.plot.y.text = "",
           ncensor.plot.title.font = c(1),
           ncensor.plot.height = c(0.24),
           # Change ggplot2 theme
           theme = theme_bw(), 
           #theme_classic(),
           #theme = minimal(),
           #palette = c(rgb(149,88,57, max=255), rgb(105,100,123, max=255)))
           palette = "nejm")

ggsurvplot(survfit(Surv(as.numeric(survival_index_new$OS),as.logical(as.numeric(survival_index_new$OS.E)))~as.factor(survival_index_new$index_group),data = survival_index_new),
           #title and legend of survival plot
           title = "",
           #font.title = c(14, "bold", "darkred"),
           subtitle = "",
           #font.subtitle = c(14, "plain", "darkred"),
           #legend
           legend.title = "Tumor Progressive Index",
           legend.labs = c("High", "Low"),
           #legend = c(0.84,0.92),
           font.legend = c(12,"plain","black"),
           #x and y
           xlab = "Time (months)",
           font.x = c(14,"bold","black"),
           ylab = "Recurrence-free survival",
           font.y = c(14,"bold","black"),
           conf.int.fill = "yellow", 
           axes.offset = T,
           surv.scale = "percent",
           break.time.by = c(12),
           font.tickslab = c(12, "bold", "black"),
           #p value and confidence interval
           pval = TRUE,
           pval.size = c(5),
           #conf.int.fill = c("#9F0B16", "#0C689F"),
           conf.int = F,
           conf.int.style = "step",
           conf.int.alpha = 0.7,
           #censor and survival curve
           censor = F,
           censor.shape = 124,
           censor.size = c(3),
           linetype = "solid", 
           surv.median.line = "none", 
           #risk plot
           risk.table = T,
           #risk.table.fontsize = c(3.7,"bold","black"),
           font.table = c(3.2,"bold","black"),
           tables.y.text = F,
           tables.y.text.col = T,
           tables.height = 0.27,
           risk.table.title = "",
           risk.table.pos = "in",
           risk.table.col = "darkblue",
           #censor plot
           ncensor.plot = F,
           ncensor.plot.title = "",
           ncensor.plot.y.text = "",
           ncensor.plot.title.font = c(1),
           ncensor.plot.height = c(0.24),
           # Change ggplot2 theme
           theme = theme_bw(), 
           #theme_classic(),
           #theme = minimal(),
           #palette = c(rgb(149,88,57, max=255), rgb(105,100,123, max=255)))
           palette = "nejm")

### Validation using TCGA-LUAD data:
mean_tcga_final <- read.csv("index/index_new/index_TCGA_LUAD.csv",row.names=1)
clin_tcga <- read.delim("TCGA_LUAD/data_clinical_patient.txt", sep="\t", header=T)
surv_tcga <- clin_tcga[,c("PATIENT_ID","PFS_MONTHS","PFS_STATUS","OS_MONTHS","OS_STATUS")]
surv_tcga$SampleID <- gsub("-",".",surv_tcga$PATIENT_ID)
library(stringr)
mean_tcga_final$SampleID <- str_match(mean_tcga_final$SampleID, "(TCGA\\..{2}\\..{4})\\.[0-9]{2}")[,2]
mean_tcga_final_tumor <- mean_tcga_final[mean_tcga_final$Group=="Tumor",]
merged_tcga <- merge(mean_tcga_final_tumor, surv_tcga, by="SampleID", all.x=F, all.y=F)
m_tcga <- median(merged_tcga$index) #2.907834
for (i in 1:nrow(merged_tcga)){
  if (merged_tcga[i,"index"] >= m_tcga){
    merged_tcga[i,"index_group"] <- "high"
  } else {merged_tcga[i,"index_group"] <- "low"}
  if (merged_tcga[i,"index"] >= 0){
    merged_tcga[i,"index_group_0"] <- "high"
  } else {merged_tcga[i,"index_group_0"] <- "low"}
}

# Survival plots:
ggsurvplot(survfit(Surv(as.numeric(merged_tcga$PFS_MONTHS),as.logical(as.numeric(merged_tcga$PFS_STATUS)))~as.factor(merged_tcga$index_group),data = merged_tcga),
           #title and legend of survival plot
           title = "",
           #font.title = c(14, "bold", "darkred"),
           subtitle = "",
           #font.subtitle = c(14, "plain", "darkred"),
           #legend
           legend.title = "Tumor Progressive Index",
           legend.labs = c("High", "Low"),
           #legend = c(0.84,0.92),
           font.legend = c(12,"plain","black"),
           #x and y
           xlab = "Time (months)",
           font.x = c(14,"bold","black"),
           ylab = "Progression-free survival",
           font.y = c(14,"bold","black"),
           conf.int.fill = "yellow", 
           axes.offset = T,
           surv.scale = "percent",
           break.time.by = c(12),
           font.tickslab = c(12, "bold", "black"),
           #p value and confidence interval
           pval = TRUE,
           pval.size = c(5),
           #conf.int.fill = c("#9F0B16", "#0C689F"),
           conf.int = F,
           conf.int.style = "step",
           conf.int.alpha = 0.7,
           #censor and survival curve
           censor = F,
           censor.shape = 124,
           censor.size = c(3),
           linetype = "solid", 
           surv.median.line = "none", 
           #risk plot
           risk.table = T,
           #risk.table.fontsize = c(3.7,"bold","black"),
           font.table = c(3.2,"bold","black"),
           tables.y.text = F,
           tables.y.text.col = T,
           tables.height = 0.27,
           risk.table.title = "",
           risk.table.pos = "in",
           risk.table.col = "darkblue",
           #censor plot
           ncensor.plot = F,
           ncensor.plot.title = "",
           ncensor.plot.y.text = "",
           ncensor.plot.title.font = c(1),
           ncensor.plot.height = c(0.24),
           # Change ggplot2 theme
           theme = theme_bw(), 
           #theme_classic(),
           #theme = minimal(),
           #palette = c(rgb(149,88,57, max=255), rgb(105,100,123, max=255)))
           palette = "nejm")

ggsurvplot(survfit(Surv(as.numeric(merged_tcga$OS_MONTHS),as.logical(as.numeric(merged_tcga$OS_STATUS)))~as.factor(merged_tcga$index_group),data = merged_tcga),
           #title and legend of survival plot
           title = "",
           #font.title = c(14, "bold", "darkred"),
           subtitle = "",
           #font.subtitle = c(14, "plain", "darkred"),
           #legend
           legend.title = "Tumor Progressive Index",
           legend.labs = c("High", "Low"),
           #legend = c(0.84,0.92),
           font.legend = c(12,"plain","black"),
           #x and y
           xlab = "Time (months)",
           font.x = c(14,"bold","black"),
           ylab = "Recurrence-free survival",
           font.y = c(14,"bold","black"),
           conf.int.fill = "yellow", 
           axes.offset = T,
           surv.scale = "percent",
           break.time.by = c(12),
           font.tickslab = c(12, "bold", "black"),
           #p value and confidence interval
           pval = TRUE,
           pval.size = c(5),
           #conf.int.fill = c("#9F0B16", "#0C689F"),
           conf.int = F,
           conf.int.style = "step",
           conf.int.alpha = 0.7,
           #censor and survival curve
           censor = F,
           censor.shape = 124,
           censor.size = c(3),
           linetype = "solid", 
           surv.median.line = "none", 
           #risk plot
           risk.table = T,
           #risk.table.fontsize = c(3.7,"bold","black"),
           font.table = c(3.2,"bold","black"),
           tables.y.text = F,
           tables.y.text.col = T,
           tables.height = 0.27,
           risk.table.title = "",
           risk.table.pos = "in",
           risk.table.col = "darkblue",
           #censor plot
           ncensor.plot = F,
           ncensor.plot.title = "",
           ncensor.plot.y.text = "",
           ncensor.plot.title.font = c(1),
           ncensor.plot.height = c(0.24),
           # Change ggplot2 theme
           theme = theme_bw(), 
           #theme_classic(),
           #theme = minimal(),
           #palette = c(rgb(149,88,57, max=255), rgb(105,100,123, max=255)))
           palette = "nejm")
