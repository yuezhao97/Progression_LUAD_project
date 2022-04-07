setwd("/Users/zhaoy2/Desktop/Copied_20220309/Added_after_copy/Source_data_and_codes")

library(stringr)
library(lme4)
metadata <- read.delim("group_info_all_3groups_pathology.txt", sep="\t", header=T)
head(metadata)
samples <- metadata$SampleID

exp <- read.delim("TPM_RIN5_coding_genes.txt", sep="\t", header=T, row.names=1)
exp[1:5,1:5]

up <- c("BCL2L15","COMP","CST1","FAM83A")
down8 <- c("ITLN2","MARCO","C8B","MASP1","CD36","TAL1","PPBP","CDH5")
exp_up <- exp[rownames(exp) %in% up,1:300]
exp_d8 <- exp[rownames(exp) %in% down8,1:300]

head(exp_up)
exp_up[1:4,1:5]
exp_up_log <- log2(exp_up+1)
exp_up_log[1:4,1:5]
exp_d8_log <- log2(exp_d8+1)
mean_up <- c()
mean_down <- c()
for (i in 1:300){
  mean_sample_up <- mean(exp_up_log[,i])
  mean_up <- rbind(mean_up, mean_sample_up)
  mean_sample_down <- mean(exp_d8_log[,i])
  mean_down <- rbind(mean_down, mean_sample_down)
}
mean_final <- cbind(mean_up, mean_down)
rownames(mean_final) <- colnames(exp_up_log)
colnames(mean_final) <- c("mean_up","mean_down")
mean_final <- as.data.frame(mean_final)
mean_final$index <- mean_final$mean_up - mean_final$mean_down
mean_final$SampleID <- rownames(mean_final)
head(mean_final)
write.csv(mean_final, file="index_new.csv")

# Boxplot:
library(ggplot2)
library(ggpubr)
metadata_x <- read.delim("group_info_all_3groups_pathology.txt", sep="\t", header=T)
head(metadata_x)
merged <- merge(mean_final, metadata_x, by="SampleID", all.x=F, all.y=F)
merged$Group <- factor(merged$Group, levels=c("Normal","AIS","MIA","LUAD"))

ggplot(data=merged, aes(x=Group, y=index, fill=Group)) +
  geom_boxplot() +
  labs(x="Group", y = "Tumor Progressive Index", fill="Group") +
  #scale_fill_brewer(palette="Blues") + 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  stat_compare_means()+
  scale_fill_brewer(palette="Set1")
ggsave("index/index_new/up_down8_4groups_pathology_log.pdf", plot=last_plot(), width=8, height=6, dpi=600)

# Validation:
### Kadara et al. Cancer Research 2017:
tpm_k <- read.delim("Kadara_CR_data/GSE102511_Smruthy_etal_allsamples_TPM.txt", header=T, sep="\t")
tpm_k[1:5,1:5]
tpm_k_up <- tpm_k[tpm_k$Unique.ID %in% up,]
tpm_k8 <- tpm_k[tpm_k$Unique.ID %in% down8,]
head(tpm_k_up)
rownames(tpm_k_up) <- tpm_k_up$Unique.ID
tpm_k_up <- tpm_k_up[,-1]
rownames(tpm_k8) <- tpm_k8$Unique.ID
tpm_k8 <- tpm_k8[,-1]

tpm_k_up[1:4,1:5]
tpm_k_up_log <- log2(tpm_k_up+1)
tpm_k_up_log[1:4,1:5]
tpm_k8_log <- log2(tpm_k8+1)
mean_k_up <- c()
mean_k_down <- c()
for (i in 1:ncol(tpm_k_up_log)){
  mean_k_sample_up <- mean(tpm_k_up_log[,i])
  mean_k_up <- rbind(mean_k_up, mean_k_sample_up)
  mean_k_sample_down <- mean(tpm_k8_log[,i])
  mean_k_down <- rbind(mean_k_down, mean_k_sample_down)
}
mean_k_final <- cbind(mean_k_up, mean_k_down)
rownames(mean_k_final) <- colnames(tpm_k_up_log)
colnames(mean_k_final) <- c("mean_up","mean_down")
mean_k_final <- as.data.frame(mean_k_final)
mean_k_final$index <- mean_k_final$mean_up - mean_k_final$mean_down
mean_k_final$SampleID <- rownames(mean_k_final)
head(mean_k_final)
write.csv(mean_k_final, file="index/index_new/index_Kadara_new.csv")

# Boxplot:
metadata_k <- read.delim("Kadara_CR_data/metadata_Kadara.txt", sep="\t", header=T)
head(metadata_k)
metadata_k$SampleID <- gsub("-",".",metadata_k$SampleID)
mean_k_final$SampleID <- gsub("X","",mean_k_final$SampleID)
merged_k <- merge(mean_k_final, metadata_k, by="SampleID", all.x=F, all.y=F)
merged_k$Group <- factor(merged_k$Group, levels=c("NL","AAH","LUAD"))

ggplot(data=merged_k, aes(x=Group, y=index, fill=Group)) +
  geom_boxplot() +
  labs(x="Group", y = "Tumor Progressive Index", fill="Group") +
  #scale_fill_brewer(palette="Blues") + 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  stat_compare_means()+
  scale_fill_brewer(palette="Set1")
ggsave("index/index_new/up_down8_Kadara_validation_log.pdf", plot=last_plot(), width=8, height=6, dpi=600)

### Validation using TCGA data:
rna_t <- read.delim("TCGA_LUAD/data_mrna_seq_v2_rsem.txt", sep="\t", header=T)
rna_n <- read.delim("TCGA_LUAD/data_mrna_seq_v2_rsem_normal_samples.txt", sep="\t", header=T)
rna_t_up <- rna_t[rna_t$Hugo_Symbol %in% up,]
rna_t_down <- rna_t[rna_t$Hugo_Symbol %in% down8,]
rna_n_up <- rna_n[rna_n$Hugo_Symbol %in% up,]
rna_n_down <- rna_n[rna_n$Hugo_Symbol %in% down8,]
rownames(rna_t_up) <- rna_t_up$Hugo_Symbol
rna_t_up <- rna_t_up[,-c(1,2)]
rownames(rna_t_down) <- rna_t_down$Hugo_Symbol
rna_t_down <- rna_t_down[,-c(1,2)]
rownames(rna_n_up) <- rna_n_up$Hugo_Symbol
rna_n_up <- rna_n_up[,-c(1,2)]
rownames(rna_n_down) <- rna_n_down$Hugo_Symbol
rna_n_down <- rna_n_down[,-c(1,2)]
rna_up <- cbind(rna_t_up, rna_n_up)
rna_down <- cbind(rna_t_down, rna_n_down)

rna_up_log <- log2(rna_up+1)
rna_up_log[1:4,1:5]
rna_down_log <- log2(rna_down+1)
mean_tcga_up <- c()
mean_tcga_down <- c()
for (i in 1:ncol(rna_up_log)){
  mean_tcga_sample_up <- mean(rna_up_log[,i])
  mean_tcga_up <- rbind(mean_tcga_up, mean_tcga_sample_up)
  mean_tcga_sample_down <- mean(rna_down_log[,i])
  mean_tcga_down <- rbind(mean_tcga_down, mean_tcga_sample_down)
}
mean_tcga_final <- cbind(mean_tcga_up, mean_tcga_down)
rownames(mean_tcga_final) <- colnames(rna_up_log)
colnames(mean_tcga_final) <- c("mean_up","mean_down")
mean_tcga_final <- as.data.frame(mean_tcga_final)
mean_tcga_final$index <- mean_tcga_final$mean_up - mean_tcga_final$mean_down
mean_tcga_final$SampleID <- rownames(mean_tcga_final)
tumor <- colnames(rna_t_up)
normal <- colnames(rna_n_up)
for (i in 1:nrow(mean_tcga_final)){
  if (mean_tcga_final[i,"SampleID"] %in% tumor){
    mean_tcga_final[i,"Group"] <- "Tumor"
  } else {mean_tcga_final[i,"Group"] <- "Normal"}
}
write.csv(mean_tcga_final, file="index/index_new/index_TCGA_LUAD.csv")

# Boxplot:
library(ggplot2)
library(ggpubr)

ggplot(data=mean_tcga_final, aes(x=Group, y=index, fill=Group)) +
  geom_boxplot() +
  labs(x="Group", y = "index", fill="Group") +
  #scale_fill_brewer(palette="Blues") + 
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  stat_compare_means()+
  scale_fill_brewer(palette="Set1")
ggsave("TCGA_validation_2groups_TN_log.pdf", plot=last_plot(), width=8, height=6, dpi=600)
