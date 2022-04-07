setwd("/Users/zhaoy2/Desktop/Copied_20220309/Added_after_copy/Source_data_and_codes")

library(ggplot2)
library(ggpubr)
library(reshape2)
library(Rmisc)
library(RColorBrewer)

metadata_x <- read.delim("group_info_all_3groups_pathology.txt", sep="\t", header=T)
head(metadata_x)
metadata_x$Group
metadata_x$SampleID

cibersort_300 <- read.delim("cibersort_300.txt", sep="\t", header=T)
colnames(cibersort_300)
merged_300 <- merge(cibersort_300, metadata_x[,c(1,2)],by.x="Mixture",by.y="SampleID",all.x=F, all.y=F)
head(merged_300)

merged_300_m <- melt(merged_300)
cells <- unique(merged_300_m$variable)[1:22]
cells
merged_300_m <- merged_300_m[merged_300_m$variable %in% cells,]
merged_300_m$Group
merged_300_m$Group <- factor(merged_300_m$Group, levels=c("Normal","AIS","MIA","LUAD"))
merged_300_m$variable <- factor(merged_300_m$variable, levels=c("Dendritic.cells.activated", "Dendritic.cells.resting",
                                                                "Mast.cells.activated","Mast.cells.resting",
                                                                "T.cells.regulatory..Tregs.","T.cells.gamma.delta",
                                                                "T.cells.follicular.helper","T.cells.CD4.memory.resting",
                                                                "T.cells.CD4.memory.activated","T.cells.CD4.naive",
                                                                "T.cells.CD8","B.cells.naive","B.cells.memory",
                                                                "Plasma.cells","NK.cells.activated","NK.cells.resting",
                                                                "Eosinophils","Neutrophils","Monocytes",
                                                                "Macrophages.M2","Macrophages.M1","Macrophages.M0"),
                                labels=c("DC activated","DC resting","Mast cells activated","Mast cells resting",
                                         "Tregs","Tgd","Tfh","CD4 memory resting","CD4 memory activated",
                                         "T cells CD4 naive","T cells CD8","B cells naive","B cells memory",
                                         "Plasma cells","NK cells activated","NK cells resting","Eosinophils",
                                         "Neutrophils","Monocytes","Macrophages M2","Macrophages M1","Macrophages M0"))
merged_300_ms <- summarySE(merged_300_m, measurevar="value", groupvars=c("Group","variable"))
head(merged_300_ms)

### File preparation:
DC <- rbind(merged_300_ms[merged_300_ms$variable=="DC activated",], merged_300_ms[merged_300_ms$variable=="DC resting",])
mast <- rbind(merged_300_ms[merged_300_ms$variable=="Mast cells activated",], merged_300_ms[merged_300_ms$variable=="Mast cells resting"])
tcells <- rbind(merged_300_ms[merged_300_ms$variable=="Tregs",],
                merged_300_ms[merged_300_ms$variable=="Tfh",],merged_300_ms[merged_300_ms$variable=="CD4 memory resting",],
                merged_300_ms[merged_300_ms$variable=="CD4 memory activated",],merged_300_ms[merged_300_ms$variable=="T cells CD4 naive",])
tcells_part <- rbind(merged_300_ms[merged_300_ms$variable=="CD4 memory activated",],merged_300_ms[merged_300_ms$variable=="CD4 memory resting",],
                     merged_300_ms[merged_300_ms$variable=="T cells CD4 naive",])
tcells_cd8 <- merged_300_ms[merged_300_ms$variable=="T cells CD8",]
bcells <- rbind(merged_300_ms[merged_300_ms$variable=="B cells naive",],merged_300_ms[merged_300_ms$variable=="B cells memory",])
NK <- rbind(merged_300_ms[merged_300_ms$variable=="NK cells activated",],merged_300_ms[merged_300_ms$variable=="NK cells resting",])
macro <- rbind(merged_300_ms[merged_300_ms$variable=="Macrophages M2",],merged_300_ms[merged_300_ms$variable=="Macrophages M1",],
               merged_300_ms[merged_300_ms$variable=="Macrophages M0",])
treg <- merged_300_ms[merged_300_ms$variable=="Tregs",]
# Plotting:
ggplot(treg)+
  labs(x="Developmental stage",y="Absolute immune estimate",color="Group")+
  geom_line(aes(x=Group,y=value,group=variable,color=variable))+
  geom_errorbar(aes(x=Group,ymin=value-se, ymax=value+se,color=variable), width=0.1,size=1.2) +
  geom_point(mapping=aes(Group,value,color=variable),size=5)+
  scale_color_manual(values=c("royalblue"))+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(hjust=0.5))
ggsave("Immune/Pathology/treg.pdf", height=6, width=8, dpi=600)

ggplot(tcells_cd8)+
  labs(x="Developmental stage",y="Absolute immune estimate",color="Group")+
  geom_line(aes(x=Group,y=value,group=variable,color=variable))+
  geom_errorbar(aes(x=Group,ymin=value-se, ymax=value+se,color=variable), width=0.1,size=1.2) +
  geom_point(mapping=aes(Group,value,color=variable),size=5)+
  scale_color_manual(values=c("firebrick"))+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(hjust=0.5))
ggsave("Immune/Pathology/cd8_t_cell.pdf", height=6, width=8, dpi=600)

ggplot(NK)+
  labs(x="Developmental stage",y="Absolute immune estimate",color="Group")+
  geom_line(aes(x=Group,y=value,group=variable,color=variable))+
  geom_errorbar(aes(x=Group,ymin=value-se, ymax=value+se,color=variable), width=0.1,size=1.2) +
  geom_point(mapping=aes(Group,value,color=variable),size=5)+
  scale_color_manual(values=c(rgb(17,119,70,max=255),rgb(208,152,184,max=255)))+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(hjust=0.5))
ggsave("Immune/Pathology/NK.pdf", height=6, width=8, dpi=600)

# Validation:
# Jianjun_NC_2021:
jianjun <- read.delim("ZhangJianjun_NC_data/CIBERSORTx_Jianjun_NC_Results.txt",sep="\t",header=T)
jianjun_m <- melt(jianjun)
head(jianjun_m)
unique(jianjun_m$variable)
cells <- unique(jianjun_m$variable)[1:22]
jianjun_m <- jianjun_m[jianjun_m$variable %in% cells,]
jianjun_m$Group <- factor(jianjun_m$Group, levels=c("NL","AAH","AIS","MIA","LUAD"), labels=c("Normal","AAH","AIS","MIA","LUAD"))
jianjun_m$variable <- factor(jianjun_m$variable, levels=c("Dendritic.cells.activated", "Dendritic.cells.resting",
                                                          "Mast.cells.activated","Mast.cells.resting",
                                                          "T.cells.regulatory..Tregs.","T.cells.gamma.delta",
                                                          "T.cells.follicular.helper","T.cells.CD4.memory.resting",
                                                          "T.cells.CD4.memory.activated","T.cells.CD4.naive",
                                                          "T.cells.CD8","B.cells.naive","B.cells.memory",
                                                          "Plasma.cells","NK.cells.activated","NK.cells.resting",
                                                          "Eosinophils","Neutrophils","Monocytes",
                                                          "Macrophages.M2","Macrophages.M1","Macrophages.M0"),
                             labels=c("DC activated","DC resting","Mast cells activated","Mast cells resting",
                                      "Tregs","Tgd","Tfh","CD4 memory resting","CD4 memory activated",
                                      "T cells CD4 naive","T cells CD8","B cells naive","B cells memory",
                                      "Plasma cells","NK cells activated","NK cells resting","Eosinophils",
                                      "Neutrophils","Monocytes","Macrophages M2","Macrophages M1","Macrophages M0"))
jianjun_ms <- summarySE(jianjun_m, measurevar="value", groupvars=c("Group","variable"))

### File preparation:
DC <- rbind(jianjun_ms[jianjun_ms$variable=="DC activated",], jianjun_ms[jianjun_ms$variable=="DC resting",])
mast <- rbind(jianjun_ms[jianjun_ms$variable=="Mast cells activated",], jianjun_ms[jianjun_ms$variable=="Mast cells resting"])
tcells <- rbind(jianjun_ms[jianjun_ms$variable=="Tregs",],
                jianjun_ms[jianjun_ms$variable=="Tfh",],jianjun_ms[jianjun_ms$variable=="CD4 memory resting",],
                jianjun_ms[jianjun_ms$variable=="CD4 memory activated",],jianjun_ms[jianjun_ms$variable=="T cells CD4 naive",])
tcells_part <- rbind(jianjun_ms[jianjun_ms$variable=="CD4 memory activated",],jianjun_ms[jianjun_ms$variable=="CD4 memory resting",],
                     jianjun_ms[jianjun_ms$variable=="T cells CD4 naive",])
tcells_cd8 <- jianjun_ms[jianjun_ms$variable=="T cells CD8",]
bcells <- rbind(jianjun_ms[jianjun_ms$variable=="B cells naive",],jianjun_ms[jianjun_ms$variable=="B cells memory",])
NK <- rbind(jianjun_ms[jianjun_ms$variable=="NK cells activated",],jianjun_ms[jianjun_ms$variable=="NK cells resting",])
macro <- rbind(jianjun_ms[jianjun_ms$variable=="Macrophages M2",],jianjun_ms[jianjun_ms$variable=="Macrophages M1",],
               jianjun_ms[jianjun_ms$variable=="Macrophages M0",])
treg <- jianjun_ms[jianjun_ms$variable=="Tregs",]

ggplot(treg)+
  labs(x="Developmental stage",y="Absolute immune estimate",color="Group")+
  geom_line(aes(x=Group,y=value,group=variable,color=variable))+
  geom_errorbar(aes(x=Group,ymin=value-se, ymax=value+se,color=variable), width=0.1,size=1.2) +
  geom_point(mapping=aes(Group,value,color=variable),size=5)+
  scale_color_manual(values=c(rgb(124,201,205,max=255),rgb(69,168,175,max=255)))+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(hjust=0.5))
ggsave("Immune/Validation/treg.pdf", height=6, width=8, dpi=600)

ggplot(NK)+
  labs(x="Developmental stage",y="Absolute immune estimate",color="Group")+
  geom_line(aes(x=Group,y=value,group=variable,color=variable))+
  geom_errorbar(aes(x=Group,ymin=value-se, ymax=value+se,color=variable), width=0.1,size=1.2) +
  geom_point(mapping=aes(Group,value,color=variable),size=5)+
  scale_color_manual(values=c(rgb(17,119,70,max=255),rgb(208,152,184,max=255)))+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="top",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(hjust=0.5))
ggsave("Immune/Validation/NK.pdf", height=6, width=8, dpi=600)

