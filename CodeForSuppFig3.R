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

## Merged_300_detailed_group_info:
merged_300_d <- merge(cibersort_300, metadata_x[c(1,3)],by.x="Mixture",by.y="SampleID",all.x=F, all.y=F)
merged_300_d_m <- melt(merged_300_d)
head(merged_300_d_m)
unique(merged_300_d_m$variable)
colnames(merged_300_d_m)[2] <- "Group"
cells <- unique(merged_300_d_m$variable)[1:22]
cells
merged_300_d_m <- merged_300_d_m[merged_300_d_m$variable %in% cells,]
merged_300_d_m$Group
merged_300_d_m$Group <- factor(merged_300_d_m$Group, levels=c("Normal","AIS","MIA","IA","IB","IIIA"))
merged_300_d_m$variable <- factor(merged_300_d_m$variable, levels=c("Dendritic.cells.activated", "Dendritic.cells.resting",
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
merged_300_d_ms <- summarySE(merged_300_d_m, measurevar="value", groupvars=c("Group","variable"))
head(merged_300_d_ms)

### File preparation:
DC <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="DC activated",], merged_300_d_ms[merged_300_d_ms$variable=="DC resting",])
mast <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="Mast cells activated",], merged_300_d_ms[merged_300_d_ms$variable=="Mast cells resting"])
tcells <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="Tregs",],
                merged_300_d_ms[merged_300_d_ms$variable=="Tfh",],merged_300_d_ms[merged_300_d_ms$variable=="CD4 memory resting",],
                merged_300_d_ms[merged_300_d_ms$variable=="CD4 memory activated",],merged_300_d_ms[merged_300_d_ms$variable=="T cells CD4 naive",])
tcells_part <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="CD4 memory activated",],merged_300_d_ms[merged_300_d_ms$variable=="CD4 memory resting",],
                     merged_300_d_ms[merged_300_d_ms$variable=="T cells CD4 naive",])
tcells_cd8 <- merged_300_d_ms[merged_300_d_ms$variable=="T cells CD8",]
bcells <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="B cells naive",],merged_300_d_ms[merged_300_d_ms$variable=="B cells memory",])
NK <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="NK cells activated",],merged_300_d_ms[merged_300_d_ms$variable=="NK cells resting",])
macro <- rbind(merged_300_d_ms[merged_300_d_ms$variable=="Macrophages M2",],merged_300_d_ms[merged_300_d_ms$variable=="Macrophages M1",],
               merged_300_d_ms[merged_300_d_ms$variable=="Macrophages M0",])
treg <- merged_300_d_ms[merged_300_d_ms$variable=="Tregs",]
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
ggsave("Immune/Pathology/treg_detailed.pdf", height=6, width=8, dpi=600)

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
ggsave("Immune/Pathology/cd8_t_cell_detailed.pdf", height=6, width=8, dpi=600)

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
ggsave("Immune/Pathology/NK_detailed.pdf", height=6, width=8, dpi=600)
