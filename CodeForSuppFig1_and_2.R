setwd("/Users/zhaoy2/Desktop/Copied_20220309/Added_after_copy/Source_data_and_codes")

library(ggplot2)
file1 <- read.delim("TMB.txt",sep="\t",header=T)
file1$GGO_group <- factor(file1$GGO_group, levels=c("Pure_GGO","Subsolid","Solid"))
head(file1)
file1$group <- factor(file1$group, levels=c("AIS","MIA","LUAD"))

# Calculate p values:
anova1 <- aov(file1$TMB~file1$GGO_group)
summary(anova1)
TukeyHSD(anova1)

anova2 <- aov(file1$TMB~file1$group)
summary(anova2)
TukeyHSD(anova2)

TMB_GGO <- ggplot(data=file1, aes(x=GGO_group, y=TMB, fill=GGO_group)) +
  geom_boxplot(color="black",size=0.2) +
  labs(x="Group", y = "No. of mutations / Mb", fill="Radiology") +
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  annotate("segment", x=c(1,1,2),xend=c(1,2,2), y= c(9,9.5,9), yend=c(9.5,9.5,9.5))+
  annotate("text", x=1.5, y=10.5, label="p=0.657", size=5)+
  annotate("segment", x=c(2,2,3),xend=c(2,3,3), y= c(33.5,34,33.5), yend=c(34,34,34))+
  annotate("text", x=2.5, y=35, label="p<0.001", size=5)+
  annotate("segment", x=c(1,1,3),xend=c(1,3,3), y= c(37,37.5,37), yend=c(37.5,37.5,37.5))+
  annotate("text", x=2, y=38.5, label="p<0.001", size=5)+
  scale_x_discrete(labels=c("Pure GGO","Mixed GGO","Solid"))+
  scale_fill_brewer(palette="Accent")
TMB_GGO
ggsave("FUSCC_TMB_GGO.pdf",plot=TMB_GGO,width = 8, height = 6, dpi = 600)

TMB_path <- ggplot(data=file1, aes(x=group, y=TMB, fill=group)) +
  geom_boxplot(color="black",size=0.2) +
  labs(x="Group", y = "No. of mutations / Mb", fill="Pathology") +
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold", size=13),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=14),
        panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=12))+
  annotate("segment", x=c(1,1,2),xend=c(1,2,2), y= c(9,9.5,9), yend=c(9.5,9.5,9.5))+
  annotate("text", x=1.5, y=10.5, label="p=0.998", size=5)+
  annotate("segment", x=c(2,2,3),xend=c(2,3,3), y= c(33.5,34,33.5), yend=c(34,34,34))+
  annotate("text", x=2.5, y=35, label="p<0.001", size=5)+
  annotate("segment", x=c(1,1,3),xend=c(1,3,3), y= c(37,37.5,37), yend=c(37.5,37.5,37.5))+
  annotate("text", x=2, y=38.5, label="p=0.003", size=5)+
  scale_x_discrete(labels=c("Pure GGO","Mixed GGO","Solid"))+
  scale_fill_brewer(palette="Set1")
TMB_path
ggsave("C:/Users/terry/Desktop/Projects/GGO_WES/results_RIN5_new/FUSCC_TMB_pathlogy.pdf",plot=TMB_path,width = 8, height = 6, dpi = 600)

