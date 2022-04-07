setwd("/Users/zhaoy2/Desktop/Copied_20220309/Added_after_copy/Source_data_and_codes")

library(stringr)
library(lme4)

tpm_hugo_coding <- read.delim("TPM_RIN5_coding_genes.txt", sep="\t", header=T, row.names=1)

library(emmeans)

# Read metadata with X:
metadata_x <- read.delim("group_info_all_3groups_pathology.txt", sep="\t", header=T)
head(metadata_x)
tpm_hugo_coding <- read.delim("TPM_RIN5_coding_genes.txt", sep="\t", header=T, row.names=1)
tpm_hugo_coding[1:5,1:5]
dim(tpm_hugo_coding)
tpm_hugo_coding_t <- as.data.frame(t(tpm_hugo_coding))
tpm_hugo_coding_t[1:5,1:5]
tpm_hugo_coding_t$SampleID <- rownames(tpm_hugo_coding_t)
merged <- merge(tpm_hugo_coding_t, metadata_x, by="SampleID", all.x=F, all.y=F)
dim(merged)
merged[1:3,c(1:3,18071:18073)]
merged$group2 <- factor(merged$Group, levels=c("LUAD","MIA","AIS","Normal"))

### 4 groups:
# Set sig = 0.0001
start_time_00001 <- Sys.time()
genes_anova_00001 <- c()
a_00001s <- c()
for (i in 2:18071){
  print(i)
  print(colnames(merged)[i])
  tryCatch({
    fit <- aov(as.numeric(as.character(merged[,i])) ~ merged$group2)
    emm <- emmeans(fit, ~group2)
    means <- summary(emm)$emmean
    a <- pairs(emm)
    a <- as.data.frame(a)
    if ((a$p.value[1] <= 0.0001 | a$p.value[2] <= 0.0001 | a$p.value[3] <= 0.0001 | a$p.value[4] <= 0.0001 | a$p.value[5] <= 0.0001 | a$p.value[6] <= 0.0001) & 
        (means[1]/means[2]>=2 | means[1]/means[2]<=0.5 | means[1]/means[3]>=2 | means[1]/means[3]<= 0.5 |
         means[1]/means[4]>=2 | means[1]/means[4]<=0.5 | means[2]/means[3]>=2 | means[2]/means[3]<=0.5 | 
         means[2]/means[4]>=2 | means[2]/means[4]<=0.5 | means[3]/means[4]>=2 | means[3]/means[4]<=0.5)){
      genes_anova_00001 <- rbind(genes_anova_00001, colnames(merged)[i])
      a$gene <- colnames(merged)[i]
      a_00001s <- rbind(a_00001s,a)
    }
  }, error=function(e){cat("ERROR: ", conditionMessage(e), "\n")})
}
end_time_00001 <- Sys.time()
time_diff_00001 <- end_time_00001 - start_time_00001
time_diff_00001

# Use 0.0001 for the next steps:
genes_anova_00001 <- as.data.frame(genes_anova_00001)
genes_00001 <- genes_anova_00001$V1
merged[1:5,1:5]
merged[1:5,18072:18073]
nrow(genes_anova_00001)

trend <- data.frame(Gene=as.character,AIS_vs_Normal=as.numeric(),MIA_vs_Normal=as.numeric(),LUAD_vs_Normal=as.numeric(),
                    MIA_vs_AIS=as.numeric(),LUAD_vs_AIS=as.numeric(),LUAD_vs_MIA=as.numeric())
for (g in genes_00001){
  print(g)
  fit <- aov(as.numeric(as.character(merged[,colnames(merged)==g])) ~ merged$group2)
  emm <- emmeans(fit, ~group2)
  means <- summary(emm)$emmean
  luad <- means[1]
  mia <- means[2]
  ais <- means[3]
  normal <- means[4]
  LUAD_vs_Normal <- luad/normal
  MIA_vs_Normal <- mia/normal
  AIS_vs_Normal <- ais/normal
  LUAD_vs_AIS <- luad/ais
  MIA_vs_AIS <- mia/ais
  LUAD_vs_MIA <- luad/mia
  line <- cbind(g,AIS_vs_Normal,MIA_vs_Normal,LUAD_vs_Normal,MIA_vs_AIS,LUAD_vs_AIS,LUAD_vs_MIA)
  trend <- rbind(trend, line)
}
colnames(trend)[1] <- "Gene"
trend_backup <- trend
head(trend,30)

for (i in 2:7){
  trend[,i] <- as.numeric(as.character(trend[,i]))
}

trend_log2 <- cbind(trend$Gene,log2(trend[,2:7]))
colnames(trend_log2)[1] <- "Gene"
class(trend[,2])
head(trend_log2)
head(trend)
write.table(trend_log2, file="ANOVA_0.0001.txt", sep="\t", row.names=F)

# Get rid of NAs:
trend_log2 <- read.delim("ANOVA_0.0001.txt", sep="\t", header=T)
trend_log2 <- na.omit(trend_log2)
nrow(trend_log2) # 3034
head(trend_log2)

# Find patterns:
## Pattern 1: up-up-up
p1 <- trend_log2[trend_log2[,2]>=1 & trend_log2[,5]>=1 & trend_log2[,7]>=1,]
p1$Pattern <- "p1"
nrow(p1)
p1
write.csv(p1, file="ANOVA/4groups/p1.csv", row.names=F)
## Pattern 2: up-same-up
p2 <- trend_log2[trend_log2[,2]>=1 & trend_log2[,5]<1 & trend_log2[,5]>-1 & trend_log2[,7]>=1,]
nrow(p2)
p2$Pattern <- "p2"
p2
write.csv(p2, file="ANOVA/4groups/p2.csv", row.names=F)
## Pattern 3: up-same-same
p3 <- trend_log2[trend_log2[,2]>=1 & trend_log2[,5]<1 & trend_log2[,5]>-1 & trend_log2[,7]<1 & trend_log2[,7]>-1,]
nrow(p3)
p3$Pattern <- "p3"
p3
write.csv(p3, file="ANOVA/4groups/p3.csv", row.names=F)
## Pattern 4: same-up-up
p4 <- trend_log2[trend_log2[,2]<1 & trend_log2[,2]>-1 & trend_log2[,5]>=1 & trend_log2[,7]>=1,]
nrow(p4)
p4$Pattern <- "p4"
p4
write.csv(p4, file=".ANOVA/4groups/p4.csv", row.names=F)
## Pattern 5: same-up-same
p5 <- trend_log2[trend_log2[,2]<1 & trend_log2[,2]>-1 & trend_log2[,5]>=1 & trend_log2[,7]>-1 & trend_log2[,7]<1,]
nrow(p5)
p5$Pattern <- "p5"
p5
write.csv(p5, file="ANOVA/4groups/p5.csv", row.names=F)
## Pattern 6: same-same-up:
p6 <- trend_log2[trend_log2[,2]<1 & trend_log2[,2]>-1 & trend_log2[,5]<1 & trend_log2[,5]>-1 & trend_log2[,7]>=1,]
nrow(p6)
p6$Pattern <- "p6"
p6
write.csv(p6, file="ANOVA/4groups/p6.csv", row.names=F)
## Pattern 7: down-down-down:
p7 <- trend_log2[trend_log2[,2]<=-1 & trend_log2[,5]<=-1 & trend_log2[,7]<=-1,]
nrow(p7)
p7$Pattern <- "p7"
p7
write.csv(p7, file="ANOVA/4groups/p7.csv", row.names=F)
## Pattern 8: down-same-down:
p8 <- trend_log2[trend_log2[,2]<=-1 & trend_log2[,5]>-1 & trend_log2[,5]<1 & trend_log2[,7]<=-1,]
nrow(p8)
p8$Pattern <- "p8"
p8
write.csv(p8, file="ANOVA/4groups/p8.csv", row.names=F)
## Pattern 9: down-same-same:
p9 <- trend_log2[trend_log2[,2]<=-1 & trend_log2[,5]<1 & trend_log2[,5]>-1 & trend_log2[,7]<1 & trend_log2[,7]>-1,]
nrow(p9)
p9$Pattern <- "p9"
p9
write.csv(p9, file="ANOVA/4groups/p9.csv", row.names=F)
## Pattern 10: same-down-down
p10 <- trend_log2[trend_log2[,2]>-1 & trend_log2[,2] < 1 & trend_log2[,5]<=-1 & trend_log2[,7]<=-1,]
nrow(p10)
p10$Pattern <- "p10"
p10
write.csv(p10, file="ANOVA/4groups/p10.csv", row.names=F)
## Pattern 11: same-down-same:
p11 <- trend_log2[trend_log2[,2]>-1 & trend_log2[,2] < 1 & trend_log2[,5]<=-1 & trend_log2[,7]>-1 & trend_log2[,7]<1,]
nrow(p11)
p11$Pattern <- "p11"
p11
write.csv(p11, file="ANOVA/4groups/p11.csv", row.names=F)
## Pattern 12: same-same-down:
p12 <- trend_log2[trend_log2[,2]>-1 & trend_log2[,2] < 1 & trend_log2[,5]>-1 & trend_log2[,5]<1 & trend_log2[,7]<=-1,]
nrow(p12)
p12$Pattern <- "p12"
p12
write.csv(p12, file="ANOVA/4groups/p12.csv", row.names=F)
all <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
write.csv(all, file="all_patterns_combined.csv", row.names=F)
# Get full names of genes:
patterns <- read.csv("all_patterns_combined.csv")
head(patterns)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
symbols <- as.character(patterns$Gene)
symbols
out <- select(org.Hs.eg.db,keys = symbols,columns = "GENENAME",keytype = "SYMBOL")
head(out)
patterns_fullnames <- merge(out,patterns,by.x="SYMBOL",by.y="Gene",all.x=F,all.y=F)
patterns_fullnames <- patterns_fullnames[order(patterns_fullnames$Pattern),]
head(patterns_fullnames)

colnames(patterns_fullnames)[1] <- "Gene"
write.csv(patterns_fullnames, file="all_patterns_combined_full_genenames.csv",row.names=F)
patterns_fullnames <- read.csv("all_patterns_combined_full_genenames.csv") ## Added a manual ordering step.
head(patterns_fullnames)

patterns_to_plot <- tpm_hugo_coding[rownames(tpm_hugo_coding) %in% patterns_fullnames$Gene,]
patterns_to_plot[1:5,1:5]

## Re-order by patterns:
patterns_to_plot$Gene <- rownames(patterns_to_plot)
patterns_to_plot1 <- merge(patterns_to_plot, patterns_fullnames[,c("Gene","Pattern")],by="Gene", all.x=F, all.y=F)
patterns_to_plot1[1:5,]
patterns_to_plot2 <- patterns_to_plot1[order(patterns_to_plot1$Pattern),]
rownames(patterns_to_plot2) <- patterns_to_plot2$Gene
write.csv(patterns_to_plot2, file="all_patterns_to_plot.csv", row.names=F) ## Added a manual ordering step.
patterns_to_plot_final <- read.csv("all_patterns_to_plot.csv", row.names=1)
patterns_to_plot_final <- log2(patterns_to_plot_final+0.1)
head(patterns_to_plot_final)
patterns_to_plot_final[1:5,1:5]

head(patterns)
patterns_anno <- as.data.frame(patterns$Pattern)
head(patterns_anno)
rownames(patterns_anno) <- patterns$Gene
colnames(patterns_anno) <- "Pattern"
head(patterns_anno)
patterns_anno$Pattern <- factor(patterns_anno$Pattern, levels=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12"))

ann_colors <- list(
  Pattern=c(p1="#A6CEE3",p2="#1F78B4",p3="#B2DF8A",p4="#33A02C",p5="#FB9A99",p6="#E31A1C",p7="#FDBF6F",p8="#FF7F00",p9="#CAB2D6",p10="#6A3D9A",p11="#FFFF99",p12="#B15928"),
  Group=c(Normal="#E41A1C",AIS="#377EB8",MIA="#4DAF4A",LUAD="#984EA3")
)
ann_colors

pheatmap(patterns_to_plot_final,
         cluster_cols=F,
         cluster_rows=F,
         show_rownames = F,
         show_colnames = F,
         annotation_col=anno1,
         annotation_row=patterns_anno,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         #breaks=bk,
         scale="row")