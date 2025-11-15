library(GeomxTools)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(writexl)
library(readxl)
library(ggpubr)
library(rstatix)

load("DataSubsets_fordownstream.RData")

# only including patients that have both A and I on the slide. 
iall_aall <- intersect(targetDataSubset_pos$`slide name`[targetDataSubset_pos$Epithelia %in% c("I")], targetDataSubset_pos$`slide name`[targetDataSubset_pos$Epithelia %in% c("A")]) #11
A_I_pos_both <- targetDataSubset_pos[,targetDataSubset_pos$Epithelial %in% c("I","A") & targetDataSubset_pos$`slide name` %in% iall_aall] 

classic_genes <- c("C1QA","C1QB", "C1QC", "C1R","C1S", "C2","C3", "C4B")
b10_genes <- rownames(targetDataSubset_pos)

geomx_df <- A_I_pos_both
geneList <- intersect(classic_genes,b10_genes)
i = 1
box2 <- list()

######### Average classic complement gene expression in groups of interest (Figure 5c-d) ######### 
logNormCounts <- assayDataElement(geomx_df[geneList,], elt = "log_q")
colnames(logNormCounts) <- geomx_df$ROI_REL
logNormCounts <- melt(logNormCounts)
colnames(logNormCounts) <- c("gene","ROI_REL","log_expression")
meta <- data.frame(ROI_REL = pData(geomx_df)$ROI_REL,
                   Epithelial = factor(pData(geomx_df)$Epithelial,levels = c("A","I")),
                   Prognosis = pData(geomx_df)$Prognosis,
                   Patient = pData(geomx_df)$`slide name`)
test <- merge(logNormCounts,meta,by = "ROI_REL")
test$Prog_EP <- paste0(test$Prognosis,"_",test$Epithelial) %>% factor()

mean_expression <- aggregate(log_expression ~ ROI_REL + Prognosis + Epithelial, data = test, FUN = mean)
mean_expression$Prog_EP <- paste0(mean_expression$Prognosis,"_",mean_expression$Epithelial) %>% factor()

# PPx A vs PPx I
subset <- mean_expression[mean_expression$Prognosis == "PPx",]
box_avg1 <- ggplot(subset, aes(x = Prog_EP, y = log_expression, fill = Epithelial)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Epithelial), size = 1) + 
  scale_fill_manual(values = c("A" = "steelblue1", "I" = "lightskyblue2"), name = "Epithelial") +
  scale_color_manual(values = c("A" = "steelblue1", "I" = "lightskyblue2")) +
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
        panel.background = element_blank(), legend.position = "none")+
  labs(x ="Prognosis", y = "Average expression (log)", title = "Classic Complement Genes") +
  stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("PPx_A","PPx_I")))

# PPx A vs GPx A
subset <- mean_expression[mean_expression$Epithelial == "A",]
box_avg2 <- ggplot(subset, aes(x = Prog_EP, y = log_expression, fill = Prognosis)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Prognosis), size = 1) + 
  scale_fill_manual(values = c("GPx" = "pink4", "PPx" = "steelblue1"), name = "Prognosis") +
  scale_color_manual(values = c("GPx" = "pink4", "PPx" = "steelblue1")) +
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
        panel.background = element_blank(), legend.position = "none")+
  labs(x ="Prognosis", y = "Average expression (log)", title = "Classic Complement Genes") +
  stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("GPx_A","PPx_A")))

pdf(file="",width = 3,height = 3.5)
box_avg1
dev.off()

pdf(file="",width = 3,height = 3.5)
box_avg2
dev.off()

######### Individual boxplots for all classic complement genes (Supplementary Figure 6a) ######### 
for (gene in geneList) {
  logNormCounts <- assayDataElement(geomx_df[gene,], elt = "log_q")
  colnames(logNormCounts) <- geomx_df$ROI_REL
  logNormCounts <- melt(logNormCounts)
  colnames(logNormCounts) <- c("gene","ROI_REL","log_expression")
  meta <- data.frame(ROI_REL = pData(geomx_df)$ROI_REL,
                     Epithelial = factor(pData(geomx_df)$Epithelial,levels = c("A","I")),
                     Prognosis = pData(geomx_df)$Prognosis,
                     Patient = pData(geomx_df)$`slide name`)
  test <- merge(logNormCounts,meta,by = "ROI_REL")
  test$Prog_EP <- paste0(test$Prognosis,"_",test$Epithelial) %>% factor()
  
  # Create a box plot
  box2[[i]] <- ggplot(test, aes(x = Prog_EP, y = log_expression, fill = Prog_EP)) + 
    geom_boxplot() +
    geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Prog_EP), size = 1) + 
    scale_fill_manual(values = c("GPx_A" = "pink4", "GPx_I" = "lightpink","PPx_A" = "steelblue1", "PPx_I" = "lightskyblue2"), name = "Prog_EP") +
    scale_color_manual(values = c("GPx_A" = "pink4", "GPx_I" = "lightpink","PPx_A" = "steelblue1", "PPx_I" = "lightskyblue2")) +
    theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
          axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
          panel.background = element_blank(), legend.position = "none")+
    labs(x ="Prognosis", y = "Expression (log)", title = paste(gene)) +
    stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("GPx_A","GPx_I"),c("PPx_A","PPx_I"),c("GPx_A","PPx_A"),c("GPx_I","PPx_I")))
  
  i = i + 1
}

pdf(file = "",width = 5*4, height = 4*4)
grid.arrange(box2[[1]],box2[[2]],box2[[3]],box2[[4]],box2[[5]],box2[[6]],box2[[7]],box2[[8]],nrow = 3,ncol = 3)
dev.off()

######### Average C3 expression across PanCK+/PanCK- AOIs in Atypia (Supplementary Figure 6b) ###########
iall_aall_neg <- intersect(targetDataSubset_neg$`slide name`[targetDataSubset_neg$Epithelia %in% c("I")], targetDataSubset_neg$`slide name`[targetDataSubset_neg$Epithelia %in% c("A")]) #11
A_I_neg_both <- targetDataSubset_neg[,targetDataSubset_neg$Epithelial %in% c("I","A") & targetDataSubset_neg$`slide name` %in% iall_aall_neg] # only including patients that have both A and I on the slide. 

iall_aall_pos <- intersect(targetDataSubset_pos$`slide name`[targetDataSubset_pos$Epithelia %in% c("I")], targetDataSubset_pos$`slide name`[targetDataSubset_pos$Epithelia %in% c("A")]) #11
A_I_pos_both <- targetDataSubset_pos[,targetDataSubset_pos$Epithelial %in% c("I","A") & targetDataSubset_pos$`slide name` %in% iall_aall_pos] # only including patients that have both A and I on the slide. 

ROIs <- intersect(A_I_neg_both$ROI_REL, A_I_pos_both$ROI_REL) # 166 ROIs

# extracting these ROIs from pos and neg 
pos <- A_I_pos_both[,A_I_pos_both$ROI_REL %in% ROIs]
neg <- A_I_neg_both[,A_I_neg_both$ROI_REL %in% ROIs]

# just the complement genes
logNormCounts1 <- assayDataElement(pos[geneList,], elt = "log_q")
colnames(logNormCounts1) <- paste0(pos$ROI_REL,"_",pos$segment)

logNormCounts2 <- assayDataElement(neg[geneList,], elt = "log_q")
colnames(logNormCounts2) <- paste0(neg$ROI_REL,"_",neg$segment)

logNormCounts <- cbind(logNormCounts1,logNormCounts2)
logNormCounts <- melt(logNormCounts)
colnames(logNormCounts) <- c("gene","ROI_REL","log_expression")
logNormCounts$segment <- sub(".*_(PanCK[+-])$", "\\1", logNormCounts$ROI_REL)
logNormCounts$ROI <- sub("_(?!.*_).*", "", logNormCounts$ROI_REL, perl = TRUE)

# group by ROI and gene, and take the average of the two segments 
test <- logNormCounts %>% group_by(ROI,gene) %>% summarise(avg = mean(log_expression))
meta <- pData(targetDataSubset_pos) %>% filter(ROI_REL %in% ROIs) %>% arrange(ROI_REL)

# extracting just C3
c3 <- test[test$gene == "C3",]
all(c3$ROI == meta$ROI_REL) # true
c3 <- c3 %>% merge(meta, by.x = "ROI", by.y = "ROI_REL")
c3$Prog_EP <- paste0(c3$Prognosis,"_",c3$Epithelial) %>% factor()

# PPx A vs GPx A
subset <- c3[c3$Epithelial == "A",]
box_avg3 <- ggplot(subset, aes(x = Prog_EP, y = avg, fill = Prognosis)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Prognosis), size = 1) + 
  scale_fill_manual(values = c("GPx" = "pink4", "PPx" = "steelblue1"), name = "Prognosis") +
  scale_color_manual(values = c("GPx" = "pink4", "PPx" = "steelblue1")) +
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
        panel.background = element_blank(), legend.position = "none")+
  labs(x ="Prognosis", y = "Average expression (log)", title = "C3 expression per ROI") +
  stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("GPx_A","PPx_A")))

pdf(file="",width = 3,height = 3.5)
box_avg3
dev.off()
