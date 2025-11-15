#make sure Rtools is installed(Rtools)
library(R.utils)
library(Seurat)
library(GeomxTools)
library(dplyr)
library(knitr)
library(ggforce)
library(writexl)
library(RColorBrewer)
library(gridExtra)
library(tidygraph)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggplotify)
library(ComplexHeatmap)
library(ggpubr)
library(enrichR)
library(enrichplot)
options(java.parameters = "-Xmx12000m")



load(file="Consolidated_AllPlotting.RData") #from 1_ProgData_QC_forGitHub.R
set1_blue <- brewer.pal(9, "Blues")[2:9]
set2_orange <- brewer.pal(9, "Oranges")[2:9]
# Combine and interpolate the two sets to create a divergent palette for Pre-QC
divergent_blueOrange_palette <- colorRampPalette(c(rev(set1_blue), set2_orange))(32)

set1_green <- brewer.pal(9, "Greens")[2:9]
set2_purple <- brewer.pal(9, "Purples")[2:9]
# Combine and interpolate the two sets to create a divergent palette for Post-QC
divergent_GreenPurple_palette <- colorRampPalette(c(rev(set1_green), set2_purple))(32)

my_color_palette <- c("palegreen", "goldenrod", "lightgoldenrod1", 
                      "blue4", 
                      "mediumpurple1",
                      "dodgerblue2", "lightskyblue1" )


distinct_colors = (c( "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4",
						"#46f0f0", "pink4", "#d2f53c", "#fabebe", "#008080", "#e6beff",
						"#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",
						"#bcf60c", "#808080", "violet", "#000000", "#a9a9a9", "#C68799",
						"#4363d8", "#dcbeff", "#9a6324", "#fffac8", "#469990", "#ff0000",
						"#000075", "#90B6FE" ))
						
						
#-----------------------------------------------------------------------------#
###Plotting Figure 1 and Supplementary Figure 1
#-----------------------------------------------------------------------------#

metData<- pData(ProgData_data)
test00<- as.data.frame(metData %>% group_by(metData$Prognosis, metData$segment) %>% summarise(count=n()))
colnames(test00) <- c("Prognosis", "Segment", "#AOIs")

figs1a=  ggplot(test00, aes(x= Prognosis,y= `#AOIs`, fill= Segment)) +geom_bar(stat="identity", width=.75)+theme_minimal()+
  theme(legend.text= element_text(face="bold",size=8),legend.title=element_text(face="bold",size=8), legend.position = "bottom",
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank(), aspect.ratio = 3/1, axis.line = element_line(color="black") )+
  labs(x = "", y = "#AOIs")+
  scale_fill_manual(values= c("lightsalmon2","lightskyblue" )) +scale_y_continuous(limits=c(0,750))
rm(test00)

metData<- pData(ProgData_data)
test00<- as.data.frame(metData %>% group_by(metData$`slide name`, metData$Prognosis) %>% summarise(count=n()))
colnames(test00) <- c("Slide Name", "Prognosis", "#AOI")

figs1b=  ggplot(test00, aes(x= `Prognosis`,y= `#AOI`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge", color="gray")+
  theme_minimal()+
  theme(legend.position = "none", 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank()) + scale_fill_manual(values= divergent_blueOrange_palette)
rm(test00)

metData<- pData(ProgData_downstream)
test00<- as.data.frame(metData %>% group_by(metData$`slide name`, metData$Prognosis) %>% summarise(count=n()))
colnames(test00) <- c("Slide Name", "Prognosis", "#AOI")

figs1b.1=  ggplot(test00, aes(x= `Prognosis`,y= `#AOI`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge", color="gray")+theme_minimal()+
  theme(legend.position = "none", 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_manual(values= divergent_GreenPurple_palette)
rm(test00)

#Number of  ROIs across regions (might go to supplementary see consolidated_allploting)
metData<- sData(ProgData_data)
test00<- as.data.frame(metData %>% group_by(metData$`slide name`, metData$Region) %>%filter(!is.na(Region)) %>% summarise(Num_ROIs=n_distinct(roi)))
colnames(test00) <- c("Slide Name", "Region","#ROIs") 
fig1c=  ggplot(test00, aes(x= `Region`,y= `#ROIs`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge", color="gray", linewidth=0.3)+theme_minimal()+
  theme(axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color="black"), legend.position = "none")+
  scale_fill_manual(values= divergent_blueOrange_palette)

rm(test00)

metData<- sData(ProgData_downstream)
test00<- as.data.frame(metData %>%  distinct() %>% group_by(metData$`slide name`, metData$Region) %>% summarise(Num_ROIs=n_distinct(roi)))
colnames(test00) <- c("Slide Name", "Region","#ROIs")

fig1c.1=  ggplot(test00, aes(x= `Region`,y= `#ROIs`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge",color="gray", linewidth=0.3)+theme_minimal()+
  theme(axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color="black"), legend.position = "none")+ scale_fill_manual(values= divergent_GreenPurple_palette)

rm(test00)

#Number of epithelia AOIS 
ProgData_datapos = subset(ProgData_data, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK+")
metData <-  sData(ProgData_datapos)
test00<- as.data.frame(metData %>% distinct()%>% group_by(metData$`slide name`, metData$Epithelial) %>% filter(!is.na(Region)) %>% summarise(Num_ROIs=n_distinct(roi)))
colnames(test00) <- c("Slide Name", "Epithelia","#AOIs")


fig1d=  ggplot(test00, aes(x= `Epithelia`,y= `#AOIs`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge", color="gray", linewidth=0.3)+theme_minimal()+
  theme(axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color="black"), legend.position = "none") + scale_fill_manual(values= divergent_blueOrange_palette)
rm(test00)

metData <-  sData(targetDataSubset_pos)
test00<- as.data.frame(metData %>% distinct() %>% group_by(metData$`slide name`, metData$Epithelial) %>% summarise(Num_ROIs=n_distinct(roi)))
colnames(test00) <- c("Slide Name", "Epithelia","#AOIs")


fig1d.1=  ggplot(test00, aes(x= Epithelia,y= `#AOIs`, fill=`Slide Name`)) +geom_bar(stat="identity", position="dodge",color="gray", linewidth=0.3)+theme_minimal()+
  theme(axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color="black"), legend.position = "none")+ scale_fill_manual(values= divergent_GreenPurple_palette)+
  scale_y_continuous(limits = c(0, 25))
rm(test00)

pdf("Figs S1a-1d.pdf", width=6, height=8)
fig1a
fig1b+fig1b.1
fig1c+fig1c.1
fig1d+fig1d.1
dev.off()

###### ------------------Figures 1d, and S1e-------------------------------- ########
library(umap)
library(ggplot2)
compute_umap_data <- function(geomx_df, assay_type, n_neighbors = 15, n_epochs = 200, metric = "cosine") {
  logNormCounts <- assayDataElement(geomx_df, elt = assay_type)
  colnames(logNormCounts) <- pData(geomx_df)$ROI_REL
  logNormCounts <- t(logNormCounts)
  
  prog_label <- data.frame(
    batch      = pData(geomx_df)$batch,
    ROI_REL    = pData(geomx_df)$ROI_REL, 
    Patient    = pData(geomx_df)$`slide name`,
    Epithelial = pData(geomx_df)$Epithelial,
    Region     = pData(geomx_df)$Region,
    segment    = pData(geomx_df)$segment,
    ROI        = as.character(sData(geomx_df)$roi),
    Prognosis  = pData(geomx_df)$Prognosis
  )
  
  config <- umap::umap.defaults
  config$n_neighbors <- n_neighbors
  config$n_epochs <- n_epochs
  config$random_state <- 40
  config$metric <- metric
  
  umap_out <- umap::umap(logNormCounts, config = config)
  umap_coords <- as.data.frame(umap_out$layout)
  colnames(umap_coords) <- c("V1", "V2")
  
  umap_df <- cbind(umap_coords, prog_label)
  return(umap_df)
}
plot_umap_from_df <- function(umap_df, label_type, col = NULL, main = "", legend.position = "right") {
  p <- ggplot(umap_df, aes_string(x = "V1", y = "V2", color = label_type)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = col) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.position = legend.position,
      legend.title = element_text(size=15, face="bold"),
      legend.text = element_text(size=15, face ="bold"),
      text = element_text(face = "bold"),
      title = element_text(size = 30)
    ) +
    labs(x = "UMAP1", y = "UMAP2", title = main)
  return(p)
}

# Step 1: Compute UMAP once
umap_data_ProgData <- compute_umap_data(ProgData_downstream, "log_q", n_neighbors = 50, n_epochs = 500, metric = "euclidean")

# Step 2: Generate multiple plots with different coloring
umap1_umap_all<- plot_umap_from_df(umap_data_ProgData, label_type = "Patient", main = "ProgData", legend.position = "none", col = distinct_colors)
umap1_umap_seg<- plot_umap_from_df(umap_data_ProgData, label_type = "segment", main = "ProgData", legend.position = "right", 
                                   col = c("PanCK-"= "mediumpurple3", "PanCK+"= "aquamarine3"))
umap1_umap_prog <- plot_umap_from_df(umap_data_ProgData, label_type = "Prognosis", main = "ProgData", legend.position = "right", 
                                     col = c("GPx" = "lightpink2", "PPx" = "lightskyblue2"))
umap1_umap_reg <- plot_umap_from_df(umap_data_ProgData, label_type = "Region", main = "ProgData", legend.position = "right",
                                     col =  c("C" = "khaki3", "E" = "wheat", "F"="lightsalmon"))


umap_data_gpx <- compute_umap_data(targetDataSubset_gpx, "log_q", n_neighbors = 50, n_epochs = 500, metric = "euclidean")
umap1_umap_gpx<- plot_umap_from_df(umap_data_gpx, label_type = "segment", main = "GPx", legend.position = "right", 
                                   col = c("PanCK-"= "mediumpurple3", "PanCK+"= "aquamarine3"))
umap1_gpx_patient <- plot_umap_from_df(umap_data_gpx, label_type = "Patient", main = "GPx", legend.position = "none", 
                                       col = distinct_colors[1:17])



umap_data_ppx <- compute_umap_data(targetDataSubset_ppx, "log_q", n_neighbors = 50, n_epochs = 500, metric = "euclidean")
umap1_umap_ppx<- plot_umap_from_df(umap_data_ppx, label_type = "segment", main = "PPx", legend.position = "right", 
                                   col = c("PanCK-"= "mediumpurple3", "PanCK+"= "aquamarine3"))
umap1_ppx_patient <- plot_umap_from_df(umap_data_ppx, label_type = "Patient", main = "PPx", legend.position = "none", 
                                       col = distinct_colors[18:32])

pdf(file="UMAPs_Fig1dAndFigs1e.pdf", width= 15*2, height=6*2 )
umap1_umap_all+umap1_umap_prog+umap1_umap_reg+umap1_umap_seg
umap1_umap_gpx+umap1_gpx_patient
umap1_umap_ppx+umap1_ppx_patient
dev.off()

###### ------------------Figure S1f-------------------------------- ########
#Correlation dendrogram of samples for GPx and PPx 
#Sample correlation using mean of PanCK+ and PanCK- to estimate sample heterogeneity, we do only with those ROIs that have both channels present

load(file="DataSubsets_fordownstream.RData")#From 1_Batch10_QC_forGitHub.R

colors <- c("#CB83F2","#F4899B","#F6AAD1","#C68799","#E0ADF8",
            "#D18399","#D3D6CB","#9E82F9","#C4E8A9",
            "#8FC2C6","#E9FC96","#EDA7AF","#9EA28F","#8FDCD4","#EEB7D6","#CAA2B2","#A4E0FF","#C3DEA7","#D1A1AD",
            "#CEB5DE","#8DA4E7","#B9FAA7","#90B6FE","#84E4D2","#8B85A2","#C0B7CD","#C2E58B","#B0D492","#8B8480",
            "#E298B2","#C594AE","#C0A2FB")

{
###For GPx
protocolData(targetDataSubset_gpx)[["roi_seg"]] <- paste0(sData(targetDataSubset_gpx)$`slide name`,".",sData(targetDataSubset_gpx)$roi,".", sData(targetDataSubset_gpx)$segment)
pData(targetDataSubset_gpx)$roi_seg <- sData(targetDataSubset_gpx)$roi_seg

colnames(targetDataSubset_gpx)<- targetDataSubset_gpx$roi_seg
df= assayDataElement(targetDataSubset_gpx, elt="log_q")

# Extract column names until the second dot
column_names <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", colnames(df)) #604 AOIs
# Get unique prefixes in column names
unique_prefixes <- unique(column_names) #334 ROIs


# Filter out column names with a pair
paired_columns <- unique_prefixes[sapply(unique_prefixes, function(prefix) sum(grepl(paste0("^", prefix, "\\."), colnames(df))) == 2) > 0]
#270/334 ROIs have PanCK+ and PanCK-

# Calculate row-wise mean for paired columns
n <- 2
dfpaired <- df[, colnames(df) %in% c(paste0(paired_columns, ".PanCK-"), paste0(paired_columns, ".PanCK+"))]
#unique colnames in teh order of the dfpaired matrix
col2 <- unique ( sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", colnames(dfpaired))) #270

MeanRens <- t(rowsum(t(dfpaired), as.integer(gl(ncol(dfpaired), n, ncol(dfpaired))))) / n 
#This gives the mean of PanCK- and PanCK+ per ROI

colnames(MeanRens) <- col2

#Clustering of the ROIs
# Calculate correlation matrix
cor_matrix <- cor(MeanRens)

# Perform hierarchical clustering using correlation-based distance
dist_matrix <- as.dist(1 - cor_matrix)  # Convert correlation to distance
hclust_result <- hclust(dist_matrix)
patients_gpx<- unique(sub("\\..*", "", hclust_result$labels))

gpx_col_pos <- list()
i=1
for (patient in patients_gpx) {
  ind <- grep(patients_gpx[[i]], hclust_result$labels)
  gpx_col_pos[ind] <- colors[[i]]
  i=i+1
}
gpx_col_pos <-unlist(gpx_col_pos) #Colors for each of the samples identified in hclust_result$labels
dend <- as.dendrogram(hclust_result)%>% set("labels_cex", c(0.4))
labels_colors(dend) <- gpx_col_pos[order.dendrogram(dend)]

# Plot the dendrogram as a circular phylogenetic tree using circlize
# Adjust margin if needed 
par(mar = c(1, 2, 2, 2))
circlize_dendrogram(dend, main = "GPx", labels_track_height = 0.3) 

####doing the same for PPx####

protocolData(targetDataSubset_ppx)[["roi_seg"]] <- paste0(sData(targetDataSubset_ppx)$`slide name`,".",sData(targetDataSubset_ppx)$roi,".", sData(targetDataSubset_ppx)$segment)
pData(targetDataSubset_ppx)$roi_seg <- sData(targetDataSubset_ppx)$roi_seg

colnames(targetDataSubset_ppx)<- targetDataSubset_ppx$roi_seg
df= assayDataElement(targetDataSubset_ppx, elt="log_q")

column_names <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", colnames(df)) #547 AOIs
unique_prefixes <- unique(column_names) #294 ROIs

# Filter out column names with a pair
paired_columns <- unique_prefixes[sapply(unique_prefixes, function(prefix) sum(grepl(paste0("^", prefix, "\\."), colnames(df))) == 2) > 0]
#253/294 ROIs have PanCK+ and PanCK-

# Calculate row-wise mean for paired columns
n <- 2
dfpaired <- df[, colnames(df) %in% c(paste0(paired_columns, ".PanCK-"), paste0(paired_columns, ".PanCK+"))]
col2 <- unique ( sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1", colnames(dfpaired))) #270

MeanRens <- t(rowsum(t(dfpaired), as.integer(gl(ncol(dfpaired), n, ncol(dfpaired))))) / n 
colnames(MeanRens) <- col2

#Clustering of the ROIs
# Calculate correlation matrix
cor_matrix <- cor(MeanRens)

# Perform hierarchical clustering using correlation-based distance
dist_matrix <- as.dist(1 - cor_matrix)  # Convert correlation to distance
hclust_result <- hclust(dist_matrix)
patients_ppx<- unique(sub("\\..*", "", hclust_result$labels))

ppx_col_pos <- list()
i=1
for (patient in patients_ppx) {
  ind <- grep(patients_ppx[[i]], hclust_result$labels)
  ppx_col_pos[ind] <- colors[[i]]
  i=i+1
}
ppx_col_pos <-unlist(ppx_col_pos) #Colors for each of the samples identified in hclust_result$labels
dend_ppx <- as.dendrogram(hclust_result)%>% set("labels_cex", c(0.4))
labels_colors(dend_ppx) <- ppx_col_pos[order.dendrogram(dend_ppx)]

# Plot the dend_ppxrogram as a circular phylogenetic tree using circlize
# Adjust margin if needed 
par(mar = c(1, 2, 2, 2))
circlize_dendrogram(dend_ppx, main = "PPx", labels_track_height = 0.3) 

par(mar = c(1, 2, 2, 2))
pdf("SupplementaryFigure1f.pdf", width=3*1.5, height=3)
circlize_dendrogram(dend, main = "GPx", labels_track_height = 0.4) 
circlize_dendrogram(dend_ppx, main = "PPx", labels_track_height = 0.4) 
dev.off()
}

#-----------------------------------------------------------------------------#
###End of Figure 1 and Figure S1  plotting in R###
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
###Plotting Figure 2, and supplementary Figures 2 and 3
#-----------------------------------------------------------------------------#

###### ------------------Figure 2d-------------------------------- ########

load(file="ProgData_DiffEX7_Neg_CvsE.RData")
#In the Neg space (Fig 2b) DiffEx of E vs C
results.90$results$LabelColor <- ifelse(results.90$results$Estimate > 0, "khaki4", "wheat3")
c90= ggplot(results.90$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                    color = Color, label = Gene)) +
  geom_vline(xintercept = c(0), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast + Subset, scales = "free_y") +
  geom_label_repel(
    data = subset(results.90$results, FDR < 0.05),
    aes(label = Gene),
    color = subset(results.90$results, FDR < 0.05)$LabelColor,
    fill = "white", 
    box.padding = 0.3,
    point.padding = 0.2,
   fontface = "bold",
    size = 5,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
    theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )
pdf("Fig2d.pdf", width=7, height=7)
c90
dev.off()

###### ------------------Figure S2a-------------------------------- ######
results.99$results$LabelColor <- ifelse(results.90$results$Estimate > 0, "khaki4", "wheat3")
c99= ggplot(results.99$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                    color = Color, label = Gene)) +
  geom_vline(xintercept = c(0), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast + Subset, scales = "free_y") +
  geom_label_repel(
    data = subset(results.99$results, FDR < 0.05),
    aes(label = Gene), 
    color = subset(results.99$results, FDR < 0.05)$LabelColor,
    fill = "white",
    box.padding = 0.2,
    point.padding = 0.15,
    fontface = "bold",
    size = 5,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )
pdf("FigS2a_GPxNeg_CvsE.pdf", width=7, height=7)
c99
dev.off()
}

###### ------------------Figure 2e-------------------------------- ######

##-------Plot top 100 genes with the highest average expression in all samples of PanCK-, in GPX and PPx 
# Calculate mean expression for GPx and PPx datasets
mean_expr_GPx <- rowMeans(assayDataElement(GPx_subset_neg, elt="log_q"))
mean_expr_PPx <- rowMeans(assayDataElement(PPx_subset_neg, elt="log_q"))

# Calculate log-fold change (M) and average expression (A) for each gene
M.r <- mean_expr_GPx - (mean_expr_PPx)#since it is already in log2
A.r <- 0.5 * (mean_expr_GPx + mean_expr_PPx)

#remove ribosomal genes
M <- M.r[-c(grep(pattern="^RP", names(M.r)))] 
A <- A.r[-c(grep(pattern="^RP", names(A.r)))]


# Create a dataframe for MA plot
ma_df <- data.frame(M = M, A = A)

# Highlight top 100 genes based on absolute mean expression
top_genes <- head(order(-abs(A)), 100)

#label the top  genes based on absolute mean expression
top_genes_lab <- head(order(-abs(A)), 30)

fig2e<- ggplot(ma_df, aes(x = A, y = M)) +
  geom_point(color = "grey85", size = 1) +  # Background points
  geom_point(data = ma_df[top_genes, ], aes(color = Direction), size = 1.5) +
  scale_color_manual(values = c("GPx" = "lightpink3", "PPx" = "lightblue3")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel(
    data = ma_df[top_genes_lab, ],
    aes(label = rownames(ma_df)[top_genes_lab]),
    fill = "white",
    color = ma_df[top_genes_lab, ]$Color,
    size = 3.5,
    fontface = "bold",
    label.size=0.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE, max.overlaps = 100) +
   labs(
    x = expression("Average Expression (A)"),
    y = expression("Log"[2] * " Fold Change (M)"),
    title = "Differential Expression: GPx vs PPx (PanCK-)") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none")

pdf("Fig2e.pdf", width=6, height=6)
fig2e
dev.off()

###### ------------------Figure 2f,2g S3c, S3d and S3e-------------------------------- ######

#Please see 4_PlottingSpatialDecon_onlyMajorCellTypes_forGitHub.R
#Please see 4_PlottingSpatialDecon_byRegion_forGitHub.R

#-----------------------------------------------------------------------------#
###End of Figure 2 and Figures S2 and S3  plotting in R###
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
###Plotting Figure 3 and supplementary figures 4
#-----------------------------------------------------------------------------#

##Figure 3b
load(file="4.2_TEnr_DiffEx.RData") #from 3_DiffExForGitHub.R

  fig3b= ggplot(results$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                    color = Color, label = Gene)) +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = "log2(FC)",
         y = "Significance, ‐log10(P)",
         color = "Significance") +
    scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                  `FDR < 0.05` = "dodgerblue",
                                  `P < 0.05` = "#F6A958",
                                  `NS` = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_text_repel(data = subset(results$results, Gene %in% results$top_g & FDR < 0.05),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,
                    max.overlaps = 50) +
    theme_bw(base_size = 16) +
    theme(legend.text= element_text(face="bold",size=8),legend.title= element_blank(), 
          axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12),
          panel.background = element_blank(), panel.grid = element_blank(), legend.position = "bottom") +
    facet_wrap(~Contrast+Subset, scales = "free_y")

  
  
#Figure 3c see 5_PlottingEcotyperResultsForGitHub.R 

#Figure 3d, 4b see xxxxxxxx- Dasha

#-----------------------------------------------------------------------------#
###End of Figures 3 and supplementary figures 3 & 4  plotting in R###
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
###Plotting Figure 4 and supplementary figure 5
#-----------------------------------------------------------------------------#

#-------Does the region influence epithelial state?
metData4 <- pData(GPx_subset_pos)
test4<- as.data.frame(metData4 %>% group_by(metData4$Epithelial,metData4$Region) %>% summarise(count=n()))
colnames(test4) <- c("Epithelia", "Region" , "#AOIs in GPx")
figa=ggplot(test4, aes(x= Epithelia,y= `#AOIs in GPx`, fill= Region)) +geom_bar(stat="identity", position="dodge")+ theme_minimal()+
  theme(legend.position = "none", 
        axis.text = element_text(angle = 90, size=10, face="bold"), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank()) + scale_fill_manual(values=  c("C" = "khaki3", "E" = "wheat", "F"="lightsalmon"))

metData5 <- pData(PPx_subset_pos)
test5<- as.data.frame(metData5 %>% group_by(metData5$Epithelial,metData5$Region) %>% summarise(count=n()))
colnames(test5) <- c("Epithelia", "Region" , "#AOIs in PPx")
figa.1=ggplot(test5, aes(x= Epithelia, y= `#AOIs in PPx`, fill= Region)) +geom_bar(stat="identity", position="dodge")+ theme_minimal()+
  theme(
        axis.text = element_text(angle = 90, size=10, face="bold"), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank()) + scale_fill_manual(values= c("C" = "khaki3", "E" = "wheat", "F"="lightsalmon"))

pdf("FigS5a_Supp.pdf", width=3*2, height=3)
figa+figa.1
dev.off()

#Main figure just plot distribution of region in I AOis
metData6 <- pData(targetDataSubset_pos)%>% filter(Epithelial=="I")
test6<- as.data.frame(metData6 %>% group_by(metData6$Prognosis,metData6$Region) %>% summarise(count=n())) 
colnames(test6) <- c("Prognosis", "Region" , "#AOIs with invasive epithelia")


figa.2=ggplot(test6, aes(x= Prognosis,y= `#AOIs with invasive epithelia`, fill= Region)) +geom_bar(stat="identity", position="dodge")+ theme_minimal()+
  theme(axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank(),
        axis.line = element_line(color="black"))+
  scale_fill_manual(values= c("C" = "khaki3", "E" = "wheat", "F"="lightsalmon"))

pdf("Fig4a.pdf", width=3, height=3)
figa.2
dev.off()

###Figure 4b Plotting C vs E volcano plots

load(file="DiffEx9_NMF_GPX_PPX_CvsE_Invasive.RData") #From 3_DiffExForGitHub.R

figS3b= ggplot(results.9$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                     color = Color, label = Gene)) +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title = "GPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "#F6A958",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast+Subset, scales = "free_y")+theme_bw(base_size = 16) +
   theme(
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )


figS3b.1= ggplot(results.10$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                      color = Color, label = Gene)) +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title= "PPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast+Subset, scales = "free_y")+theme_bw(base_size = 16) +
  theme(legend.text= element_text(face="bold",size=8),legend.title= element_blank(), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank(), legend.position = "bottom") +
  facet_wrap(~Contrast+Subset, scales = "free_y")


pdf("FigS3b.pdf", width=4*2, height=6)
figS3b+figS3b.1
dev.off()


#-------------------Plotting 3d DiffEx results------_#

IvsJ_GPX= loadToEnv(file="9_DiffEx_nmfGPX.RData")[["results.7"]]
AvsI_GPX= loadToEnv(file="9_DiffEx_nmfGPX.RData")[["results.8"]]

IvsJ_PPX= loadToEnv(file="9_DiffEx_nmfPPX.RData")[["results.7"]]
AvsI_PPX= loadToEnv(file="9_DiffEx_nmfPPX.RData")[["results.8"]]


IvsJ_PPX.sig <- subset(IvsJ_PPX$results, `FDR` <0.05 & abs(Estimate)>=0.5)
AvsI_PPX.sig <- subset(AvsI_PPX$results, `FDR` <0.05 & abs(Estimate)>=0.5)

IvsJ_GPX.sig <- subset(IvsJ_GPX$results, `FDR` <0.05 & abs(Estimate)>=0.5)
AvsI_GPX.sig <- subset(AvsI_GPX$results, `FDR` <0.05 & abs(Estimate)>=0.5)


# For paper, plotting the top genes identified in IvsJ comparision in GPX and PPX
Gpx_up_ij <- IvsJ_GPX.sig %>% arrange(desc(Estimate))%>% slice_head(n=(20))
Gpx_down_ij <- IvsJ_GPX.sig %>% arrange(Estimate)%>% slice_head(n=(20))

top_g_gpxij <- union(Gpx_up_ij$Gene, Gpx_down_ij$Gene)


Ppx_up_ij <- IvsJ_PPX.sig %>% arrange(desc(Estimate))%>% slice_head(n=(20))
Ppx_down_ij <- IvsJ_PPX.sig %>% arrange(Estimate)%>% slice_head(n=(20))

top_g_ppxij <- union(Ppx_up_ij$Gene, Ppx_down_ij$Gene)

IvsJ_GPX$results$LabelColor <- ifelse(IvsJ_GPX$results$Estimate > 0, "lightpink", "lightpink3")
fig4b=  ggplot(IvsJ_GPX$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                     color = Color, label = Gene)) +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size=1) +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title = "GPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "#F6A958",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast+Subset, scales = "free_y")+
  geom_label_repel(
    data = subset(IvsJ_GPX$results, Gene %in% top_g_gpxij ),
    aes(label = Gene),
    color = subset(IvsJ_GPX$results, Gene %in% top_g_gpxij)$LabelColor,
    fill = "white", 
    box.padding = 0.3,
    point.padding = 0.2,
    fontface = "bold",
    size = 5,
    max.overlaps = 40,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

IvsJ_PPX$results$LabelColor <- ifelse(IvsJ_PPX$results$Estimate > 0, "lightblue2", "steelblue3")
fig4b.1= ggplot(IvsJ_PPX$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                        color = Color, label = Gene)) +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size=1) +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title= "PPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(~Contrast+Subset, scales = "free_y")+  geom_label_repel(
    data = subset(IvsJ_PPX$results, Gene %in% top_g_ppxij),
    aes(label = Gene),
    color = subset(IvsJ_PPX$results,  Gene %in% top_g_ppxij)$LabelColor,
    fill = "white", 
    box.padding = 0.3,
    point.padding = 0.2,
    fontface = "bold",
    size = 5,
    max.overlaps = 40,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )


pdf("Fig4b.pdf", width=5.5*2, height=7.5)
fig4d+fig4d.1
dev.off()


#venn of overlap plots of IvsJ, GPx and PPx
dim(IvsJ_GPX.sig%>% filter(Estimate<0))
IvsJ_PPX.sigup <- 322
IvsJ_PPX.sigdown <- 42
IvsJ_GPX.sigup <- 346
IvsJ_GPX.sigdown <- 64

#Enrichment results supplementary figure 4a

load(file="EnrichmentResults/9_nmf_DiffExEnrich.RData")#Enrichment done using erichR library per vignette
save(pdf="fig4d.pdf")
IvsJ_GPX.enrich$plots$plh_up[[1]]+ IvsJ_PPX.enrich$plots$plh_up[[1]]
dev.off()

#Unique to PPx_I

gpx_unique = setdiff(subset(IvsJ_GPX.sig, Estimate>0)$Gene, subset(IvsJ_PPX.sig, Estimate>0)$Gene)
ppx_unique = setdiff(subset(IvsJ_PPX.sig, Estimate>0)$Gene, subset(IvsJ_GPX.sig, Estimate>0)$Gene)

gpxunique.enr <- enrichr(gpx_unique,databases = "MSigDB_Hallmark_2020")
ppxunique.enr <- enrichr(ppx_unique, databases = "MSigDB_Hallmark_2020")

setwd(figures)
pdf(file="Figs4c.pdf", width=7, height= 3)
plotEnrich(gpxunique.enr$MSigDB_Hallmark_2020, showTerms = 10, title= "GPx-I unique")+ theme(text=element_text(face="bold"), panel.grid = element_blank())
dev.off()

pdf(file="Figs4c.pdf", width=5, height= 3)
plotEnrich(ppxunique.enr$MSigDB_Hallmark_2020, showTerms = 10, title= "PPx-I unique")+ theme(text=element_text(face="bold"), panel.grid = element_blank())
dev.off()


#####Figure 4f- Boxplots
# 54 gene signature, obtained by immune enrichment from PPI, available in GitHub in excel.
PPxpos_nmfsub <- PPx_subset_pos[rownames(ppxPos_meta),] #NMF genes subsetted
GPxpos_nmfsub <- GPx_subset_pos[rownames(gpxPos_meta),]

imm_genes <- unlist(strsplit("ACTR3,ADAR,ANXA2,ARPC2,ARPC4,C1S,CLTA,CLTC,CSTB,CTSB,CTSD,CTSZ,DDOST,FSCN1,GRB2,GSTP1,HEBP2,HLA-B,HLA-C,HNRNPA2B1,HSP90AA1,HSPA1B,IFI27,IFI30,IFI6,IL32,ITGB1,KPNB1,LAMP1,LAMP2,LAMTOR2,MSN,MX1,MYH9,NME2,PADI2,PDXK,PFKL,PKM,PLAUR,PRDX4,PSMB1,PSMC2,PSMC4,RBX1,RNF213,S100A9,SDC1,SLC44A2,TAP1,TAPBP,VIM,XRCC6,YWHAB",","))

#ForGPX and PPx using the same invasive 
#(that is it needs to habve both I and J within the slide of interest) we used for calculation of I vs J comparision
metData4 <- pData(GPxpos_nmfsub)
test4<- as.data.frame(metData4 %>% group_by(metData4$Epithelial, metData4$`slide name`) %>% summarise(count=n()))
colnames(test4) <- c("Epithelia", "Slide Name", "#AOI in GPx pos")
iall_jall <- intersect(test4$`Slide Name`[test4$Epithelia %in% c("I")], test4$`Slide Name`[test4$Epithelia %in% c("J")])
#Basically only slides that have both I and J in them

metData4 <- pData(PPxpos_nmfsub)
test4<- as.data.frame(metData4 %>% group_by(metData4$Epithelial, metData4$`slide name`) %>% summarise(count=n()))
colnames(test4) <- c("Epithelia", "Slide Name", "#AOI in PPx pos")
iall_jall.ppx <- intersect(test4$`Slide Name`[test4$Epithelia %in% c("I")], test4$`Slide Name`[test4$Epithelia %in% c("J")]) 
#Basically only slides that have both I and J in them

#subset targetDataPos to only these above slide names and extract the AOIs with invaisve for the boxplots
newData = subset(targetDataSubset_pos,  CodeClass=="Endogenous" | CodeClass == "Negative", `slide name` %in% c(iall_jall,iall_jall.ppx) & Epithelial %in% c("I")) 
imm_epi_log <- assayDataElement(newData[imm_genes,], elt="log_q")
dim(imm_epi_log) #54 282

colnames(imm_epi_log) <- pData(newData)$ROI_REL
imm_epi_log <- as.data.frame(imm_epi_log) %>%rownames_to_column(., var="GeneID") 
imm_epi_log<- melt(imm_epi_log)
colnames(imm_epi_log) <- c("GeneID","ROI_REL","log_expression")

meta <- data.frame(ROI_REL= pData(newData)$ROI_REL,Prognosis = pData(newData)$Prognosis)
test <- merge(imm_epi_log,meta,by = "ROI_REL")

library(ggpubr)

# Calculate the mean expression per sample
mean_expression <- aggregate(log_expression ~ Prognosis + ROI_REL, data = test, FUN = mean)
# Specify the factor levels for Prognosis to match the levels in your data
mean_expression$Prognosis <- factor(mean_expression$Prognosis, levels = c("GPx", "PPx"))

# Plotting
fig4f <- ggplot(mean_expression, aes(x = Prognosis, y = log_expression, fill = Prognosis)) +
  geom_boxplot() +  # Box plot for GPx and PPx
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Prognosis), size = 1) +  # Mean expression points per sample
  scale_fill_manual(values = c("GPx" = "lightpink2", "PPx" = "lightskyblue2"), name = "Prognosis") +
  scale_color_manual(values = c("GPx" = "lightpink2", "PPx" = "lightskyblue2")) +
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
        panel.background = element_blank(), legend.position = "none")+
  labs(x ="", y = "Mean expression (log)", title ="54 immune genes")+
  stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("GPx","PPx")))

pdf(file="fig4f_ImmuneGene_boxplot.pdf", width=3, height=3.5)
fig4f
dev.off()


#ORIS gene signature Figure 4g

#Plotting for the ORIS signature as well
imm_genes1 <- c("BATF","CCR7","CD28","CD37","CD40LG","CD74","CLEC2D","DOCK2","FLT3","FOXP3","FUT7","GATA1","IKZF1","IL10","IL12B","ITGB2","KLRK1","LYN","NCOR1","NLRP12","PYCARD","SLC11A1","SPI1","SPIB","STAT5A","STAT5B","STX11","TCF3","TLR1","TNFSF4","TNFSF8","TSPAN32")
GOI= intersect(imm_genes1, rownames(newData))

newData = subset(targetDataSubset_pos,  CodeClass=="Endogenous" | CodeClass == "Negative", `slide name` %in% c(iall_jall,iall_jall.ppx) & Epithelial %in% c("I")) 
imm_epi_log1 <- assayDataElement(newData[GOI,], elt="log_q")
dim(imm_epi_log1) #23 282

colnames(imm_epi_log1) <- pData(newData)$ROI_REL
imm_epi_log1 <- as.data.frame(imm_epi_log1) %>%rownames_to_column(., var="GeneID") 
imm_epi_log1<- melt(imm_epi_log1)
colnames(imm_epi_log1) <- c("GeneID","ROI_REL","log_expression")

meta <- data.frame(ROI_REL= pData(newData)$ROI_REL,Prognosis = pData(newData)$Prognosis)
test <- merge(imm_epi_log1,meta,by = "ROI_REL")

library(ggpubr)

# Calculate the mean expression per sample
mean_expression <- aggregate(log_expression ~ Prognosis + ROI_REL, data = test, FUN = mean)
# Specify the factor levels for Prognosis to match the levels in your data
mean_expression$Prognosis <- factor(mean_expression$Prognosis, levels = c("GPx", "PPx"))

# Plotting
fig4g <- ggplot(mean_expression, aes(x = Prognosis, y = log_expression, fill = Prognosis)) +
  geom_boxplot() +  # Box plot for GPx and PPx
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(color = Prognosis), size = 1) +  # Mean expression points per sample
  scale_fill_manual(values = c("GPx" = "lightpink2", "PPx" = "lightskyblue2"), name = "Prognosis") +
  scale_color_manual(values = c("GPx" = "lightpink2", "PPx" = "lightskyblue2")) +
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12), axis.line = element_line(color="black"),
        panel.background = element_blank(), legend.position = "none")+
  labs(x ="", y = "Mean expression (log)")+
  stat_compare_means(method = "t.test", label = "p.signif", position="identity", comparisons = list(c("GPx","PPx")))

pdf(file="fig4g_ImmuneGene_boxplot.pdf", width=3, height=3.5)
fig4g
dev.off()


#-----------------------------------------------------------------------------#
###End of Plotting Figure 4 and supplementary figure 5
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
###Plotting Figure 5 and supplementary figure 6
#-----------------------------------------------------------------------------#


#Figure 5 A vs I
# For paper, plotting the top genes identified in IvsJ comparision in GPX and PPX

Gpx_up_ia <- AvsI_GPX.sig %>% filter(Estimate > 0)
Gpx_down_ia <- AvsI_GPX.sig %>% filter(Estimate < 0)

top_g_gpxia <- union(Gpx_up_ia$Gene, Gpx_down_ia$Gene)

Ppx_up_ia <- AvsI_PPX.sig %>% arrange(desc(Estimate))%>% slice_head(n=(20))
Ppx_down_ia <- AvsI_PPX.sig %>% arrange(Estimate)%>% slice_head(n=(20))

top_g_ppxia <- union(Ppx_up_ia$Gene, Ppx_down_ia$Gene)

AvsI_GPX$results$LabelColor <- ifelse(AvsI_GPX$results$Estimate > 0,"pink4", "lightpink")

figs5a.1= ggplot(AvsI_PPX$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                   color = Color, label = Gene)) +
  #geom_hline(yintercept = -log10(0.05), lty = "dashed")+ geom_vline(xintercept = 0, lty = "dashed") +
  geom_point(size=1) +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title= "PPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  geom_label_repel(
    data = subset(AvsI_PPX$results, FDR < 0.05 ),
    aes(label = Gene),
    color = subset(AvsI_PPX$results, FDR < 0.05)$LabelColor, 
    fill = "white", 
    box.padding = 0.3,
    point.padding = 0.2,
    fontface = "bold",
    size = 5,
    max.overlaps = 80,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  ) 


fig5a= ggplot(AvsI_PPX$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                      color = Color, label = Gene)) +
  #geom_hline(yintercept = -log10(0.05), lty = "dashed")+ geom_vline(xintercept = 0, lty = "dashed") +
  geom_point(size=1) +
  labs(x = "log2(FC)",
       y = "Significance, ‐log10(P)",
       color = "Significance", title= "PPx") +
  scale_color_manual(values = c(`FDR < 0.001` = "red4",
                                `FDR < 0.05` = "dodgerblue",
                                `P < 0.05` = "orange2",
                                `NS` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  geom_label_repel(
    data = subset(AvsI_PPX$results, FDR < 0.05 & abs(Estimate)>=0.5),
    aes(label = Gene),
    color = subset(AvsI_PPX$results, FDR < 0.05 & abs(Estimate)>=0.5)$LabelColor,
    fill = "white", 
    box.padding = 0.3,
    point.padding = 0.2,
    fontface = "bold",
    size = 5,
    max.overlaps = 80,
    show.legend = FALSE
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  ) 

setwd(figures)
pdf("Fig5.pdf", width=7.5, height=7.5)
fig5a
figs5a.1
dev.off()

#Just want the enrichment of AvsI upregulated genes 
g= AvsI_PPX.sig %>% filter(Estimate>0)
ppxxunique_avsi.enr <- enrichr(g$Gene,databases = "MSigDB_Hallmark_2020")

setwd(figures)
pdf(file="Fig5b.pdf", width=5.5, height= 4.5)
plotEnrich(ppxxunique_avsi.enr$MSigDB_Hallmark_2020, showTerms = 10, title= "PPx-A ")+ theme(text=element_text(face="bold", size=16), panel.grid = element_blank())
dev.off()


#####Figure 5c and supplementary figures 6b,c  please see 5_C3_boxplots_forGitHub.R####
#####Figures 5f and 3e please see 6_BulkSignalR_forGitHub.R######


#----------------------------------------------------------------------------------------#
#End of all Figures#
#----------------------------------------------------------------------------------------#