library(NanoStringNCTools)
library(GeomxTools)
library(dplyr)
library(knitr)
library(ggforce)
library(writexl)
library(R.utils)

setwd(analysis)
GPx_deconPropMat <- loadToEnv(file="4_GPx_Quantiles.RData") [["reordered_matrix0"]]
PPx_deconPropMat <- loadToEnv(file="4_PPx_Quantiles.RData") [["reordered_matrix0"]]

my_color_palette <- c("palegreen", "goldenrod", "lightgoldenrod1", 
                      "blue4", 
                      "mediumpurple1",
                      "dodgerblue2", "lightskyblue1" )
# Define the desired row name order
desired_order <- c("Endothelial","Fibroblast/Mesenchymal", "PVL", 
                   "T-cells",
                   "Myeloid",
                   "B-cells","Plasmablasts")


#Also, re-setting the order of the GPx_deconPropMat, based on center, edge and F (Post Reviewer suggestion)

region_labelsG <- sapply(strsplit(colnames(GPx_deconPropMat), "\\."), function(x) x[3])
region_orderG <- factor(region_labelsG, levels = c("E", "C", "F"))  # or c("E", "C", "F")
# Get order of columns by region
oG <- order(region_orderG)

# Transpose and attach region labels
prop_dfG <- as.data.frame(GPx_deconPropMat)
prop_df_tG <- as.data.frame(t(prop_dfG))
prop_df_tG$Region <- region_labelsG
prop_df_tG$Group <- "GPx"
prop_df_tG$Group_Region <- paste0(prop_df_tG$Group,"_", prop_df_tG$Region)

# Melt into long format
long_dfG <- melt(prop_df_tG, id.vars = c("Region","Group", "Group_Region" ), variable.name = "CellType", value.name = "Proportion")

# Step 3: Define region color mapping
region_colors <- c("C" = "khaki", "F" = "lightsalmon", "E" = "wheat")

#-----------------Doing the same for PPx-----------#####
#Also, re-setting the order of the PPx_deconPropMat, based on center, edge and F (Post Reviewer suggestion)

region_labels <- sapply(strsplit(colnames(PPx_deconPropMat), "\\."), function(x) x[3])
region_order <- factor(region_labels, levels = c("E", "C", "F"))  # or c("E", "C", "F")
# Get order of columns by region
o <- order(region_order)

# Transpose and attach region labels
prop_df <- as.data.frame(PPx_deconPropMat)
prop_df_t <- as.data.frame(t(prop_df))
prop_df_t$Region <- region_labels
prop_df_t$Group <- "PPx"
prop_df_t$Group_Region <- paste0(prop_df_t$Group,"_", prop_df_t$Region)


# Melt into long format
long_df <- melt(prop_df_t, id.vars = c("Region","Group","Group_Region"), variable.name = "CellType", value.name = "Proportion")

# Step 3: Define region color mapping
region_colors <- c("C" = "khaki3", "E" = "wheat2", "F" = "lightsalmon")

######-----------Plotting side by side-------------#######

combined_df = rbind(long_df, long_dfG)

group_region_colors <- c(
  "GPx_C" = "khaki",        # Center - GPx
  "PPx_C" = "khaki3",       # Center - PPx
  "GPx_E" = "wheat",        # Edge - GPx
  "PPx_E" = "wheat3",       # Edge - PPx
  "GPx_F" = "lightsalmon",  # Foci - GPx
  "PPx_F" = "lightsalmon3"  # Foci - PPx
)

library(ggplot2)

pdf("Fig2g_s3d.pdf", width=3*3, height=3*3)

ggplot(combined_df, aes(x = Region, y = Proportion, fill = Group_Region)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.size = 0.5, alpha = 0.9) +
  facet_wrap(~ CellType, scales = "free_y", ncol=4) +
  scale_fill_manual(values = group_region_colors) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom") +
  labs(title = "Cell Type Proportions by Region and Group",
       x = "Region", y = "Proportion", fill = "Group_Region") +
  stat_compare_means(
         aes(group = Group_Region), 
         method = "wilcox.test",
         label = "p.signif",color="red",
         hide.ns = TRUE
       )

dev.off()


