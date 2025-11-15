library(NanoStringNCTools)
library(GeomxTools)
library(dplyr)
library(knitr)
library(ggforce)
library(writexl)
library(R.utils)

my_color_palette <- c("palegreen", "goldenrod", "lightgoldenrod1", 
                      "blue4", 
                      "mediumpurple1",
                      "dodgerblue2", "lightskyblue1" )
# Define the desired row name order
desired_order <- c("Endothelial","Fibroblast/Mesenchymal", "PVL", 
                   "T-cells",
                   "Myeloid",
                   "B-cells","Plasmablasts")

#Plotting the frequencies GPx
setwd(analysis)
load(file="4_GPx_Quantiles.RData") #From 4_SpatialDecon_onlyMajorCellTypes_forGitHub.R
t= ggplot(data = t6.gpx, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "blue4", color = "grey") +
  labs( title= "Frequency of cell proportions, GPx", x = "T-cells", y = "Frequency") +
  theme_minimal() + ylim(0,100)+xlim(0,0.75) + 
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=12, face="bold"), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  ggplot(data = b6.gpx, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "dodgerblue", color = "grey") +
  labs( x = "B-cells", y = "Frequency") +
  theme_minimal()+ylim(0,100)+xlim(0,.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  ggplot(data = pb6.gpx, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "lightskyblue1", color = "grey") +
  labs( x = "Plasmablasts", y = "Frequency") +
  theme_minimal() +ylim(0,100)+xlim(0,0.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  ggplot(data = myl6.gpx, aes(x = value)) + geom_histogram(binwidth = 0.05, fill = "slateblue", color = "grey") +
  labs( x = "Myeloid", y = "Frequency") +
  theme_minimal() +ylim(0,100)+xlim(0,0.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  ggplot(data = newmes6.gpx, aes(x = value)) +geom_histogram(binwidth = 0.05, fill = "goldenrod", color = "grey") +
  labs( x = "Mesenchymal", y = "Frequency") +
  theme_minimal() +ylim(0,100)+xlim(0,1)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  plot_layout(nrow = 1, ncol = 5, guides = "collect")

setwd(figures)
pdf(file="figs3c_suppQuantiles_gpx.pdf", width= 3*4, height=3)
t
dev.off()

GPx_deconPropMat <- reordered_matrix0

setwd(figures)
o = hclust(dist(t(GPx_deconPropMat)),method="average")$order
par(mar=c(10, 3, 1, 1) + 0.1)
TIL_barplot.km(GPx_deconPropMat[,o], cex.names = 0.75,draw_legend = F)
TIL_barplot.km(GPx_deconPropMat[,o], cex.names = 0.75,draw_legend = T)


pdf(file="Fig2f_gpx.pdf", width=3*6, height = 3*2)
par(mar=c(10, 3, 1, 1) + 0.1)
TIL_barplot.km(GPx_deconPropMat[,o], cex.names = 0.75,draw_legend = F)
TIL_barplot.km(GPx_deconPropMat[,o], cex.names = 0.75,draw_legend = T)
dev.off()

#Plotting only the Tenr subset plots
subset_gpx_tenr$ROI_REL
gpx_tenr_plot = GPx_deconPropMat[,subset_gpx_tenr$ROI_REL]

pdf(file="figs3e_gpx_tenr.pdf", width=3*4, height = 3*2)
TIL_barplot.km(gpx_tenr_plot, cex.names = 0.5,draw_legend = F)
dev.off()


#-----------------------------#
setwd(analysis)
load(file="4_PPx_Quantiles.RData")

t2= ggplot(data = t6, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "blue4", color = "grey") +
  labs( title= "Frequency of cell proportions, PPx", x = "T-cells", y = "Frequency") +
  theme_minimal() + ylim(0,75)+xlim(0,0.75) + 
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  ggplot(data = b6, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "dodgerblue", color = "grey") +
  labs( x = "B-cells", y = "Frequency") +
  theme_minimal()+ylim(0,75)+xlim(0,0.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  ggplot(data = pb6, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "lightskyblue1", color = "grey") +
  labs( x = "Plasmablasts", y = "Frequency") +
  theme_minimal() +ylim(0,75)+xlim(0,0.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  ggplot(data = myl6, aes(x = value)) + geom_histogram(binwidth = 0.05, fill = "slateblue", color = "grey") +
  labs( x = "Myeloid", y = "Frequency") +
  theme_minimal() +ylim(0,75)+xlim(0,0.75)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  ggplot(data = newmes6, aes(x = value)) +geom_histogram(binwidth = 0.05, fill = "goldenrod", color = "grey") +
  labs( x = "Mesenchymal", y = "Frequency") +
  theme_minimal() +ylim(0,75)+xlim(0,1)+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), 
        axis.text = element_text(size=10, face="bold"), axis.title = element_text(face="bold",size=10),
        panel.background = element_blank(), panel.grid = element_blank())+
  
  
  plot_layout(nrow = 1, ncol = 5, guides = "collect")

PPx_deconPropMat <- reordered_matrix0

o1 = hclust(dist(t(PPx_deconPropMat)),method="average")$order
par(mar=c(10, 3, 1, 1) + 0.1)

setwd(figures)
pdf(file="fig2c_ppx.pdf", width=3*6, height = 3*2)
TIL_barplot.km(PPx_deconPropMat[,rev(o1)], cex.names = 0.75,draw_legend = F)
TIL_barplot.km(PPx_deconPropMat[,rev(o1)], cex.names = 0.75,draw_legend = T)
dev.off()

#Plotting only the Tenr subset plots
subset_ppx_tenr$ROI_REL
ppx_tenr_plot = PPx_deconPropMat[,subset_ppx_tenr$ROI_REL]

pdf(file="figs3e_ppx_tenr.pdf", width=3*3, height = 3*2)
TIL_barplot.km(ppx_tenr_plot, cex.names = 0.5,draw_legend = F)
dev.off()


pdf(file="figs3c_ppx.pdf", width= 3*4, height=3)
t2
dev.off()



