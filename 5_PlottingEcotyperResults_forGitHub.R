library(NanoStringNCTools)
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

my_color_palette <- c( "steelblue", "hotpink","paleturquoise4","gold2","turquoise","plum2","lawngreen","#8B8480","tomato")

####------------Now consolidating the States computed for GPX by cell type------_############


pathToFile<- "StateAssigned_GPx.xlsx" #Please use the path where you placed the state assignments after running Ecotyper
library(rio)
data<- import_list(pathToFile)

GPX_tenr_cellState <- do.call(rbind, data)

final_data <- GPX_tenr_cellState%>%
  filter(!(cellType %in% c("Dendritic.cells")))%>% mutate(patient_ID = str_extract(...1, "^[^.]+"))#removing cell types that we don't have annotation for in PanCK-

#Plotting for GPx
final_data$cellType <- factor(final_data$cellType, levels = c("B.cells", "PCs", "CD4.T.cells","CD8.T.cells",
                                                              "Monocytes.and.Macrophages", 
                                                              "Fibroblasts","Endothelial.cells"))



final_data <- final_data %>% mutate(cellType = recode(cellType,
                           "B.cells" = "B-cells",
                           "Monocytes.and.Macrophages" = "Myeloid",
                           "Endothelial.cells"="Endothelial",
                           "PCs" ="Plasmablasts"))
                      
fig3c= ggplot(final_data, aes(x=cellType, fill=State)) + geom_bar(position="fill", width=0.75) + theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper, GPx-Tenr", y = "Proportion")+
  theme(legend.text= element_text(size=12, face="bold"),legend.title=element_text(size=12, face="bold"), 
        axis.text.x = element_text(size=14, angle=90, face="bold"), axis.title = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"),
                panel.background = element_blank(), panel.grid = element_blank())

# Rotating x-axis text by 90 degrees

a=ggplot(final_data, aes(x=cellType, fill=State)) + geom_bar() + theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper", y = "Proportion", title = "GPx-Tenr")+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), axis.text= element_text(face="bold"),
        axis.text.x = element_text(size=12, angle=45, vjust =0.75), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank())

b=ggplot(final_data, aes(x=cellType, fill=State)) + geom_bar()+ facet_grid (~patient_ID)+  theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper", y = "Proportion", title = "GPx-Tenr")+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), axis.text= element_text(face="bold"),
        axis.text.x = element_text(size=12, angle=45, vjust =0.75), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank())


#-------------Doing the same for PPx---------#

pathToFile<- "StateAssigned_PPx.xlsx"#Please use the path where you placed the state assignments after running Ecotyper
library(rio)
data<- import_list(pathToFile)
PPX_tenr_cellState <- do.call(rbind, data)

final_data1 <- PPX_tenr_cellState%>%
  filter(!(cellType %in% c("Dendritic.cells")))%>% mutate(patient_ID = str_extract(...1, "^[^.]+"))#removing cell types that we don't have annotation for in PanCK-

#Plotting for PPx
final_data1$cellType <- factor(final_data1$cellType, levels = c("B.cells", "PCs", "CD4.T.cells","CD8.T.cells",
                                                              "Monocytes.and.Macrophages", 
                                                              "Fibroblasts","Endothelial.cells"))


final_data1 <- final_data1 %>% mutate(cellType = recode(cellType,
                                                      "B.cells" = "B-cells",
                                                      "Monocytes.and.Macrophages" = "Myeloid",
                                                      "Endothelial.cells"="Endothelial",
                                                      "PCs" ="Plasmablasts"))

fig3c.1= ggplot(final_data1, aes(x=cellType, fill=State)) + geom_bar(position="fill", width=0.75) + theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper, PPx-Tenr", y = "Proportion")+
  theme(legend.text= element_text(size=12, face="bold"),legend.title=element_text(size=12, face="bold"), 
        axis.text.x = element_text(size=14, angle=90, face="bold"), axis.title = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=14, face="bold"),
        panel.background = element_blank(), panel.grid = element_blank())

# Rotating x-axis text by 90 degrees

c=ggplot(final_data1, aes(x=cellType, fill=State)) + geom_bar() + theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper", y = "Proportion", title = "PPx-Tenr")+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), axis.text= element_text(face="bold"),
        axis.text.x = element_text(size=12, angle=45, vjust =0.75), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank())

d=ggplot(final_data1, aes(x=cellType, fill=State)) + geom_bar()+ facet_grid (~patient_ID)+  theme_minimal() +
  scale_fill_manual(values=my_color_palette)+ labs(x = "Cell type identified in Ecotyper", y = "Proportion", title = "PPx-Tenr")+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), axis.text= element_text(face="bold"),
        axis.text.x = element_text(size=12, angle=45, vjust =0.75), axis.title = element_text(face="bold",size=12),
        panel.background = element_blank(), panel.grid = element_blank())

