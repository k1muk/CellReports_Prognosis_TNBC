library(GeomxTools)
library(reshape2)
library(dplyr)
library(scran)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(writexl)
library(readxl)
library(ComplexHeatmap)
library(BulkSignalR)

load("5_DataSubsets_fordownstream.RData")

# Set seed to get reproducible results
set.seed(123)
############## Important functions for plotting chord diagrams #########
chordDiagramLR_me <- function(dataframe.bsrinf,filename, pw.id.filter,limit,width,height, pair.to.highlight) {
  
  dataframe.bsrinf <- dataframe.bsrinf[order(dataframe.bsrinf$corr.abs,decreasing=TRUE),]
  if (limit < nrow(dataframe.bsrinf)) {
    dataframe.bsrinf <- dataframe.bsrinf[1:limit,] }
  dataframe.bsrinf <- unique(dataframe.bsrinf[,c("ligands","receptors","corr","pair","cell.type")])
  
  print(dataframe.bsrinf)
  
  cr <- circlize::colorRamp2(c(0.3,0.6,0.9),c("#FDE725FF","#21908CFF","#440154FF"))
  # function that plots correlation range onto the chosen color range
  
  myList.ligands <- rep("gray25",times=length(dataframe.bsrinf$ligands))
  names(myList.ligands)  <- as.list(dataframe.bsrinf$ligands)
  
  myList.receptors <- rep("indianred2",times=length(dataframe.bsrinf$receptors))
  names(myList.receptors)  <- as.list(dataframe.bsrinf$receptors)
  
  myList <- c(myList.receptors,myList.ligands)
  
  link.col <- rep("dodgerblue3",nrow(dataframe.bsrinf))
  
  link.lwd <- rep(1,nrow(dataframe.bsrinf))
  
  link.width <- rep(0.12,nrow(dataframe.bsrinf))
  
  if (!is.null(pair.to.highlight)){
    index.filter <- which(dataframe.bsrinf$pair %in% pair.to.highlight)
    link.col[index.filter] <- "#d40000"
    link.lwd[index.filter] <- 3
    link.width[index.filter] <- 0.15
  }
  
  
  if (!is.null(filename)){
    grDevices::pdf(paste0("BulkSignalR","/",filename,".pdf"),
                   width=width, height=height)
    par(mar = c(0, 0, 0, 0)) # set pdf margins to be 0 all around (maximize space)
  }
  
  interactions <- data.frame(from=dataframe.bsrinf$ligands,
                             to=dataframe.bsrinf$receptors,
                             value=dataframe.bsrinf$corr,
                             celltype = dataframe.bsrinf$cell.type)
  
  cdm_res <- circlize::chordDiagramFromDataFrame(interactions,
                                                 col=cr,
                                                 annotationTrack = c("grid"),
                                                 grid.col = myList,
                                                 transparency = 0, # whether the background is transparent or not (0 is not)
                                                 preAllocateTracks = 1,
                                                 directional = 1,
                                                 direction.type = "arrows",
                                                 link.arr.length = link.width,
                                                 link.arr.width  = link.width,
                                                 link.arr.type  = "big.arrow", # instead of overlaying triangle arrow, big.arrow turns the link itself into a thick arrow shape
                                                 link.arr.lty = "solid",
                                                 link.arr.lwd = link.lwd,
                                                 link.arr.col = cr(dataframe.bsrinf$corr), # coloring the arrows using the defined cr function 
                                                 big.gap = 2,
                                                 small.gap = 1)
  
  cr_cells <- c("Endothelial"= "palegreen","Fibroblast/Mesenchymal"="goldenrod","PVL"="lightgoldenrod1",
                "T-cells"="blue4", "Myeloid"="mediumpurple","B-cells"="dodgerblue2","Plasmablasts"="lightskyblue1",
                "NA"="white","Epithelial" = "black")
  
  # Define the desired row name order
  desired_order <- c("Endothelial","Fibroblast/Mesenchymal", "PVL",
                     "T-cells",
                     "Myeloid",
                     "B-cells","Plasmablasts")
  
  
  ylim = circlize::get.cell.meta.data("ylim")
  y1 = ylim[1]
  y2 = ylim[2] -circlize::mm_y(1)
  for(i in seq_len(nrow(cdm_res))) { # adding cell type as an additional ring around the plot
    if(cdm_res$value1[i] > 0) {
      circlize::circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45, 
                            col = cr_cells[interactions$celltype[i]],
                            border = cr_cells[interactions$celltype[i]],
                            sector.index = cdm_res$rn[i], track.index = 1)
      circlize::circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value1"]), y1 + (y2-y1)*0.45,
                            col = cr_cells[interactions$celltype[i]],
                            border = cr_cells[interactions$celltype[i]],
                            sector.index = cdm_res$cn[i], track.index = 1)
    }}
  
  circlize::circos.trackPlotRegion(track.index = 2,
                                   panel.fun = function(x, y) {
                                     xlim = circlize::get.cell.meta.data("xlim")
                                     ylim = circlize::get.cell.meta.data("ylim")
                                     sector.name = circlize::get.cell.meta.data("sector.index")
                                     
                                     circlize::circos.text(mean(xlim), ylim[1] + 2.4, sector.name, # ylim[] + __ adjusts how close the text label is to the ring (decrease __ for closer)
                                                           facing = "clockwise", niceFacing = TRUE,
                                                           adj = c(0, 0.5), cex = 0.5) # cex is text size
                                     
                                     circlize::circos.axis(h="top",labels=FALSE,minor.ticks=FALSE,
                                                           major.tick.length=1.2,
                                                           major.at=c(xlim),
                                                           sector.index=sector.name,
                                                           track.index=2)
                                     circlize::circos.par(circle.margin = 0.00001) # these are the margins on the side of the plot, have to be greater than 0
                                   })
  
  
  # LEGEND 
  
  lgd_points = Legend(labels = c("Ligands", "Receptors"),
                      type = "points", pch = 16,
                      legend_gp = grid::gpar(col = c("gray25","indianred2")),
                      title_position = "topleft",
                      labels_gp = grid::gpar( font = 6),
                      title = "LR")
  lgd_points2 = Legend(labels = c("Endothelial","Fibroblast", "PVL/Mesenchymal",
                                  "T-cells",
                                  "Myeloid",
                                  "B-cells","Plasmablasts"), #"Epithelial"),
                       type = "points", pch = 16,
                       legend_gp = grid::gpar(col = c("palegreen", "goldenrod", "lightgoldenrod1",
                                                      "blue4",
                                                      "mediumpurple1",
                                                      "dodgerblue2", "lightskyblue1" )),#"black")),
                       title_position = "topleft",
                       labels_gp = grid::gpar( font = 6),
                       title = "Cell Type")
  
  lgd_links = Legend(at = c(0.3, 0.9), col_fun = cr,
                     title = "Correlation",direction ="horizontal"   ,
                     grid_width = unit(0.9, "mm") ,
                     grid_height = unit(1.3, "mm") ,
                     labels_gp = grid::gpar( font = 6),
                     title_position = "topcenter",
  )
  
  lgd_list_vertical = packLegend(lgd_points)
  
  ComplexHeatmap::draw(lgd_list_vertical, x = unit(0.3, "inch"),
                       y = unit(16, "mm"),
                       just = c("left", "bottom"))
  ComplexHeatmap::draw(lgd_links, x = unit(0.3, "inch"),
                       y = unit(2, "mm"),
                       just = c("left", "bottom"))
  ComplexHeatmap::draw(lgd_points2, x = unit(0.3, "inch"),
                       y = unit(100, "mm"),
                       just = c("left", "bottom"))
  
  circlize::circos.clear()
  
  if (!is.null(filename))
    grDevices::dev.off()
  
}

plot_chord <- function(BSR_inf,L_R,celltypes,pair.to.highlight=NULL,filename,num_genes) {
  dataframe.bsrinf<- data.frame(
    ligands=unlist(ligands(BSR_inf)),
    receptors=unlist(receptors(BSR_inf)),
    corr = LRinter(BSR_inf)$LR.corr,
    corr.abs = abs(LRinter(BSR_inf)$LR.corr),
    pw.id=LRinter(BSR_inf)$pw.id,
    pathways=LRinter(BSR_inf)$pw.name,
    pval=LRinter(BSR_inf)$pval,
    qval=LRinter(BSR_inf)$qval,
    genes=unlist(lapply(tGenes(BSR_inf), function(x) paste0(x,collapse = ",")))
  )
  dataframe.bsrinf$pair <- paste(dataframe.bsrinf$ligands,
                                 dataframe.bsrinf$receptors,sep="-")
  dataframe.bsrinf <- dataframe.bsrinf  %>% filter(pair %in% L_R)
  dataframe.bsrinf <- merge(dataframe.bsrinf,celltypes,by.x="pair",by.y="L-R",all.x=TRUE)
  dataframe.bsrinf[is.na(dataframe.bsrinf)] <- "NA"
  pathways <- unique(dataframe.bsrinf$pw.id)
  
  # visualize with custom chord diagram plotting function
  chordDiagramLR_me(dataframe.bsrinf,
                     filename = filename,
                     pw.id.filter = "",
                     limit = num_genes,
                     pair.to.highlight,
                     width = 2*2,
                     height = 2*2
  )
  return(dataframe.bsrinf)
}

########## Identifying RL interactions in the T-enr subset ###########

# Subsetting T-enr AOIs
gpx_topT <- # 61 AOIs that were annotated as high T content
ppx_topT <- # 40 AOIs that were annotated as high T content

gpx_subset_T <- targetDataSubset_gpx[,targetDataSubset_gpx$ROI_REL %in% gpx_topT] # 104 
ppx_subset_T <- targetDataSubset_ppx[,targetDataSubset_ppx$ROI_REL %in% ppx_topT] # 73 

gpx_T_neg <- gpx_subset_T[,gpx_subset_T$segment == "PanCK-"]
ppx_T_neg <- ppx_subset_T[,ppx_subset_T$segment == "PanCK-"]

# Using the Q3 pre normalized data
gpx_counts <- assayDataElement(gpx_T_neg,elt = "q_norm")
colnames(gpx_counts) <- gpx_T_neg$ROI_REL
ppx_counts <- assayDataElement(ppx_T_neg,elt = "q_norm")
colnames(ppx_counts) <- ppx_T_neg$ROI_REL

# Preparing BulkSignalR objects
bs_gpx <- prepareDataset(counts = gpx_counts, min.count = 0, prop = 0, normalize = FALSE, method = 'UQ') 
bs_gpx <- learnParameters(bs_gpx, quick=FALSE, verbose=TRUE, # default min.positive = 4: min # target genes to be found in a given pathway
                          plot.folder = "BulkSignalR/", filename = "gpx_counts")
bs_ppx <- prepareDataset(counts = ppx_counts, min.count = 0, prop = 0, normalize = FALSE, method = 'UQ')
bs_ppx <- learnParameters(bs_ppx, quick=FALSE, verbose=TRUE,
                          plot.folder = "BulkSignalR/", filename = "ppx_counts")

# Generate inferences -- min cor = 0.25
gpx_inf <- initialInference(bs_gpx, min.cor = 0.25) 
LRinter.df_gpx <- LRinter(gpx_inf) 
ppx_inf <- initialInference(bs_ppx, min.cor = 0.25)
LRinter.df_ppx <- LRinter(ppx_inf)

# Save results
write_xlsx(list(GPx = LRinter.df_gpx,PPx = LRinter.df_ppx),"BulkSignalR/1_neg_initialInference_corr0.25.xlsx")

# Reduce to best pathway - GPx
gpx_inf.redBP <- reduceToBestPathway(gpx_inf)
gpx_sig.redBP <- getLRGeneSignatures(gpx_inf.redBP, qval.thres=0.05) 
scoresLR_gpx <- scoreLRGeneSignatures(bs_gpx, gpx_sig.redBP,
                                      name.by.pathway=FALSE)
# Reduce to best pathway - PPx
ppx_inf.redBP <- reduceToBestPathway(ppx_inf)
ppx_sig.redBP <- getLRGeneSignatures(ppx_inf.redBP, qval.thres=0.05)
scoresLR_ppx <- scoreLRGeneSignatures(bs_ppx, ppx_sig.redBP,
                                      name.by.pathway=FALSE)

# Save results and all data objects
write_xlsx(list(GPx = LRinter(gpx_inf.redBP),PPx = LRinter(ppx_inf.redBP)),"BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx")
save(file="BulkSignalR/BSR_Tenr_neg_corr0.25.RData",bs_gpx, bs_ppx, gpx_inf, ppx_inf, gpx_inf.redBP,ppx_inf.redBP,gpx_sig.redBP,ppx_sig.redBP,scoresLR_gpx,scoresLR_ppx)

# Assigning cell types to interactions by inputting a proportion score for each AOI
load("Tenr_spatial_proportions_all.RData")

# assignCellTypesToInteractions
# qval = 0.05
gpx_ct_0.05 <- assignCellTypesToInteractions(bs_gpx,gpx_inf,gpx_mat,qval.thres=0.05)
ppx_ct_0.05 <- assignCellTypesToInteractions(bs_ppx,ppx_inf,ppx_mat,qval.thres=0.05)
gpx_ct_0.01 <- assignCellTypesToInteractions(bs_gpx,gpx_inf,gpx_mat,qval.thres=0.01)
ppx_ct_0.01 <- assignCellTypesToInteractions(bs_ppx,ppx_inf,ppx_mat,qval.thres=0.01)
all_ct <- list(GPx_0.01 = data.frame(gpx_ct_0.01), PPx_0.01 = data.frame(ppx_ct_0.01), GPx_0.05 = data.frame(gpx_ct_0.05), PPx_0.05 = data.frame(ppx_ct_0.05))

write_xlsx(all_ct,"BulkSignalR/3_Tenr_cell_types_int_neg_corr0.25.xlsx")

############## Unique L-R in Gpx, PPx Tenr neg and plotting chord diagrams (Figure 3e) #########
load("BulkSignalR/BSR_Tenr_neg_corr0.25.RData")

# Loading in the L-R unique list from the spreadsheets, after doing the Venns FOR CORR 0.25 (no inhibitory interactions, and with q < 0.05)
gpx_unique_genes <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "gpx_unique_0.05")
ppx_unique_genes <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "ppx_unique_0.05")
common_genes <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "common")

gpx_unique <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "GPx") %>% filter(`L-R` %in% gpx_unique_genes$`L-R`)
ppx_unique <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "PPx") %>% filter(`L-R` %in% ppx_unique_genes$`L-R`)
common_gpx <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "GPx") %>% filter(`L-R` %in% common_genes$`L-R`)
common_ppx <- read_excel("BulkSignalR/2_neg_uniqueLR_corr0.25.xlsx",sheet = "PPx") %>% filter(`L-R` %in% common_genes$`L-R`)

# Incorporating cell types
gpx_cell_0.05 <- read_excel("BulkSignalR/3_Tenr_cell_types_int_neg_corr0.25.xlsx",sheet = "GPx_0.05") %>% mutate(`L-R`=paste0(L,"-",R))
ppx_cell_0.05 <- read_excel("BulkSignalR/3_Tenr_cell_types_int_neg_corr0.25.xlsx",sheet = "PPx_0.05") %>% mutate(`L-R`=paste0(L,"-",R))

# Taking the top weighted cell type if there are multiple hits for one L-R pair
gpx_cell_0.05 <- merge(aggregate(weight ~ `L-R`, gpx_cell_0.05 , max), gpx_cell_0.05 )
ppx_cell_0.05 <- merge(aggregate(weight ~ `L-R`, ppx_cell_0.05 , max), ppx_cell_0.05 )

# Chord diagrams -- corr = 0.25, q < 0.05 -- FOR THE PAPER
dataframe.bsrinf_gpx <- plot_chord(gpx_inf.redBP,gpx_unique_genes$`L-R`,gpx_cell_0.05,filename = "gpx_neg_Tenr_plot_50_corr0.25_new",num_genes = 50)
dataframe.bsrinf_ppx <- plot_chord(ppx_inf.redBP,ppx_unique_genes$`L-R`, ppx_cell_0.05,filename ="ppx_neg_Tenr_plot_50_corr0.25_new",num_genes =50)

############## Identifying RL interactions in Atypia ###############

# Subsetting Atypic AOIs
atypia_gpx <- targetDataSubset_gpx[,targetDataSubset_gpx$Epithelial == "A"] # 43
pData(atypia_gpx)$ROI_REL_seg <-  paste0(atypia_gpx$ROI_REL,"_",atypia_gpx$segment)
atypia_ppx <- targetDataSubset_ppx[,targetDataSubset_ppx$Epithelial == "A"] # 53
pData(atypia_ppx)$ROI_REL_seg <-  paste0(atypia_ppx$ROI_REL,"_",atypia_ppx$segment)

# Using the Q3 pre normalized data
gpx_counts <- assayDataElement(atypia_gpx,elt = "q_norm")
colnames(gpx_counts) <- atypia_gpx$ROI_REL_seg
ppx_counts <- assayDataElement(atypia_ppx,elt = "q_norm")
colnames(ppx_counts) <- atypia_ppx$ROI_REL_seg

# Preparing BulkSignalR objects
bs_gpx <- prepareDataset(counts = gpx_counts, min.count = 0, prop = 0, normalize = FALSE, method = 'UQ') 
bs_gpx <- learnParameters(bs_gpx, quick=FALSE, verbose=TRUE, # default min.positive is 4: min # target genes to be found in a given pathway
                          plot.folder = "BulkSignalR/", filename = "atypia_gpx_counts")
bs_ppx <- prepareDataset(counts = ppx_counts, min.count = 0, prop = 0, normalize = FALSE, method = 'UQ')
bs_ppx <- learnParameters(bs_ppx, quick=FALSE, verbose=TRUE,
                          plot.folder = "BulkSignalR/", filename = "atypia_ppx_counts")

# Generate inferences with min cor = 0.25
gpx_inf <- initialInference(bs_gpx, min.cor = 0.25) 
LRinter.df_gpx <- LRinter(gpx_inf) 
LRinter.df_gpx$Genes <- unlist(lapply(tGenes(gpx_inf), function(x) paste0(x,collapse = ",")))

ppx_inf <- initialInference(bs_ppx, min.cor = 0.25)
LRinter.df_ppx <- LRinter(ppx_inf)
LRinter.df_ppx$Genes <- unlist(lapply(tGenes(ppx_inf), function(x) paste0(x,collapse = ",")))

# Save results
write_xlsx(list(GPx = LRinter.df_gpx,PPx = LRinter.df_ppx),"BulkSignalR/1_atypia_initialInference_corr0.25.xlsx")

# Reduce to best pathway - GPx
gpx_inf.redBP <- reduceToBestPathway(gpx_inf) # one of each LR interaction
gpx_sig.redBP <- getLRGeneSignatures(gpx_inf.redBP, qval.thres=0.01)
scoresLR_gpx <- scoreLRGeneSignatures(bs_gpx, gpx_sig.redBP,
                                      name.by.pathway=FALSE)
# Reduce to best pathway - PPx
ppx_inf.redBP <- reduceToBestPathway(ppx_inf)
ppx_sig.redBP <- getLRGeneSignatures(ppx_inf.redBP, qval.thres=0.01)
scoresLR_ppx <- scoreLRGeneSignatures(bs_ppx, ppx_sig.redBP,
                                      name.by.pathway=FALSE)

# Saving results and data objects
write_xlsx(list(GPx = LRinter(gpx_inf.redBP),PPx = LRinter(ppx_inf.redBP)),"BulkSignalR/2_atypia_uniqueLR_corr0.25.xlsx")
save(file="BulkSignalR/BSR_atypia_corr0.25.RData",bs_gpx, bs_ppx, gpx_inf, ppx_inf,gpx_inf.redBP, ppx_inf.redBP)

# Further, can assign cell types to interactions by inputting a proportion score for each AOI (from SpatialDecon)
load("4_SpatialDecon_major_GPx.RData") 
atypia_gpx_neg <- atypia_gpx[,atypia_gpx$segment=="PanCK-"]
gpx_mat <- reordered_matrix0[,atypia_gpx_neg$ROI_REL]

load("4_SpatialDecon_major_PPx.RData") 
atypia_ppx_neg <- atypia_ppx[,atypia_ppx$segment=="PanCK-"]
ppx_mat <- reordered_matrix0[,atypia_ppx_neg$ROI_REL]

atypia_spatialprop <- cbind(gpx_mat,ppx_mat)

# Adding epithelial cell proportions to gpx_mat, ppx_mat
gpx_mat <- rbind(gpx_mat, rep(0,ncol(gpx_mat)))
rownames(gpx_mat)[8] <- "Epithelial" 
colnames(gpx_mat) <- atypia_gpx_neg$ROI_REL_seg

ppx_mat <- rbind(ppx_mat, rep(0,ncol(ppx_mat)))
rownames(ppx_mat)[8] <- "Epithelial"
colnames(ppx_mat) <- atypia_ppx_neg$ROI_REL_seg

cols <- atypia_gpx[,atypia_gpx$segment == "PanCK+"]$ROI_REL_seg
gpx_mat_pos <- matrix(0,7,length(cols))
gpx_mat_pos <- rbind(gpx_mat_pos, rep(1,ncol(gpx_mat_pos)))
rownames(gpx_mat_pos) <- rownames(gpx_mat)
colnames(gpx_mat_pos) <- cols

cols <- atypia_ppx[,atypia_ppx$segment == "PanCK+"]$ROI_REL_seg
ppx_mat_pos <- matrix(0,7,length(cols))
ppx_mat_pos <- rbind(ppx_mat_pos, rep(1,ncol(ppx_mat_pos)))
rownames(ppx_mat_pos) <- rownames(ppx_mat)
colnames(ppx_mat_pos) <- cols

gpx_mat_all <- cbind(gpx_mat,gpx_mat_pos)
ppx_mat_all <- cbind(ppx_mat,ppx_mat_pos)

setequal(colnames(gpx_mat_all),colnames(gpx_counts)) # TRUE
setequal(colnames(ppx_mat_all),colnames(ppx_counts)) # TRUE

# assignCellTypesToInteractions
gpx_ct_0.05 <- assignCellTypesToInteractions(bs_gpx,gpx_inf,gpx_mat_all,qval.thres=0.05)
ppx_ct_0.05 <- assignCellTypesToInteractions(bs_ppx,ppx_inf,ppx_mat_all,qval.thres=0.05)
gpx_ct_0.01 <- assignCellTypesToInteractions(bs_gpx,gpx_inf,gpx_mat_all,qval.thres=0.01)
ppx_ct_0.01 <- assignCellTypesToInteractions(bs_ppx,ppx_inf,ppx_mat_all,qval.thres=0.01)
all_ct <- list(GPx_0.01 = data.frame(gpx_ct_0.01), PPx_0.01 = data.frame(ppx_ct_0.01), GPx_0.05 = data.frame(gpx_ct_0.05), PPx_0.05 = data.frame(ppx_ct_0.05))

write_xlsx(all_ct,"BulkSignalR/3_atypia_cell_types_int_all.xlsx")

##############  Unique L-R in GPx, PPx Atypia and plotting chord diagrams (Figure 5f) #########
# Loading in the L-R unique list from the spreadsheets, after doing a Venn diagram of gene intersections between Atypia-GPx and Atypia-PPx
load("BulkSignalR/BSR_atypia_corr0.25.RData")

gpx_unique_genes <- read_excel("BulkSignalR/2_atypia_uniqueLR_corr0.25_final.xlsx",sheet = "gpx_unique")
ppx_unique_genes <- read_excel("BulkSignalR/2_atypia_uniqueLR_corr0.25_final.xlsx",sheet = "ppx_unique")

gpx_unique <- read_excel("BulkSignalR/2_atypia_uniqueLR_corr0.25_final.xlsx",sheet = "GPx") %>% filter(`L-R` %in% gpx_unique_genes$`L-R`) # 115
ppx_unique <- read_excel("BulkSignalR/2_atypia_uniqueLR_corr0.25_final.xlsx",sheet = "PPx") %>% filter(`L-R` %in% ppx_unique_genes$`L-R`) # 140

gpx_L <- gpx_unique$L
gpx_R <- gpx_unique$R

ppx_L <- ppx_unique$L
ppx_R <- ppx_unique$R

# Reading in cell type information
gpx_cell_0.05 <- read_excel("BulkSignalR/3_atypia_cell_types_int_all_final.xlsx",sheet = "GPx_0.05") %>% mutate(`L-R`=paste0(L,"-",R))
ppx_cell_0.05 <- read_excel("BulkSignalR/3_atypia_cell_types_int_all_final.xlsx",sheet = "PPx_0.05") %>% mutate(`L-R`=paste0(L,"-",R))

# Taking the top weighted cell type if there are multiple hits for one L-R pair
gpx_cell_0.05 <- merge(aggregate(weight ~ `L-R`, gpx_cell_0.05 , max), gpx_cell_0.05 )
ppx_cell_0.05 <- merge(aggregate(weight ~ `L-R`, ppx_cell_0.05 , max), ppx_cell_0.05 )

# Plotting chord diagrams -- corr = 0.25, q < 0.05 -- FOR THE PAPER
dataframe.bsrinf_gpx <- plot_chord(BSR_inf = gpx_inf.redBP,L_R=gpx_unique_genes$`L-R`,gpx_cell_0.05,filename ="gpx_atypia_all_plot_25",num_genes=25)
plot_chord(BSR_inf = gpx_inf.redBP,L_R=gpx_unique_genes$`L-R`,gpx_cell_0.05,filename ="gpx_atypia_all_plot_50",num_genes=50)
dataframe.bsrinf_ppx <- plot_chord(BSR_inf = ppx_inf.redBP,L_R=ppx_unique_genes$`L-R`,ppx_cell_0.05,filename ="ppx_atypia_all_plot_25",num_genes=25)
plot_chord(BSR_inf = ppx_inf.redBP,L_R=ppx_unique_genes$`L-R`,ppx_cell_0.05,filename ="ppx_atypia_all_plot_50",num_genes=50)

write_xlsx(list(GPx = dataframe.bsrinf_gpx,PPx = dataframe.bsrinf_ppx),"BulkSignalR/4_atypia_all_genes.xlsx")

