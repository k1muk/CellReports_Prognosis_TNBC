library(NanoStringNCTools)
library(GeomxTools)
library(dplyr)
library(knitr)
library(ggforce)
library(writexl)

load("/BreastCancer_ImmuneProfile_majorCellTypes.RData") #In Github
#This has the major annotations 

my_color_palette <- c("palegreen", "goldenrod", "lightgoldenrod1", 
                      "blue4", 
                      "mediumpurple1",
                      "dodgerblue2", "lightskyblue1" )

# Define the desired row name order
desired_order <- c("Endothelial","Fibroblast/Mesenchymal", "PVL", 
                   "T-cells",
                   "Myeloid",
                   "B-cells","Plasmablasts")


setwd(analysis)
load(file="5_DataSubsets_fordownstream.RData")


########## Running Spatial Decon for GPX###################
res1 = runspatialdecon(object= GPx_subset_neg,
                       norm_elt = "q_norm",
                       raw_elt = "exprs",
                       X=  tnbc.profile_matrix.major,
                       is_pure_tumor = NULL,
                       align_genes = TRUE,
                       cell_counts= GPx_subset_neg$nuclei)

#the output for runspatial decon will be like this
#https://rdrr.io/github/Nanostring-Biostats/SpatialDecon/man/runspatialdecon.html
res00= res1

all.equal(rownames(pData(GPx_subset_neg)), rownames(pData(res00)$cell.counts$cell.counts))
rownames(pData(res00)$cell.counts$cell.counts)<- GPx_subset_neg$ROI_REL
all.equal(rownames(pData(GPx_subset_neg)), rownames(res00$prop_of_all))
rownames(res00$prop_of_all)<- GPx_subset_neg$ROI_REL


dim(res00$beta)
dim(pData(res00)$cell.counts$cell.counts)
dim(pData(GPx_subset_neg))

dim(res00$prop_of_all)

#check if the order is identical, 
newcell <- t(pData(res00)$cell.counts$cell.counts)
newprop <- t(res00$prop_of_all)


# Reorder the rows based on the desired order
reordered_matrix0 <- newprop[match(desired_order,row.names(newprop)),] #proportions of counts
reordered_matrix <- newcell[match(desired_order,row.names(newcell)),] 

t3= as.data.frame(GPx_subset_neg$nuclei, rownames(pData(GPx_subset_neg)))
t4 = as.data.frame(rowSums(pData(res00)$cell.counts$cell.counts)) #cell.counts: beta rescaled to estimate cell numbers, based on prop_of_all and nuclei count
t5 = reshape2::melt(cbind(t3, t4) %>% rownames_to_column("AOI"))

#% of T cells alone in each of the samples 
t6.gpx <- data.frame(AOI = colnames(reordered_matrix0),
                     Tcell_pct = reordered_matrix0["T-cells", ]) %>%
  mutate(sampleName = sub("\\..*", "", AOI),
         region = sub("^[^.]*\\.[^.]*\\.([A-Z]).*", "\\1", AOI),
         epithelial = sub("^[^.]*\\.[^.]*\\.[^.]*\\.([A-Z]).*", "\\1", AOI))


#% of B cells and plasmablasts in the sample
b6.gpx= reshape2::melt(reordered_matrix0["B-cells",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 
pb6.gpx= reshape2::melt(reordered_matrix0["Plasmablasts",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 

#% of Myeloid cells
myl6.gpx = reshape2::melt(reordered_matrix0["Myeloid",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 

#% endo cells
endo6.gpx = reshape2::melt(reordered_matrix0["Endothelial",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 

#% of fibroblast and PVL cells combined?
mes6.gpx =  reshape2::melt(reordered_matrix0[c("Fibroblast/Mesenchymal", "PVL") ,])  %>% mutate(sampleName = sub("\\..*", "", Var2)) 

newmes6.gpx <- mes6.gpx %>%
  group_by(Var2) %>%
  summarise(value = sum(value), sampleName = first(sampleName)) %>%
  ungroup() #merging the fibroblast and PVL to the same class called mesenchymal


gpxtQuantiles = quantile(t6.gpx$Tcell_pct)
q= quantile(t6.gpx$Tcell_pct, 0.25) 
length(which(t6.gpx$Tcell_pct>q)) 
dim(t6.gpx) #234

#B cell quantiles
gpxbQuantiles= quantile(b6.gpx$value)

#PB cell quantiles
gpxpbQuantiles= quantile(pb6.gpx$value)

#Myeloid cell quantiles
gpxmylQuantiles = quantile(myl6.gpx$value)

#mes values
gpxmesQuantiles = quantile(newmes6.gpx$value)

#endo values
gpxendoQuantiles = quantile(endo6.gpx$value)


#Identifying a T-enriched subset
q1= quantile(t6.gpx$Tcell_pct, 0.75) 
t6.gpx[which(t6.gpx$Tcell_pct>q1),] %>% group_by(sampleName) %>%summarize(count = n())
#Getting only the AOIs that have a good T as seen in the above table, from only the * patients
gpx_tenr <-  t6.gpx[which(t6.gpx$Tcell_pct>q1),] %>% group_by(sampleName) %>% 
  filter(sampleName %in% c("xx"))

#Extract only these AOIs from the original gpx_subset_neg object
subset_gpx_tenr <- subset(GPx_subset_neg, CodeClass=="Endogenous" | CodeClass == "Negative", 
                          pData(GPx_subset_neg)$ROI_REL %in% gpx_tenr$AOI )  


#sanity check
all.equal(pData(subset_gpx_tenr)$ROI_REL, gpx_tenr$AOI) #TRUE

#######
save(file="4_SpatialDecon_major_GPx_updated.RData",t6.gpx, t4, t5, t3,  res1, reordered_matrix0,reordered_matrix)
########

#######
save(file="4_GPx_Quantiles_updated.RData", gpx_tenr, subset_gpx_tenr, 
     gpxbQuantiles, gpxpbQuantiles, gpxtQuantiles, gpxmylQuantiles, gpxmesQuantiles,gpxendoQuantiles,
     t6.gpx, b6.gpx, pb6.gpx, myl6.gpx, newmes6.gpx, mes6.gpx, endo6.gpx, reordered_matrix0) #Used to generate Fig S3c

#####


################ Doing the deconvolution for the PPX############################
res2 = runspatialdecon(object= PPx_subset_neg,
                       norm_elt = "q_norm",
                       raw_elt = "exprs",
                       X=  tnbc.profile_matrix.major,
                       is_pure_tumor = NULL,
                       align_genes = TRUE,
                       cell_counts= PPx_subset_neg$nuclei)
res02= res2
all.equal(rownames(pData(PPx_subset_neg)), rownames(pData(res02)$cell.counts$cell.counts))
rownames(pData(res02)$cell.counts$cell.counts)<- PPx_subset_neg$ROI_REL
all.equal(rownames(pData(PPx_subset_neg)), rownames(res02$prop_of_all))
rownames(res02$prop_of_all)<- PPx_subset_neg$ROI_REL
dim(res02$beta)
dim(pData(PPx_subset_neg))
dim(res02$prop_of_all)
#check if the order is identical
newprop <- t(res02$prop_of_all)
newcell <- t(pData(res02)$cell.counts$cell.counts)


# Define the desired row name order
# Reorder the rows based on the desired order
reordered_matrix <- newcell[match(desired_order,row.names(newcell)),]
reordered_matrix0 <- newprop[match(desired_order,row.names(newprop)),]
t3= as.data.frame(PPx_subset_neg$nuclei, rownames(pData(PPx_subset_neg)))
t4 = as.data.frame(rowSums(pData(res02)$cell.counts$cell.counts))
t5 = melt(cbind(t3, t4) %>% rownames_to_column("AOI"))

#% proportion of T cells alone in each of the samples 
t6.ppx <- data.frame(AOI = colnames(reordered_matrix0),
                     Tcell_pct = reordered_matrix0["T-cells", ]) %>%
  mutate(sampleName = sub("\\..*", "", AOI),
         region = sub("^[^.]*\\.[^.]*\\.([A-Z]).*", "\\1", AOI),
         epithelial = sub("^[^.]*\\.[^.]*\\.[^.]*\\.([A-Z]).*", "\\1", AOI))


#% proportion of B cells and plasmablasts in the sample
b6= reshape2::melt(reordered_matrix0["B-cells",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 
pb6= reshape2::melt(reordered_matrix0["Plasmablasts",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 


#% of Myeloid cells
myl6 = reshape2::melt(reordered_matrix0["Myeloid",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 

#% endo cells
endo6 = reshape2::melt(reordered_matrix0["Endothelial",])  %>% rownames_to_column("AOI")  %>% mutate(sampleName = sub("\\..*", "", AOI)) 

#% of fibroblast and PVL cells combined?
mes6 =  reshape2::melt(reordered_matrix0[c("Fibroblast/Mesenchymal", "PVL") ,])  %>% mutate(sampleName = sub("\\..*", "", Var2)) 

newmes6 <- mes6 %>%
  group_by(Var2) %>%
  summarise(value = sum(value), sampleName = first(sampleName)) %>%
  ungroup() #merging the fibroblast and PVL to the same class called mesenchymal


#--------quantile calculations---------#

ppxtQuantiles = quantile(t6.ppx$Tcell_pct)
q= quantile(t6.ppx$Tcell_pct, 0.25) 

length(which(t6.ppx$Tcell_pct>q)) #192

#B cell quantiles
ppxbQuantiles= quantile(b6$value)

#PB cell quantiles
ppxpbQuantiles= quantile(pb6$value)

#Myeloid cell quantiles
ppxmylQuantiles = quantile(myl6$value)

#mes values
ppxmesQuantiles = quantile(newmes6$value)

#endo values
ppxendoQuantiles = quantile(endo6$value)


#Finding the Tenr subset
q1= quantile(t6.ppx$Tcell_pct, 0.75) 
t6.ppx[which(t6.ppx$Tcell_pct>q1),] %>% group_by(sampleName) %>% summarize(count = n())

#Getting only the AOIs that have a good T as seen in the above table, from only the * patients, and atleast >7-8 ROIS from each patient
ppx_tenr <-  t6.ppx[which(t6.ppx$Tcell_pct>q1),] %>% group_by(sampleName) %>% 
  filter(sampleName %in% c("xx")) #only T high patients from PPx

#Extract only these AOIs from the original PPx_subset_neg object
subset_ppx_tenr <- subset(PPx_subset_neg, CodeClass=="Endogenous" | CodeClass == "Negative", 
                      pData(PPx_subset_neg)$ROI_REL %in% ppx_tenr$AOI )  


#sanity check
all.equal(pData(subset_ppx_tenr)$ROI_REL, ppx_tenr$AOI)

######
save(file="4_SpatialDecon_major_PPx_updatd.RData",t3, t4, t5, t6.ppx, res02, reordered_matrix0,reordered_matrix)
#####


#######
save(file="4_PPx_Quantiles_updated.RData", ppx_tenr, subset_ppx_tenr, 
     ppxbQuantiles, ppxpbQuantiles, ppxtQuantiles, ppxmylQuantiles, ppxmesQuantiles,ppxendoQuantiles,
     t6.ppx, b6, pb6, myl6, newmes6, mes6, endo6, reordered_matrix0) #Used to generate Fig S3c
#####

