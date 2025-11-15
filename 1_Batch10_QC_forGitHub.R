library(NanoStringNCTools)
library(GeomxTools)
library(dplyr)
library(knitr)
library(ggforce)
library(writexl)
library(purrr)
library(tidyverse)

setwd("YourPath")
options(java.parameters = "-Xmx12000m")

jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
} 

#### IMPORTANT: ####
#This code can be adapted and utilized to process the data uploaded to GEO.

#The annotation file was adjusted as per the following notes and placed with the DCC folder allCombined.
# Make sure that the order of the DCC files in SampletypeFile (metadata excel) is correct -- No Template Control is the first DCC for each plate 
# (aka Sample_ID is sorted alphabetically)
# Also do not modify the headers -- "sample name" NOT "sample_name" or "sample.name"

# --- Parameters changed from recommended ---
# segment QC: minNegativeCount = 1, maxNTCCount = 100000, minNuclei = 10, minArea = 1000 (lost 24 segments) 
#                   -- 23 are panck+, brief check seems these are ROI in immune regions with very little pancytokeratin staining 
# probe QC: removeLocalOutliers = F (lost 0 segments)


#### LOADING DATA ####----------------------------------------------------------------------------------------------
data.dir <- file.path("YourPath") 
DCCFiles <- dir(file.path(data.dir), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
PKCFiles <- "YourPath to PKC"
SampletypeFile <-dir(file.path(data.dir), pattern = ".xlsx$",full.names = TRUE, recursive = TRUE)


ProgData_data <- readNanoStringGeoMxSet(dccFiles = DCCFiles, pkcFiles = PKCFiles,phenoDataFile = SampletypeFile,
                                       phenoDataDccColName = "Sample_ID", 
                                       protocolDataColNames = c("roi", "aoi"), 
                                       experimentDataColNames =c("panel"),
                                       phenoDataSheet= "Sheet1")


# shift counts by one
ProgData_data_shifted <- shiftCountsOne(ProgData_data, useDALogic=TRUE)
# 18815 features, 1423 samples

# checking modules used 
library(knitr)
pkcs <- annotation(ProgData_data_shifted)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))


#### PRE QC PLOTS ####----------------------------------------------------------------------------------------------
# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(ProgData_data_shifted)$NTC),col.names = c("NTC Count", "# of Segments"))
'''
|NTC Count | # of Segments|
|:---------|-------------:|
|0         |            87|
|52        |            95|
|56        |            95|
|99        |            95|
|228       |            95|
|480       |            95|
|486       |            81|
|1134      |            95|
|1252      |            95|
|1319      |            95|
|1619      |            90|
|3997      |            64|
|15878     |            70|
|17479     |            85|
|41205     |            92|
|54776     |            95|
'''

f1 = QC_histogram(sData(ProgData_data_shifted), "NTC", col_by, 1000) +
  QC_histogram(sData(ProgData_data_shifted), "area", col_by, 5000, scale_trans = "log10") + 
  QC_histogram(sData(ProgData_data_shifted), "nuclei", col_by, 100) 

f2 = QC_histogram(sData(ProgData_data_shifted), "Raw", col_by, 1000) + 
  QC_histogram(sData(ProgData_data_shifted), "Trimmed", col_by, 80) +
  QC_histogram(sData(ProgData_data_shifted), "Stitched", col_by, 80)+
  QC_histogram(sData(ProgData_data_shifted), "Aligned", col_by, 80)

#### SEGMENT QC ####----------------------------------------------------------------------------------------------
paste("## of Negative Probes:", sum(fData(ProgData_data_shifted)$Negative))
#139

min(exprs(ProgData_data_shifted[fData(ProgData_data_shifted)$Negative,]))
max(exprs(ProgData_data_shifted[fData(ProgData_data_shifted)$Negative,])) #1615, not an unreasonable negative count, is probe RTS0039347

#hist(exprs(ProgData_data_shifted[fData(ProgData_data_shifted)$Negative,])[100,])

metData <- sData(ProgData_data_shifted)
 f3 = ggplot(metData, aes(y=log10(metData$area), x=log10(metData$nuclei))) +
  geom_point()+ labs(y = "log10(Area)",x = "log10(Nuclei)") + geom_smooth(method = "lm", formula = y ~ x, geom = "smooth")+
  ggplot(metData, aes(x=metData$Raw, y=metData$Aligned, color=metData$`slide name`, shape=metData$PP)) +
  geom_point()+ labs(x = "Raw Reads",y = "AlignedReads")+stat_smooth(method = "lm", formula = y ~ x, geom = "smooth")+
  ggplot(metData, aes(x=metData$DeduplicatedReads, y=metData$Aligned, color=metData$`slide name`, shape=metData$Prognosis)) +
  geom_point()+ labs(x = "Deduplicated Reads",y = "AlignedReads")


 # checking that slides are labelled properly, also tells us how many segments in each
 table(Slide = ProgData_data_shifted$`slide name`,Prognosis = ProgData_data_shifted$Prognosis)
 

QC_params <- # adjusted params for this batch
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 80,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%).
       minNegativeCount = 1,   # Minimum negative control counts (10).
       maxNTCCount = 60000,     # Maximum counts observed in NTC well (1000).
       minNuclei = 40,         # Minimum # of nuclei estimated (100).
       minArea = 1000)         # Minimum segment area (5000).

ProgData <- setSegmentQCFlags(ProgData_data_shifted, qcCutoffs = QC_params)

# Collate QC Results
QCResults <- protocolData(ProgData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")

'''
Table: QC Summary Table for each Segment

|              | Pass| Warning|
  |:-------------|----:|-------:|
  |LowReads      | 1423|       1|
  |LowTrimmed    | 1424|       0|
  |LowStitched   | 1424|       0|
  |LowAligned    | 1424|       0|
  |LowSaturation | 1421|       3|
  |LowNegatives  | 1424|       0|
  |HighNTC       | 1424|       0|
  |LowNuclei     | 1382|      42|
  |LowArea       | 1405|      19|
  |TOTAL FLAGS   | 1380|      45|
  
  '''
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(ProgData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(ProgData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(ProgData)[, negCols] <- sData(ProgData)[["NegGeoMean"]]
for(ann in negCols) { # negCols is 1 
  plt <- QC_histogram(pData(ProgData), ann, col_by, 2, scale_trans = "log10")
}

pdf(file="ProgData_preQC_plots.pdf",width = 15, height=8)
f1
f2
f3
plt
dev.off()

# detatch neg_geomean columns ahead of aggregateCounts call
pData(ProgData) <- pData(ProgData)[, !colnames(pData(ProgData)) %in% negCols]

# remove flagged QC segments
ProgData_pass <- ProgData[, QCResults$QCStatus == "PASS"]
ProgData_fail <- ProgData[, QCResults$QCStatus == "WARNING"]
# ProgData_failed_info <- data.frame(subset(sData(ProgData_fail), select = -c(FileVersion,SoftwareVersion)))

protocolData(ProgData_fail)$ROI_REL <- paste0(sData(ProgData_fail)$`slide name`,".",sData(ProgData_fail)$roi,".", sData(ProgData_fail)$Region,".",sData(ProgData_fail)$Epithelial, ".", sData(ProgData_fail)$ImmuneCellLocalization)
qc_info <- sData(ProgData_fail)$QCFlags
rownames(qc_info)<-sData(ProgData_fail)$ROI_REL

xlsx::write.xlsx(qc_info,file="ProgData_SegmentQC_Fail_new.xlsx")

metData1 <- sData(ProgData_pass)
f4 = ggplot(metData1, aes(y=log10(metData1$area), x=log10(metData1$nuclei), color= metData1$`slide name`, shape=metData1$Prognosis)) +
  geom_point()+ labs(y = "log10(Area)",x = "log10(Nuclei)")+
  ggplot(metData1, aes(x=metData1$Raw, y=metData1$Aligned, color=metData1$`slide name`, shape=metData1$PP)) +
  geom_point()+ labs(x = "Raw Reads",y = "AlignedReads")+
  ggplot(metData1, aes(x=metData1$DeduplicatedReads, y=metData1$Aligned, color=metData1$`slide name`, shape=metData1$Prognosis)) +
  geom_point()+ labs(x = "Deduplicated Reads",y = "AlignedReads")+
  ggplot(metData1, aes(x=as.numeric(unlist(metData1$`Saturated (%)`)), y= metData1$Raw, color=metData1$`slide name`, shape=metData1$Prognosis)) +
  geom_point()+ labs(x = "Sequencing Saturation %",y = "#Raw Reads (log10)")

pdf(file="ProgData_postQC_plots.pdf",width = 15, height=8)
f4
dev.off()

metData2 <- sData(ProgData)
ggplot(metData2, aes(x=as.numeric(unlist(metData2$`Saturated (%)`)), y= metData2$Raw, color=metData2$`slide name`, shape=metData2$Prognosis)) +
  geom_point()+ labs(x = "Sequencing Saturation %",y = "#Raw Reads (log10)")

dim(ProgData)
# 18815 features, 1424 samples
dim(ProgData_pass)
# 18815 features, 1379 samples

#### PROBE QC ####----------------------------------------------------------------------------------------------
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
ProgData_pass <- setBioProbeQCFlags(ProgData_pass, 
                                 qcCutoffs = list(minProbeRatio = 0.1,
                                                  percentFailGrubbs = 20), 
                                 removeLocalOutliers = F)

ProbeQCResults <- fData(ProgData_pass)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
kable(qc_df, caption = "Probes flagged or passed as outliers")
'''
| Passed| Global| Local|
|------:|------:|-----:|
|  18738|      0|    77|
'''
# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(ProgData_pass, 
         fData(ProgData_pass)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(ProgData_pass)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
# Features  Samples 
# 18815     1379
ProgData_passed_QC <- ProbeQCPassed 

# checking how many 0 probes there are. 
probeidx <- rowSums(exprs (ProgData_passed_QC))
all_lowprobeidx <- which(probeidx==0)  # no 0 probes

#### GENE-LEVEL COUNT DATA ####----------------------------------------------------------------------------------------------
# Check how many unique targets the object has
length(unique(featureData(ProgData_passed_QC)[["TargetName"]]))
# 18677

# collapse to targets
target_ProgData <- aggregateCounts(ProgData_passed_QC)
dim(target_ProgData)
# Features  Samples 
# 18677     1379 

#### LOQ ####-----------------------------------------------------------------------------------------------------------------
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_ProgData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_ProgData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_ProgData)[, vars[1]] * 
             pData(target_ProgData)[, vars[2]] ^ cutoff)
  }
}
pData(target_ProgData)$LOQ <- LOQ

# Filtering based on LOQ
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_ProgData)$Module == module
  Mat_i <- t(esApply(target_ProgData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_ProgData)$TargetName, ]

# Segment gene detection rate: filter out segments with very low signal

# Save detection rate information to pheno data
pData(target_ProgData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_ProgData)$GeneDetectionRate <-
  pData(target_ProgData)$GenesDetected / nrow(target_ProgData)

metData3= pData(target_ProgData)
f5 = ggplot(metData3, aes(x=as.numeric(unlist(metData3$LOQ)), y=metData3$GenesDetected, color= metData3$`slide name`, shape=metData3$segment)) +
  geom_point()+ labs(x = "LOQ",y = "genesDetected") 

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_ProgData)$DetectionThreshold <- 
  cut(pData(target_ProgData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-2%", "2-3%", "3-4%", "4-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
f6 = ggplot(pData(target_ProgData),
            aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = `slide name`)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Sample Name")


f7 = ggplot(pData(target_ProgData),
            aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

f8= ggplot(pData(target_ProgData),
           aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Epithelial)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Epithelial")

f6+f7+f8

table(pData(target_ProgData)$DetectionThreshold,
            pData(target_ProgData)$`slide name`)

#For the broad questions of GpX vs PPx 
table(pData(target_ProgData)$DetectionThreshold, pData(target_ProgData)$Prognosis)

'        GPx PPx
<1%      0   0
1-2%    25  14
2-3%    26  37
3-4%    28  35
4-5%    26  37
5-10%  125 149
10-15%  86  86
>15%   394 311
'

############ SEGMENT FILTERING For BROAD QUESTIONS OF PPX vs GPX what is an ideal threshold for Segment filtering
# the distribution of a reasonable number of segments >5% in the PPx for pseudobulk processsing, it would be okay to filter at 5%
target_ProgData.sf <-
  target_ProgData[, pData(target_ProgData)$GeneDetectionRate >= .05] 

dim(target_ProgData.sf)
# Features  Samples
# 18677     1151

# gene detection rate
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_ProgData.sf)]
fData(target_ProgData.sf)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_ProgData.sf)$DetectionRate <-
  fData(target_ProgData.sf)$DetectedSegments / nrow(pData(target_ProgData.sf))


#### GENE Filtering ####---------------------------------------------------------------------------------------
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_ProgData.sf)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_ProgData.sf))
rownames(plot_detect) <- plot_detect$Freq

f8.1 = ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

negativeProbefData <- subset(fData(target_ProgData.sf), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# Retaining all genes seen in at least 10% of the segments
target_ProgData_GF_10 <- target_ProgData.sf[fData(target_ProgData.sf)$DetectionRate >= 0.1 |
                                         fData(target_ProgData.sf)$TargetName %in% neg_probes, ]

# Retaining all genes seen in at least 5% of the segments
target_ProgData_GF_5 <- target_ProgData.sf[fData(target_ProgData.sf)$DetectionRate >= 0.05 |
                                         fData(target_ProgData.sf)$TargetName %in% neg_probes, ]
dim(target_ProgData_GF_10)
# Features  Samples 
#8365 1151

dim(target_ProgData_GF_5)
#Features  Samples 
#10150    1151

table(pData(target_ProgData_GF_10)$DetectionThreshold,
            pData(target_ProgData_GF_10)$segment)

'
       PanCK- PanCK+
  <1%         0      0
  1-2%        0      0
  2-3%        0      0
  3-4%        0      0
  4-5%        0      0
  5-10%     133    141
  10-15%     83     89
  >15%      369    336
  
  
'
#For Pre-normalization UMAPs, Convert q_norm to log2 values (called logcounts) and save in 
assayDataElement(object = target_ProgData_GF_5, elt = "log_exprs") <- assayDataApply(target_ProgData_GF_5, 2, FUN = log, base = 2, elt = "exprs")

#### NORMALIZATION ####----------------------------------------------------------------------------------------------

ProgData_downstream <- GeomxTools::normalize(target_ProgData_GF_5 ,norm_method = "quant", desiredQuantile = .75,toElt = "q_norm")

boxplot(exprs(ProgData_downstream[,1:10]),col = "#9EDAE5", main = "Raw Counts",
             log = "y", names = 1:10, xlab = "Segment",ylab = "Counts, Raw")

boxplot(assayDataElement(ProgData_downstream[,1:10], elt = "q_norm"),
             col = "#2CA02C", main = "Q3 Norm Counts",
             log = "y", names = 1:10, xlab = "Segment",
             ylab = "Counts, Q3 Normalized")

#Convert q_norm to log2 values (called logcounts) for all downstream analysis
assayDataElement(object = ProgData_downstream, elt = "log_q") <- assayDataApply(ProgData_downstream, 2, FUN = log, base = 2, elt = "q_norm")


#### -----------------------------------------------------------------------
save(file="ProgData_downstream_GF5.RData",ProgData_downstream)
#####--------------------------------------------------------------------------

se <- SummarizedExperiment(
  assays = list(counts = exprs(ProgData_downstream), 
                log_counts= assayDataElement(ProgData_downstream, elt="log_q")), #log of the quantile normalized data
  rowData = fData(ProgData_downstream),
  colData = pData(ProgData_downstream)
)


base_palette <- "Paired"
my_colors <- colorRampPalette(brewer.pal(12, base_palette))(32)
groups <- unique(colData(se)$`slide name`)
my_colors_named <- setNames(my_colors, groups)

standR::drawPCA(se, assay= 2, col= colData(se)$`slide name`)+ theme(legend.position = "none")+
  scale_color_manual(values = my_colors_named) +
standR::drawPCA(se, assay= 2, col= colData(se)$Epithelial)+
standR::drawPCA(se, assay= 2, col= colData(se)$Region)+
standR::drawPCA(se, assay= 2, col= colData(se)$segment)


#####---------------------------------------------------------------------------
save(file="Consolidated_AllPlotting.RData", ProgData_data, ProgData_downstream, targetDataSubset_gpx, targetDataSubset_neg,
     targetDataSubset_pos, targetDataSubset_ppx, GPx_subset_pos, GPx_subset_neg, PPx_subset_pos, PPx_subset_neg, )
#####_-------------------------------------------------------------------------


#For ssGSEA slightly reformatting the data
{ #For ssGSEA 
#10150 features, 1151 samples 


protocolData(ProgData_downstream)[["ROI_L"]] <- paste0(sData(ProgData_downstream)$`slide name`,".",sData(ProgData_downstream)$roi,".", sData(ProgData_downstream)$ImmuneCellLocalization)
pData(ProgData_downstream)$ROI_L <- sData(ProgData_downstream)$ROI_L


protocolData(ProgData_downstream)[["R_E"]] <- paste0(sData(ProgData_downstream)$Region,"_", sData(ProgData_downstream)$Epithelial)
pData(ProgData_downstream)$R_E <- sData(ProgData_downstream)$R_E


protocolData(ProgData_downstream)[["ROI_REL"]] <- paste0(sData(ProgData_downstream)$`slide name`,".",sData(ProgData_downstream)$roi,".", sData(ProgData_downstream)$Region,".",sData(ProgData_downstream)$Epithelial, ".", sData(ProgData_downstream)$ImmuneCellLocalization)
pData(ProgData_downstream)$ROI_REL <- sData(ProgData_downstream)$ROI_REL



targetDataSubset_pos= subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK+")
targetDataSubset_neg= subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK-")


GPx_subset_neg = subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK-" & Prognosis %in% "GPx")
PPx_subset_neg = subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK-" & Prognosis %in% "PPx")

targetDataSubset_gpx= subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", Prognosis=="GPx")
targetDataSubset_ppx= subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", Prognosis=="PPx")

GPx_subset_pos = subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK+" & Prognosis %in% "GPx")
PPx_subset_pos = subset(ProgData_downstream, CodeClass=="Endogenous" | CodeClass == "Negative", segment=="PanCK+" & Prognosis %in% "PPx")


save(file="DataSubsets_fordownstreamAnalysis.RData", targetDataSubset_pos, targetDataSubset_neg, GPx_subset_neg, PPx_subset_neg, targetDataSubset_gpx, targetDataSubset_ppx, GPx_subset_pos,PPx_subset_pos)



x1= assayDataElement(GPx_subset_pos , elt = "log_q")
colnames(x1)<- sData(GPx_subset_pos)$ROI_REL

x2= assayDataElement(GPx_subset_neg , elt = "log_q")
colnames(x2)<- sData(GPx_subset_neg)$ROI_REL


x3= assayDataElement(PPx_subset_pos , elt = "log_q")
colnames(x3)<- sData(PPx_subset_pos)$ROI_REL

x4= assayDataElement(PPx_subset_neg , elt = "log_q")
colnames(x4)<- sData(PPx_subset_neg)$ROI_REL

save(file="For_ssGSEA.RData", x1, x2, x3, x4) #This served as input to the ssGSEA software available through Broad.

