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



load(file="Consolidated_AllPlotting.RData")
load(file="DataSubsets_fordownstream.RData")#From 1_Batch10_QC_forGitHub.R
#------------------------------------------------------------------------#

# formula follows conventions defined by the lme4 package
runLMM.enmass.across <- function (targetData, testClass, interceptTerm ){#runs all possible testClass combinations across slides
  results2 <- c()
  newD = targetData
  pData(newD)$testClass <-factor(pData(newD)[[testClass]])
  pData(newD)[["slide"]] <- factor(pData(newD)[[interceptTerm]])
  mixedOutmc <-
    mixedModelDE(newD,
                 elt = "log_q",
                 modelFormula = ~ testClass+ (1 | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <-
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate",
                       "Pr(>|t|)", "FDR")]
  results2 <- rbind(results2, r_test)
  
  # Categorize results2 based on P‐value & FDR for plotting
  results2$Color <- "NS"
  results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
  results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
  #results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
  results2$Color <- factor(results2$Color,
                           levels = c("NS", "P < 0.05","FDR < 0.05", "FDR < 0.001"))
  
  condNames<- names(table(results2$Contrast))
  
  results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
  top_g <- c()
  for(cond in condNames) {
    ind <- results2$Contrast == cond
    top_g <- c(top_g,
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]])
  }
  top_g <- unique(top_g)
  results2 <- results2[,-1*ncol(results2)] # remove invert_P from matrix
  
  listR <- list(results= results2, top_g= top_g)
  return(listR)
}

#-----------------------------------------------------------------------#

runLMM.enmass.within <- function (targetData, testClass, interceptTerm){ #This basically runs all possible combinations of the chose testclass
  results2 <- c()
  newD= targetData
  pData(newD)$testClass <- factor(pData(newD)[[testClass]])
  pData(newD)[["slide"]] <- factor(pData(newD)[[interceptTerm]])
  
  mixedOutmc <-
    mixedModelDE(newD,
                 elt = "log_q",
                 modelFormula = ~ testClass+ (1+testClass | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <-
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate",
                       "Pr(>|t|)", "FDR")]
  results2 <- rbind(results2, r_test)
  
  # Categorize results2 based on P‐value & FDR for plotting
  results2$Color <- "NS"
  results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
  results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
  #results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
  results2$Color <- factor(results2$Color,
                           levels = c("NS", "P < 0.05","FDR < 0.05", "FDR < 0.001"))
  
  condNames<- names(table(results2$Contrast))
  
  results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
  top_g <- c()
  for(cond in condNames) {
    ind <- results2$Contrast == cond
    top_g <- c(top_g,
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]])
  }
  top_g <- unique(top_g)
  results2 <- results2[,-1*ncol(results2)] # remove invert_P from matrix
  
  listR <- list(results= results2, top_g= top_g)
  return(listR)
}

#-----------------------------------------------------------------------#

runLMM.subset.across <-function(subset= subset, subsetName= subsetName, targetData= targetData, interceptTerm= interceptTerm, testClass= testClass){
  
  results2 <- c()
  for (rely in subset){ 
    ind <- pData(targetData)[,subsetName] == rely
    newD= targetData[,ind]
    pData(newD)$testClass <-factor(pData(newD)[[testClass]])
    pData(newD)[["slide"]] <- factor(pData(newD)[[interceptTerm]])
    mixedOutmc <- mixedModelDE(newD,
                               elt = "log_q",
                               modelFormula = ~ testClass + (1 | slide),
                               groupVar = "testClass",
                               nCores = parallel::detectCores(),
                               multiCore = FALSE)
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Gene <-unlist(lapply(colnames(mixedOutmc),
                                rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- rely
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results2 <- rbind(results2, r_test)
  }
  
  
  results2$Color <- "NS"
  results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
  results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
  #results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
  results2$Color <- factor(results2$Color,
                           levels = c("NS", "P < 0.05","FDR < 0.05", "FDR < 0.001"))
  
  # Top_g for Volcano plotting
  results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
  top_g <- c()
  for(cond in c(subset)) {
    ind <- results2$Subset == subset
    top_g <- c(top_g,
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]])
  }
  top_g <- unique(top_g)
  results2 <- results2[,-1*ncol(results2)]
  
  listR <- list(results= results2, top_g= top_g)
  return (listR)
  
}

#---------------------------------------------------------------------------------#

runLMM.subset.within <-function(subset= subset, categoryName= categoryName, targetData= targetData, interceptTerm= interceptTerm, testClass=testClass){
  
  results2 <- c()
  for (rely in subset){ 
    ind <- pData(targetData)[,categoryName] == rely
    newD= targetData[,ind]
    pData(newD)$testClass <-factor(pData(newD)[[testClass]])
    pData(newD)[["slide"]] <- factor(pData(newD)[[interceptTerm]])
    mixedOutmc <- mixedModelDE(newD,
                               elt = "log_q",
                               modelFormula = ~ testClass + (1+testClass | slide),
                               groupVar = "testClass",
                               nCores = parallel::detectCores(),
                               multiCore = T)
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Gene <-unlist(lapply(colnames(mixedOutmc),
                                rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- rely
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results2 <- rbind(results2, r_test)
  }
  
  
  results2$Color <- "NS"
  results2$Color[results2$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results2$Color[results2$FDR < 0.05] <- "FDR < 0.05"
  results2$Color[results2$FDR < 0.001] <- "FDR < 0.001"
  #results2$Color[abs(results2$Estimate) < 0.5] <- "NS or FC < 0.5"
  results2$Color <- factor(results2$Color,
                           levels = c("NS", "P < 0.05","FDR < 0.05", "FDR < 0.001"))
  
  # Top_g for Volcano plotting
  results2$invert_P <- (-log10(results2$`Pr(>|t|)`)) * sign(results2$Estimate)
  top_g <- c()
  for(cond in c(subset)) {
    ind <- results2$Subset == subset
    top_g <- c(top_g,
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results2[ind, 'Gene'][
                 order(results2[ind, 'invert_P'], decreasing = FALSE)[1:15]])
  }
  top_g <- unique(top_g)
  results2 <- results2[,-1*ncol(results2)]
  
  listR <- list(results= results2, top_g= top_g)
  return (listR)
  
}

#-----------------------------------------------------------------------#


############## DiffEx 7: In the TME- how does edge vs center differ####################

#In GPX
{
  
  #Finding all slides which have both edge and center from GPX
  metData66 <- pData(GPx_subset_neg)
  test66<- as.data.frame(metData66 %>% group_by(metData66$`slide name`,metData66$Region) %>% summarise(count=n()))
  colnames(test66) <- c( "Slide Name", "Region" ,"#AOI in GPx Neg")
  
  
  i_ce <- intersect(test66$`Slide Name`[test66$Region %in% c("E")],
                    test66$`Slide Name`[test66$Region %in% c("C")])
  
  subset_ice1 <- subset(GPx_subset_neg, CodeClass=="Endogenous" | CodeClass == "Negative", 
                        `slide name`%in% i_ce & Region %in% c("C","E") )  , 
  
  
  testClass = "Region"
  subset= c("PanCK-")
  categoryName= "segment"
  interceptTerm = "slide name"
  targetData= subset_ice1
  
  results.99 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)
  
}


#In PPX
{
  #Finding all slides which have both edge and center from GPX
  metData76 <- pData(PPx_subset_neg)
  test76<- as.data.frame(metData76 %>% group_by(metData76$`slide name`,metData76$Region) %>% summarise(count=n()))
  colnames(test76) <- c( "Slide Name", "Region" ,"#AOI in GPx Neg")
  
  
  i_ce <- intersect(test76$`Slide Name`[test76$Region %in% c("E")],
                    test76$`Slide Name`[test76$Region %in% c("C")])
  
  
  subset_ice11 <- subset(PPx_subset_neg, CodeClass=="Endogenous" | CodeClass == "Negative", 
                         `slide name`%in% i_ce & Region %in% c("C","E") )  , 
  
  
  testClass = "Region"
  subset= c("PanCK-")
  categoryName= "segment"
  interceptTerm = "slide name"
  targetData= subset_ice11
  
  results.90 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)
  
  }
  
  
}
  save(file="ProgData_DiffEX7_Neg_CvsE.RData", results.90, results.99, c90, c99)


############ DiffEx of Tenr subsets for Fig 3b###################################

gpx_topT <-  # 61 AOIs that were annotated as high T content
ppx_topT <-  # 40 AOIs that were annotated as high T content

#Extract only these AOIs from the original PPx_subset_neg object
subset_tenr <- subset(targetDataSubset_neg, CodeClass=="Endogenous" | CodeClass == "Negative", 
                          pData(targetDataSubset_neg)$ROI_REL %in% c(ppx_topT, gpx_topT) )  

#Performing across sample differential expression Tenr       

subsetName = "segment"
subset <- c( "PanCK-")
interceptTerm= "slide name"
testClass = "Prognosis"
targetData = subset_tenr

# run LMM: by subset 
results<- runLMM.subset.across(targetData = targetData, subset= subset, subsetName = subsetName, testClass = testClass, interceptTerm = interceptTerm)

{
  results$results$LabelColor <- ifelse(results$results$Estimate > 0, "lightpink3", "lightblue3")
  c111= ggplot(results$results, aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                                        color = Color, label = Gene)) +
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
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
    geom_label_repel(
      data = subset(results$results, FDR < 0.05),
      aes(label = Gene),
      color = subset(results$results, FDR < 0.05)$LabelColor,
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
  
  
}

save(file="4.2_TEnr_DiffEx.RData", results, c111)

##############  Invasive Epithelia- how does edge vs center differ- So exploring only in the pos region####################

#For PanCK+ we utilize an NMF reduced gene sets for PPx and GPx which can be generated as outlined in teh methods.

#In GPX

#Finding all slides which have both edge and center from GPX
metData6 <- pData(GPxpos_nmfsub)
test6<- as.data.frame(metData6 %>% group_by(metData6$Epithelial, metData6$`slide name`,metData6$Region, metData6$R_E) %>% summarise(count=n()))
colnames(test6) <- c("Epithelia", "Slide Name", "Region" , "R_E" ,"#AOI in gpx pos")

i_ce <- intersect(test6$`Slide Name`[test6$R_E %in% c("E_I")],
                  test6$`Slide Name`[test6$R_E %in% c("C_I")])

subset_ice1 <- subset(GPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                      `slide name`%in% i_ce & R_E %in% c("C_I","E_I") ) 


testClass = "Region"
subset= c("I")
categoryName= "Epithelial"
interceptTerm = "slide name"
targetData= subset_ice1

results.9 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)

###################-----------doing the same fro PPX--------------########################
#In PPX
metData7 <- pData(PPxpos_nmfsub)
test7<- as.data.frame(metData7 %>% group_by(metData7$Epithelial, metData7$`slide name`,metData7$Region, metData7$R_E) %>% summarise(count=n()))
colnames(test7) <- c("Epithelia", "Slide Name", "Region" , "R_E", "#AOI in ppx pos")

ggplot(test7, aes(x= Epithelia,y= `#AOI in ppx pos`, fill= Region)) +geom_bar(stat="identity", position="dodge")+ theme_minimal()+
  theme(legend.text= element_text(face="bold",size=10),legend.title=element_text(face="bold",size=10), axis.text = element_text(angle = 90, size=10, face="bold"), 
        axis.title = element_text(face="bold",size=12),panel.background = element_blank())


i_ce2 <- intersect(test7$`Slide Name`[test7$R_E %in% c("E_I")],
                   test7$`Slide Name`[test7$R_E %in% c("C_I")])
 subset_ice2 <- subset(PPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                      `slide name`%in% i_ce2 & R_E %in% c("C_I","E_I") )


testClass = "Region"
subset= c("I")
categoryName= "Epithelial"
interceptTerm = "slide name"
targetData= subset_ice2

results.10 <- runLMM.subset.within(subset= subset, categoryName= categoryName,targetData= targetData, interceptTerm=interceptTerm, testClass= testClass)

save(file="DiffEx9_NMF_GPX_PPX_CvsE_Invasive.RData",results.9, results.10)


####### DiffEx comparing all invasive epithelial  vs all normal adjacent epithelia J in PPx
subsetIJ_PPx3 <- subset(PPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                        `slide name` %in% iall_jall & Epithelial %in% c("I","J") ) 


testClass = "Epithelial"
subset= c("PanCK+")
categoryName= "segment"
interceptTerm = "slide name"
targetData= subsetIJ_PPx3

results.7 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)

########DiffEx comparing all invasive epithelial  vs all normal adjacent epithelia J in PPx

subsetAI_PPx4 <- subset(PPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                        `slide name` %in% iall_aall & Epithelial %in% c("I","A") ) 


testClass = "Epithelial"
subset= c("PanCK+")
categoryName= "segment"
interceptTerm = "slide name"
targetData= subsetAI_PPx4

results.8 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)

########
save(file="9_DiffEx_nmfPPX.RData",  results.7,results.8 )
########

########DiffEx comparing all invasive epithelial  vs all normal adjacent epithelia J in GPx

subsetIJ_GPx3 <- subset(GPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                        `slide name` %in% iall_jall & Epithelial %in% c("I","J") ) 


testClass = "Epithelial"
subset= c("PanCK+")
categoryName= "segment"
interceptTerm = "slide name"
targetData= subsetIJ_GPx3

results.7 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)

########DiffEx comparing all invasive epithelial  vs all normal adjacent epithelia J in GPx

subsetIJ_GPx4 <- subset(GPxpos_nmfsub, CodeClass=="Endogenous" | CodeClass == "Negative", 
                        `slide name` %in% iall_aall & Epithelial %in% c("I","A") ) 


testClass = "Epithelial"
subset= c("PanCK+")
categoryName= "segment"
interceptTerm = "slide name"
targetData= subsetIJ_GPx4

results.8 <- runLMM.subset.within(subset, categoryName,targetData= targetData, interceptTerm, testClass)

save(file="9_DiffEx_nmfGPX.RData",  results.7,results.8 )
#