############ Plotting ssGSEA as a dotplot with the average score ##########
# gct files have to be generated using GenePattern ssGSEA

gpx_pos = mat(parse_gctx("gpx_pos_results.gct"))
gpx_neg = mat(parse_gctx("gpx_neg_results.gct"))
ppx_pos = mat(parse_gctx("ppx_pos_results.gct"))
ppx_neg = mat(parse_gctx("ppx_neg_results.gct"))

rownames(gpx_pos) <- gsub("HALLMARK_","", rownames(gpx_pos))
rownames(gpx_neg) <- gsub("HALLMARK_","", rownames(gpx_neg))
rownames(ppx_pos) <- gsub("HALLMARK_","", rownames(ppx_pos))
rownames(ppx_neg) <- gsub("HALLMARK_","", rownames(ppx_neg))

# find average score of gene sets in invasive and non-invasive for each matrix
gpx_pos_I <- gpx_pos[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(gpx_pos)) == "I" ] # 211
gpx_pos_NI <- gpx_pos[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(gpx_pos)) != "I" ] # 80

gpx_neg_I <- gpx_neg[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(gpx_neg)) == "I" ] # 222
gpx_neg_NI <- gpx_neg[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(gpx_neg)) != "I" ] # 91

ppx_pos_I <- ppx_pos[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(ppx_pos)) == "I" ] # 190
ppx_pos_NI <- ppx_pos[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(ppx_pos)) != "I" ] # 85

ppx_neg_I <- ppx_neg[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(ppx_neg)) == "I" ] # 183
ppx_neg_NI <- ppx_neg[,gsub("^(?:[^.]+\\.){3}([^.])(?:\\..*)$", "\\1", colnames(ppx_neg)) != "I" ] # 89

a <- rowMeans(gpx_pos_I)
b <- rowMeans(gpx_pos_NI)
c <- rowMeans(ppx_pos_I)
d <- rowMeans(ppx_pos_NI)

df <- data.frame(
  Pathway = names(a),
  `GPx-NI` = b,
  `GPx-I` = a,
  `PPx-NI` = d,
  `PPx-I` = c,
  check.names = FALSE
)

library(tidyr)
df_long <- pivot_longer(df, cols = -Pathway, names_to = "Condition", values_to = "Value")
df_long$Condition <- factor(df_long$Condition, levels = c("GPx-NI","GPx-I","PPx-NI","PPx-I"))
df_long <- df_long %>%
  mutate(Pathway = fct_reorder(Pathway, 
                               Value * (Condition == "PPx-I"), 
                               .fun = sum))
# Plot
library(ggplot2)
ggplot(df_long, aes(x = Condition, y = Pathway, color = Value, size = abs(Value))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal()

p1 <- ggplot(df_long, aes(x = Value, y = Pathway, color = Condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("GPx-NI" = "pink4", "GPx-I" = "lightpink","PPx-NI" = "steelblue1", "PPx-I" = "lightskyblue2")) +
  theme_bw() +
  theme(legend.text= element_text(size=12),legend.title=element_text(face="bold",size=12), 
        axis.text = element_text(size=12), axis.title = element_text(face="bold",size=14), axis.line = element_line(color="black"),
        panel.background = element_blank())+  
  labs(x = "Mean ssGSEA Score", y = "Pathway")

# add your path
pdf(file = "",width = 11, height = 12)
p1
dev.off()