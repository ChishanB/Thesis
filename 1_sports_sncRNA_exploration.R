
# Written: 26/3/24
#Purpose: Explore pre-analysed sncRNA dysregulation data from SPORTS1.0 output

#Author: Chishan Burch
#Contact: chishanburch@gmail.com

##################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)


#Threshold for FDR which is the same as p-adjusted
FDR_threshold <- 0.05


#Start with miRNA
miRNA <- fread("./1_data/miRNA_SCH_glmLRT_FC.csv")


#Adding column for significance and up/down regulation (NS means not significant)
miRNA <- miRNA %>%
  mutate(Significance = if_else(FDR < FDR_threshold & logFC < 0 & logFC < -1.5, "Downreg",
                           if_else(FDR < FDR_threshold & logFC > 0 & logFC > 1.5, "Upreg", "NS")
                          
        )) 

#Adding column for -log10FDR
miRNA <- miRNA %>% mutate(log_FDR = -log10(FDR))

#Declare graph colours and factor values
palette <- c("Upreg" = "#D55E00", "Downreg" = "#0072B2", "NS" = "#808080")

miRNA$Significance <- factor(miRNA$Significance, levels = c("Upreg","Downreg", "NS"))
 
#plot miRNAs
miRNA_plot <- ggplot(data = miRNA, aes(x =logFC, y= log_FDR)) + geom_point() + theme_minimal() +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = -1.5, linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = 1.5, linetype = "dashed", col = "#999999") +
  aes(x = logFC, y = log_FDR, col = Significance, fill = Significance) +
  scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
  xlab("log2FoldChange") + ylab("-log10(FDR)") + theme(legend.title = element_blank(), 
  axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
        legend.text = element_text(size = 14))
  
ggsave(plot = miRNA_plot, filename = "./2_figures/miRNA_volcano_plot.png", width = 6, height = 6, dpi = 600)



#tRNA
tRNA <- read.xlsx("./1_data/Sch_edgeR_glmRT.xlsx", sheet = "tRNA")

#Adding column for significance and up/down regulation (NS means not significant)
tRNA <- tRNA %>%
  mutate(Significance = if_else(FDR < FDR_threshold & logFC < 0, "Downreg",
                                if_else(FDR < FDR_threshold & logFC > 0, "Upreg", "NS")
  )) 
#Adding column for -log10FDR
tRNA <- tRNA %>% mutate(log_FDR = -log10(FDR))

tRNA$Significance <- factor(tRNA$Significance, levels = c("Upreg","Downreg", "NS"))

#plot tRNAs
tRNA_plot <- ggplot(data = tRNA, aes(x =logFC, y= log_FDR)) + geom_point() + theme_minimal() +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = -1.5, linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = 1.5, linetype = "dashed", col = "#999999") +
  aes(x = logFC, y = log_FDR, col = Significance, fill = Significance) +
  scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
  xlab("log2FoldChange") + ylab("-log10(FDR)") + theme(legend.title = element_blank(), 
                                                   axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
                                                       legend.text = element_text(size = 14))

ggsave(plot = tRNA_plot, filename = "./2_figures/tRNA_volcano_plot.png", width = 6, height = 6, dpi = 600)

#rRNA
rRNA <- read.xlsx("./1_data/Sch_edgeR_glmRT.xlsx", sheet = "rRNA")

#Adding column for significance and up/down regulation (NS means not significant)
rRNA <- rRNA %>%
  mutate(Significance = if_else(FDR < FDR_threshold & logFC < 0, "Downreg",
                                if_else(FDR < FDR_threshold & logFC > 0, "Upreg", "NS")
  )) 
#Adding column for -log10FDR
rRNA <- rRNA %>% mutate(log_FDR = -log10(FDR))

rRNA$Significance <- factor(tRNA$Significance, levels = c("Upreg","Downreg", "NS"))

#plot rRNAs
rRNA_plot <- ggplot(data = rRNA, aes(x =logFC, y= log_FDR)) + geom_point() + theme_minimal() +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = -1.5, linetype = "dashed", col = "#999999") +
  geom_vline(xintercept = 1.5, linetype = "dashed", col = "#999999") +
  aes(x = logFC, y = log_FDR, col = Significance, fill = Significance) +
  scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
  xlab("log2FoldChange") + ylab("-log10(FDR)") + theme(legend.title = element_blank(), 
                                                       axis.text = element_text(size = 14), axis.title = element_text(size = 18), 
                                                       legend.text = element_text(size = 14))

ggsave(plot = rRNA_plot, filename = "./2_figures/rRNA_volcano_plot.png", width = 6, height = 6, dpi = 600)








"#388ECC", "#F68B33", "#009E73", 
"#CC79A7", "#F0E442", "black", "#D55E00", "#0072B2", 
"#999999"



