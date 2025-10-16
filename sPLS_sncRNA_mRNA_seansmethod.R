# Purpose : sPLS analysis on TMM standardised sncRNA and mRNA data.

# I'm using the TMM results to get the 100 sncRNAs and 100 genes with the most variance

# Fixed the issue where sPLS plots were not displaying the colours I set as level_colors

# Based on https://mixomics.org/case-studies/spls-liver-toxicity-case-study-2/

# Date : 23/8/25
# Author: Chishan

# try next:
# next try the 1682 TMM sncRNAs x the top 2000 genes
# keep in novel miRNAs and piRNAs cause the info is still valuable
# omit arrow plot
# don't forget to send the tuning plots, Q2 plots etc

##########################################

# Load packages
library(mixOmics) # import the mixOmics library
library(dplyr)
library(tibble)
#install.packages('matrixStats')
#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
library(matrixStats)
library(data.table)

# 1. Data preparation. 

sncRNA <- read.csv("./3_results/sncRNA_july_TMM_matrix.csv", row.names = 1, check.names = FALSE)
# 1682 rows
genes  <- read.csv("3_results/transcriptomic_rerun_TMM_matrix_genes.csv",  row.names = 1, check.names = FALSE)
# 17 862 rows

# I wanted get the 100 sncRNAs and 100 genes with the most variance

ntop <- 100

# calculate the variance for each gene
rv_genes <- rowVars(as.matrix(genes))
rv_sncRNA<- rowVars(as.matrix(sncRNA))

# select the ntop genes by variance
select_genes <- order(rv_genes, decreasing = TRUE)[seq_len(min(ntop, length(rv_genes)))]

select_sncRNA <- order(rv_sncRNA, decreasing = TRUE)[seq_len(min(ntop, length(rv_sncRNA)))]



top_genes <- genes[select_genes, , drop = FALSE]

top_sncRNA <- sncRNA[select_sncRNA, , drop = FALSE]

# check dimensions
dim(top_genes)
dim(top_sncRNA)

#############


# Y (outcome) = gene
# X (predictor) = _RNA


# In this script, we're exploring sncRNA x mRNA interactions. Hence,
# Transposable element (TE) and other sncRNA data (piRNA, rRNA etc) will be excluded.



# In case I only want to use mature miRNAs
#miRNA <- miRNA %>%
#  filter(grepl(" mature", mirna_id))



# Use these if I wanna remove the 'mature' and 'hairpin' suffixes
#miRNA_raw$mirna_id <- gsub(" mature$", "", miRNA_raw$mirna_id)
#miRNA_raw$mirna_id <- gsub(" hairpin$", "", miRNA_raw$mirna_id)


#miRNA <- column_to_rownames(miRNA, var = "mirna_id")

# 2. Currently, samples are in columns and features are in rows.
# Transpose them so they can go into the PCA and sPLS correctly.

# This is so that the function interprets the samples as 'observations' and the 
# We do this because spls() expects the format of 8 rows (samples) and many columns (variables).

X <- t(as.matrix(top_sncRNA))   # _RNA → predictor block
Y <- t(as.matrix(top_genes))     # gene   → outcome block

# 3. Preliminary Analysis with PCA
pca.sncRNA <- pca(X, ncomp = 6, center = TRUE, scale = TRUE)
pca.genes <- pca(Y, ncomp = 6, center = TRUE, scale = TRUE)

# Barplot of the variance each principal component explains of the seminal vesicle data.
plot(pca.sncRNA) 

#png("./2_figures/aug25th_plots/pca_sncRNA_august.png", width = 1200, height = 900, res = 150)
#plot(pca.sncRNA)
#dev.off()

plot(pca.genes)

#png("./2_figures/aug25th_plots/pca_genes_august.png", width = 1200, height = 900, res = 150)
#plot(pca.genes)
#dev.off()


# 3. Make PCA plots

# Read in the treatment info

# A = Acrylamide, P = Control
sample_metadata <- read.csv("./1_data/sncRNA_sample_metadata.csv")

sample_metadata$Treatment_group <- as.factor(sample_metadata$Treatment_group)
levels(sample_metadata$Treatment_group)

level_colors <- c("Acr" = "#994455", "Control" = "#0077BB")

# Not sure why I'm not getting the colours I want, but that can be fixed later

# X
plotIndiv(pca.sncRNA, comp = c(1, 2), 
          group = sample_metadata$Treatment_group, 
          ind.names = sample_metadata$SampleID, 
          legend = TRUE, title = 'SV sncRNA, PCA comp 1 - 2')



png("./2_figures/aug25th_plots/mixomics_pca_X_sncRNA.png", width = 1600, height = 1200, res = 150)
plotIndiv(pca.sncRNA,
          comp = c(1, 2),
          group = sample_metadata$Treatment_group,
          ind.names = sample_metadata$SampleID,
          legend = TRUE,
          title = 'SV sncRNA, PCA comp 1 - 2')
dev.off()


# Y
plotIndiv(pca.genes, comp = c(1, 2), 
          group = sample_metadata$Treatment_group, 
          ind.names = sample_metadata$SampleID, 
          legend = TRUE, title = 'SV genes, PCA comp 1 - 2')

png("./2_figures/aug25th_plots/mixomics_pca_Y_genes.png", width = 1600, height = 1200, res = 150)
plotIndiv(pca.genes,
          comp = c(1, 2),
          group = sample_metadata$Treatment_group,
          ind.names = sample_metadata$SampleID,
          legend = TRUE,
          title = 'SV genes, PCA comp 1 - 2')
dev.off()

# 4. Initial sPLS model
spls.SV <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression') #########################trying canonical
#spls.SV <- spls(X = X, Y = Y, ncomp = 5, mode = 'canonical')

# 5. Tuning sPLS model by altering number of components

# repeated CV tuning of component count
perf.spls.SV  <- perf(spls.SV, validation = 'Mfold',
                      folds = 2, nrepeat = 30) # Sean reccomended Mfold = 2, repeats = 30. I increased folds to 5 to prevent warnings

plot(perf.spls.SV, criterion = 'Q2.total') 



png("./2_figures/aug25th_plots/mixomics_initial_spls.png", width = 1600, height = 1200, res = 150)
plot(perf.spls.SV, criterion = 'Q2.total')
dev.off()


# 6. Selecting the number of variables

# set range of test values for number of variables to use from X dataframe
list.keepX <- seq(5, 15, 5)# edited from 20, 50, 5
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 


tune.spls.SV <- tune.spls(X, Y, ncomp = 2,
                          test.keepX = list.keepX,
                          test.keepY = list.keepY,
                          nrepeat = 10, folds = 2, # Sean recommended Mfold = 2, repeats = 10 (adjusted from 30)
                          mode = 'canonical', measure = 'cor') 

plot(tune.spls.SV)         # use the correlation measure for tuning

# I got a bunch of this warning at this stage:
#In internal_wrapper.mint(X = X, Y = Y, ncomp = ncomp,  ... :
#                           At least one study has less than 5 samples, mean centering might
#                         not do as expected

png("./2_figures/aug25th_plots/mixomics_tuning_spls.png", width = 1600, height = 1200, res = 150)
plot(tune.spls.SV)  
dev.off()



# The optimal number of features to use for both datasets can be extracted through the below calls.
tune.spls.SV$choice.keepX
tune.spls.SV$choice.keepY

# These values will be stored to form the final model

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.SV$choice.keepX

# extract optimal number of variables for Y dataframe
optimal.keepY <- tune.spls.SV$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

# 7. Final model. Using the tuned parameters generated above, the final sPLS model can be constructed.
# use all tuned values from above
final.spls.SV <- spls(X, Y, ncomp = optimal.ncomp, 
                      keepX = optimal.keepX,
                      keepY = optimal.keepY,
                      mode = "regression") # explanatory approach being used, 
# hence use regression mode

# 8. Plots!

# Sample plots (must edit to fit my dataset)

# Plot the X variate space
plotIndiv(
  final.spls.SV,
  comp = c(1, 2),
  ind.names = sample_metadata$SampleID,
  rep.space = "X-variate",
  group = sample_metadata$Treatment_group,
  col.per.group = level_colors,
  legend = TRUE,
  legend.title = 'Treatment group',
  main = "SV sPLS – X variate space",
  title = "SV sPLS – X variate space"
)


png("./2_figures/aug25th_plots/mixomics_final_spls_X.png", width = 1200, height = 800, res = 150)
plotIndiv(
  final.spls.SV,
  comp = c(1, 2),
  ind.names = sample_metadata$SampleID,
  rep.space = "X-variate",
  group = sample_metadata$Treatment_group,
  col.per.group = level_colors,
  legend = TRUE,
  legend.title = 'Treatment group',
  main = "SV sPLS – X variate space",
  title = "SV sPLS – X variate space"
)
dev.off()


# Plot the Y variate space


plotIndiv(
  final.spls.SV,
  comp = c(1, 2),
  ind.names = sample_metadata$SampleID,
  rep.space = "Y-variate",
  group = sample_metadata$Treatment_group,
  col.per.group = level_colors,
  legend = TRUE,
  legend.title = 'Treatment group',
  title = "SV sPLS – Y variate space"
)


png("./2_figures/aug25th_plots/mixomics_final_spls_Y.png", width = 1200, height = 800, res = 150)
plotIndiv(
  final.spls.SV,
  comp = c(1, 2),
  ind.names = sample_metadata$SampleID,
  rep.space = "Y-variate",
  group = sample_metadata$Treatment_group,
  col.per.group = level_colors,
  legend = TRUE,
  legend.title = 'Treatment group',
  title = "SV sPLS – Y variate space"  
)
dev.off()


# Plot in XY variate space
plotIndiv(
  final.spls.SV,
  ind.names = sample_metadata$SampleID,     
  rep.space = "XY-variate",                 
  group = sample_metadata$Treatment_group, 
  col.per.group = level_colors,          
  legend = TRUE,
  legend.title = "Treatment",
  title = "SV sPLS – XY variate space"
)


png("./2_figures/aug25th_plots/mixomics_final_spls_XY.png", width = 1200, height = 800, res = 150)
plotIndiv(
  final.spls.SV,
  ind.names = sample_metadata$SampleID,     
  rep.space = "XY-variate",                 
  group = sample_metadata$Treatment_group, 
  col.per.group = level_colors,          
  legend = TRUE,
  legend.title = "Treatment",
  title = "SV sPLS – XY variate space"
)
dev.off()


# Arrow plot

# level_colors was set wayyy earlier

# Plot arrows
plotArrow(
  final.spls.SV,
  ind.names = TRUE,                     # show sample labels if desired
  group = sample_metadata$Treatment_group,  # colour by treatment
  col.per.group = level_colors,
  legend = TRUE,
  legend.title = 'Treatment Group',
  title = "SV sPLS – Arrow Plot"
)

png("./2_figures/aug25th_plots/mixomics_final_spls_arrow.png", width = 1200, height = 800, res = 150)
plotArrow(
  final.spls.SV,
  ind.names = TRUE,                     # show sample labels if desired
  group = sample_metadata$Treatment_group,  # colour by treatment
  col.per.group = col.treat,
  legend = TRUE,
  legend.title = 'Treatment Group',
  main = "SV sPLS – Arrow Plot"
)
dev.off()



# 1. Performance of final model using repeated M-fold CV
perf.spls.SV <- perf(final.spls.SV, 
                     folds = 2, nrepeat = 30,  # fewer folds due to small sample size
                     validation = "Mfold", 
                     dist = "max.dist",       # use max.dist measure
                     progressBar = FALSE)

# 2. Plot feature stability for the first two components
par(mfrow = c(1, 2))

# Component 1
plot(perf.spls.SV$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las = 2,
     xlim = c(0, ncol(X)))  # X has 19 variables


png("./2_figures/aug25th_plots/mixomics_final_spls_featurestability1.png", width = 1200, height = 800, res = 150)

plot(perf.spls.SV$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las = 2,
     xlim = c(0, ncol(X)))  # X has 19 variables

dev.off()

# Component 2
plot(perf.spls.SV$features$stability.X[[2]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las = 2,
     xlim = c(0, ncol(X)))  # same as above

png("./2_figures/aug25th_plots/mixomics_final_spls_featurestability2.png", width = 1200, height = 800, res = 150)

# Component 2
plot(perf.spls.SV$features$stability.X[[2]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las = 2,
     xlim = c(0, ncol(X)))  # same as above

dev.off()

# 3. Variable plot
plotVar(final.spls.SV, cex = c(3, 4), var.names = c(FALSE, TRUE))

png("./2_figures/aug25th_plots/mixomics_final_spls_correlationcircle.png", width = 1200, height = 800, res = 150)

plotVar(final.spls.SV, cex = c(3, 4), var.names = c(FALSE, TRUE))

dev.off()

# 4. Correlation circle plot (optional)
#color.edge <- color.GreenRed(50)  # set edge colors for connecting lines ## ngl I have no idea what this line does


# X11() # To open a new window for Rstudio

network(final.spls.SV, comp = 1:2,       # use the 2 components of model
        cutoff = 0.7,                    # only show connections with correlation > 0.7
        shape.node = c("rectangle", "circle"), 
        color.node = c("cyan", "pink"),  # X vs Y variables
        color.edge = color.edge,         # previously defined color.GreenRed(50)
        save = 'png',                    # save as PNG
        name.save = 'sPLS_SV_Network_Plot')  # filename


# Network representation 
par(mar = c(20, 10, 4, 20))  
cim(final.spls.SV, comp = 1:2, xlab = "Genes", ylab = "sncRNAs",
    cex.lab = 1,      # axis labels
    cex.axis = 0.8)   # tick labels

png("./2_figures/aug25th_plots/mixomics_final_spls_network.png", width = 1600, height = 1200, res = 150)
par(mar = c(20, 10, 4, 20))  
cim(final.spls.SV, comp = 1:2, xlab = "Genes", ylab = "sncRNAs")
dev.off()


################ Troubleshooting to find out how to make labels readable
png("./2_figures/mixomics_final_network2spls_hairpinmaturemirna100.png", 
    width = 2000, height = 1400, res = 150)

# Increase margins and give more room for labels
par(mar = c(16, 14, 4, 6))  

# CIM heatmap with smaller cells (via smaller labels)
cim(final.spls.SV, comp = 1:2, 
    xlab = "Genes", ylab = "sncRNAs",
    cex.row = 0.4,    # smaller y-axis labels
    cex.col = 0.6,    # smaller x-axis labels
    margins = c(12, 12) # adds extra space around heatmap
)

dev.off()


