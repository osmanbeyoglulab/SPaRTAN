#########################################################################################################################
# Simple example for visualizing SPaRTAN outputs
#########################################################################################################################
rm(list=ls())
options(stringsAsFactors = FALSE)
## We load the required packages
require(ggplot2)
require(ggrepel)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(pvclust)
require(circlize)
require(reshape)
dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
##read input files
P = read.csv("../data/inputs/Ppbmc5kn_CD16.csv", row.names = 1) ##cell by surface protein based on scADT-seq
Y = read.csv("../data/inputs/Ypbmc5kn_CD16.csv", row.names = 1) ##genes by cells based on scRNA-seq
D = read.csv("../data/inputs/Dpbmc.csv", row.names = 1) ##genes by transcription factors
##read SPaRTAN output
#cell-type SPaRTAN-predicted TF activity
TFa <- read.csv("../data/outputs/projD.csv", row.names = 1)
scaledTFa <- scale(t(TFa))
scaledP <- scale(P)


# adt.pv <- pvclust(scaledP, nboot=100, method.hclust = "ward.D2", method.dist = "correlation")
# cell.pv <- pvclust(t(scaledP), nboot=100, method.hclust = "ward.D2", method.dist = "correlation")
# tf.pv <- pvclust(scaledTFa, nboot=100, method.hclust = "ward.D2", method.dist = "correlation")

ctype = "cd16+ pbmc cells "
# Plot the Heatmap of surface protein and SPaRTAN-predicted TF activities
my_palette = colorRamp2(seq(-0.7, 0.7, length = 3), c("#377eb8", "#EEEEEE", "#e41a1c"))
hp = Heatmap(scaledP, show_column_names = T, show_row_names = F, column_names_gp = gpar(fontsize = 3.5), 
             name="surface protein", col=my_palette)  
my_palette =  colorRamp2(seq(-2, 2, length = 3), c("#377eb8", "#EEEEEE", "#4daf4a"))
ht = Heatmap(scaledTFa, show_column_names = F, show_row_names = F, name="TF activity",col=my_palette)
hp+ht
##Calculate Pearson’s Correlation between surface protein expression and SPaRTAN-predicted TF activities
pcc <- cor(P,t(TFa))

## Plot heatmaps for  Pearson`s correlations
Heatmap(pcc, col = brewer.pal(n=9,name="PRGn"),  show_column_names = F, show_row_names = T, row_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(color_bar = "continuous", grid_height = unit(2, "mm"), grid_width = unit(2, "mm"), 
                                    title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
        name = "correlation", width = 1.2,  row_gap = unit(0.4, "mm"), column_gap = unit(0.2, "mm"), border = F)

## Plot heatmaps for filtered Pearson`s correlations
cutoff <- 0.7
pcc_filtered <- pcc[names(which(apply(abs(pcc),1,max) > cutoff)),names(which(apply(abs(pcc),2,max) > cutoff))]
Heatmap(pcc_filtered, col = brewer.pal(n=9,name="PRGn"),  show_column_names = F, show_row_names = T, row_names_gp = gpar(fontsize = 7),
        heatmap_legend_param = list(color_bar = "continuous", grid_height = unit(2, "mm"), grid_width = unit(2, "mm"), 
                                    title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
        name = "correlation", width = 1.2,  row_gap = unit(0.4, "mm"), column_gap = unit(0.2, "mm"), border = F)

## Plot correlations between a particular surface protein and SPaRTAN-predicted TF activities
datLong <- melt(data = pcc)
colnames(datLong) <- c("adt","tf","PCC")
surfaceProtein <- "CD278"
dat <- datLong[datLong$adt %in% surfaceProtein,]
dat <- dat[order(dat$PCC),]
dat <- cbind(dat, 1:nrow(dat))
colnames(dat) <- c("adt","tf","PCC","rank")

th=quantile(abs(dat$PCC),probs = seq(0, 1, 0.01))[91]
p <- ggplot(dat, aes(x=rank, y=PCC, label=tf)) +
  geom_point()  +  ylim(-1, 1) + xlim(-20, 280) + #ggtitle(paste0(type," cell")) +
  xlab("Transcription factors") +
  ylab(paste0("Correlation betw ", surfaceProtein, " expr and TF activity")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_label_repel(data=subset(dat, abs(PCC) > 0.8), aes(label=tf), size = 3, max.overlaps = 30 ) +
  theme(legend.position='none') + theme(axis.title.y =element_text(size=15)) +
  theme(axis.title.x =element_text(size=15)) +
  theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"))
print(p)

