
################################################################################
###################################################################################################
###################################################################################################
###Script for ATAC-Seq: using gplots to reproduce ArchR markers heatmap using top genes per cluster 	     ##
###################################################################################################
###################################################################################################
# README
# Maintained by Mar Muniz.  Last update July 2022.
# This scrip have a serie of functions that will allow to reproduce the markers heatmap 
# implemented in archR using the function  plotMarkerHeatmap in ggplot to be able to 
# customize the data used in this case we will do the plot selecting the top 50 Differentially
# expressed genes per cluster
# for more info in ArchR:
# ArchR doc https://github.com/GreenleafLab/ArchR/
# color palette archR: https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
###########################################################################################
# NOTES: data object from T.Griswold : Save-ArchR-Project.rds
###########################################################################################
#  R version 4.2.1 (2022-06-23) 
############################################################################################

library(gplots);
library("colorspace");library("Seurat");library("ArchR");library("dplyr");library("cowplot");
library("ggplot2");library("xlsx");library("grid");library("gridBase");library("gridExtra");
library("patchwork");library("ggrepel");library("Matrix");library("tidyr");library("RColorBrewer");
library("parallel");library("ggpubr");library("openxlsx");

#####################
#set seeed and wd
#####################

set.seed(22); # need to set it to get always the same random results and plots
wd <- "C:/Users/scATAC_snRNASeq/"
setwd(wd);

atacArch <- readRDS(paste0(wd,f_input,"ArchR-Project.rds"));
loadArchRProject(path = paste0(wd,f_input), force = TRUE, showLogo = TRUE);


#imputation with MAGIC
atacArch <- addImputeWeights(atacArch)

markersIntegrate <- getMarkerFeatures(
    ArchRProj = atacArch, 
    useMatrix = "GeneIntegrationMatrix",  
    groupBy = "integrated_snn_res.0.5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersIntegrateList <- getMarkers(markersIntegrate, cutOff ="FDR <= 0.05 & Log2FC >= 1 | FDR <= 0.05 & Log2FC <= 1");
valueLog2FC <- max(markersIntegrate@assays@data$Log2FC,abs(min(markersIntegrate@assays@data$Log2FC)));



#extract info for ggplot heatmap
log2FCmatrix<- as.matrix(markersIntegrate@assays@data$Log2FC)
colnames(log2FCmatrix)<-  markersIntegrate@colData@rownames
#extracting genomci location info
Gencoordinates<- as.data.frame(markersIntegrate@elementMetadata)
rownames(log2FCmatrix) <- Gencoordinates$name

 
#to group per expression inside of the features : https://jokergoo.github.io/ComplexHeatmap/articles/most_probably_asked_questions.html
 
for (i in 1:as.numeric(dim(log2FCmatrix)[2]) ){
	if (i == 1){
		selGenesF <- as.data.frame(log2FCmatrix[which(log2FCmatrix[,i]<=0.05 & abs(log2FCmatrix[,i]) >=1),]);
		selGenesF$Gene <- rownames(selGenesF);
		selGenesF <- selGenesF[order(selGenesF[,i]),];

	} else{
		selGenes.tmp <- as.data.frame(log2FCmatrix[which(log2FCmatrix[,i]<=0.05 & abs(log2FCmatrix[,i]) >=1),]);
		selGenes.tmp$Gene <- rownames(selGenes.tmp);
		selGenes.tmp <- selGenes.tmp[order(selGenes.tmp[,i]),];
		selGenesF <- unique(bind_rows(selGenesF,selGenes.tmp));
	};
};

#1) matrix for heatmap using all the genes selected based to log2fc and FDR
selGenesF.m <- as.matrix(selGenesF[,-12])
rownames(selGenesF.m) <- selGenesF$Gene
t_selGenesF.m <- t(selGenesF.m)
hc <- hclust(dist(t_selGenesF.m))
group = cutree(hc, k = dim(t_selGenesF.m)[1])

#identifying the top 50 genes per cluster
selGenesF_changed <- gather(selGenesF,key="cluster",Log2FC,1:11)
selGenesF_changed$Gene <-selGenesF$Gene
top50_per_cluster <- unique(as.data.frame(selGenesF_changed %>%
        group_by(cluster) %>%
        top_n(n = 50,
              wt = Log2FC))[,'Gene'])

#2) matrix for heatmap using only the top 50  genes per cluster of those selected based to log2fc and FDR

selGenesF_top50.m <- selGenesF.m[which(rownames(selGenesF.m) %in% top50_per_cluster),]
t_selGenesF_top50.m <- t(selGenesF_top50.m)
hc_top50 <- hclust(dist(t_selGenesF_top50.m))
group_top50 = cutree(hc_top50, k = dim(t_selGenesF_top50.m)[1])

#heatmap right now#
mycol3 <-sequential_hcl(7, palette = "Viridis");  
#horizonExtra PALLETE FROM ARCHr 9-colors as defined in their MAN
horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC",
	"3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A")




pdf(paste0(wd, f_results,f_plots,f_heatmaps,"groupedPaths_heatmap_integrated_vtop50.pdf"), width =8.6, height =6.2 )
nameExp <- "clusters_vtop50";
f_heatmaps <- "heatmap/"
space <-"               "
rownames(t_selGenesF_top50.m) <- gsub("_"," ",rownames(t_selGenesF_top50.m))

heatmapTonycols <- gplots::heatmap.2(t_selGenesF_top50.m, symm=F,symkey=T,symbreaks=T, scale="column",col = horizonExtra,dendrogram="row", Colv = F, Rowv=T,density.info="none",labCol = FALSE,
 lhei = c(2,9), keysize =1.5,trace="none", margin=c(10, 14),cexCol=1.6, cexRow=1.2, cex.main=0.5, main=paste0( space,"cell Types markers"),);
heatmapDef <- gplots::heatmap.2(t_selGenesF_top50.m, symm=F,symkey=T,symbreaks=T, scale="column",col = mycol3,dendrogram="row", Colv = F, Rowv=T,density.info="none",labCol = FALSE,
 lhei = c(2,9), keysize =1.5,trace="none", margin=c(10, 14),cexCol=1.6, cexRow=1.2, cex.main=0.5, main=paste0( space,"cell Types markers"),);
plot(heatmapTonycols)
plot(heatmapDef)
dev.off();

