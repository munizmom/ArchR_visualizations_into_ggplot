library(colorspace); library(ArchR);library(ggplot2);library(dplyr);

#########################################################################################################
# How to extract from the output of getMarkerFeatures in ArchR 
#  (stored in the R object called markersPeaks) the info needed to plot the heatmap using ggplot2
# instead of using plotMarkerHeatmap from ArchR
#########################################################################################################
markersPeaks <- getMarkerFeatures(
    ArchRProj = atacArch, 
    useMatrix = "PeakMatrix", 
    groupBy = "clustera",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#selecting the 2 clusters we want to plot from all the clusters
markersPeaks4astro <- markersPeaks[,c("AF_1","AF_2")]

#Storing the data needed for the creation fo the plot in a dataframe 
log2Df <- as.data.frame(markersPeaks4astro@assays@data$Log2FC)
colnames(log2Df)<-  markersPeaks4astro@colData@rownames

#extracting genomic location info, is going to be the rownames and labels of the plot
Gencoordinates<- as.data.frame(markersPeaks4astro@elementMetadata)
rownames(log2Df) <- paste0(gsub("chr", "", Gencoordinates$seqnames), ": ",
  Gencoordinates$start, "-",Gencoordinates$end);

#selecting chr19 specific region noted as "LA"
log2Df_chr19 <- log2Df[grep("^19: ",rownames(log2Df)),];
#max/min value for the colour scale 
max_Chr19 <- max(log2Df_chr19,abs(min(log2Df_chr19)));


# Trying to annotate beginning end of LA in heatmap
##########################################################
log2Df_chr19$nrow <- 1:length(log2Df_chr19[,1]);
log2Df_chr19$start <- as.numeric(gsub("-.*","",gsub("^19: ","",rownames(log2Df_chr19))));
log2Df_chr19$end <- as.numeric(gsub(".*-","",gsub("^19: ","",rownames(log2Df_chr19))));
# Defining beginning and end positions to add a vertical line in the heatmap
LA_init <- log2Df_chr19[which(log2Df_chr19$start>44000000),"nrow"][1];
LA_end <- tail(log2Df_chr19[which(log2Df_chr19$end<46000000),"nrow"])[6];

#plotting only the peaks in chr19
log2Df_chr19$Label <- rownames(log2Df_chr19);
log2Df_chr19 <- gather(log2Df_chr19, "LA",log2FC,1:2);


# Plotting a heatmap using ggplot 
###########################################################

#colour palette for the heatmap
###############################
mycol2 <-diverging_hcl(length(mycol),
 palette = "Tropic");  
################################

pheatmap <- ggplot(log2Df_chr19, aes(as.factor(Label), LA, fill= log2FC)) + 
  geom_tile() + scale_fill_gradientn(colours = c(mycol2)) + labs(x= "chr19 peaks") +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size=2))+
  geom_vline(xintercept = LA_init, linetype="dotdash") + geom_vline(xintercept = LA_end, linetype="dotdash");

#remove x axis labels
pheatmap2 <- pheatmap + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank());

pdf(file=paste0(wd,f_results,f_plots,"heatmap_Chr19Peaks_v_ggplot.pdf"),onefile=TRUE,width = 250, height =4)
print(pheatmap);
dev.off();
