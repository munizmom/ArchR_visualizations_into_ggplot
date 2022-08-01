
# proportion cells per sample atac
atacArch$propCells_per_cluster_per_sample <-paste0(atacArch$condition, ";",atacArch$Clusters_number ,"; ",atacArch$orig.ident);

cells <- as.data.frame(table(atacArch$propCells_per_cluster_per_sample));
colnames(cells)[2]<- "nbCells"
cells_plot<- separate(cells, Var1, c("condition","CellCluster","SampleID"  ),sep=";");
cells_plot$SampleID <- gsub(" ","",cells_plot$SampleID);
cells_plot[is.na(cells_plot$nbCells),'nbCells'] <-0

totals <- as.data.frame(cells_plot %>% dplyr::group_by(LA,SampleID) %>% dplyr::summarise(totalPerLASample = sum(nbCells)))

cells_plot <- left_join(cells_plot,totals,by=c("condition","SampleID"))
cells_plot <- cells_plot[,c(1,3,2,4,5)];
cells_plot<- cells_plot[order(cells_plot$LA, cells_plot$SampleID, cells_plot$CellCluster),]
cells_plot$percentage <- round((cells_plot$nbCells/cells_plot$totalPerLASample)*100,digits=1)

cells_table<- separate(cells, Var1, c("LA_Cellcluster","SampleID"),sep="; ");
cells_table<-  spread(cells_table, "SampleID","nbCells" );
cells_table_prop <- cells_table
for (i in 1:ncol(cells_table_prop)){
	cells_table_prop[is.na(cells_table_prop[,i]),i] <-0
	i <- i+1
}


res <- list(cellsProportions=cells_plot, cellNb=cells_table,peaksNbs=table)
openxlsx::write.xlsx(res, file = paste0(wd,f_results,f_supFig1, "cells_assignation.xlsx"));



#cells_plot[130,]<-c("AF","0220","Astrocytes_12",0,1116 ,0 )
cells_plot.df <- left_join(cells_plot,ClusterColors, by=c("CellCluster") )
cells_plot.df <- cells_plot.df[order(cells_plot.df$condition, cells_plot.df$SampleID,cells_plot.df$labels),]
cells_plot.df <- cells_plot.df[order(cells_plot.df$condition, cells_plot.df$SampleID,cells_plot.df$labels),]

rownames(cells_plot.df) <- 1:nrow(cells_plot.df)
dfF <-cells_plot.df


sampleID <- unique(data.frame(SampleID =dfF$SampleID,LA=dfF$condition))
sampleID$SampleLabel <- 1
sampleID[which(sampleID$condition=="cond1"),"SampleLabel"] <- 1:length(sampleID[which(sampleID$condition=="cond1"),"SampleLabel"] )
sampleID[which(sampleID$condition=="cond2"),"SampleLabel"] <- 1:length(sampleID[which(sampleID$condition=="cond2"),"SampleLabel"] )


dfF <- left_join(dfF,sampleID, by=c("SampleID","condition") )
dfF$SampleLabel <- factor(dfF$SampleLabel, levels = unique(dfF$SampleLabel))
dfF$CellCluster  <- gsub("_"," ", dfF$CellCluster )
dfF$CellCluster<- factor(dfF$CellCluster, levels = unique(dfF$CellCluster))


##################################################################
##################### plot stacked histograme for each sample all together ########
##################################################################

ClusterColors <- data.frame(CellCluster=unique(cells_plot[,'CellCluster']))

ClusterColors$Col <- ifelse(ClusterColors[,1]=="Astrocytes_4", "#140bb8",
ifelse(ClusterColors[,1]=="Astrocytes_12", "#47a4d6",
	ifelse(ClusterColors[,1]=="VLMC_21", "#cb09ed",
		ifelse(ClusterColors[,1]=="OPC_8", "#FFEC8B",
			ifelse(ClusterColors[,1]=="Oligodendrocytes_0", "#66CDAA",
				ifelse(ClusterColors[,1]=="Microglia_7", "#ed096c",
					ifelse(ClusterColors[,1]=="Inhibitory-neurons_11", "#ca8fe3",
						ifelse(ClusterColors[,1]=="Excitatory-neurons_1", "#FF0000",
						ifelse(ClusterColors[,1]=="Excitatory-neurons_5", "#d96868",
						ifelse(ClusterColors[,1]=="Excitatory-neurons_14", "#bf3232",
						ifelse(ClusterColors[,1]=="Excitatory-neurons_17", "#F08080", NA)))))))))))

unique(ClusterColors$Col) #ook no NAs
ClusterColors$labels <- paste0(gsub("_.*","",ClusterColors$CellCluster)," ", c(1,1,3,4,2,1,1,1,1,1,2))
ClusterColors <- ClusterColors[order(ClusterColors$labels),]



postscript(file=paste0(wd,f_results,f_supFig1,"_cellsPer_cluster_per_sample_perCondition.ps"),width=5, height=4);
print(p_stackHisto);
dev.off();

#print(p)

pdf(file=paste0(wd,f_results,f_supFig1,"_cellsPer_cluster_per_sample_perCondition.pdf"),width=4, height=3);
print(p_stackHisto);
dev.off();

jpeg(file=paste0(wd,f_results,f_supFig1,"_cellsPer_cluster_per_sample_perCondition.jpeg"),units = "in",width=6.5, height=4,res = 600);
print(p_stackHisto);
dev.off();



stack_histo_clusterSample.function <- function(input,nameOutput,ClusterColors){
	
    input<- input[gtools::mixedorder(as.character(input$SampleID)),]
	input$SampleID <- factor(input$SampleID, levels = unique(input$SampleID))
	input$CellCluster<- factor(input$CellCluster, levels = unique(input$CellCluster))

	p_stackHisto <- ggplot(input, aes(x = SampleID, y= Percentage , fill = CellCluster)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(labels= unique(input$CellCluster), values = unique(input$Col) ) +
    theme_classic()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 13, face="bold"), 
	 legend.title = element_text(size=12, face="bold"), legend.text=element_text(size=9, face="bold"),
	 axis.text=element_text(size=9, face="bold"), 
	 legend.key=element_rect(fill=NA),axis.title=element_text(size=12,face="bold"),
	 legend.key.size = unit(1.2,"line")) + labs(fill = "peak genomic region");
	
	p_stackHisto_percond <- ggplot(input, aes(x = SampleLabel, y= percentage , fill = CellCluster))+
	  facet_grid(. ~ condition, switch = "x") + geom_bar(stat="identity") + 
	scale_fill_manual(labels= levels(input$CellCluster), values = unique(input$Col) ) +
	theme_classic()+ ylab("Percentage of cells (%)")+xlab("") + theme(axis.text.x = element_blank(),
	        axis.ticks.x = element_blank(), plot.title = element_text(size = 13, face="bold"), 
	 legend.title = element_text(size=12, face="bold"), legend.text=element_text(size=9, face="bold"),
	 axis.text.y=element_text(size=9, face="bold"),  legend.key=element_rect(fill=NA),axis.title=element_text(size=12,face="bold"),
	 legend.key.size = unit(1.2,"line")) + labs(fill = "Cell cluster");

	postscript(file=paste0(wd,f_results,f_plots,nameOutput,"_peaks_binding_Regions_per_SampleID.ps"),width=5, height=4);
	print(p_stackHisto);
	dev.off();


	pdf(file=paste0(wd,f_results,f_plots,nameOutput,"_peaks_binding_Regions_per_SampleID.pdf"),width=6.5, height=3);
	print(p_stackHisto);
	dev.off();

	jpeg(file=paste0(wd,f_results,f_plots,nameOutput,"_peaks_binding_Regions_per_SampleID.jpeg"),units = "in",width=6.5, height=4,res = 600);
	print(p_stackHisto);
	dev.off();

	#per condition facet grid
	postscript(file=paste0(wd,f_results,f_plots,nameOutput,"_cellsPer_cluster_per_sample_perCondition.ps"),width=5, height=4);
	print(p_stackHisto_percond);
	dev.off();

	pdf(file=paste0(wd,f_results,f_plots,nameOutput,"_cellsPer_cluster_per_sample_perCondition.pdf"),width=4, height=3);
	print(p_stackHisto_percond);
	dev.off();

	jpeg(file=paste0(wd,f_results,f_plots,nameOutput,"_cellsPer_cluster_per_sample_perCondition.jpeg"),units = "in",width=6.5, height=4,res = 600);
	print(p_stackHisto_percond);
	dev.off();


	assign(paste0("p_stackHisto_",nameOutput),input,.GlobalEnv)
};

stack_histo_clusterSampleion.function(dfF,"rnaseq",ClusterColors)


