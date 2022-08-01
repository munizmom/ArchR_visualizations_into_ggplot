# ArchR_visualizations_into_ggplot

snATACseq analyses performed using [ArchR](https://www.archrproject.com) are really easy to perform thanks to the big effort of Corces team, however in the visualization area if you dont have computational biology background can be quite challenging to try to modify the plots. 
This collection of scripts will aid to the community, providing an easy way to extract the raw data from archR snATACSeq objects to produce in house visualizations of the results using ggplot or ggplot2 that are widely used by the biological community.

## Visualizations provided:
- plotMarkerHeatmap of all the peaks/genes altered per cluster using ggplot:heatmap.2 function. 
- plotMarkerHeatmap of the top number of altered peaks/genes based in log2FC per cluster using ggplot:heatmap.2 function. The user cna define the top number of features or peaks they want to plot, currently set to the top 50.
- stacked barplot to show the percentage of cells per sample per cluster.

Last update 2022.
