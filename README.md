# scRNA_wokflow
Personal notes and scripts for single cell analysis 


<img width="994" height="605" alt="chord" src="https://github.com/user-attachments/assets/c87d02cc-d94e-40c9-94ba-482d1bb49e53" />









#### The main workflow is : 
##### 1- QC  : checking MT percentage, batch effect ?,filtering base on number of transcripts (nFeature_RNA)
###### hint : you can find better thresholds by looking at the vln plot of "nFeature_RNA" , "nCount_RNA","percent.mt"


##### 2- stat : normalization: scale, dimensionality reduction (pca and umap)
###### hint : pca is linear but umap is better for visualization 


##### 3- we get to ELbowplot : choose the number of pcs based on the plot for down stream analysis. use FindNeighbors() and FindClusters() and RunUmap() for the clusters umap by DimPlot()

##### 4- cell type anotation : mainly conducted in 2 ways: 1- manual 2- by SingleR 


###### this workflow is manual , by checking significant cell markers.


###### 4-1 : choosing markers : we inteparate FeaturePlots for different cell markers while we have an eye on our initial umap plot we manully align brighter parts of a feature plot to a cluster in our umap. we assing our seurat object levels to the cell tyes that we found.

##### 5- want DEGs ? FindAllMarkers(pbmc,) then go for ORA or GSEA

##### 6- cell cell comunication :/cellchat GetAssayData()/ choose species / subset, identifyoverexpression/compute common probe/aggregate/ different plots


##### 7- trajectory : slingshot/tradeSeq/SingleCellExprement/
