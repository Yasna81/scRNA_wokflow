install.packages("BiocManager")
library(BiocManager)
BiocManager::install("Seurat")
library(Seurat)
#-----------------------------------------------
if (!requireNamespace("remotes",quietly = TRUE))
    install.packages("remotes")
library(remotes)
remotes::install_github("satijalab/seurat-data")
library(SeuratData)
#InstallData("pbmc3k")
#options(timeout = 600)
pbmc.data <- Read10X(data.dir = "~/practice/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(count = pbmc.data , project ="PBMC3K",min.cells = 3 ,min.features = 200)
#fixing dashes for mit genes :
rownames(pbmc) <- gsub("_","-",rownames(pbmc))
#QC :
#first mitochondria 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^MT-")
VlnPlot(pbmc,features = c ("nFeature_RNA","nCount_RNA","percent.mt"))
# a basic filtering :
pbmc <- subset(pbmc , subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#norm :
pbmc <- NormalizeData(pbmc)
#featureselection
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc)
# 15 pc sounds good for downstream analysis
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc , dim = 1:15 )
DimPlot(pbmc , reduction = "umap")
#finding markers :
library(dplyr)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE , min.pct = 0.25 ,logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)
Idents(pbmc) <- pbmc$seurat_clusters

#library(SingleR)
#library(SingleCellExperiment)
library(dplyr)

DimPlot(pbmc, reduction = "umap", label = TRUE,pt.size = 0.5)
FeaturePlot(pbmc,features = c("CD3D","CD8A","MS4A1","MKG7","CD14","FCER1A","PPBP"), reduction = "umap",pt.size = 0.5)
#
FeaturePlot(pbmc,features = c("IL7R","CCR7","S100A4","MS4A1","CD79A","NKG7","GNLY"), reduction = "umap",pt.size = 0.5)
FeaturePlot(pbmc,features = c("GZMB","NKG7","MKI67","HLA-DRA","PDCD1","CD4"), reduction = "umap",pt.size = 0.5)
new.cluster.ids <- c(
    "Memory CD4 T",#0
    "Naive CD4 T",#1
    "CD14+ Monocytes",#2
    "B Cell",#3
    "NK Cell",#4
    "Activated Memory CD4 T",#5
    "CD8 T",#6
    "DC",#7
    "Platalets" #8
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction =  "umap", label = TRUE , pt.size = 0.5,repel = TRUE) + ggplot2::ggtitle("PBMC anotation")
#GSEA :
#install.packages("clusterProfiler")
library(clusterProfiler)
#install.packages("org.Hs.eg.db")
library(org.Hs.eg.db)
library(enrichplot)
# lets look for DEGs 
library(Seurat)
#Idents(pbmc) <- pbmc$seurat_clusters
markers <- FindAllMarkers(pbmc ,only.ops = TRUE , min.pct = 0.25,logfc.threshold = 0.25)
###lets check the pathways for activated cd4 t : 
cluster_name <- "Activated Memory CD4 T"
cluster_genes <- subset(markers, cluster == cluster_name)
gene_ranks <- cluster_genes %>%
    arrange(desc(avg_log2FC))
gene_ranks$gene_sym <- rownames(gene_ranks)
gene_names <- bitr(gene_ranks$gene,
                   fromType =  "SYMBOL",
                   toType = "ENTREZID",
                   org.Hs.eg.db)
final_gene <- merge(gene_ranks,gene_names, by.x = "gene",by.y = "SYMBOL")
#oRA 
library(dplyr)
gene_list<- final_gene %>%
    arrange(desc(avg_log2FC))

sig_genes <- final_gene$ENTREZID[abs(final_gene$avg_log2FC) > 1]
ego <- enrichGO(
    gene = sig_genes,
    OrgDb = org.Hs.eg.db ,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05 , 
    qvalueCutoff = 0.2
)
dotplot(ego, showCategory = 20)

#### B cells :
cluster_name <- "B Cell"
cluster_genes <- subset(markers, cluster == cluster_name)
gene_ranks <- cluster_genes %>%
    arrange(desc(avg_log2FC))
gene_ranks$gene_sym <- rownames(gene_ranks)
gene_names <- bitr(gene_ranks$gene,
                   fromType =  "SYMBOL",
                   toType = "ENTREZID",
                   org.Hs.eg.db)
final_gene <- merge(gene_ranks,gene_names, by.x = "gene",by.y = "SYMBOL")
#oRA 
library(dplyr)
gene_list<- final_gene %>%
    arrange(desc(avg_log2FC))

sig_genes <- final_gene$ENTREZID[abs(final_gene$avg_log2FC) > 1]
ego <- enrichGO(
    gene = sig_genes,
    OrgDb = org.Hs.eg.db ,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05 , 
    qvalueCutoff = 0.2
)
dotplot(ego, showCategory = 10)
#computing cell cell comunication :
#BiocManager::install("CellChat")
library(remotes)
remotes::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
data.input <- GetAssayData(pbmc, assay = "RNA",slot = "data")
labels <- Idents(pbmc)
meta <- data.frame(group = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input,meta= meta ,group.by = "group")
#choose species :
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
#subset to relvent signaling :
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#overall network strength
netVisual_circle(cellchat@net$count,vertex.weight = TRUE, weight.scale = TRUE)
#pathway specific network : which clusters communicate via the MHC-I pathway 
netVisual_aggregate(cellchat,signaling = "MHC-I",layout = "circle")
#bubble plot : top ligand receptor pairs
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(3,4))
#which ligands and receptors contribute to a given pathway
netVisual_chord_gene(cellchat,signaling = "MHC-I")
#the most influential ligand-receptors
netAnalysis_contribution(cellchat,signaling = "MHC-I")
####----------------------------------------------------------------------------------
#trajectory_analysis :
#remotes::install_github("kstreet13/slingshot")
#remotes::install_github("statOmics/tradeSeq")
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(Seurat)
sce <- as.SingleCellExperiment(pbmc)
sce$clusters <- Idents(pbmc)
sce <- slingshot(
    sce,
    clusterLabels = "clusters",
    reducedDim = "UMAP"
)
summary(sce$slingPseudotime_1)
library(RColorBrewer)
colors <- brewer.pal(length(unique(sce$seurat_clusters)), "Set1")
plot(reducedDims(sce)$UMAP, col = sce$clusters,pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
library(matrixStats)
# genes changing along pseudotime :
top_genes <- head(order(rowVars(log1p(counts(sce))),decreasing = TRUE),200)

#fitting GAM 
sce_subset <- sce[top_genes,]
set.seed(1)
sce_subset <- fitGAM(sce_subset,nknots = 4 , parallel = FALSE)
plotSmoothers(sce_subset,assays(sce_subset)$counts, gene = "CD3D")
plotSmoothers(sce_subset,assays(sce_subset)$counts, gene = "CD3D", alpha = 0.6)
#association test :
assoRes <- associationTest(sce_subset)
head(assoRes)
#lets get to the significant onse: gene changes significantly along pseudo time in at least one linage
topGenes <- rownames(assoRes)[assoRes$pvalue <0.05]
head(top_genes)
sig <- assoRes %>% subset(assoRes$pvalue != 0.000000e+00)
# which genes differ between linages (genes with different expression at the endpoint of linage ):
diffEndRes <- diffEndTest(sce_subset)
head(diffEndRes)
sig_end <- diffEndRes %>% subset(diffEndRes$pvalue <=0.05 & diffEndRes$pvalue != 0.000000e+00)
top5_lowest_p <- sig_end %>%
    arrange(pvalue) %>%
    head(5)
ggplot(top5_lowest_p, aes(x=rownames(top5_lowest_p), y = pvalue))+ 
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    labs(title = "Top 5 Most Significant Variables",
         x = "Variables",
         y = "P-value") +
    theme_minimal() +


#LETS CHECK FR A SIGNIFICANT GENE "ACTB"
plotSmoothers(sce_subset,assays(sce_subset)$counts, gene = "ACTB")
#
lineages <- slingLineages(sce)
# we want to get dominant cells in each linage
sce_1 <- slingshot(sce,clusterLabels = sce$clusters,reducedDim = "UMAP")
pseudo <- as.data.frame(slingPseudotime(sce_1))
pseudo$cell <- rownames(pseudo)
meta <- data.frame(cell = rownames(colData(sce_1)),
                   celltype = sce_1$clusters
                   )
df <- merge(meta,pseudo, by = "cell")
head(df)
#identify lineage :
pseudo$lineage <- apply(pseudo[,grep("Lineage",colnames(pseudo))],1,function(x) {
    if (all(is.na(x))) return(NA)
    colnames(pseudo)[grep("Lineage",colnames(pseudo))][which.min(x)]
})
head(pseudo)

merged <- merge(meta,pseudo[,c("cell","lineage")], by = "cell")
#summary :
library(dplyr)
dominant <- merged %>% 
    group_by(lineage,celltype) %>%
    summarise(count = n()) %>%
    mutate(precent = count / sum(count * 100)) %>%
    arrange(lineage, desc(precent))
