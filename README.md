# scRNA_wokflow
Personal notes and scripts for single cell analysis 

<img width="1470" height="1667" alt="image" src="https://github.com/user-attachments/assets/172c764d-8ac3-435a-b63e-89a0965208c4" />





### keep these in mind :


#### creating a suerat object by CreateSeuratObject(),loading data by Read10x()

the data structure in a single cell analysis is as follows >>



<img width="859" height="373" alt="image" src="https://github.com/user-attachments/assets/99b697e1-030a-41e3-986b-f685ed8c3945" />


When we want to add new columns to the metadata, we use [[]].like as 



```
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```


#### want to normalize ? NormalizeData()


#### Identification of highly variable features (feature selection) >> FindVariableFeatures ()


#### scaling ? ScaleData()


#### dimentional reduction :  > Runpca()- its linear-


#### clustering : FindNeighbors(),FindClusters(),

#### RunUMAP/RunTSNE : reducing dimentionality but non-linear

