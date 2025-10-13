# scRNA_wokflow
Personal notes and scripts for single cell analysis 


<img width="994" height="605" alt="chord" src="https://github.com/user-attachments/assets/c87d02cc-d94e-40c9-94ba-482d1bb49e53" />




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

