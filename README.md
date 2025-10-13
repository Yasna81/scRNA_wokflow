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




#### Resulting Plots :



#### umap :

<img width="341" height="297" alt="umap" src="https://github.com/user-attachments/assets/cb93ed1e-68e4-477b-93cb-2bc1e0904004" />


#### marker annotation :

<img width="382" height="297" alt="marker_anot_2" src="https://github.com/user-attachments/assets/918c24dc-fcd0-4bcb-b845-b863ffceb42a" />

<img width="382" height="297" alt="marker_anot_3" src="https://github.com/user-attachments/assets/84a8ceee-4d17-4fa4-ad8b-7b1ab719a9ab" />






#### network of cells :


<img width="404" height="305" alt="networks" src="https://github.com/user-attachments/assets/a51350ea-88c0-4fb5-b0a7-ece8e4a598b3" />



#### significant genes in trajectory analysis :


<img width="462" height="382" alt="trajectory_pvalue" src="https://github.com/user-attachments/assets/4370f70d-19ac-475b-b14a-0377989b2631" />




