#########
#########
#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

#########
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R


library(Seurat)
library(Matrix)

######


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
all_seurat_RNA_big_cl_detail_new = readRDS('all_seurat_RNA_big_cl_detail_new_202204')

table(all_seurat_RNA_big_cl_detail_new$new_celltypes)

###### separate Cone cells to M and S ####
######

LGS13_RNA_Cone = all_seurat_RNA_big_cl_detail_new[,which(all_seurat_RNA_big_cl_detail_new$new_celltypes=="Cone")]

###### redo UMAPs ####
LGS13_RNA_Cone <- FindVariableFeatures(LGS13_RNA_Cone)
LGS13_RNA_Cone <- RunPCA(LGS13_RNA_Cone)

LGS13_RNA_Cone <- RunUMAP(
  LGS13_RNA_Cone,
  reduction = "pca",       # which reduction to use
  dims = 1:50,             # which PCs to use
  reduction.name = "umap", # name of the slot (default “umap”)
  reduction.key  = "UMAP_" # prefix for cell embeddings (default “UMAP_”)
)

LGS13_RNA_Cone <- FindNeighbors(LGS13_RNA_Cone,reduction = "pca")
LGS13_RNA_Cone <- FindClusters(LGS13_RNA_Cone,resolution=0.8)

Idents(LGS13_RNA_Cone) <- LGS13_RNA_Cone$Cone_type

rownames(LGS13_RNA_Cone)[grep("Arr3",rownames(LGS13_RNA_Cone))]

png(
  filename = "umap_plot.png",
  width    = 1800,
  height   = 1600,
  res      = 150  # resolution in DPI
)

# Draw the UMAP
print(FeaturePlot(
  LGS13_RNA_Cone,
  reduction = "umap",
  features = c("Opn1sw","OPN1MW","Ccdc136","Sag","Isl2","Arr3")
))

# Close the device
dev.off()

####
png(
  filename = "umap_plot1.png",
  width    = 700,
  height   = 600,
  res      = 150  # resolution in DPI
)

# Draw the UMAP
print(DimPlot(
  LGS13_RNA_Cone,
  reduction = "umap",
  label=T
))

# Close the device
dev.off()

######
######
######

LGS13_RNA_Cone$Cone_type = "M"

k = which(LGS13_RNA_Cone$seurat_clusters == "7")

LGS13_RNA_Cone$Cone_type[k] = "S"


####### see DEGs #####

Idents(LGS13_RNA_Cone) = LGS13_RNA_Cone$Cone_type

DefaultAssay(LGS13_RNA_Cone) <- "RNA"
LGS13_RNA_Cone <- NormalizeData(LGS13_RNA_Cone)

Res = FindMarkers(LGS13_RNA_Cone,ident.1="S",ident.2="M",only.pos=F)


head(Res)

Res[which(rownames(Res) == "OPN1MW"),]

library(writexl)

Res = cbind(data.frame(Genes=rownames(Res)),Res)
Res = Res[which(Res$p_val_adj < 0.05),]
Res = Res[order(Res$avg_log2FC,decreasing=T),]

library(writexl)
write_xlsx(Res,path = "LGS13_SvsM_Cone_DEGs .xlsx")

getwd()

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")

save(LGS13_RNA_Cone,file="LGS13_RNA_Cone_2025")

######
load("LGS13_RNA_Cone_2025")


######
table(LGS13_RNA_Cone$timepoint)

######
LGS13_RNA_Cone_P21 = LGS13_RNA_Cone[,which(LGS13_RNA_Cone$timepoint == "P21")]

######
DefaultAssay(LGS13_RNA_Cone_P21) = "RNA"
LGS13_RNA_Cone_P21 <- NormalizeData(LGS13_RNA_Cone_P21)
LGS13_RNA_Cone_P21 <- FindVariableFeatures(LGS13_RNA_Cone_P21)
LGS13_RNA_Cone_P21 <- ScaleData(LGS13_RNA_Cone_P21)
LGS13_RNA_Cone_P21 <- RunPCA(LGS13_RNA_Cone_P21)

LGS13_RNA_Cone_P21 <- RunUMAP(
  LGS13_RNA_Cone_P21,
  reduction = "pca",       # which reduction to use
  dims = 1:50,             # which PCs to use
  reduction.name = "umap", # name of the slot (default “umap”)
  reduction.key  = "UMAP_" # prefix for cell embeddings (default “UMAP_”)
)

#####
#####

png(
  filename = "umap_plotP21.png",
  width    = 1800,
  height   = 1600,
  res      = 150  # resolution in DPI
)

# Draw the UMAP
print(FeaturePlot(
  LGS13_RNA_Cone_P21,
  reduction = "umap",
  features = c("Opn1sw","OPN1MW","Ccdc136","Sag","Isl2","Arr3")
))

# Close the device
dev.off()

####
png(
  filename = "umap_plot1.png",
  width    = 700,
  height   = 600,
  res      = 150  # resolution in DPI
)

# Draw the UMAP
print(DimPlot(
  LGS13_RNA_Cone,
  reduction = "umap",
  label=T
))

# Close the device
dev.off()

#######----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


