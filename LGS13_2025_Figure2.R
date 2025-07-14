##########
########## This file is the code for the figures in 13LGS ########
##########
########## first we will load the mouse seurat files ########
########## we want the gene expression in all the mouse cells #####
##########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")

#### saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_cl_Nov18.rds")
#### saveRDS(Mouse_RNA_seurat_merge_cl,file="Mouse_RNA_seurat_merge_cl_Nov18.rds")
#### Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_cl_Nov18.rds")

Mouse_RNA_seurat_merge_cl <- readRDS("Mouse_RNA_seurat_merge_cl_Nov18.rds")

library(ArchR)
library(Seurat)

head(Mouse_RNA_seurat_merge_cl@meta.data)

table(Mouse_RNA_seurat_merge_cl$time)
table(Mouse_RNA_seurat_merge_cl$celltype)

###### AC/HC     BC   Cone Cone_p    E_N    L_N     MG    RGC    Rod  Rod_p RPC_S1 
######  8170   2276   2477   1022   3081   7861   1391   6561   7743   7797   7853 
###### RPC_S2 RPC_S3 
###### 16646  22969 

###### Next we will load the total LGS Gene expression values ###########
###### for the UMAP in the figure1 ###########
######


##### see the scRNA-seq project #####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
all_seurat_RNA_big_cl_detail_new = readRDS('all_seurat_RNA_big_cl_detail_new_202204')

head(all_seurat_RNA_big_cl_detail_new@meta.data)
table(all_seurat_RNA_big_cl_detail_new$new_celltypes)



#####
##### convert to gene names  #############
#####

LGS to before
Mouse to before

#####
##### 
#####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

library(Seurat)
library(Matrix)


#####---------------------------------------------------------------------------------------------------------------------------------------
##### mat = Matrix2 ####
Remove_dup_rownames <- function(mat){
    ##########
    k = which(duplicated(rownames(mat)) == T)
    dup_G = rownames(mat)[k]
    ##########
    mat_1 = mat[which(rownames(mat) %in% dup_G == F),]
    ########## which(duplicated(rownames(mat_1)) == T)
    mat_2 = mat[which(rownames(mat) %in% dup_G == T),]
    ##########
    dup_G = dup_G[!duplicated(dup_G)]
    ##########
    mat_2_new = list()
    for(i in 1:length(dup_G)){
        print(i)
        ######
        tmp_mat = mat_2[which(rownames(mat_2) == dup_G[i]),]
        ######
        tmp_mat_sum = colSums(tmp_mat)
        ######
        mat_2_new = c(mat_2_new,tmp_mat_sum)
    }
    ########
    mat_2_new = as.numeric(mat_2_new)
    mat_2_new2 = matrix(mat_2_new,ncol=length(dup_G))
    ########
    colnames(mat_2_new2) <- dup_G
    rownames(mat_2_new2) <- colnames(mat)
    mat_2_new2 = t(mat_2_new2)
    #######
    ####### head(mat_2_new)
    #######
    ###mat_2_new = as.matrix(mat_2_new)
    mat_2_new2 <- as(mat_2_new2, "dgCMatrix")
    #######
    mat_combined = rbind(mat_2_new2,mat_1)
    ######
    all.equal(colnames(mat_2_new2),colnames(mat_1))
    ######
    return(mat_combined)
}
#####---------------------------------------------------------------------------------------------------------------------------------------
#####----------- 这个function 读入原始的 Mouse 的 RNAseq 文件  ---------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development")
Before_Mouse_all_Gene_Names = read.table("10x_mouse_retina_development_feature.csv",sep=",",header=T)
Before_Mouse_all_Gene_Names$gene_short_name
Cells = read.table("10x_Mouse_retina_pData_umap2_CellType_annot_w_horiz.csv",sep=",",header=T)
#####
Matrix2 = readMM("10x_mouse_retina_development.mtx")
rownames(Matrix2) <- Before_Mouse_all_Gene_Names$gene_short_name
colnames(Matrix2) <- Cells$barcode
#####
dim(Matrix2)
k = which(duplicated(rownames(Matrix2)) == T)
##### rownames(Matrix2)[k]
Matrix_cl = Remove_dup_rownames(Matrix2)
#####
Raw_seurat = CreateSeuratObject(Matrix_cl)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development")
saveRDS(Raw_seurat,file="Mouse_RNA_dev_raw_seurat_Jan29_2025.rds")

#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
##### 这个function 读入我们自己的 Mouse 的 celltype annotation ###################
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
Mouse_RNA_seurat_merge_cl <- readRDS("Mouse_RNA_seurat_merge_cl_Nov18.rds")
Mouse_RNA_seurat_merge_cl$cell_id = colnames(Mouse_RNA_seurat_merge_cl)
Mouse_RNA_seurat_merge_cl$cell_id = gsub("#",".",Mouse_RNA_seurat_merge_cl$cell_id)


k = which(colnames(Raw_seurat) %in% Mouse_RNA_seurat_merge_cl$cell_id == T)
Raw_seurat_cl = Raw_seurat[,k]
m = match(colnames(Raw_seurat_cl),Mouse_RNA_seurat_merge_cl$cell_id)

table(as.character(Mouse_RNA_seurat_merge_cl$celltype))
length(which(is.na(Mouse_RNA_seurat_merge_cl$celltype) == T))

Raw_seurat_cl$celltype = as.character(Mouse_RNA_seurat_merge_cl$celltype[m])
Raw_seurat_cl$timepoint = as.character(Mouse_RNA_seurat_merge_cl$time[m])

table(as.character(Raw_seurat_cl$celltype))
table(as.character(Raw_seurat_cl$timepoint))
length(which(is.na(Raw_seurat_cl$celltype) == T))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development")
saveRDS(Raw_seurat_cl,file="Mouse_RNA_dev_raw_seurat_Anno_Jan29_2025.rds")

#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
##### 这个function 读入 Mouse 和 13LGS 所有对应的基因 ##########
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure2')
squirrel_to_mouse <- read.table('squirrel_to_mouse_gene_names.csv',sep=',',header=T)

squirrel_to_mouse$M_INDEX = 'NO'
squirrel_to_mouse$LGS_INDEX = 'NO'

##### load the Mouse seurat #####
Mouse_RNA_dev_raw_seurat_Anno <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development/Mouse_RNA_dev_raw_seurat_Anno_Jan29_2025.rds")
Mouse_total_genes = rownames(Mouse_RNA_dev_raw_seurat_Anno)
k1 = which(squirrel_to_mouse$Gene.name %in% Mouse_total_genes == T)
squirrel_to_mouse$M_INDEX[k1] = 'YES'

##### load the LGS seurat ########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
all_seurat_RNA_big_cl_detail_new = readRDS('all_seurat_RNA_big_cl_detail_new_202204')
LGS_total_genes = rownames(all_seurat_RNA_big_cl_detail_new)
k2 = which(squirrel_to_mouse$Squirrel.name %in% LGS_total_genes == T)
squirrel_to_mouse$LGS_INDEX[k2] = 'YES'

##### clean the squirrel_to_mouse ######
k3 = which(squirrel_to_mouse$M_INDEX == "YES" & squirrel_to_mouse$LGS_INDEX == "YES")
squirrel_to_mouse_cl = squirrel_to_mouse[k3,]

######## save the squirrel_to_mouse_cl ###########
########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(squirrel_to_mouse_cl,file="squirrel_to_mouse_cl_2025.rds")
head(squirrel_to_mouse_cl,n=1000)


#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
##### 这个function 读入 Mouse 和 13LGS 所有对应的基因 ##########
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
######## Next we will generate comparable Seurat files #######
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

##### for LGS13 ##########
LGS_count_mat = all_seurat_RNA_big_cl_detail_new[["RNA"]]@counts
m0 = which(rownames(LGS_count_mat) %in% squirrel_to_mouse_cl$Squirrel.name == T)
LGS_count_mat_cl = LGS_count_mat[m0,]
##### convert LGS13 genes to Mouse ########
LGS_count_mat_toM = LGS_count_mat_cl
m = match(rownames(LGS_count_mat_toM),squirrel_to_mouse_cl$Squirrel.name)
rownames(LGS_count_mat_toM) = squirrel_to_mouse_cl$Gene.name[m]
LGS_count_mat_toM_Reduce = Remove_dup_rownames(LGS_count_mat_toM)
##### 14011 !!!! #########
LGS_seurat_compare = CreateSeuratObject(LGS_count_mat_toM_Reduce)
##### add cell types and save #########
m = match(colnames(LGS_seurat_compare),colnames(all_seurat_RNA_big_cl_detail_new))
LGS_seurat_compare$celltype = all_seurat_RNA_big_cl_detail_new$new_celltypes[m]
LGS_seurat_compare$timepoint = all_seurat_RNA_big_cl_detail_new$timepoint[m]
table(LGS_seurat_compare$celltype)
#####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS_seurat_compare,file="LGS_seurat_compare_Jan29_2025.rds")
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

##### one to one list ####
LGStoM_datatable = data.frame(LGS_G =rownames(LGS_count_mat_cl) ,M_G=rownames(LGS_count_mat_toM))
saveRDS(LGStoM_datatable,file="LGStoM_datatable_2025.rds")
res = table(LGStoM_datatable[,2])
res = sort(res)
LGStoM_datatable[which(LGStoM_datatable$M_G == "Tmem132b"),]
LGStoM_datatable[which(LGStoM_datatable$M_G == "Yy1"),]

#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------


##### for Mouse ##########
Mouse_count_mat = Raw_seurat_cl[["RNA"]]@counts
m0 = which(rownames(Mouse_count_mat) %in% squirrel_to_mouse_cl$Gene.name == T)
Mouse_count_mat_cl = Mouse_count_mat[m0,]
##### 14011 !!! ##########
Mouse_seurat_compare = CreateSeuratObject(Mouse_count_mat_cl)
##### add celltype anno and save ########
m = match(colnames(Mouse_seurat_compare),colnames(Raw_seurat_cl))
Mouse_seurat_compare$celltype = Raw_seurat_cl$celltype[m]
Mouse_seurat_compare$timepoint = Raw_seurat_cl$timepoint[m]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(Mouse_seurat_compare,file="Mouse_seurat_compare_Jan29_2025.rds")


#####
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
##### Function Average Expression ######
#####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_seurat_compare <- readRDS(file="LGS_seurat_compare_Jan29_2025.rds")
Mouse_seurat_compare <- readRDS(file="Mouse_seurat_compare_Jan29_2025.rds")

#####
#####
Get_Avg_Exp <- function(Mat,Meta,Celltypes_Need){
    #########
    all_avg = c()
    for(i in 1:length(Celltypes_Need)){
        #########
        tmp_ct = Celltypes_Need[i]
        print(tmp_ct)
        k = which(Meta$index == tmp_ct)
        tmp_mat = Mat[,k]
        ##########
        tmp_sum = apply(tmp_mat,1,sum)
        ##########
        all_avg <- c(all_avg,tmp_sum)
    }
    ####
    all_avg_mat = matrix(all_avg,ncol=length(Celltypes_Need))
    rownames(all_avg_mat) = rownames(Mat)
    colnames(all_avg_mat) = Celltypes_Need
    ####
    factor = colSums(all_avg_mat) / 1e5
    all_avg_mat = sweep(all_avg_mat,2,factor,FUN="/")
    all_avg_mat_log = log2(all_avg_mat+1)
    ####
    return(all_avg_mat_log)
}



#####
##### 接下来计算一下每个 celltype 的基因表达水平 ##########
#####
##### Mouse_seurat_compare
##### LGS_seurat_compare
##### 
##### We need early RPC late RPC early NG late NG ########
#####

table(Mouse_seurat_compare$celltype)
Early_NG = "E_N"
Late_NG = "L_N"
Early_RPC = c("RPC_S1","RPC_S2")
Late_RPC = c("RPC_S3")
Photo_pre = c("Cone_p","Rod_p")

M_list = list(Early_NG,Late_NG,Early_RPC,Late_RPC,Photo_pre)
names(M_list) = c("Early_NG","Late_NG","Early_RPC","Late_RPC","Photo_pre")

k = which(Mouse_seurat_compare$celltype %in% unlist(M_list) == T)
Mouse_seurat_compare_Sub = Mouse_seurat_compare[,k]
table(Mouse_seurat_compare_Sub$celltype)
Mouse_seurat_compare_Sub$newcelltype = Mouse_seurat_compare_Sub$celltype
for(i in 1:length(M_list)){
    k = which(Mouse_seurat_compare_Sub$newcelltype %in% M_list[[i]] == T)
    Mouse_seurat_compare_Sub$newcelltype[k] = names(M_list)[i]
}
table(Mouse_seurat_compare_Sub$newcelltype)
Mouse_seurat_compare_Sub$index = Mouse_seurat_compare_Sub$newcelltype
Mouse_Avg = Get_Avg_Exp(Mouse_seurat_compare_Sub[['RNA']]@counts,Mouse_seurat_compare_Sub@meta.data,c('Early_RPC',"Late_RPC","Early_NG","Late_NG","Photo_pre"))


#####---------------------------------------------------------------------------------------------------------------------
#####----- Next for LGS13 Avg expression ----------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------

table(LGS_seurat_compare$celltype)
Early_NG = "Early_NG"
Late_NG = "Late_NG"
Early_RPC = c("Early_RPC")
Late_RPC = c("Late_RPC")
Photo_pre = c("Photoprecursor")
LGS_list = list(Early_NG,Late_NG,Early_RPC,Late_RPC,Photo_pre)
names(LGS_list) = c("Early_NG","Late_NG","Early_RPC","Late_RPC","Photo_pre")

k = which(LGS_seurat_compare$celltype %in% unlist(LGS_list) == T)
LGS_seurat_compare_Sub = LGS_seurat_compare[,k]
table(LGS_seurat_compare_Sub$celltype)
LGS_seurat_compare_Sub$newcelltype = LGS_seurat_compare_Sub$celltype
for(i in 1:length(LGS_list)){
    k = which(LGS_seurat_compare_Sub$newcelltype %in% LGS_list[[i]] == T)
    LGS_seurat_compare_Sub$newcelltype[k] = names(LGS_list)[i]
}
table(LGS_seurat_compare_Sub$newcelltype)
LGS_seurat_compare_Sub$index = LGS_seurat_compare_Sub$newcelltype

####
LGS_Avg = Get_Avg_Exp(LGS_seurat_compare_Sub[['RNA']]@counts,LGS_seurat_compare_Sub@meta.data,c('Early_RPC',"Late_RPC","Early_NG","Late_NG","Photo_pre"))
LGS_Avg['Onecut2',]

#####
#####---------------------------------------------------------------------------------------------------------------------
#####----- Next for LGS13 Figure2H ----------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####----- see the expression levels on these TFs #######
#####
Perform Qnorm in the LGS_Avg and Mouse_Avg ####
#####
colnames(LGS_Avg) = paste0("LGS13:",colnames(LGS_Avg))
colnames(Mouse_Avg) = paste0("Mouse:",colnames(Mouse_Avg))

all.equal(rownames(LGS_Avg),rownames(Mouse_Avg))
m = match(rownames(LGS_Avg),rownames(Mouse_Avg))
Mouse_Avg = Mouse_Avg[m,]
all.equal(rownames(LGS_Avg),rownames(Mouse_Avg))

LGS_M_Avg = cbind(Mouse_Avg,LGS_Avg)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS_M_Avg,file="LGS_M_Avg_2025.rds")

######
######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_M_Avg <- readRDS("LGS_M_Avg_2025.rds")

library(limma)
LGS_M_Avg_mat = as.matrix(LGS_M_Avg)
LGS_M_Avg_Norm <- normalizeBetweenArrays(LGS_M_Avg_mat, method = "quantile")

summary(LGS_M_Avg_Norm[,1])
summary(LGS_M_Avg_Norm[,3])

######
#####---------------------------------------------------------------------------------------------------------------------
#####----- Next read the enriched TFs result !!!!! ----------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
######

###### 

Perform_TF_enrich <- function(GRNs,Class="pos",Cone_Genes,TAG2="Cone"){
    ######
    all_TFs = unique(GRNs$TF)
    ######
    GRNs = data.frame(GRNs)
    GRNs = GRNs[which(GRNs$Class==Class),]
    Total_res = list()
    ######
    for(i in 1:length(all_TFs)){
        #######
        Tmp_TFs = all_TFs[i]
        #######
        total_background = dim(GRNs)[1]
        #######
        total_success = GRNs[which(GRNs$gene %in% Cone_Genes == T & GRNs$Class == Class),]
        total_success = dim(total_success)[1]
        #######
        total_TFs = GRNs[which(GRNs$TF %in% Tmp_TFs == T),]
        total_TFs = dim(total_TFs)[1]
        #######
        sample_success = GRNs[which(GRNs$TF %in% Tmp_TFs == T & GRNs$gene %in% Cone_Genes == T & GRNs$Class == Class),]
        sample_success = dim(sample_success)[1]
        ######
        Inset_Insample = sample_success
        Notset_Insample = total_TFs-sample_success
        Inset_Notsample = total_success-sample_success
        Notset_Nosample = total_background-total_success-sample_success+sample_success
        ######
        mytable <- matrix(
            c(Inset_Insample, Notset_Insample,  # DE行: In Set=0, Not in Set=k
            Inset_Notsample, Notset_Nosample),  # Not DE行: In Set=m, Not in Set=???
            nrow=2, byrow=TRUE
        )
        res_p = fisher.test(mytable, alternative = "greater")$p
        coverage = sample_success / total_TFs 
        ######
        res = data.frame(TF=Tmp_TFs,pvalue=res_p,coverage=coverage,Class=Class,Num_of_target=sample_success,GRN_CT=TAG2)
        Total_res <- c(Total_res,list(res))
    }
    ####
    Total_res_tab = do.call("rbind",Total_res)
    ####
    Total_res_tab = Total_res_tab[order(Total_res_tab$pvalue),]
    return(Total_res_tab)
}


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR2
R


###### for LGS Cone !!!!! ########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_seurat_Cone_markers <- readRDS("LGS13_seurat_Cone_markers_Jan2025")

rownames(LGS13_seurat_Cone_markers)[grep("-like",rownames(LGS13_seurat_Cone_markers))]
rownames(LGS13_seurat_Cone_markers)[grep("-containing",rownames(LGS13_seurat_Cone_markers))]

Cone_Genes = rownames(LGS13_seurat_Cone_markers)
Class = "pos"

GRNs = readRDS("LGS13_Cone_GRNs_Mar2025")[[3]]
LGS13_Cone_promote_Cone = Perform_TF_enrich(GRNs,Class="pos",Cone_Genes,TAG2="Cone")
LGS13_Cone_promote_Cone = LGS13_Cone_promote_Cone[order(LGS13_Cone_promote_Cone$pvalue),]

LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Zic3"),]
LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Pou2f1"),]
LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Thrb"),]
LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Rxrg"),]
LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Sall3"),]
LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Mef2c"),]

GRNs = readRDS("LGS13_NG_GRNs_Jan2025")[[3]]
LGS13_NG_promote_Cone = Perform_TF_enrich(GRNs,Class="pos",Cone_Genes,TAG2="NG")

###### Method: Cone specific genes #####################
###### Method: 1e-6 rank by number of targets ##########
GRNs = readRDS("LGS13_Photopre_GRNs_Mar2025")[[3]]
LGS13_Pre_promote_Cone = Perform_TF_enrich(GRNs,Class="pos",Cone_Genes,TAG2="PhotoPre")
LGS13_Pre_promote_Cone = LGS13_Pre_promote_Cone[order(LGS13_Pre_promote_Cone$Num_of_target,decreasing=T),]
LGS13_Pre_promote_Cone = LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$pvalue < 1e-6),]

####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Pre_promote_Cone,file="LGS13_TFs_Pre_promote_Cone_2025.rds")
####
#### write to table: ######
####

library(openxlsx)
write.xlsx(LGS13_Pre_promote_Cone, file = "TableS6_LGS13_Pre_promote_Cone.xlsx", sheetName = "TFs", rowNames = FALSE)

####
####

LGS13_Pre_promote_Cone[which(LGS13_Pre_promote_Cone$TF == "Mef2c"),]
LGS13_Pre_promote_Cone_Top50 = LGS13_Pre_promote_Cone[1:50,]
LGS13_Pre_promote_Cone_Top30 = LGS13_Pre_promote_Cone[1:30,]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Pre_promote_Cone_Top30,file="LGS13_Pre_promote_Cone_Top30.rds")
saveRDS(LGS13_Pre_promote_Cone_Top50,file="LGS13_Pre_promote_Cone_Top50.rds")

#######
#######
#######

Zic3 Pou2f1
Thrb Rxrg Sall3 Mef2c

#######
#######
#######

GRNs = readRDS("LGS13_RPC_GRNs_Jan2025")[[3]]
LGS13_RPC_promote_Cone = Perform_TF_enrich(GRNs,Class="pos",Cone_Genes,TAG2="RPC")

####
#### PlotTFs =  LGS13_Pre_promote_Cone$TF[1:20]
####
####

which(res$Gene == "Zic3")
which(res$Gene == "Pou2f1")
which(res$Gene == "Thrb")
which(res$Gene == "Rxrg")
which(res$Gene == "Sall3")
which(res$Gene == "Mef2c")


##########
##########

######
#####---------------------------------------------------------------------------------------------------------------------
#####----- combine the results and plot the heatmap !!!! ----------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
######
###### Figure2H ########
######

summary(LGS_M_Avg_Norm[,1])
summary(LGS_M_Avg_Norm[,3])


which(res$Gene == "Zic3")
which(res$Gene == "Pou2f1")
which(res$Gene == "Thrb")
which(res$Gene == "Rxrg")
which(res$Gene == "Sall3")
which(res$Gene == "Mef2c")


LGS_M_Avg_Norm["Zic3",]
LGS_M_Avg_Norm["Pou2f1",]
LGS_M_Avg_Norm["Thrb",]
LGS_M_Avg_Norm["Rxrg",]
LGS_M_Avg_Norm["Sall3",]
LGS_M_Avg_Norm["Mef2c",]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Pre_promote_Cone_Top30 <- readRDS(file="LGS13_Pre_promote_Cone_Top30.rds")


library('ComplexHeatmap')
library('circlize')

PlotTFs = LGS13_Pre_promote_Cone_Top30$TF

PlotTFs = gsub("CRX-like","Crx",PlotTFs)
PlotTFs = gsub("MEF2D-containing","Mef2d",PlotTFs)

Avg_Plot = LGS_M_Avg_Norm[PlotTFs,]
Avg_Plot_S <- t(apply(Avg_Plot,1,scale))
colnames(Avg_Plot_S) = colnames(Avg_Plot)
rownames(Avg_Plot_S) = rownames(Avg_Plot)

Avg_Plot_S = t(Avg_Plot_S)

library(ComplexHeatmap)
library(ggplot2)
library('circlize')

col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), c('#026AB1','#61BFB9',"#9EA220",'white','#EF9000','red',"#936DAD"))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
png('F2H_expreson.png',height=3000,width = 8000,res=72*12)
Heatmap(Avg_Plot_S, name = "XX", border = T,use_raster=TRUE,show_row_names=T,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 1),cluster_rows = F,cluster_columns = T,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-3,-2,-1,0,1,2,3)))
dev.off()

######
###### 下一步我们需要输出这些 networks ！！！！ ################
######
######
###### 我们看一下 Thrb #####
######

LGS_Cone_GRNs = readRDS("LGS13_Cone_GRNs_Jan2025")[[2]]
Mouse_Cone_GRNs = readRDS("Mouse_Cone_GRNs_Jan2025")[[2]]

LGS_Cone_GRNs_Thrb = LGS_Cone_GRNs[which(LGS_Cone_GRNs$gene == "Thrb"),]
Mouse_Cone_GRNs_Thrb = Mouse_Cone_GRNs[which(Mouse_Cone_GRNs$gene == "Thrb"),]

######
######
LGS_Cone_GRNs_Thrb[which(LGS_Cone_GRNs_Thrb$TF == "Otx2"),]
LGS_Cone_GRNs_Thrb[which(LGS_Cone_GRNs_Thrb$TF == "Neurod1"),]


######
######
######
######

GRNs = LGS_Cone_GRNs_Thrb

Output_GRNs_to_cyto(LGS_Cone_GRNs_Thrb,TAG="Thrb_LGS13_cyto")
Output_GRNs_to_cyto(Mouse_Cone_GRNs_Thrb,TAG="Thrb_Mouse_cyto")

#########
#########
#########

Output_GRNs_to_cyto <- function(GRNs,TAG="Thrb_Mouse_cyto"){
    ############
    GRNs$triple = paste(GRNs$TF,GRNs$peaks,GRNs$gene)
    k = which(duplicated(GRNs$triple) == T)
    ###########
    ############
    start1 = GRNs$TF
    start2 = GRNs$peaks
    end1 = GRNs$peaks
    end2 = GRNs$gene
    #####
    tab1 = data.frame(s=start1,e=end1)
    tab1_index = paste0(tab1$s,tab1$e)
    tab1 = tab1[!duplicated(tab1_index),]
    tab2 = data.frame(s=start2,e=end2)
    tab2_index = paste0(tab2$s,tab2$e)
    tab2 = tab2[!duplicated(tab2_index),]
    out_GRNS = rbind(tab1,tab2)
    ##### Next we will output the features #####
    #####
    f1 = data.frame(f=start1[!duplicated(start1)],class="TF")
    f2 = data.frame(f=start2[!duplicated(start2)],class="Peak")
    f3 = data.frame(f=end2[1],class="Target")
    ####
    out_Feaures = rbind(f1,f2,f3)
    ####
    write.table(out_GRNS,file=paste0(TAG,'_GRNs.txt'),sep="\t",row.names=F,quote=F)
    write.table(out_Feaures,file=paste0(TAG,'_Features.txt'),sep="\t",row.names=F,quote=F)
    ####
}

######
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
######


###### 接下来我们打开一下 AI 文件 ########
###### 我们检测一下 Figure2 ############
###### c.	Differentially expressed genes by species/cell type ###################
###### 

Add additional genes here and remove Onecut2, Otx2, and Neurod1.  Include the following:
C1: Tbx2, Rb1.
C5: Casz1, Sox11
C6: Nr2f1, Dlx2, Sox8.
C7: Nr2e3, Nrl, Insm2, Samd7.
It may be helpful to rearrange the figure somewhat so that the main cone (C3) and rod (C7)-associated clusters – and their associated TFs – can be highlighted.

####### we will convert 13LGS genes to Mouse genes in comparision ################
#######
####### 我们首先需要重复一下 Figure2 的 Heatmaps #################----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

library(Seurat)
library(Matrix)


setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')
head(LGS13_Mouse.combined@meta.data)
tail(LGS13_Mouse.combined@meta.data)
############ 这个 features 的数目好像不太对 ###############
############ 我们重新合并一下这个 seuart objects ##########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")

Mouse_seurat_compare = readRDS("Mouse_seurat_compare_Jan29_2025.rds")
LGS13_seurat_compare = readRDS("LGS_seurat_compare_Jan29_2025.rds")

Mouse_seurat_compare[["RNA"]]@counts
LGS13_seurat_compare[["RNA"]]@counts

############ Next we will add the annotations to the seurat and merge them #######
############ 直接merge两个seurat objects ####

LGS13_Mouse_Merged = merge(LGS13_seurat_compare, y = Mouse_seurat_compare,add.cell.ids = c("LGS", "Mouse"))
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Mouse_Merged,file="LGS13_Mouse_Merged_Jan2025")

############

k = which(colnames(LGS13_Mouse_Merged) %in% colnames(LGS13_Mouse.combined))
length(k)

dim(LGS13_Mouse.combined)
LGS13_Mouse_Merged_cl = LGS13_Mouse_Merged[,k]

############
table(LGS13_Mouse_Merged_cl$timepoint)
table(LGS13_Mouse_Merged_cl$celltype)

############ Next we will add others ###
m = match(colnames(LGS13_Mouse_Merged_cl),colnames(LGS13_Mouse.combined))

LGS13_Mouse_Merged_cl$celltype = LGS13_Mouse.combined$celltype[m]
LGS13_Mouse_Merged_cl$sp = LGS13_Mouse.combined$sp[m]

table(LGS13_Mouse.combined$celltype)

LGS13_Mouse_Merged_cl$index <- paste0(LGS13_Mouse_Merged_cl$sp,":",LGS13_Mouse_Merged_cl$celltype)
table(LGS13_Mouse_Merged_cl$index)
Idents(LGS13_Mouse_Merged_cl) <- "index"

#######
#######
#######
LGS13_Mouse_Merged_cl <- NormalizeData(LGS13_Mouse_Merged_cl)
saveRDS(LGS13_Mouse_Merged_cl,file="LGS13_Mouse_Merged_cl_Jan2025")

Markers1 = FindMarkers(LGS13_Mouse_Merged_cl,ident.1='LGS:RPC',ident.2='Mouse:RPC',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers1[which(rownames(Markers1) == "Otx2"),]
Markers1[which(rownames(Markers1) == "Neurod1"),]
Markers1[which(rownames(Markers1) == "Pou2f1"),]

Markers2 = FindMarkers(LGS13_Mouse_Merged_cl,ident.1='LGS:NG',ident.2='Mouse:NG',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers2[which(rownames(Markers2) == "Otx2"),]
Markers2[which(rownames(Markers2) == "Neurod1"),]
Markers2[which(rownames(Markers2) == "Pou2f1"),]

Markers3 = FindMarkers(LGS13_Mouse_Merged_cl,ident.1='LGS:BC_Photo_pre',ident.2='Mouse:BC_Photo_pre',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers3[which(rownames(Markers3) == "Onecut1"),]
Markers3[which(rownames(Markers3) == "Onecut2"),]
Markers3[which(rownames(Markers3) == "Zic3"),]

Markers4 = FindMarkers(LGS13_Mouse_Merged_cl,ident.1='LGS:Cone',ident.2='Mouse:Cone',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers4[which(rownames(Markers4) == "Onecut1"),]
Markers4[which(rownames(Markers4) == "Onecut2"),]
Markers4[which(rownames(Markers4) == "Zic3"),]
Markers4[which(rownames(Markers4) == "Thrb"),]
Markers4[which(rownames(Markers4) == "Rxrg"),]
Markers4[which(rownames(Markers4) == "Sall3"),]
Markers4[which(rownames(Markers4) == "Mef2c"),]

Markers5 = FindMarkers(LGS13_Mouse_Merged_cl,ident.1='LGS:Cone',ident.2='Mouse:Rod',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers5[which(rownames(Markers5) == "Onecut1"),]
Markers5[which(rownames(Markers5) == "Onecut2"),]
Markers5[which(rownames(Markers5) == "Zic3"),]
Markers5[which(rownames(Markers5) == "Thrb"),]
Markers5[which(rownames(Markers5) == "Rxrg"),]
Markers5[which(rownames(Markers5) == "Sall3"),]
Markers5[which(rownames(Markers5) == "Mef2c"),]

Markers_All = list(Markers1,Markers2,Markers3,Markers4,Markers5)

for(i in 1:length(Markers_All)){
    tmp = Markers_All[[i]]
    tmp$Gene = rownames(tmp)
    Markers_All[[i]] = tmp
}

###########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
saveRDS(Markers_All,file="Markers_All_Jan2025.rds")


####### --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### in the next section, we will get the genes on the heatmaps perform cluster analysis ###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####### we first will get the average expression for all these cells !!!!! #####
#######

Celltypes = names(table(LGS13_Mouse_Merged_cl$index))
Celltypes_Need = Celltypes[c(10,7,3,4,20,17,13,14,19)]

############
RNA_Counts_Mat = LGS13_Mouse_Merged_cl[['RNA']]@counts
RNA_Counts_Mat_Meta = LGS13_Mouse_Merged_cl@meta.data

###Mat = RNA_Counts_Mat
###Meta = RNA_Counts_Mat_Meta

Get_Avg_Exp <- function(Mat,Meta,Celltypes_Need){
    #########
    all_avg = c()
    for(i in 1:length(Celltypes_Need)){
        #########
        tmp_ct = Celltypes_Need[i]
        print(tmp_ct)
        k = which(RNA_Counts_Mat_Meta$index == tmp_ct)
        tmp_mat = Mat[,k]
        ##########
        tmp_sum = apply(tmp_mat,1,sum)
        ##########
        all_avg <- c(all_avg,tmp_sum)
    }
    ####
    all_avg_mat = matrix(all_avg,ncol=length(Celltypes_Need))
    rownames(all_avg_mat) = rownames(Mat)
    colnames(all_avg_mat) = Celltypes_Need
    ####
    factor = colSums(all_avg_mat) / 1e5
    all_avg_mat = sweep(all_avg_mat,2,factor,FUN="/")
    all_avg_mat_log = log2(all_avg_mat+1)
    ####
    return(all_avg_mat_log)
}

Avg_Exp_for_Figure2 = Get_Avg_Exp(Mat = RNA_Counts_Mat,Meta = RNA_Counts_Mat_Meta,Celltypes_Need=Celltypes_Need)

Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Thrb"),]
Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Zic3"),]
Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Rxrg"),]
Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Onecut1"),]
Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Onecut2"),]
Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) == "Mef2c"),]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
saveRDS(Avg_Exp_for_Figure2,file="Avg_Exp_for_Figure2_Jan2025")

######
###### filter the Exp Mat by DEGs ######
######

Markers_All_Merge = do.call("rbind",Markers_All)

Markers_All_Merge[which(Markers_All_Merge$Gene == "Onecut2"),]
Markers_All_Merge[which(Markers_All_Merge$Gene == "Onecut1"),]

Markers_All_Merge_Cl = Markers_All_Merge[which(Markers_All_Merge$p_val_adj < 1e-6 & abs(Markers_All_Merge$avg_log2FC) > 0.35),]
Average_Exp_Mat_Cl = Avg_Exp_for_Figure2[which(rownames(Avg_Exp_for_Figure2) %in% Markers_All_Merge_Cl$Gene == T),]

Average_Exp_Mat_Cl_Scale = t(apply(Average_Exp_Mat_Cl,1,scale))
colnames(Average_Exp_Mat_Cl_Scale) = colnames(Average_Exp_Mat_Cl)
rownames(Average_Exp_Mat_Cl_Scale) = rownames(Average_Exp_Mat_Cl) 


####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### in the next section, we will perform the K-means on these DEGs !!!!!!!! ###------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

k = which(Average_Exp_Mat_Cl_Scale > 2)
Average_Exp_Mat_Cl_Scale[k] = 2
k = which(Average_Exp_Mat_Cl_Scale < -2)
Average_Exp_Mat_Cl_Scale[k] = -2

kc = kmeans(Average_Exp_Mat_Cl_Scale, 8, iter.max = 100)
kc_dat = data.frame(genes = names(kc$cluster),cluster=as.numeric(kc$cluster))
kc_dat[which(kc_dat$genes == 'Zic3'),]

strsp_index = sapply(strsplit(colnames(Average_Exp_Mat_Cl_Scale),split=":"),function(x) x[[1]])
######
col_sp = as.factor(strsp_index)
row_sp = as.factor(kc_dat$cluster)

library('ComplexHeatmap')
library('circlize')

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")

col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))
#col_fun = colorRamp2(c(-2,-1,0,1,2), c('#9EA220','#86BF38','white','#EF9000','#D11536'))
#col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c('#026AB1','#61BFB9','white','#F6BA00','#EF9000'))

labels = c("Tb2","Rb1","Casz1","Sox11","Nr2f1","Dlx2","Sox8","Nr2e3","Nrl","Insm2","Samd7","Mef2c","Onecut1","Zic3","Sall3","Thrb","Rxrg")
at = match(labels,rownames(Average_Exp_Mat_Cl_Scale))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('RNA_test1.png',height=4000,width=3000,res=72*12)
Heatmap(Average_Exp_Mat_Cl_Scale, name = "XX", border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#####
##### we will add the Genes annotations ######
#####

order = c(7,2,5,4,8,6,1,3)

Reorder_Kmeans <- function(x,order){
	m = match(x,order)
	m = factor(m,levels=c(1:max(m)))
	return(m)
}

kc_dat2 = kc_dat
kc_dat2$cluster = Reorder_Kmeans(kc_dat2$cluster,order)

col_sp = as.factor(strsp_index)
row_sp = as.factor(kc_dat2$cluster)

library('ComplexHeatmap')
library('circlize')

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")

col_fun = colorRamp2(c(-2,-1,0,1,2), c('#026AB1','lightblue','white','#EF9000','#D11536'))
#col_fun = colorRamp2(c(-2,-1,0,1,2), c('#9EA220','#86BF38','white','#EF9000','#D11536'))
#col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c('#026AB1','#61BFB9','white','#F6BA00','#EF9000'))

labels = c("Tb2","Rb1","Casz1","Sox11","Nr2f1","Dlx2","Sox8","Nr2e3","Nrl","Insm2","Samd7","Mef2c","Onecut1","Zic3","Sall3","Thrb","Rxrg")
at = match(labels,rownames(Average_Exp_Mat_Cl_Scale))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('RNA_test2.png',height=4000,width=3000,res=72*12)
Heatmap(Average_Exp_Mat_Cl_Scale, name = "XX", border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#### kc_dat2_tab = kc_dat2_tab[,c(1,2)]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
kc_dat2_tab = kc_dat2[order(kc_dat2$cluster,decreasing=F),]
saveRDS(kc_dat2_tab,file="kc_dat2_tab_Jan_2025.rds")

m = match(kc_dat2_tab$genes,rownames(Average_Exp_Mat_Cl))
Average_Exp_Mat_Cl_re = Average_Exp_Mat_Cl[m,]
all.equal(rownames(Average_Exp_Mat_Cl_re),kc_dat2_tab$genes)
kc_dat2_tab = cbind(kc_dat2_tab,Average_Exp_Mat_Cl_re)
saveRDS(kc_dat2_tab,file="kc_dat2_tab_Jan_2025.rds")

table(kc_dat2_tab$cluster)

install.packages("writexl")
library(writexl)
write_xlsx(kc_dat2_tab, "Figure2c_MouseVS13LGS_DEGs.xlsx")

##############
############## OK!!!! Good!! ###
##############

############## 接下来我们需要处理一下 PtoG correlations ###########

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

############### PtoG 的问题在于 Genes 不能用 covert 过的，只能用原来的基因 ##################
###############
############### Peak 不需要重新 call ！！！！ ############

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

######
###### OK!!! Let integrate each time points #######
######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

###########-------------
conda activate ArchR
R

library(ArchR)
library(Seurat)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
#### saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_cl_Nov18.rds")
#### saveRDS(Mouse_RNA_seurat_merge_cl,file="Mouse_RNA_seurat_merge_cl_Nov18.rds")

Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_cl_Nov18.rds")

###### 这个gene 少了！！！ ##### Mouse_RNA_seurat_merge_cl <- readRDS("Mouse_RNA_seurat_merge_cl_Nov18.rds")
###### 需要更换 ！！！！ #######

table(Mouse_Project_new_server_cl$Sample)
Mouse_RNA_seurat_merge_cl <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development/Mouse_RNA_dev_raw_seurat_Anno_Jan29_2025.rds")
Mouse_RNA_seurat_merge_cl$time = Mouse_RNA_seurat_merge_cl$timepoint

###### for E11 !!! ####-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
E11_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E11")]
E11_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E11") == T)]
DefaultAssay(E11_rna) <- "RNA"
E11_rna <- NormalizeData(E11_rna)
E11_project <- addGeneScoreMatrix(E11_project,force = TRUE)

E11_project <- addIterativeLSI(
    ArchRProj = E11_project,
    useMatrix = "TileMatrix", 
    name = "E11_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

E11_project <- addGeneIntegrationMatrix(
    ArchRProj = E11_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "E11_IterativeLSI",
    seRNA = E11_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
############
############
table(Mouse_Project_new_server_cl$Sample)
E12_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E12")]
table(Mouse_RNA_seurat_merge_cl$time)
E12_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E12_rep1") == T)]
DefaultAssay(E12_rna) <- "RNA"
E12_rna <- NormalizeData(E12_rna)

E12_project <- addGeneScoreMatrix(E12_project,force = TRUE)

E12_project <- addIterativeLSI(
    ArchRProj = E12_project,
    useMatrix = "TileMatrix", 
    name = "E12_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

E12_project <- addGeneIntegrationMatrix(
    ArchRProj = E12_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "E12_IterativeLSI",
    seRNA = E12_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table(Mouse_Project_new_server_cl$Sample)
E14_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E14")]
table(Mouse_RNA_seurat_merge_cl$time)
E14_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E14_rep1","E14_rep2") == T)]
DefaultAssay(E14_rna) <- "RNA"
E14_rna <- NormalizeData(E14_rna)

E14_project <- addGeneScoreMatrix(E14_project,force = TRUE)

E14_project <- addIterativeLSI(
    ArchRProj = E14_project,
    useMatrix = "TileMatrix", 
    name = "E14_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

E14_project <- addGeneIntegrationMatrix(
    ArchRProj = E14_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "E14_IterativeLSI",
    seRNA = E14_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(Mouse_Project_new_server_cl$Sample)
E16_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E16")]
table(Mouse_RNA_seurat_merge_cl$time)
E16_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E16") == T)]
DefaultAssay(E16_rna) <- "RNA"
E16_rna <- NormalizeData(E16_rna)

E16_project <- addGeneScoreMatrix(E16_project,force = TRUE)

E16_project <- addIterativeLSI(
    ArchRProj = E16_project,
    useMatrix = "TileMatrix", 
    name = "E16_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

E16_project <- addGeneIntegrationMatrix(
    ArchRProj = E16_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "E16_IterativeLSI",
    seRNA = E16_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(Mouse_Project_new_server_cl$Sample)
E18_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E18")]
table(Mouse_RNA_seurat_merge_cl$time)
E18_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E18_rep2","E18_rep3") == T)]
DefaultAssay(E18_rna) <- "RNA"
E18_rna <- NormalizeData(E18_rna)

E18_project <- addGeneScoreMatrix(E18_project,force = TRUE)

E18_project <- addIterativeLSI(
    ArchRProj = E18_project,
    useMatrix = "TileMatrix", 
    name = "E18_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)



E18_project <- addGeneIntegrationMatrix(
    ArchRProj = E18_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "E18_IterativeLSI",
    seRNA = E18_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table(Mouse_Project_new_server_cl$Sample)
P0_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "P0")]
table(Mouse_RNA_seurat_merge_cl$time)
P0_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("P0") == T)]
DefaultAssay(P0_rna) <- "RNA"
P0_rna <- NormalizeData(P0_rna)

P0_project <- addGeneScoreMatrix(P0_project,force = TRUE)

P0_project <- addIterativeLSI(
    ArchRProj = P0_project,
    useMatrix = "TileMatrix", 
    name = "P0_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)


P0_project <- addGeneIntegrationMatrix(
    ArchRProj = P0_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "P0_IterativeLSI",
    seRNA = P0_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table(Mouse_Project_new_server_cl$Sample)
P2_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "P2")]
table(Mouse_RNA_seurat_merge_cl$time)
P2_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("P2_rep2","P2_rep3") == T)]
DefaultAssay(P2_rna) <- "RNA"
P2_rna <- NormalizeData(P2_rna)

P2_project <- addGeneScoreMatrix(P2_project,force = TRUE)

P2_project <- addIterativeLSI(
    ArchRProj = P2_project,
    useMatrix = "TileMatrix", 
    name = "P2_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)


P2_project <- addGeneIntegrationMatrix(
    ArchRProj = P2_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "P2_IterativeLSI",
    seRNA = P2_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)
############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(Mouse_Project_new_server_cl$Sample)
P5_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "P5")]
table(Mouse_RNA_seurat_merge_cl$time)
P5_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("P5") == T)]
DefaultAssay(P5_rna) <- "RNA"
P5_rna <- NormalizeData(P5_rna)

P5_project <- addGeneScoreMatrix(P5_project,force = TRUE)

P5_project <- addIterativeLSI(
    ArchRProj = P5_project,
    useMatrix = "TileMatrix", 
    name = "P5_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

P5_project <- addGeneIntegrationMatrix(
    ArchRProj = P5_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "P5_IterativeLSI",
    seRNA = P5_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table(Mouse_Project_new_server_cl$Sample)
P8_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "P8")]
table(Mouse_RNA_seurat_merge_cl$time)
P8_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("P8_rep1","P8_rep2") == T)]
DefaultAssay(P8_rna) <- "RNA"
P8_rna <- NormalizeData(P8_rna)

P8_project <- addGeneScoreMatrix(P8_project,force = TRUE)

P8_project <- addIterativeLSI(
    ArchRProj = P8_project,
    useMatrix = "TileMatrix", 
    name = "P8_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)


P8_project <- addGeneIntegrationMatrix(
    ArchRProj = P8_project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "P8_IterativeLSI",
    seRNA = P8_rna,
    addToArrow = TRUE,
    force=TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

############----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
############ Next will be for the LGS integration!!!!! #######
############------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")
all_seurat_RNA_big_cl_detail_new <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge/all_seurat_RNA_big_cl_detail_new_202204")

LGS13_Project_new_server_cl@geneAnnotation$exons$symbol
LGS13_Project_new_server_cl@geneAnnotation$genes
LGS13_Project_new_server_cl@geneAnnotation$TSS

LGS13_Project_new_server_cl@geneAnnotation$exons$gene_id = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$exons$gene_id)
LGS13_Project_new_server_cl@geneAnnotation$exons$gene_id = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$exons$gene_id)
LGS13_Project_new_server_cl@geneAnnotation$exons$symbol = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$exons$symbol)
LGS13_Project_new_server_cl@geneAnnotation$exons$symbol = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$exons$symbol)

LGS13_Project_new_server_cl@geneAnnotation$genes$gene_id = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$genes$gene_id)
LGS13_Project_new_server_cl@geneAnnotation$genes$gene_id = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$genes$gene_id)
LGS13_Project_new_server_cl@geneAnnotation$genes$symbol = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$genes$symbol)
LGS13_Project_new_server_cl@geneAnnotation$genes$symbol = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$genes$symbol)

LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_id = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_id)
LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_id = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_id)
LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_name = gsub("_like","-like",LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_name)
LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_name = gsub("_containing","-containing",LGS13_Project_new_server_cl@geneAnnotation$TSS$tx_name)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_Nov2024")
######## DefaultAssay(all_seurat_RNA_big_cl_detail_new) <- "RNA"
######## -------------------------------------------------------------------------------
Integrate_Fun <- function(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples){
    ########
    ########
    tmp_project = ArchR_project[which(ArchR_project$Sample %in% ATAC_samples == T)]
    tmp_rna = RNA_seurat[,which(RNA_seurat$timepoint %in% RNA_samples == T)]
    DefaultAssay(tmp_rna) <- "RNA"
    tmp_rna <- NormalizeData(tmp_rna)
    print(dim(tmp_rna))
    ########
    ########
    ###getAvailableMatrices(tmp_project)
    tmp_project <- addGeneScoreMatrix(tmp_project,force = TRUE)
    #######
    tmp_project <- addIterativeLSI(
        ArchRProj = tmp_project,
        useMatrix = "TileMatrix", 
        name = "tmp_IterativeLSI", 
        iterations = 2, 
        clusterParams = list( #See Seurat::FindClusters
            resolution = c(0.2), 
            sampleCells = 10000, 
            n.start = 10
        ),  
        varFeatures = 25000, 
        dimsToUse = 1:30,
        force = TRUE
    )
    ####
    tmp_project <- addGeneIntegrationMatrix(
        ArchRProj = tmp_project, 
        useMatrix = "GeneScoreMatrix",
        matrixName = "GeneIntegrationMatrix",
        reducedDims = "tmp_IterativeLSI",
        seRNA = tmp_rna,
        addToArrow = TRUE,
        force=TRUE,
        groupRNA = "new_celltypes",
        nameCell = "predictedCell_Un",
        nameGroup = "predictedGroup_Un",
        nameScore = "predictedScore_Un"
    )
    #####
    print("Done!!!!")
    ########
}
#####-------------------------------------------------------------------------------------------
 
##### Found 30646 overlapping gene names ####
##### Found 31372 overlapping gene names 
##### Found 32357 overlapping gene names

##### E18 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E18_1","E18_2A","E18_2B")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E18")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

##### E21 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E21_1","E21_2","E21_3B")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E21")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

##### E24 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E24_1","E24_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E24")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

##### E26 ------------------------------------------------------------------------------------------------------------------------------------------------------
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E26_1","E26_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E26")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

##### P1 ######
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P1_1","P1_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P1")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

##### P4 #######
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P4_1","P4_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P4")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

###### P8 #########
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P8_1","P8_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P8")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P12  #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P12_1","P12_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P12")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P17  #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P17_1","P17_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P17")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)


####### P21 #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P21_1","P21_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P21")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)


####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### 下一步我们计算一下 Mouse 的PtoG correlation #########
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Get_p2g_fun <- function(x){
	corCutOff = 0.25
	FDRCutOff = 1e-6
	varCutOffATAC = 0.7
	varCutOffRNA = 0.3
	p2g <- metadata(x@peakSet)$Peak2GeneLinks
	p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
	if(!is.null(varCutOffATAC)){
    	p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
	}
	if(!is.null(varCutOffRNA)){
    	p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
	}
	mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
	mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
	p2g$peak <- paste0(rowRanges(mATAC))
	p2g$gene <- rowData(mRNA)$name
	return(p2g)
}
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#### Nov18 ###

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_cl_Nov18.rds")

length(Mouse_Project_new_server_cl@peakSet)
#### 301739 !!! ###
####

Mouse_Project_new_server_cl <- addIterativeLSI(
    ArchRProj = Mouse_Project_new_server_cl,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)


Mouse_Project_new_server_cl <- addHarmony(
    ArchRProj = Mouse_Project_new_server_cl,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)


Mouse_Project_new_server_cl <- addUMAP(
    ArchRProj = Mouse_Project_new_server_cl, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)


p1 <- plotEmbedding(ArchRProj = Mouse_Project_new_server_cl, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-Sample-Clusters_2025.pdf", ArchRProj = Mouse_Project_new_server_cl, addDOC = FALSE, width = 5, height = 5)


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_Jan_2025")

##########-------------------

cell_index_meta = Mouse_Project_new_server_cl@cellColData
table(cell_index_meta$Sample)
cell_index_meta_cl = cell_index_meta[which(cell_index_meta$Sample %in% c("E18","P0","P2","P5","P8") == T),]
table(cell_index_meta_cl$Sample)
table(cell_index_meta_cl$celltype2)

#### 我们不加 MG cells ------------------------------------------------------------------------------------
cell_index_meta_cl = cell_index_meta_cl[which(cell_index_meta_cl$celltype2 %in% c("10_Cone","11_Early_Rod","12_Rod","2_RPCs_S2","3_RPCs_S3","5_Early_NG","6_Late_NG","9_Early_Cone") == T),]
table(cell_index_meta_cl$celltype2)
getAvailableMatrices(Mouse_Project_new_server_cl)
cell_list = rownames(cell_index_meta_cl)

Mouse_Project_new_server_cl_PtoG_Sub = Mouse_Project_new_server_cl[cell_list]

Mouse_Project_new_server_cl_PtoG_Sub <- addIterativeLSI(
    ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

######

Mouse_Project_new_server_cl_PtoG_Sub <- addHarmony(
    ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub,
    reducedDims = "IterativeLSI",
    name = "Harmony_Late",
    groupBy = "Sample",
    force = TRUE
)


Mouse_Project_new_server_cl_PtoG_Sub <- addUMAP(
    ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub, 
    reducedDims = "Harmony_Late", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

#######

p1 <- plotEmbedding(ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1, name = "Mouse_sub_CT.pdf", ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub, addDOC = FALSE, width = 5, height = 5)

#######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
saveRDS(Mouse_Project_new_server_cl_PtoG_Sub,file="Mouse_Project_new_server_cl_PtoG_Sub_Jan_2025")


Mouse_Project_new_server_cl_PtoG_Sub <- addPeak2GeneLinks(
    ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub,
    reducedDims = "Harmony_Late",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)

#########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
saveRDS(Mouse_Project_new_server_cl_PtoG_Sub,file="Mouse_Project_new_server_cl_PtoG_Sub_Jan_2025")

Mouse_Project_new_server_cl_PtoG_Sub_RES = Get_p2g_fun(Mouse_Project_new_server_cl_PtoG_Sub)
saveRDS(Mouse_Project_new_server_cl_PtoG_Sub_RES,file="Mouse_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### 下一步我们计算一下 LGS13 的PtoG correlation #########
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

####### 读入之前的LGS13 project ########
conda activate ArchR 
R
library(ArchR)
library(Seurat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")

####### check the peaks ######
length(LGS13_Project_new_server_cl@peakSet)
429178 ranges  OK!!! ###
##### Finished Creating Union Peak Set (428909)!, 4.561 mins elapsed. ######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Avg_PeakMat_Late <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3/LGS13_Avg_PeakMat_Late_2024Dec9")
#####
LGS13_Project_new_server_cl = addPeakSet(LGS13_Project_new_server_cl,peakSet=GRanges(rownames(LGS13_Avg_PeakMat_Late)),force=TRUE)
length(LGS13_Project_new_server_cl@peakSet)
LGS13_Project_new_server_cl = addPeakMatrix(LGS13_Project_new_server_cl)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_Nov2024")
#####
##### Only use the late time points ######
#####

cell_index_meta = LGS13_Project_new_server_cl@cellColData
table(cell_index_meta$Sample)
cell_index_meta_cl = cell_index_meta[which(cell_index_meta$Sample %in% c("P1_1","P1_2","P12_1","P12_2","P17_1","P17_2","P21_1","P21_2","P4_1","P4_2","P8_1","P8_2") == T),]
table(cell_index_meta_cl$Sample)
table(cell_index_meta_cl$celltype2)
cell_index_meta_cl = cell_index_meta_cl[which(cell_index_meta_cl$celltype2 %in% c("1_Early_RPC","10_Cone","11_Rod","2_Late_RPC","4_Early_NG","5_Late_NG","9_Photoprecursor") == T),]
table(cell_index_meta_cl$celltype2)
getAvailableMatrices(LGS13_Project_new_server_cl)
cell_list = rownames(cell_index_meta_cl)

LGS13_Project_new_server_cl_PtoG_Sub = LGS13_Project_new_server_cl[cell_list]

#####


LGS13_Project_new_server_cl_PtoG_Sub <- addIterativeLSI(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub,
    useMatrix = "TileMatrix", 
    name = "Late_IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = TRUE
)

######

LGS13_Project_new_server_cl_PtoG_Sub <- addHarmony(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub,
    reducedDims = "Late_IterativeLSI",
    name = "Harmony_Late",
    groupBy = "Sample",
    force = TRUE
)


LGS13_Project_new_server_cl_PtoG_Sub <- addUMAP(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub, 
    reducedDims = "Harmony_Late", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)


p1 <- plotEmbedding(ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1, name = "LGS_sub_CT2.pdf", ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub, addDOC = FALSE, width = 5, height = 5)

##########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_cl_Jan_2025")


LGS13_Project_new_server_cl_PtoG_Sub <- addPeak2GeneLinks(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub,
    reducedDims = "Harmony_Late",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)

saveRDS(LGS13_Project_new_server_cl_PtoG_Sub,file="LGS13_Project_new_server_cl_PtoG_Sub_Jan_2025")

length(LGS13_Project_new_server_cl_PtoG_Sub@peakSet)

LGS13_Project_new_server_cl_PtoG_Sub_RES = Get_p2g_fun(LGS13_Project_new_server_cl_PtoG_Sub)
saveRDS(LGS13_Project_new_server_cl_PtoG_Sub_RES,file="LGS13_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")

dim(table(LGS13_Project_new_server_cl_PtoG_Sub_RES$gene))
grep("-containing",LGS13_Project_new_server_cl_PtoG_Sub_RES$gene)

##### 我知道为什么了 _like _containing #########
##### 所以得重新 integrate ！！！ ##############

LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation$genes$symbol

#####
############
##### ######
############

############
#####
############
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### 下一步我们需要把 peaks 对应到 gene 上面 ！！！！ ###########
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Find_one_genes_peak_fun = function(genename,Gtf,allpeaks){
    ############
    library(GenomicRanges)
    ############
    Target_G = Gtf[which(Gtf$symbol == genename)][1]
    Other_G = Gtf[which(Gtf$symbol != genename)]
    ############
    Target_G_tss_GR = promoters(Target_G, upstream = 1000, downstream = 1000)
    Target_G_ext_GR = promoters(Target_G, upstream = 500000, downstream = 500000)
    ############
    ############
    TSS_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G_tss_GR) > 0)]
    if(length(TSS_peaks) == 0){
        TSS_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        TSS_peaks_tab = data.frame(peaks=TSS_peaks,class="TSS")
    }
    Body_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G) > 0)]
    if(length(Body_peaks) == 0){
        Body_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        Body_peaks_tab = data.frame(peaks=Body_peaks,class="Body")
    }
    Inter_peaks = as.character(allpeaks)[which(countOverlaps(allpeaks,Target_G_ext_GR) > 0)]
    if(length(Inter_peaks) == 0){
        Inter_peaks_tab = data.frame(peaks="NA",class="NA")
    }else{
        Inter_peaks_tab = data.frame(peaks=Inter_peaks,class="Inter")
    }
    ############
    ############
    Gene_table = rbind(TSS_peaks_tab,Body_peaks_tab,Inter_peaks_tab)
    K = which(Gene_table$peaks != "NA")
    if(length(K) == 0){
        return(data.frame(peaks="NA",class="NA",gene=genename))
    }
    Gene_table = Gene_table[K,]
    ############
    Gene_table = Gene_table[!duplicated(Gene_table$peaks),]
    ############
    ############ Next we will filterout peaks with overlap with Other genes TSS region and TSS body #####
    Gene_table_inter_peaks = GRanges(Gene_table$peaks[which(Gene_table$class == "Inter")])
    ############ get other G tss #####
    Other_G_TSS <- promoters(Other_G, upstream = 1000, downstream = 1000)
    ############ filterout with other gene body ####
    ############
    body_index = which(countOverlaps(Gene_table_inter_peaks,Other_G) > 0)
    tss_index = which(countOverlaps(Gene_table_inter_peaks,Other_G_TSS) > 0)
    ############
    remove_peaks = Gene_table_inter_peaks[c(body_index,tss_index)]
    ############
    if(length(remove_peaks) > 0){
        k = which(Gene_table$class == "Inter" & Gene_table$peaks %in% as.character(remove_peaks) == T)
        Gene_table = Gene_table[-k,]
    }
    ############
    if(dim(Gene_table)[1] == 0){
        return(data.frame(peaks="NA",class="NA",gene=genename))
    }
    Gene_table$gene=genename
    return(Gene_table)
}


##### for all the peaks in the analysis ######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
LGS13_Avg_PeakMat_Late <- readRDS("LGS13_Avg_PeakMat_Late_2024Dec9")
Mouse_Avg_PeakMat_Late <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3/Mouse_Avg_PeakMat_Late_2024Dec9")

#####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### for 13LGS #####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_13_Peak_to_Gene <- readRDS("LGS_13_Peak_to_Gene_2025")
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
Mouse_Peak_to_Gene <- readRDS("Mouse_Peak_to_Gene_2025")

LGS_13_Peak_to_Gene_merge = do.call("rbind",LGS_13_Peak_to_Gene)
Mouse_Peak_to_Gene_merge = do.call("rbind",Mouse_Peak_to_Gene)

##### 这个是下划线 ！！！ ###################
grep("-like",LGS_13_Peak_to_Gene_merge$gene)
grep("-containing",LGS_13_Peak_to_Gene_merge$gene)

LGS_13_Peak_to_Gene_merge$gene = gsub("_like","-like",LGS_13_Peak_to_Gene_merge$gene)
LGS_13_Peak_to_Gene_merge$gene = gsub("_containing","-containing",LGS_13_Peak_to_Gene_merge$gene)

#######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS_13_Peak_to_Gene_merge,file="LGS_13_Peak_to_Gene_merge_2025")
saveRDS(Mouse_Peak_to_Gene_merge,file="Mouse_Peak_to_Gene_merge_2025")

#######--------
#######--------
##### for Mouse #####

Mouse_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Mouse_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")
LGS13_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")

##### "-like" #####
##### ####### ##### This time we will add the PtoGs according to each gene, and Plot the accessibilty ########
#####
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)

########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_PeakAvg = readRDS("Mouse_Avg_PeakMat_Nov2024")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_PeakAvg = readRDS("LGS13_Avg_PeakMat_Nov2024")

######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
Mouse_Peak_to_Gene_merge <- readRDS(file="Mouse_Peak_to_Gene_merge_2025")
LGS_13_Peak_to_Gene_merge <- readRDS(file="LGS_13_Peak_to_Gene_merge_2025")

Mouse_Peak_to_Gene_merge_split = split(Mouse_Peak_to_Gene_merge,Mouse_Peak_to_Gene_merge$gene)
LGS_13_Peak_to_Gene_merge_split = split(LGS_13_Peak_to_Gene_merge,LGS_13_Peak_to_Gene_merge$gene)

#######
####### PeaktoGene = LGS_13_Peak_to_Gene
#######
Add_peak_anno_tothe_list <- function(PeaktoGene){
    #######
    library(parallel)
    #######
    ####### tmp = PeaktoGene[[1]]
    ####### add the
    One_fun <- function(x){
        tmp_list = split(x,x$class)
        for(i in 1:length(tmp_list)){
            tmp_list[[i]]$index = paste0(tmp_list[[i]]$gene,":",tmp_list[[i]]$class,":",1:length(tmp_list[[i]]$class))
        }
        tmp_list_res = do.call("rbind",tmp_list)
        return(tmp_list_res)
    }
    #######
    num_cores <- detectCores() - 1
    res_list = mclapply(PeaktoGene,One_fun,mc.cores = num_cores)
    #######
    return(res_list)
}

LGS_13_Peak_to_Gene_Anno = Add_peak_anno_tothe_list(LGS_13_Peak_to_Gene_merge_split)
Mouse_Peak_to_Gene_Anno = Add_peak_anno_tothe_list(Mouse_Peak_to_Gene_merge_split)

###########
#### check the -like -containing in the Gene_Anno ####
###########

grep("-like",names(LGS_13_Peak_to_Gene_Anno))
grep("-containing",names(LGS_13_Peak_to_Gene_Anno))

###########
###########
###########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
kc_dat2_tab = readRDS("kc_dat2_tab_Jan_2025.rds")

###########
########### for each gene in that cluster, get the high correlated peaks #######
########### and add the peak annotations 
###########

########### if one gene matched to different genes, remove that ############
###########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")
squirrel_dup_genes = squirrel_to_mouse_cl$Squirrel.name[which(duplicated(squirrel_to_mouse_cl$Squirrel.name) == T)]
mouse_dup_genes = squirrel_to_mouse_cl$Gene.name[which(duplicated(squirrel_to_mouse_cl$Gene.name) == T)]
k = which(squirrel_to_mouse_cl$Squirrel.name %in% squirrel_dup_genes == T | squirrel_to_mouse_cl$Gene.name %in% mouse_dup_genes == T)
squirrel_to_mouse_cl_onetoone = squirrel_to_mouse_cl[-k,]

saveRDS(squirrel_to_mouse_cl_onetoone,file="squirrel_to_mouse_cl_onetoone_2025.rds")

###########
###########
Mouse_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Mouse_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")
summary(abs(Mouse_PtoG$Correlation))
summary(abs(LGS13_PtoG$Correlation))


###########
###########

Mouse_Peak_to_Gene_merge

###########
###########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
Mouse_Avg_PeakMat_Late <- readRDS(file="Mouse_Avg_PeakMat_Late_2024Dec9")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
LGS13_Avg_PeakMat_Late <- readRDS("LGS13_Avg_PeakMat_Late_2024Dec9")

########### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Mouse_Peak_to_Gene_Anno_merge = do.call("rbind",Mouse_Peak_to_Gene_Anno)
LGS13_Peak_to_Gene_Anno_merge = do.call("rbind",LGS_13_Peak_to_Gene_Anno)

Peak_to_Gene_merge = Mouse_Peak_to_Gene_Anno_merge
Cluster_res = kc_dat2_tab 
Peak_Mat = Mouse_Avg_PeakMat_Late
PtoG = Mouse_PtoG

########## which(OnetoOne$Gene.name == "Sox11")

Prepare_Peak_Heatmaps_Mouse <- function(Peak_to_Gene_merge,Cluster_res,Peak_Mat,PtoG){
    ###########
    ###########
    Cluster_res$cluster = as.numeric(Cluster_res$cluster)
    all_cluster = min(Cluster_res$cluster):max(Cluster_res$cluster)
    all_Plot_res = list()
    for(i in all_cluster){
        print(i)
        #######
        tmp_cluster = Cluster_res[which(Cluster_res$cluster == i),]
        tmp_cluster$gene_order = 1:dim(tmp_cluster)[1]
        print(dim(tmp_cluster))
        tmp_cluster_cl = tmp_cluster
        ####### using the gene list to find the peaks #######
        tmp_Peak_to_Gene_merge = Peak_to_Gene_merge[which(Peak_to_Gene_merge$gene %in% tmp_cluster_cl$genes == T),]
        ####### merged with PtoG ###########
        PtoG$index2 = paste0(PtoG$gene,"@",PtoG$peak)
        tmp_Peak_to_Gene_merge$index2 = paste0(tmp_Peak_to_Gene_merge$gene,"@",tmp_Peak_to_Gene_merge$peaks)
        k = which(tmp_Peak_to_Gene_merge$index2 %in%  PtoG$index2 == T)
        tmp_Peak_to_Gene_merge_cl = tmp_Peak_to_Gene_merge[k,]
        m = match(tmp_Peak_to_Gene_merge_cl$index2,PtoG$index2)
        tmp_Peak_to_Gene_merge_cl$PtoG = PtoG$Correlation[m]
        tmp_Peak_to_Gene_merge_cl$PtoG_FDR = PtoG$FDR[m]
        ####### filtered by the posititve #################
        tmp_Peak_to_Gene_merge_clcl = tmp_Peak_to_Gene_merge_cl[which(tmp_Peak_to_Gene_merge_cl$PtoG > 0.25),]
        ### length(table(tmp_Peak_to_Gene_merge_clcl$gene))
        m = match(tmp_Peak_to_Gene_merge_clcl$gene,tmp_cluster$genes)
        tmp_Peak_to_Gene_merge_clcl$order = tmp_cluster$gene_order[m]
        tmp_Peak_to_Gene_merge_clcl = tmp_Peak_to_Gene_merge_clcl[order(tmp_Peak_to_Gene_merge_clcl$order),]
        ####
        #### Get Peak matrix for these: ######
        ####
        match_1 = which(tmp_Peak_to_Gene_merge_clcl$peaks %in% rownames(Peak_Mat) == T)
        length(match_1)
        dim(tmp_Peak_to_Gene_merge_clcl)
        ####
        match_2 = match(tmp_Peak_to_Gene_merge_clcl$peaks,rownames(Peak_Mat))
        tmp_peak_matrix = Peak_Mat[match_2,]
        tmp_peak_matrix = tmp_peak_matrix[,c("Late_RPCs","Late_NG","Photo_BC","Cone","Rod")]
        #####
        tmp_Peak_to_Gene_merge_clcl$index3 = paste0(tmp_Peak_to_Gene_merge_clcl$index,":::",tmp_Peak_to_Gene_merge_clcl$index2)
        rownames(tmp_peak_matrix) = tmp_Peak_to_Gene_merge_clcl$index3 
        #####
        all_Plot_res = c(all_Plot_res,list(tmp_peak_matrix))
        ####
    }
    ########
    all_Plot_res_merge = do.call("rbind",all_Plot_res)
    ########
    length = sapply(all_Plot_res, function(x) dim(x)[1])
    rowSp = rep(all_cluster,length)
    ########
    ######## we need scale the matrix ######
    ########
    all_Plot_res_merge_s = t(apply(all_Plot_res_merge,1,scale))
    rownames(all_Plot_res_merge_s) = rownames(all_Plot_res_merge)
    colnames(all_Plot_res_merge_s) = colnames(all_Plot_res_merge)
    ########
    return(list(Peak_mat_plot = all_Plot_res_merge_s, rowSp= rowSp))
}

Mouse_Plot_RES = Prepare_Peak_Heatmaps_Mouse(Peak_to_Gene_merge,Cluster_res,Peak_Mat,PtoG)

####### ######

library(ggplot2)
library('ComplexHeatmap')
library('circlize')

colorList = c("#EC6F64","#61BFB9","#D11536","#936DAD","#A74997","#AAA9A9","#EF9000","#9EA220","#026AB1","#804537")

#col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c('#936DAD','#A74997','white','#EF9000','#F6BA00'))
#col_fun = colorRamp2(c(-2,-1,0,1,2), c('#9EA220','#86BF38','white','#EF9000','#D11536'))
col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c('#026AB1','#61BFB9','white','#F6BA00','#EF9000'))
row_sp = Mouse_Plot_RES[[2]]

#####
##### add gene anno for Rod genes ##############
#####

Plot = Mouse_Plot_RES[[1]]
rownames(Plot) = sapply(strsplit(rownames(Plot),split=":::"),function(x) x[[1]])

Rod_Genes = c("Sox11","Sox8","Dlx2","Casz1","Nrl","Nr2e3","Samd7")

Need_show = c()
for(i in 1:length(Rod_Genes)){
    #######
    index = grep(Rod_Genes[i],rownames(Plot))
    if(length(index > 0)){
        names = rownames(Plot)[index]
        names = names[order(names,decreasing=T)]
        ######
        if(length(index) == 1){
            Need_show = c(Need_show,names)
        }
        if(length(index) > 1){
            Need_show = c(Need_show,names[c(1:2)])
        }
    }
}
########
########


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('ATAC_Mouse_test2_2025.png',height=6000,width = 4500,res=72*12)
Heatmap(Plot,name = "XX", border = T,use_raster=TRUE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-1,-0.5,0,0.5,1)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = match(Need_show,rownames(Plot)),labels = Need_show))
dev.off()

#######
####### Next we will see the LGS peak matrix ########
####### 

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")

#######


Prepare_Peak_Heatmaps_LGS13 <- function(Peak_to_Gene_merge,Cluster_res,Peak_Mat,PtoG,squirrel_to_mouse_cl){
    ###########
    ###########
    Cluster_res$cluster = as.numeric(Cluster_res$cluster)
    all_cluster = min(Cluster_res$cluster):max(Cluster_res$cluster)
    all_Plot_res = list()
    for(i in all_cluster){
        print(i)
        #######
        tmp_cluster = Cluster_res[which(Cluster_res$cluster == i),]
        tmp_cluster$gene_order = 1:dim(tmp_cluster)[1]
        print(dim(tmp_cluster))
        tmp_cluster_cl = tmp_cluster
        #######
        ####### convert to LGS_genes #########
        #######
        M_Gs = tmp_cluster_cl$genes
        LGS_Gs = squirrel_to_mouse_cl$"Squirrel.name"[which(squirrel_to_mouse_cl$Gene.name %in% M_Gs == T)]
        #######
        LGS_Gs_table = data.frame(LGS_G = LGS_Gs)
        LGS_Gs_table$M_Gs = squirrel_to_mouse_cl$Gene.name[match(LGS_Gs_table$LGS_G,squirrel_to_mouse_cl$"Squirrel.name")]
        LGS_Gs_table$order = tmp_cluster$gene_order[match(LGS_Gs_table$M_Gs,tmp_cluster$genes)]
        #######
        ####### using the gene list to find the peaks #######
        tmp_Peak_to_Gene_merge = Peak_to_Gene_merge[which(Peak_to_Gene_merge$gene %in% LGS_Gs_table$LGS_G == T),]
        ####### merged with PtoG ###########
        PtoG$index2 = paste0(PtoG$gene,"@",PtoG$peak)
        tmp_Peak_to_Gene_merge$index2 = paste0(tmp_Peak_to_Gene_merge$gene,"@",tmp_Peak_to_Gene_merge$peaks)
        k = which(tmp_Peak_to_Gene_merge$index2 %in%  PtoG$index2 == T)
        tmp_Peak_to_Gene_merge_cl = tmp_Peak_to_Gene_merge[k,]
        m = match(tmp_Peak_to_Gene_merge_cl$index2,PtoG$index2)
        tmp_Peak_to_Gene_merge_cl$PtoG = PtoG$Correlation[m]
        tmp_Peak_to_Gene_merge_cl$PtoG_FDR = PtoG$FDR[m]
        ####### filtered by the posititve #################
        tmp_Peak_to_Gene_merge_clcl = tmp_Peak_to_Gene_merge_cl[which(tmp_Peak_to_Gene_merge_cl$PtoG > 0.25),]
        ### length(table(tmp_Peak_to_Gene_merge_clcl$gene))
        m = match(tmp_Peak_to_Gene_merge_clcl$gene,LGS_Gs_table$LGS_G)
        tmp_Peak_to_Gene_merge_clcl$order = LGS_Gs_table$order[m]
        tmp_Peak_to_Gene_merge_clcl = tmp_Peak_to_Gene_merge_clcl[order(tmp_Peak_to_Gene_merge_clcl$order),]
        ####
        #### Get Peak matrix for these: ######
        ####
        match_1 = which(tmp_Peak_to_Gene_merge_clcl$peaks %in% rownames(Peak_Mat) == T)
        length(match_1)
        dim(tmp_Peak_to_Gene_merge_clcl)
        ####
        match_2 = match(tmp_Peak_to_Gene_merge_clcl$peaks,rownames(Peak_Mat))
        tmp_peak_matrix = Peak_Mat[match_2,]
        tmp_peak_matrix = tmp_peak_matrix[,c("Late_RPCs","Late_NG","Photo_BC","Cone")]
        #####
        tmp_Peak_to_Gene_merge_clcl$index3 = paste0(tmp_Peak_to_Gene_merge_clcl$index,":::",tmp_Peak_to_Gene_merge_clcl$index2)
        rownames(tmp_peak_matrix) = tmp_Peak_to_Gene_merge_clcl$index3 
        #####
        all_Plot_res = c(all_Plot_res,list(tmp_peak_matrix))
        ####
    }
    ########
    all_Plot_res_merge = do.call("rbind",all_Plot_res)
    ########
    length = sapply(all_Plot_res, function(x) dim(x)[1])
    rowSp = rep(all_cluster,length)
    ########
    ######## we need scale the matrix ######
    ########
    all_Plot_res_merge_s = t(apply(all_Plot_res_merge,1,scale))
    rownames(all_Plot_res_merge_s) = rownames(all_Plot_res_merge)
    colnames(all_Plot_res_merge_s) = colnames(all_Plot_res_merge)
    ########
    return(list(Peak_mat_plot = all_Plot_res_merge_s, rowSp= rowSp))
}


Peak_to_Gene_merge = LGS13_Peak_to_Gene_Anno_merge
Cluster_res = kc_dat2_tab 
Peak_Mat = LGS13_Avg_PeakMat_Late
PtoG = LGS13_PtoG

LGS13_Plot_RES = Prepare_Peak_Heatmaps_LGS13(Peak_to_Gene_merge,Cluster_res,Peak_Mat,PtoG,squirrel_to_mouse_cl=squirrel_to_mouse_cl)

col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), c('#026AB1','#61BFB9','white','#F6BA00','#EF9000'))
row_sp = LGS13_Plot_RES[[2]]


Plot = LGS13_Plot_RES[[1]]
rownames(Plot) = sapply(strsplit(rownames(Plot),split=":::"),function(x) x[[1]])

Cone_Genes = c("Mef2c","Onecut1","Zic3","Sall3","Thrb","Rxrg")

Need_show = c()
for(i in 1:length(Cone_Genes)){
    #######
    index = grep(Cone_Genes[i],rownames(Plot))
    if(length(index > 0)){
        names = rownames(Plot)[index]
        names = names[order(names,decreasing=T)]
        ######
        if(length(index) == 1){
            Need_show = c(Need_show,names)
        }
        if(length(index) > 1){
            Need_show = c(Need_show,names[c(1:2)])
        }
    }
}
########
########


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('ATAC_LGS13_test2_2025.png',height=6000,width = 3500,res=72*12)
Heatmap(Plot,name = "XX", border = T,use_raster=TRUE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-1,-0.5,0,0.5,1)),row_split=row_sp) +
rowAnnotation(link = anno_mark(at = match(Need_show,rownames(Plot)),labels = Need_show))
dev.off()

########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######## Next we will calculate the Peak numbers between the LGS and Mouse !!! ####
########----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)


######## we need the clusters 
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
kc_dat2_tab <- readRDS("kc_dat2_tab_Jan_2025.rds")
Cluster_res = kc_dat2_tab 
######## we need the PtoGs #########
Mouse_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Mouse_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")
LGS13_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")

######## we need the conversion between Mouse and LGS13 ########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")
########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Peak_to_Gene_Anno_merge <- readRDS("LGS_13_Peak_to_Gene_merge_2025")
Mouse_Peak_to_Gene_Anno_merge <- readRDS("Mouse_Peak_to_Gene_merge_2025")


Find_Peaks_for_LGS13_and_Mouse <- function(PtoG_cutoff = 0.25,Mouse_PtoG,LGS13_PtoG,squirrel_to_mouse_cl,Cluster_res,Mouse_Peak_to_Gene_Anno_merge,LGS13_Peak_to_Gene_Anno_merge){
    ##########
    Mouse_PtoG_cl = Mouse_PtoG[which(Mouse_PtoG$Correlation > PtoG_cutoff & Mouse_PtoG$FDR < 0.01),]
    LGS13_PtoG_cl = LGS13_PtoG[which(LGS13_PtoG$Correlation > PtoG_cutoff & LGS13_PtoG$FDR < 0.01),]
    ########## 这些PtoG 必须还在 Mouse_Peak_to_Gene_Anno_merge LGS13_Peak_to_Gene_Anno_merge 里面 ########
    Mouse_PtoG_cl$index = paste0(Mouse_PtoG_cl$gene,"@",Mouse_PtoG_cl$peak)
    LGS13_PtoG_cl$index = paste0(LGS13_PtoG_cl$gene,"@",LGS13_PtoG_cl$peak)
    Mouse_Peak_to_Gene_Anno_merge$index = paste0(Mouse_Peak_to_Gene_Anno_merge$gene,"@",Mouse_Peak_to_Gene_Anno_merge$peaks)
    LGS13_Peak_to_Gene_Anno_merge$index = paste0(LGS13_Peak_to_Gene_Anno_merge$gene,"@",LGS13_Peak_to_Gene_Anno_merge$peaks)
    ##########
    Mouse_PtoG_cl = Mouse_PtoG_cl[which(Mouse_PtoG_cl$index %in% Mouse_Peak_to_Gene_Anno_merge$index == T),]
    LGS13_PtoG_cl = LGS13_PtoG_cl[which(LGS13_PtoG_cl$index %in% LGS13_Peak_to_Gene_Anno_merge$index == T),]
    ##########
    ## first we will add the cluster results to the squirrel_to_mouse_cl ####
    ##########
    k = which(squirrel_to_mouse_cl$Gene.name %in% Cluster_res$genes == T)
    squirrel_to_mouse_clcl = squirrel_to_mouse_cl[k,]
    m = match(squirrel_to_mouse_clcl$Gene.name,Cluster_res$genes)
    squirrel_to_mouse_clcl$cluster = Cluster_res$cluster[m]
    ##########
    squirrel_to_mouse_clcl_list = split(squirrel_to_mouse_clcl,squirrel_to_mouse_clcl$cluster)
    ##########
    Peaks_res = list()
    ##########
    for(i in 1:length(squirrel_to_mouse_clcl_list)){
        ############
        tmp_res = squirrel_to_mouse_clcl_list[[i]]
        ####
        tmp_res$LGS_Peaks = 0
        tmp_res$Mouse_Peaks = 0
        ############ find the peaks ############
        for(j in 1:length(tmp_res$Squirrel.name)){
            tmp = tmp_res$Squirrel.name[j]
            k = which(LGS13_PtoG_cl$gene == tmp)
            tmp_res$LGS_Peaks[j] = length(k)
        }
        ###########
        for(j in 1:length(tmp_res$Gene.name)){
            tmp = tmp_res$Gene.name[j]
            k = which(Mouse_PtoG_cl$gene == tmp)
            tmp_res$Mouse_Peaks[j] = length(k)
        }
        #######
        Peaks_res = c(Peaks_res,list(tmp_res))
    }
    Peaks_res_merge = do.call("rbind",Peaks_res)
    #####
    ##### Peaks_res_merge[which(Peaks_res_merge$Gene.name == "Thrb"),]
    ##### Peaks_res_merge[which(Peaks_res_merge$Gene.name == "Sall3"),]
    ##### Peaks_res_merge[which(Peaks_res_merge$Gene.name == "Zic3"),]
    ##### Peaks_res_merge[which(Peaks_res_merge$Gene.name == "Rxrg"),]
    #####
    ##### Next we will see the median numbers ###
    #####
    LGS13_PtoG_cl_median_res = median(table(LGS13_PtoG_cl$gene))
    Mouse_PtoG_cl_median_res = median(table(Mouse_PtoG_cl$gene))
    LGS13_PtoG_cl_median_res = mean(table(LGS13_PtoG_cl$gene))
    Mouse_PtoG_cl_median_res = mean(table(Mouse_PtoG_cl$gene))
    #####
    #####
    factor = LGS13_PtoG_cl_median_res/Mouse_PtoG_cl_median_res
    #####
    Peaks_res_merge$LGS_Peaks_norm = Peaks_res_merge$LGS_Peaks / factor
    Peaks_res_merge$Mouse_Peaks_norm = Peaks_res_merge$Mouse_Peaks
    #####
    return(Peaks_res_merge)
}

##########
#########
##########

Mouse_LGS13_peaks_compare_res = Find_Peaks_for_LGS13_and_Mouse(PtoG_cutoff = 0.25,Mouse_PtoG,LGS13_PtoG,squirrel_to_mouse_cl,Cluster_res,Mouse_Peak_to_Gene_Anno_merge,LGS13_Peak_to_Gene_Anno_merge)
Mouse_LGS13_peaks_compare_res$cluster = paste0("C",Mouse_LGS13_peaks_compare_res$cluster)
##########
##########

library(viridis)

ggplot(Mouse_LGS13_peaks_compare_res, aes(x = Mouse_Peaks_norm, y = LGS_Peaks_norm)) +
  geom_point(color = "black",alpha = 0.5,size=0.2) +
  geom_abline(intercept = 0, slope = 1 , color ="red",alpha=0.5)   +      # Scatter points
  facet_wrap(~cluster, nrow = 1, ncol = 8) +             # Facet into a 3x2 grid
  labs(x = "Condition", y = "Value", title = "") + 
  theme_classic() +
  theme(
    legend.position = "right",                              # Remove the legend
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add black borders to subplots
  ) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))

ggsave("Test123.png",height=2.5,width=14)

#########
#########
#########

tapply(Mouse_LGS13_peaks_compare_res$LGS_Peaks_norm,Mouse_LGS13_peaks_compare_res$cluster,mean)
tapply(Mouse_LGS13_peaks_compare_res$Mouse_Peaks_norm,Mouse_LGS13_peaks_compare_res$cluster,mean)

########
######## test #########
########

sub = Mouse_LGS13_peaks_compare_res[which(Mouse_LGS13_peaks_compare_res$cluster == "C8"),]
summary(sub$LGS_Peaks_norm)
summary(sub$Mouse_Peaks_norm)
t.test(sub$LGS_Peaks_norm,sub$Mouse_Peaks_norm,paired = TRUE)$p.val


#######
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######

####### Next we will continuous to search the GRNs ###########
#######

########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------















####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### GRNs #######
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### 这一步是重新跑一下 GRNs，只研究同源基因 -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


####### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### Get the total matrix, with the existing cell type annotations ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####### for Mouse:
####### load the Normalized matrix for Mouse !!!! ###########
#######
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development")
Mouse_count_matrix_total <- readRDS(file="Mouse_RNA_dev_raw_seurat_Anno_Jan29_2025.rds")

####### for LGS13:
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
LGS13_count_matrix_total = readRDS('all_seurat_RNA_big_cl_detail_new_202204')

####### ------------------------------------------------------------------------------------------------------------------------------------------------
####### 我们读入 combined seurat see the cell annotations #########
####### run the cell types "Cone","RPC","BC_Photo_pre","Rod","NG" ########
####### ------------------------------------------------------------------------------------------------------------------------------------------------

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')
table(LGS13_Mouse.combined$celltype)

Need_celltypes = c("Cone","RPC","BC_Photo_pre","Rod","NG")
Mouse_seurat = LGS13_Mouse.combined[,LGS13_Mouse.combined$sp == "Mouse"]
LGS13_seurat = LGS13_Mouse.combined[,LGS13_Mouse.combined$sp == "LGS"]

Mouse_Index = LGS13_Mouse.combined[,which(LGS13_Mouse.combined$celltype %in% Need_celltypes == T & LGS13_Mouse.combined$sp == "Mouse")]
LGS13_Index = LGS13_Mouse.combined[,which(LGS13_Mouse.combined$celltype %in% Need_celltypes == T & LGS13_Mouse.combined$sp == "LGS")]


######## Get Mouse GRNs input #####
Mouse_Index$cellid = gsub("Mouse_","",colnames(Mouse_Index))
Mouse_count_matrix_sub = Mouse_count_matrix_total[,which(colnames(Mouse_count_matrix_total) %in% Mouse_Index$cellid == T)]
m = match(colnames(Mouse_count_matrix_sub),Mouse_Index$cellid)
Mouse_count_matrix_sub$CT = Mouse_Index$celltype[m]
table(Mouse_count_matrix_sub$CT)


######## Get LGS13 GRNs input #####
LGS13_Index$cellid = gsub("LGS_","",colnames(LGS13_Index))
LGS13_count_matrix_sub = LGS13_count_matrix_total[,which(colnames(LGS13_count_matrix_total) %in% LGS13_Index$cellid == T)]
m = match(colnames(LGS13_count_matrix_sub),LGS13_Index$cellid)
LGS13_count_matrix_sub$CT = LGS13_Index$celltype[m]
table(LGS13_count_matrix_sub$CT)

########
######## ######## ####### Next we will filter genes ######## ######## ###########
########

DefaultAssay(Mouse_count_matrix_sub) <- "RNA"
DefaultAssay(LGS13_count_matrix_sub) <- "RNA"

########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
load("PWM_list_combine_cl")

########
All_names = c()
All_ID = c()
for(i in 1:length(PWM_list_combine_cl)){
    tmp = PWM_list_combine_cl[[i]]
    tmp_res = tmp@name
    All_names = c(All_names,tmp_res)
    All_ID = c(All_ID,tmp@ID)
}

Table_list = data.frame(ID=All_ID,TF=All_names)
TFs = c()
for(i in 1:length(Table_list$TF)){
    tmp = Table_list$TF[i]
    tmp_list = unlist(strsplit(tmp,split=";"))
    TFs = c(TFs,tmp_list)
}
TFs[5566] = "Atoh7"
TFs = TFs[!duplicated(TFs)]

################
################ Next we will see the homo TFs #########
################

################ filter we will filter the matrix ######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
##### one to one list ####
LGStoM_datatable = readRDS("LGStoM_datatable_2025.rds")

#####------------- this step we will output the oneVSone gene matched table for excel #######
colnames(LGStoM_datatable) <- c("13-LGS_Genes","Mouse_Genes")

library(openxlsx)
write.xlsx(LGStoM_datatable, file = "13LGS_Mouse_Genes_conserved.xlsx")
##### add this into the methods sections ########
#####


#####------------- 
k = which(rownames(Mouse_count_matrix_sub) %in% LGStoM_datatable$M_G == T)
Mouse_count_matrix_sub_cl = Mouse_count_matrix_sub[k,]
Mouse_count_matrix_sub_cl = NormalizeData(Mouse_count_matrix_sub_cl)


k = which(rownames(LGS13_count_matrix_sub) %in% LGStoM_datatable$LGS_G == T)
LGS13_count_matrix_sub_cl = LGS13_count_matrix_sub[k,]
LGS13_count_matrix_sub_cl = NormalizeData(LGS13_count_matrix_sub_cl)

rownames(LGS13_count_matrix_sub_cl)

######
saveRDS(Mouse_count_matrix_sub_cl,file="Mouse_count_matrix_sub_cl_Jan2025")
saveRDS(LGS13_count_matrix_sub_cl,file="LGS13_count_matrix_sub_cl_Jan2025")

######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
FN = "Mouse_matrix_GRN_input_Jan2025.txt"
mat = Mouse_count_matrix_sub_cl[['RNA']]@data
rownames(mat)
mat = round(mat,5)
x = Matrix::t(mat)
write.table(as(x, "dgTMatrix"), file = FN, quote = FALSE, sep = "\t", row.names = FALSE)

Mouse_Corr = RNA_Corr_Add_cutoff(mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(Mouse_Corr,file="Mouse_Corr_LateGRNs_Jan2025.rds")

#########------------------------------------------------------------------------------------------------------------------------------------
RNA_Corr_Add_cutoff <- function(smooth_matrix){
	smooth_matrix_ori = smooth_matrix
	####
	k_index = which(smooth_matrix@x > 0)
	smooth_matrix_k = smooth_matrix
	smooth_matrix_k@x[k_index] = 1
	####
	k1 = which(Matrix::rowSums(smooth_matrix_k) < 25)
	k2 = which(is.na(Matrix::rowSums(smooth_matrix_k)) ==T)
	####
	genes_k1 = rownames(smooth_matrix_k)[k1]
	genes_k2 = rownames(smooth_matrix_k)[k2]
	####
	k_Cells = which(rownames(smooth_matrix) %in% c(genes_k1,genes_k2) == T)
	smooth_matrix = smooth_matrix[-k_Cells,]
	####
	print(dim(smooth_matrix))
	####
	Input = t(as.matrix(smooth_matrix))
	####
	Corr_res = sparse.cor3(Input)
	print(dim(Corr_res))
	####
	Corr_res[lower.tri(Corr_res,diag=T)]= 2
	####
	library(reshape2)
	Corr_res_out = melt(Corr_res)
	####
	k3 = which(Corr_res_out$value == 2)
	Corr_res_out_cl = Corr_res_out[-k3,]
	####
	#### Corr_res_out_cl[which(Corr_res_out_cl$Var1 == "Zic3"),]
	return(Corr_res_out_cl)
}


sparse.cor3 <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    cSums <- colSums(x)
    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp
    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}
#########------------------------------------------------------------------------------------------------------------------------------------

FN = "LGS13_matrix_GRN_input_Jan2025.txt"
mat = LGS13_count_matrix_sub_cl[['RNA']]@data
mat = round(mat,5)
rownames(mat)
x = Matrix::t(mat)
write.table(as(x, "dgTMatrix"), file = FN, quote = FALSE, sep = "\t", row.names = FALSE)

LGS13_Corr = RNA_Corr_Add_cutoff(mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Corr,file="LGS13_Corr_LateGRNs_Jan2025.rds")


#####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######
####### Next we will output the TFs name!!! #####----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### this step is to add the names of the TFs in GRNs steps ######
#####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

LGStoM_datatable$TF_index = "NO"
k = which(LGStoM_datatable$M_G %in% TFs == T)
LGStoM_datatable$TF_index[k] = "YES"

######

Mouse_TFs_GRN = LGStoM_datatable$M_G[which(LGStoM_datatable$TF_index == "YES")]
Mouse_TFs_GRN = Mouse_TFs_GRN[!duplicated(Mouse_TFs_GRN)]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
write.table(Mouse_TFs_GRN, file = 'Mouse_TFs_GRN.txt', quote = FALSE, sep = "\t", row.names = FALSE, col.names=F)


LGS13_TFs_GRN = LGStoM_datatable$LGS_G[which(LGStoM_datatable$TF_index == "YES")]
LGS13_TFs_GRN = LGS13_TFs_GRN[!duplicated(LGS13_TFs_GRN)]
write.table(LGS13_TFs_GRN, file = 'LGS13_TFs_GRN.txt', quote = FALSE, sep = "\t", row.names = FALSE, col.names=F)


#######-------------------------------------------------------------------------------------------------------------------------------------
#######------------------#######------------------#######------------------#######------------------#######------------------
#######-------------------------------------------------------------------------------------------------------------------------------------

####### Let us run the GRNs script !!!! #########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate /home/plyu3/anaconda3_new/envs/arboreto-env

python

import os
import pandas as pd
from distributed import Client, LocalCluster
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

import pandas as pd
import arboreto
import dask
from dask.distributed import Client
from distributed import LocalCluster, Client


os.chdir("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
ex_matrix = pd.read_csv("LGS13_matrix_GRN_input_Jan2025.txt", sep='\t')

### ex_matrix = ex_matrix.round(3)

col_sums = ex_matrix.sum(axis=0)
ex_matrix = ex_matrix.loc[:, col_sums > 10]

ex_matrix.columns = [col.replace('.', '-') for col in ex_matrix.columns]

tf_names = load_tf_names("LGS13_TFs_GRN.txt")

# instantiate a custom Dask distributed Client
client = Client(LocalCluster())
out_file = "LGS13_total_network_Jan2025.tsv"
# compute the GRN
network = grnboost2(expression_data=ex_matrix,
                        tf_names=tf_names,
                        client_or_address=client)
# write the GRN to file
network.to_csv(out_file, sep='\t', index=False, header=False)
#
print("Done!")

##### 在LGS13 的 GRNs 里面 gene 的名字 #######
##### 查看是否有. 或者是 - ###################



####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

os.chdir("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
ex_matrix = pd.read_csv("Mouse_matrix_GRN_input_Jan2025.txt", sep='\t')

### ex_matrix = ex_matrix.round(3)

col_sums = ex_matrix.sum(axis=0)
ex_matrix = ex_matrix.loc[:, col_sums > 10]
tf_names = load_tf_names("Mouse_TFs_GRN.txt")

# instantiate a custom Dask distributed Client
client = Client(LocalCluster())
out_file = "Mouse_total_network_Jan2025.tsv"
# compute the GRN
network = grnboost2(expression_data=ex_matrix,
                        tf_names=tf_names,
                        client_or_address=client)
# write the GRN to file
network.to_csv(out_file, sep='\t', index=False, header=False)
#
print("Done!")

####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####### #######
####### #######
####### #######


Mouse_total_GRNs = read.table("Mouse_total_network_Jan2025.tsv",sep="\t")
LGS13_total_GRNs = read.table("LGS13_total_network_Jan2025.tsv",sep="\t")

LGS13_total_GRNs[which(LGS13_total_GRNs$V1 == "Zic3" & LGS13_total_GRNs$V2 == "Thrb"),]
#### convert .like to -like #####
#### convert .containing to -containing #########
###LGS13_total_GRNs$V1 = gsub(".like","-like",LGS13_total_GRNs$V1)
###LGS13_total_GRNs$V1 = gsub(".containing","-containing",LGS13_total_GRNs$V1)
###LGS13_total_GRNs$V2 = gsub(".like","-like",LGS13_total_GRNs$V2)
###LGS13_total_GRNs$V2 = gsub(".containing","-containing",LGS13_total_GRNs$V2)

grep("-like",LGS13_total_GRNs$V1)
grep("-like",LGS13_total_GRNs$V2)



Mouse_total_GRNs_Corr = readRDS("Mouse_Corr_LateGRNs_Jan2025.rds")
LGS13_total_GRNs_Corr = readRDS("LGS13_Corr_LateGRNs_Jan2025.rds")

######
LGS13_total_GRNs_Corr[which(LGS13_total_GRNs_Corr$Var1 == "Zic3" & LGS13_total_GRNs_Corr$Var2 == "Thrb"),]

######
GRNs = Mouse_total_GRNs
Cor = Mouse_total_GRNs_Corr

Add_cor_to_GRNs <- function(GRNs,Cor){
    GRNs$V1 = gsub(".like","-like",GRNs$V1)
    GRNs$V2 = gsub(".like","-like",GRNs$V2)
    GRNs$V1 = gsub(".containing","-containing",GRNs$V1)
    GRNs$V2 = gsub(".containing","-containing",GRNs$V2)
    ########
    index1 = paste(GRNs$V1,GRNs$V2)
    ########
    index_Cor1 =  paste(Cor$Var1,Cor$Var2)
    index_Cor2 =  paste(Cor$Var2,Cor$Var1)
    ########
    m1 = match(index1,index_Cor1)
    m2 = match(index1,index_Cor2)
    ########
    function_1 <- function(x){
        x1 = x[1]
        x2 = x[2]
        #####
        if(is.na(x1) == T){
            return(x2)
        }
        if(is.na(x2) == T){
            return(x1)
        }
    }
    ########
    tab = data.frame(m1=m1,m2=m2)
    tab$res = apply(tab,1,function_1)
    ########
    GRNs$Corr = Cor$value[tab$res]
    ########
    return(GRNs)
    ########   
}

############
############
Mouse_total_GRNs_add = Add_cor_to_GRNs(Mouse_total_GRNs,Mouse_total_GRNs_Corr)
saveRDS(Mouse_total_GRNs_add,file="Mouse_total_GRNs_add_Jan2025")

LGS13_total_GRNs_add = Add_cor_to_GRNs(LGS13_total_GRNs,LGS13_total_GRNs_Corr)
saveRDS(LGS13_total_GRNs_add,file="LGS13_total_GRNs_add_Jan2025")

############
############

######-----#######----------------------------------------------------------------------------------------------------------------------------------------------------------
######-----#######----------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_13_Peak_to_Gene_merge <- readRDS(file="LGS_13_Peak_to_Gene_merge_2025")
Mouse_Peak_to_Gene_merge <- readRDS(file="Mouse_Peak_to_Gene_merge_2025")

table(LGS_13_Peak_to_Gene_merge$class)
table(Mouse_Peak_to_Gene_merge$class)

######
###### load the PtoG #######
######

Mouse_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Mouse_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")
LGS13_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Jan_2025")


Add_PtoG_and_filter <- function(Peak_to_Gene_merge,PtoG){
    ##########
    Peak_to_Gene_merge = Peak_to_Gene_merge
    k = which(Peak_to_Gene_merge$class == "NA")
    if(length(k) > 0){
        Peak_to_Gene_merge = Peak_to_Gene_merge[-k,]
    }
    ##########
    library(data.table)
    Peak_to_Gene_merge = data.table(Peak_to_Gene_merge)
    PtoG = data.table(data.frame(PtoG))
    PtoG = PtoG[,c("peak","gene","Correlation")]
    colnames(PtoG) <- c("peaks","gene","PtoG_Cor")
    ##########
    res <- merge(
        x = Peak_to_Gene_merge,
        y = PtoG,
        by = c("peaks", "gene"),
        all.x = TRUE
    )
    ##########
    k = which(res$class %in% c("Inter","Body") == T & is.na(res$PtoG_Cor) == T)
    res_cl = res[-k,]
    ########## res_cl[which(res_cl$gene == "Thrb"),]
    ########## res_cl[which(res_cl$gene == "Nrl"),]
    ########## res_cl[which(res_cl$gene == "Rho"),]
    return(res_cl)
}

####
####

Mouse_Peak_to_Gene_add = Add_PtoG_and_filter(Mouse_Peak_to_Gene_merge,Mouse_PtoG)
LGS13_Peak_to_Gene_add = Add_PtoG_and_filter(LGS_13_Peak_to_Gene_merge,LGS13_PtoG)

LGS13_Peak_to_Gene_add$gene

saveRDS(Mouse_Peak_to_Gene_add,file="Mouse_Peak_to_Gene_add_Jan2025.rds") 
saveRDS(LGS13_Peak_to_Gene_add,file="LGS13_Peak_to_Gene_add_Jan2025.rds") 

########
Mouse_Motif_matchtopeaks <- readRDS("Mouse_Motif_matchtopeaks_Jan13.rds")
LGS13_Motif_matchtopeaks <- readRDS("LGS13_Motif_matchtopeaks_Jan13.rds")

length(Mouse_Motif_matchtopeaks)


####### ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### process motifs #####------------------------------------------------------------------
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
load("PWM_list_combine_cl")
motifs = names(PWM_list_combine_cl)
genes = as.character(sapply(PWM_list_combine_cl,function(x) x@name))
dat = data.frame(motifs=motifs,tfs=genes)
dat$tfs[4558] = "Atoh7"
dat_list = apply(dat,1,split_fun)
########
split_fun <- function(x){
    #####
    motifs = x[1]
    tfs = x[2]
    tfs_sp = unlist(strsplit(as.character(tfs),split=";"))
    out = data.frame(motifs=motifs,tfs=tfs_sp)
    return(out)
}
dat_list = do.call("rbind",dat_list)
dat_list[4:10,]
######
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")
######
Mouse_motif_list = dat_list
Mouse_motif_list = Mouse_motif_list[which(Mouse_motif_list$tfs %in% squirrel_to_mouse_cl$Gene.name == T),]
###### let us get the LGS13 motif list #######

save(Mouse_motif_list,file="Mouse_motif_list_2025Mar")

########
########
########
########
footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Cone_footprint_res_2024")
TAG = "Cone"
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
Peak_to_Gene_add = readRDS("Mouse_Peak_to_Gene_add_Jan2025.rds")

########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
LGS_motif_tf_table = Mouse_motif_list
total_GRNs_add = readRDS("Mouse_total_GRNs_add_Jan2025")

########------------------------------------------------------------------------------------------------------------------------------------------------------------------
######## we will recalculate the Avg Exp for Mouse and LGS13 
########------------------------------------------------------------------------------------------------------------------------------------------------------------------

conda activate ArchR 
R
library(ArchR)
library(Seurat)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Downloading_datasets/Mouse_Development")
Mouse_count_matrix_total <- readRDS(file="Mouse_RNA_dev_raw_seurat_Anno_Jan29_2025.rds")

####### for LGS13:
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
LGS13_count_matrix_total = readRDS('all_seurat_RNA_big_cl_detail_new_202204')

####### ------------------------------------------------------------------------------------------------------------------------------------------------
####### 我们读入 combined seurat see the cell annotations #########
####### run the cell types "Cone","RPC","BC_Photo_pre","Rod","NG" ########
####### ------------------------------------------------------------------------------------------------------------------------------------------------

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')
table(LGS13_Mouse.combined$celltype)

Need_celltypes = c("Cone","RPC","BC_Photo_pre","Rod","NG")
Mouse_seurat = LGS13_Mouse.combined[,LGS13_Mouse.combined$sp == "Mouse"]
LGS13_seurat = LGS13_Mouse.combined[,LGS13_Mouse.combined$sp == "LGS"]

Mouse_Index = LGS13_Mouse.combined[,which(LGS13_Mouse.combined$celltype %in% Need_celltypes == T & LGS13_Mouse.combined$sp == "Mouse")]
LGS13_Index = LGS13_Mouse.combined[,which(LGS13_Mouse.combined$celltype %in% Need_celltypes == T & LGS13_Mouse.combined$sp == "LGS")]


######## Get Mouse GRNs input #####
Mouse_Index$cellid = gsub("Mouse_","",colnames(Mouse_Index))
Mouse_count_matrix_sub = Mouse_count_matrix_total[,which(colnames(Mouse_count_matrix_total) %in% Mouse_Index$cellid == T)]
m = match(colnames(Mouse_count_matrix_sub),Mouse_Index$cellid)
Mouse_count_matrix_sub$CT = as.character(Mouse_Index$celltype[m])
table(Mouse_count_matrix_sub$CT)


######## Get LGS13 GRNs input #####
LGS13_Index$cellid = gsub("LGS_","",colnames(LGS13_Index))
LGS13_count_matrix_sub = LGS13_count_matrix_total[,which(colnames(LGS13_count_matrix_total) %in% LGS13_Index$cellid == T)]
m = match(colnames(LGS13_count_matrix_sub),LGS13_Index$cellid)
LGS13_count_matrix_sub$CT = as.character(LGS13_Index$celltype[m])
table(LGS13_count_matrix_sub$CT)

########
######## ######## ####### Next we will filter genes ######## ######## ###########
########
library(Seurat)
########

DefaultAssay(Mouse_count_matrix_sub) <- "RNA"
DefaultAssay(LGS13_count_matrix_sub) <- "RNA"



Get_Avg_Exp <- function(Mat,Meta,Celltypes_Need){
    #########
    all_avg = c()
    for(i in 1:length(Celltypes_Need)){
        #########
        tmp_ct = Celltypes_Need[i]
        print(tmp_ct)
        k = which(Meta$index == tmp_ct)
        tmp_mat = Mat[,k]
        ##########
        tmp_sum = apply(tmp_mat,1,sum)
        ##########
        all_avg <- c(all_avg,tmp_sum)
    }
    ####
    all_avg_mat = matrix(all_avg,ncol=length(Celltypes_Need))
    rownames(all_avg_mat) = rownames(Mat)
    colnames(all_avg_mat) = Celltypes_Need
    ####
    factor = colSums(all_avg_mat) / 1e5
    all_avg_mat = sweep(all_avg_mat,2,factor,FUN="/")
    all_avg_mat_log = log2(all_avg_mat+1)
    ####
    return(all_avg_mat_log)
}

######
###### update the avg by counts matrix ###### ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
######
Mouse_count_matrix_sub$index = Mouse_count_matrix_sub$CT
LGS13_count_matrix_sub$index = LGS13_count_matrix_sub$CT

table(Mouse_count_matrix_sub$CT)
table(LGS13_count_matrix_sub$CT)

Celltypes_Need = names(table(Mouse_count_matrix_sub$CT))

LGS13_Late_RNA_Avg = Get_Avg_Exp(Mat = LGS13_count_matrix_sub[['RNA']]@counts,Meta = LGS13_count_matrix_sub@meta.data,Celltypes_Need=Celltypes_Need)
Mouse_Late_RNA_Avg = Get_Avg_Exp(Mat = Mouse_count_matrix_sub[['RNA']]@counts,Meta = Mouse_count_matrix_sub@meta.data,Celltypes_Need=Celltypes_Need)

head(LGS13_Late_RNA_Avg)
LGS13_Late_RNA_Avg[which(rownames(LGS13_Late_RNA_Avg) == "Zic3"),]
Mouse_Late_RNA_Avg[which(rownames(Mouse_Late_RNA_Avg) == "Zic3"),]

LGS13_Late_RNA_Avg[which(rownames(LGS13_Late_RNA_Avg) == "Thrb"),]
Mouse_Late_RNA_Avg[which(rownames(Mouse_Late_RNA_Avg) == "Thrb"),]


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Late_RNA_Avg,file="LGS13_Late_RNA_Avg_Mar2025")
saveRDS(Mouse_Late_RNA_Avg,file="Mouse_Late_RNA_Avg_Mar2025")

rownames(LGS13_Late_RNA_Avg)[grep("-like",rownames(LGS13_Late_RNA_Avg))]
rownames(LGS13_Late_RNA_Avg)[grep("-containing",rownames(LGS13_Late_RNA_Avg))]

########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")

########------------------------------------------------------------------------------------------------------------------------------------------------------------------

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Cone_footprint_res_2024")
TAG = "Cone"
Peak_to_Gene_add = readRDS("Mouse_Peak_to_Gene_add_Jan2025.rds")
LGS_motif_tf_table = Mouse_motif_list
total_GRNs_add = readRDS("Mouse_total_GRNs_add_Jan2025")
Avg_Mat = readRDS("Mouse_Late_RNA_Avg_Mar2025")

Mouse_Cone_GRNs = Add_foot_print_to_peak(footprint,"Cone",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------
########------------------------------------------------------------------------------------------------------------------------------------------------------------------

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Rod_footprint_res_2024")
TAG = "Rod"
Mouse_Rod_GRNs = Add_foot_print_to_peak(footprint,"Rod",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Photo_BC_footprint_res_2024")
TAG = "BC_Photo_pre"
Mouse_Photopre_GRNs = Add_foot_print_to_peak(footprint,"BC_Photo_pre",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Late_RPCs_footprint_res_2024")
TAG = "RPC"
Mouse_RPC_GRNs = Add_foot_print_to_peak(footprint,"RPC",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Late_NG_footprint_res_2024")
TAG = "NG"
Mouse_NG_GRNs = Add_foot_print_to_peak(footprint,"NG",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(Mouse_Cone_GRNs,"Mouse_Cone_GRNs_Mar2025")
saveRDS(Mouse_Rod_GRNs,"Mouse_Rod_GRNs_Mar2025")
saveRDS(Mouse_Photopre_GRNs,"Mouse_Photopre_GRNs_Mar2025")
saveRDS(Mouse_RPC_GRNs,"Mouse_RPC_GRNs_Mar2025")
saveRDS(Mouse_NG_GRNs,"Mouse_NG_GRNs_Mar2025")

########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------

Add_foot_print_to_peak <- function(footprint,TAG,Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat){
    ##########
    library(data.table)
    ########## filter the footprint #####
    left_delta = footprint$left - footprint$score
    right_delta = footprint$right - footprint$score
    filter_foot_index = which(left_delta > 0.05 & right_delta > 0.05)
    ##########
    footprint_cl = footprint[filter_foot_index]
    ########## first we will filter motif by corr #####
    footprint_region = as.character(footprint_cl)
    footprint_name = names(footprint_cl)
    footprint_table = data.frame(Motif=footprint_name,footprint=footprint_region)
    footprint_table = data.table(footprint_table)
    ########## merge footprint with peak ####
    findOverlaps_res = findOverlaps(GRanges(footprint_table$footprint),GRanges(Peak_to_Gene_add$peaks))
    ##########
    left_table = footprint_table[queryHits(findOverlaps_res),]
    right_table = Peak_to_Gene_add[subjectHits(findOverlaps_res),]
    merge_table = cbind(left_table,right_table)
    ########## Next merge tf motifs with footprint ####
    LGS_motif_tf_table_cl = LGS_motif_tf_table[,c("motifs","tfs")]
    colnames(LGS_motif_tf_table_cl) = c("Motif","TF")
    LGS_motif_tf_table_cl = data.table(LGS_motif_tf_table_cl)
    ##########
    merge_table2 <- merge(
        x = merge_table,
        y = LGS_motif_tf_table_cl,
        by = c("Motif"),
        allow.cartesian = TRUE
    )
    ##########
    ########## Next we will merge the gene-gene correaltion and importance ####
    total_GRNs_add = data.table(total_GRNs_add)
    ##########
    colnames(total_GRNs_add) = c("TF","gene","GG_importance","GG_Corr")
    ##########
    ##########
    merge_table3 <- merge(
        x = merge_table2,
        y = total_GRNs_add,
        by = c("TF","gene")
    )
    ########### Next we will add the average expression of the TFs ################
    Avg_mat_sub = data.frame(TF=rownames(Avg_Mat),Exp=Avg_Mat[,TAG])
    ###########
    m = match(merge_table3$TF,Avg_mat_sub$TF)
    merge_table3$TF_Exp = Avg_mat_sub$Exp[m]
    ###########
    merge_table3$celltype = TAG
    #####
    ##### first we will filter PtoG_Cor ########
    ##### only positive ####
    #####
    k1_filter = which(merge_table3$class %in% c("Inter","Body") == T & merge_table3$PtoG_Cor < 0)
    merge_table3_f1 = merge_table3[-k1_filter,]
    #####
    k2_filter = which(merge_table3_f1$GG_importance > 0.1 & abs(as.numeric(merge_table3_f1$GG_Corr)) > 0.05)
    merge_table3_f2 = merge_table3_f1[k2_filter,]
    ###### filter TF expression ######
    k3_filter = which(merge_table3_f2$TF_Exp > 0.1)
    merge_table3_f3 = merge_table3_f2[k3_filter,]
    ###### merge_table3_f3[which(merge_table3_f3$TF == "Zic3"),] ######
    ######
    ###### output the TF-peak-Gene network and TF-Gene network ####
    ###### add pos and neg ####
    merge_table3_f3$Class = "pos"
    k_neg = which(merge_table3_f3$GG_Corr < 0)
    merge_table3_f3$Class[k_neg] = 'neg'
    #######
    merge_table3_f3_3col = merge_table3_f3[,c("TF","peaks","gene","Class","celltype")]
    merge_table3_f3_2col = merge_table3_f3[,c("TF","gene","Class","celltype")]
    ###
    merge_table3_f3_3col_index = paste(merge_table3_f3_3col$TF,merge_table3_f3_3col$peaks,merge_table3_f3_3col$gene)
    merge_table3_f3_3col_cl = merge_table3_f3_3col[!duplicated(merge_table3_f3_3col_index),]
    ###
    merge_table3_f3_2col_index = paste(merge_table3_f3_2col$TF,merge_table3_f3_2col$gene)
    merge_table3_f3_2col_cl = merge_table3_f3_2col[!duplicated(merge_table3_f3_2col_index),]
    ####
    return(list(Ori=merge_table3_f3,TPG=merge_table3_f3_3col_cl,TG=merge_table3_f3_2col_cl))
}

########-----
########-------------------------------------------------------------------------------------------------------------------
########----- Next for the 13LGS !!! --------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-----

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
load("Mouse_motif_list_2025Mar")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")

########
LGS13_motif_list <- list()
for(i in 1:length(Mouse_motif_list$motifs)){
    print(i)
    ##########
    tmp_dat = Mouse_motif_list[i,]
    ##########
    tmp_motifs = tmp_dat$motifs
    tmp_tfs = tmp_dat$tfs
    ##########
    k = which(squirrel_to_mouse_cl$"Gene.name" == tmp_tfs)
    if(length(k) > 0){
        tmp_tfs_13LGS = squirrel_to_mouse_cl$"Squirrel.name"[k]
        ##########
        tmp_add = data.frame(motifs=tmp_motifs,tfs=tmp_tfs_13LGS,mouse_tfs=tmp_tfs)
        ##########
        LGS13_motif_list <- c(LGS13_motif_list,list(tmp_add))
    }
}
LGS13_motif_list = do.call("rbind",LGS13_motif_list)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_motif_list,file="LGS13_motif_list_2025Mar")

########
########-------------------- covert to the 13LGS --------------------######
########

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Cone_footprint_res_2024")
TAG = "Cone"
Peak_to_Gene_add = readRDS("LGS13_Peak_to_Gene_add_Jan2025.rds")
LGS_motif_tf_table = readRDS("LGS13_motif_list_2025Mar")
total_GRNs_add = readRDS("LGS13_total_GRNs_add_Jan2025")
Avg_Mat = readRDS("LGS13_Late_RNA_Avg_Mar2025")

head(Avg_Mat)

##-------------------------- ###### setcutoff_to_1 !!!!

Avg_Mat[which(rownames(Avg_Mat) == "OPN1MW"),]
Avg_Mat[which(rownames(Avg_Mat) == "Opn1sw"),]
Avg_Mat[which(rownames(Avg_Mat) == "Zic3"),]
Avg_Mat[which(rownames(Avg_Mat) == "Isl2"),]

########
########

Peak_to_Gene_add$gene[grep("-like",Peak_to_Gene_add$gene)]
Peak_to_Gene_add$gene[grep("-containing",Peak_to_Gene_add$gene)]

########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------

total_GRNs_add$V1[grep("-like",total_GRNs_add$V1)]
total_GRNs_add$V1[grep("-containing",total_GRNs_add$V1)]
total_GRNs_add$V2[grep("-like",total_GRNs_add$V2)]
total_GRNs_add$V2[grep("-containing",total_GRNs_add$V2)]

########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------
########-------------------------------------------------------------------------------------------------------------------

LGS13_Cone_GRNs = Add_foot_print_to_peak(footprint,"Cone",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Cone_GRNs,"LGS13_Cone_GRNs_Mar2025")

########-------------------------------------------------------------------------------------------------------------------

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_NG_footprint_res_2024")
TAG = "NG"
LGS13_NG_GRNs = Add_foot_print_to_peak(footprint,"NG",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_NG_GRNs,"LGS13_NG_GRNs_Mar2025")

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Photo_BC_footprint_res_2024")
TAG = "BC_Photo_pre"
LGS13_Photopre_GRNs = Add_foot_print_to_peak(footprint,"BC_Photo_pre",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Photopre_GRNs,"LGS13_Photopre_GRNs_Mar2025")

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_RPCs_footprint_res_2024")
TAG = "RPC"
LGS13_RPC_GRNs = Add_foot_print_to_peak(footprint,"RPC",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_RPC_GRNs,"LGS13_RPC_GRNs_Mar2025")

#########
#########
######### We will next to calcualte the GRNs !!! ################
#########
######### Next we will recalculate the markers !!!! ###################
#########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
LGS13_count_matrix_total = readRDS('all_seurat_RNA_big_cl_detail_new_202204')

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')
table(LGS13_Mouse.combined$celltype)

LGS13_Index = LGS13_Mouse.combined[,LGS13_Mouse.combined$sp == "LGS"]

LGS13_Index$cellid = gsub("LGS_","",colnames(LGS13_Index))
LGS13_count_matrix_sub = LGS13_count_matrix_total[,which(colnames(LGS13_count_matrix_total) %in% LGS13_Index$cellid == T)]
m = match(colnames(LGS13_count_matrix_sub),LGS13_Index$cellid)
LGS13_count_matrix_sub$CT = LGS13_Index$celltype[m]
table(LGS13_count_matrix_sub$CT)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
squirrel_to_mouse_cl <- readRDS(file="squirrel_to_mouse_cl_2025.rds")

k = which(rownames(LGS13_count_matrix_sub) %in% squirrel_to_mouse_cl$Squirrel.name == T)
LGS13_count_matrix_sub = LGS13_count_matrix_sub[k,]

DefaultAssay(LGS13_count_matrix_sub) = "RNA"
LGS13_count_matrix_sub = NormalizeData(LGS13_count_matrix_sub)

Idents(LGS13_count_matrix_sub) = 'CT'
LGS13_seurat_Cone_markers = FindMarkers(LGS13_count_matrix_sub,ident.1="Cone")
LGS13_seurat_Cone_markers['Zic3',]
LGS13_seurat_Cone_markers2 = LGS13_seurat_Cone_markers[which(LGS13_seurat_Cone_markers$avg_log2FC > 0.3 & LGS13_seurat_Cone_markers$p_val_adj < 0.01),]
saveRDS(LGS13_seurat_Cone_markers2,file="LGS13_seurat_Cone_markers_Jan2025")

LGS13_seurat_Cone_markers <- readRDS("LGS13_seurat_Cone_markers_Jan2025")

######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------














#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####-----------接下来我们需要 RUN 一下 GRNs --------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------
#####---------------------------------------------------------------------------------------------------------------------



















