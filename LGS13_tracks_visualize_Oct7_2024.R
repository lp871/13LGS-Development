#######
####### on the new server2 #######
####### on the new server2 #######
#######
####### 我们重新用下 小鼠的 ！！！ 在新的 server 上面 ！！！###
#######
#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

#####---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ssh plyu3@omb1.onc.jhmi.edu
njd$rft1
conda activate ArchR2
R

####### for mouse ########----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### 我们先搞一下 Mouse ########
####### setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")

setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
load("All_arrows_project_cl_2024")

library(ArchR)
getAvailableMatrices(All_arrows_project_cl)

All_arrows_project_cl$celltype <- factor(All_arrows_project_cl$celltype,levels=c("RPCs_S1","RPCs_S2","RPCs_S3","MG","Early_NG","Late_NG","RGC","AC_HC","Early_Cone","Cone","Early_Rod","Rod","BC"))
All_arrows_project_cl$celltype2 = All_arrows_project_cl$celltype
m = match(All_arrows_project_cl$celltype2,c("RPCs_S1","RPCs_S2","RPCs_S3","MG","Early_NG","Late_NG","RGC","AC_HC","Early_Cone","Cone","Early_Rod","Rod","BC"))
All_arrows_project_cl$celltype2 = paste0(m,"_",All_arrows_project_cl$celltype2)

########
######## recall peaks by celltype2 !!! #######
########

All_arrows_project_cl <- addGroupCoverages(ArchRProj = All_arrows_project_cl, groupBy = "celltype2",threads = 20)

pathToMacs2 <- findMacs2()

library(BSgenome.Mmusculus.UCSC.mm10)
All_arrows_project_cl <- addReproduciblePeakSet(
    ArchRProj = All_arrows_project_cl, 
    groupBy = "celltype2", 
    pathToMacs2 = pathToMacs2,
    genomeSize = 1.87e9
)

#####
#####
#####

All_arrows_project_cl <- addPeakMatrix(All_arrows_project_cl)

saveRDS(All_arrows_project_cl,file="All_arrows_project_cl_2024_Mouse")

All_arrows_project_cl <- readRDS("All_arrows_project_cl_2024_Mouse")

#####
##### next for the Peak matrix ######
#####
##### dazhonghua #####
#####




Get_Peak_Avg_matrix <- function(All_project_cl){
    #########
    peakMatrix <- getMatrixFromProject(All_project_cl, useMatrix = "PeakMatrix")
    Peak_matrix_count = peakMatrix@assays@data[[1]]
    Peak_matrix_row = as.character(peakMatrix@rowRanges)
    rownames(Peak_matrix_count) <- Peak_matrix_row
    #colnames(Peak_matrix_count)
    ########
    cellGroups <- All_project_cl$celltype2
    averagePeakByCellType <- sapply(unique(cellGroups), function(cluster) {
        # 找到属于当前 Cluster 的细胞索引
        cellIdx <- which(cellGroups == cluster)
        # 提取属于该 Cluster 的 Peak 值并计算平均值
        Matrix::rowSums(Peak_matrix_count[, cellIdx])
    })
    #####
    factor = colSums(averagePeakByCellType) /1e6
    averagePeakByCellType = sweep(averagePeakByCellType,2,factor,FUN="/")
    #####
    return(averagePeakByCellType)
}



Get_Gene_Avg_matrix <- function(All_project_cl){
    #########
    peakMatrix <- getMatrixFromProject(All_project_cl, useMatrix = "GeneIntegrationMatrix")
    Peak_matrix_count = peakMatrix@assays@data[[1]]
    Peak_matrix_row = as.character(peakMatrix@rowRanges)
    rownames(Peak_matrix_count) <- Peak_matrix_row
    #colnames(Peak_matrix_count)
    ########
    cellGroups <- All_project_cl$celltype2
    averagePeakByCellType <- sapply(unique(cellGroups), function(cluster) {
        # 找到属于当前 Cluster 的细胞索引
        cellIdx <- which(cellGroups == cluster)
        # 提取属于该 Cluster 的 Peak 值并计算平均值
        Matrix::rowSums(Peak_matrix_count[, cellIdx])
    })
    #####
    factor = colSums(averagePeakByCellType) /1e6
    averagePeakByCellType = sweep(averagePeakByCellType,2,factor,FUN="/")
    #####
    return(averagePeakByCellType)
}

####
####
setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Peak_Avg = Get_Peak_Avg_matrix(All_arrows_project_cl)
#### Peak_Avg = Get_Peak_Avg_matrix(All_project_cl)
saveRDS(Peak_Avg,file="Mouse_Peak_Avg_Oct14.rds")

Peak_Avg <- readRDS("Mouse_Peak_Avg_Oct14.rds")


####### OK!!! let us plot the tracks !!! ######
####### see the cell types: ######
#######
####### 目前最快的方法就是双拼了 ！！！！ ########
#######
#######
####### 我们重新 call 一下 PtoG 的links #########
#######

table(All_arrows_project_cl$celltype2)

####### 
cell_index_meta = All_arrows_project_cl@cellColData
cell_index_meta_sp = split(cell_index_meta,cell_index_meta$celltype2)

cell_list <- c()
for(i in 1:length(cell_index_meta_sp)){
    tmp = cell_index_meta_sp[[i]]
    tmp_cells = rownames(tmp)
    print(length(tmp_cells))
    ####
    if(length(tmp_cells) > 2500){
        ####
        tmp_cells = tmp_cells[sample(1:length(tmp_cells),2500)]
        print(length(tmp_cells))
    }
    cell_list <- c(cell_list,tmp_cells)
}

All_arrows_project_cl_PtoG_sub = All_arrows_project_cl[cell_list]

All_arrows_project_cl_PtoG_sub <- addPeak2GeneLinks(
    ArchRProj = All_arrows_project_cl_PtoG_sub,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)

#####----------------------------------------------------------------------------------------------------------------------------------------

saveRDS(All_arrows_project_cl_PtoG_sub,file="All_arrows_project_cl_PtoG_sub_Mouse_Oct16")

#####----------------------------------------------------------------------------------------------------------------------------------------------------------

Get_p2g_fun <- function(x,corCutOff,FDRCutOff){
	corCutOff = corCutOff
	FDRCutOff = FDRCutOff
	varCutOffATAC = 0.3
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



get_TSS <- function(project,gene = "Onecut1"){
    #####
    Genes = project@geneAnnotation$genes
    #####
    Genes_sub = Genes[which(Genes$symbol == gene)]
    #####
    if(as.character(strand(Genes_sub)) == "+"){
        end(Genes_sub) = start(Genes_sub)
    }
    if(as.character(strand(Genes_sub)) == "-"){
        start(Genes_sub) = end(Genes_sub)
    }
    strand(Genes_sub) = "*"
    return(Genes_sub)
}

construction_PtoG <- function(p2g,LGS13_p2g_tab_sub,TSS,need_peaks){
    LGS13_p2g_tab_sub = LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak %in% need_peaks == T),]
    ############
    print(dim(LGS13_p2g_tab_sub))
    GR = GRanges(LGS13_p2g_tab_sub$peak)
    Mid = round(start(GR) + end(GR))/2
    ############
    TSS_loci = start(TSS)
    chr = as.character(seqnames(TSS))
    ####
    vector <- c()
    for(i in 1:length(Mid)){
        #####
        if(Mid[i] >= TSS_loci){
            v = paste0(chr,":",TSS_loci,"-",Mid[i])
        }
        if(Mid[i] < TSS_loci){
            v = paste0(chr,":",Mid[i],"-",TSS_loci)
        }
        vector <- c(vector,v)
    }
    ######
    vector <- GRanges(vector)
    vector$value = LGS13_p2g_tab_sub$Correlation
    vector$FDR = LGS13_p2g_tab_sub$FDR
    #######
    p2g$Peak2GeneLinks = vector
    ###
    return(p2g)
}

PtoG_identify_and_Plot <- function(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat){
    ############
    ##### first call PtoG ########
    ############
    print("Get_PtoG_fun!")
    LGS13_p2g_tab = Get_p2g_fun(PtoG_ArchR,corCutOff=Corr_cutoff,FDRCutOff=FDR)
    ############
    LGS13_p2g_tab_sub = LGS13_p2g_tab[which(LGS13_p2g_tab$gene == Gene_Name),]
    length(LGS13_p2g_tab_sub$peak)
    length(which(LGS13_p2g_tab_sub$peak %in% rownames(Peak_Avg_Mat) == T))
    ############
    Need_Peaks = rownames(Peak_Avg_Mat)[which(apply(Peak_Avg_Mat,1,max) > 2)] 
    LGS13_p2g_tab_sub = LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak %in% Need_Peaks == T),]
    ############
    TSS = get_TSS(PtoG_ArchR,Gene_Name)
    ############
    p2g_EXAMPLE <- getPeak2GeneLinks(
        ArchRProj = PtoG_ArchR,
        corCutOff = -10,
        resolution = 50,
        FDRCutOff = 1e-02,
        returnLoops = TRUE,
        varCutOffATAC = 0.1,
        varCutOffRNA = 0.1
    )
    #######
    p2g = construction_PtoG(p2g_EXAMPLE,LGS13_p2g_tab_sub,TSS=TSS,need_peaks=Need_Peaks)
    #######
    print("PLOT!!!")
    p <- plotBrowserTrack(
	    ArchRProj = Ori_ArchR,
        plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
        sizes = c(8, 1.5, 5, 1.5),
	    groupBy = "celltype2", 
	    geneSymbol = Gene_Name, 
	    useMatrix = "GeneIntegrationMatrix",
	    upstream = UPregion, 
	    downstream = DOWNregion,
        loops = p2g
    )
    #######
    #######
    plotPDF(plotList = p, 
        name = FN, 
        ArchRProj = Ori_ArchR, 
        addDOC = FALSE, width = 12, height = 8
    )
    #########
    return(LGS13_p2g_tab_sub)
    #########
}

########
Ori_ArchR = All_arrows_project_cl
PtoG_ArchR = All_arrows_project_cl_PtoG_sub
Gene_Name = "Thrb"
UPregion = 50000
DOWNregion = 400000
Corr_cutoff = 0.15
FDR = 0.005
PeakAvg = 2
Peak_Avg_Mat <- readRDS("Mouse_Peak_Avg_Oct14.rds")
FN = "Mouse_Thrb_Oct15.pdf"
PtoG_Mouse_Thrb <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Thrb,file="PtoG_Mouse_Thrb_Oct16")

#########
#########

Gene_Name = "Onecut1"
UPregion = 100000
DOWNregion = 40000
FN = "Mouse_Onecut1_Oct15.pdf"
PtoG_Mouse_Onecut1 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Onecut1,file="PtoG_Mouse_Onecut1_Oct16")

#########
#########

Gene_Name = "Onecut2"
UPregion = 40000
DOWNregion = 100000
FN = "Mouse_Onecut2_Oct15.pdf"
PtoG_Mouse_Onecut2 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Onecut2,file="PtoG_Mouse_Onecut2_Oct16")

##########
##########

Gene_Name = "Zic3"
UPregion = 40000
DOWNregion = 100000
FN = "Mouse_Zic3_Oct15.pdf"
PtoG_Mouse_Zic3 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Zic3,file="PtoG_Mouse_Zic3_Oct16")


Gene_Name = "Sall3"
UPregion = 40000
DOWNregion = 100000
FN = "Mouse_Sall3_Oct15.pdf"
PtoG_Mouse_Sall3 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Sall3,file="PtoG_Mouse_Sall3_Oct16")



Gene_Name = "Mef2c"
UPregion = 40000
DOWNregion = 200000
FN = "Mouse_Mef2c_Oct15.pdf"
PtoG_Mouse_Mef2c <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Mef2c,file="PtoG_Mouse_Mef2c_Oct16")


Gene_Name = "Rxrg"
UPregion = 100000
DOWNregion = 100000
FN = "Mouse_Rxrg_Oct15.pdf"
PtoG_Mouse_Rxrg <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Rxrg,file="PtoG_Mouse_Rxrg_Oct16")



Gene_Name = "Pou2f1"
UPregion = 200000
DOWNregion = 100000
FN = "Mouse_Pou2f1_Oct15.pdf"
PtoG_Mouse_Pou2f1 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_Mouse_Pou2f1,file="PtoG_Mouse_Pou2f1_Oct16")

########
######## 我们下一步直接 saveRDS 把那个 13LGS的数据 #########
########










#####----------------------------------------------------------------------------------------------------------------------------------------------------------



####### for mouse ########----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####
#####-----------13LGS samples--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####
#### OK!!!!! Let us see the 13LGS samples !!!!!!! ######
##### for 13LGSs we will get the peaks #### !!!!! ######

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1
conda activate ArchR2
R

setwd("/zp1/data/plyu3/LGS_13_NEW")
All_project_cl <- readRDS("All_project_cl_2024_LGS")

library(ArchR)
getAvailableMatrices(All_project_cl)

table(All_project_cl$new_celltypes)

All_project_cl$celltype2 = All_project_cl$new_celltypes
m = match(All_project_cl$celltype2,c("Early_RPC","Late_RPC","MG","Early_NG","Late_NG","RGC","HC","AC","Photoprecursor","Cone","Rod","BC"))

All_project_cl$celltype2 = paste0(m,"_",All_project_cl$celltype2)

#### re-call peaks #####
####
All_project_cl <- addGroupCoverages(ArchRProj = All_project_cl, groupBy = "celltype2",threads = 20)

pathToMacs2 <- findMacs2()

All_project_cl <- addReproduciblePeakSet(
    ArchRProj = All_project_cl, 
    groupBy = "celltype2", 
    pathToMacs2 = pathToMacs2,
    genomeSize = 1.87e9
)
Creating Union Peak Set (395677) ！！！！！！！
#####
#####
All_project_cl <- addPeakMatrix(All_project_cl)

saveRDS(All_project_cl,file="All_project_cl_2024_LGS")

#####
##### for 13LGS project #####
#####


ssh plyu3@omb1.onc.jhmi.edu
njd$rft1
conda activate ArchR2
R
##### for LGS13 save RDSs ---------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/LGS_13_NEW")
All_project_cl <- readRDS("All_project_cl_2024_LGS")
All_arrows_project_cl_PtoG_sub <- readRDS("All_arrows_project_cl_PtoG_sub_Oct14.rds")
LGS13_Gene_Anno = All_project_cl@geneAnnotation$genes

k = grep("^LOC",LGS13_Gene_Anno$symbol)
LGS13_Gene_Anno = LGS13_Gene_Anno[-k]

k = grep("^ENSST",LGS13_Gene_Anno$symbol)
LGS13_Gene_Anno = LGS13_Gene_Anno[-k]

k = grep("^SNOR",LGS13_Gene_Anno$symbol)
LGS13_Gene_Anno = LGS13_Gene_Anno[-k]

head(LGS13_Gene_Anno,n=500)
saveRDS(LGS13_Gene_Anno,file="LGS13_Gene_Anno_Oct23")


#####
##### we will save the PtoGs #####
##### !!!!!! 

Get_p2g_fun_Plus <- function(x,corCutOff,FDRCutOff){
	corCutOff = corCutOff
	FDRCutOff = FDRCutOff
	varCutOffATAC = 0.3
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

library(ArchR)
PtoG_ArchR = All_arrows_project_cl_PtoG_sub
LGS13_p2g_tab = Get_p2g_fun_Plus(PtoG_ArchR,corCutOff=0.25,FDRCutOff=0.01)
saveRDS(LGS13_p2g_tab,file="LGS13_p2g_tab_Oct23")



##### for Mouse save RDSs ---------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_ArchR <- readRDS("All_arrows_project_cl_2024_Mouse")
Mouse_PtoG_ArchR <- readRDS("All_arrows_project_cl_PtoG_sub_Mouse_Oct16")
Mouse_Gene_Anno = Mouse_ArchR@geneAnnotation$genes
saveRDS(Mouse_Gene_Anno,file="Mouse_Gene_Anno_Oct23")


library(ArchR)
PtoG_ArchR = Mouse_PtoG_ArchR
Mouse_p2g_tab = Get_p2g_fun_Plus(PtoG_ArchR,corCutOff=0.25,FDRCutOff=0.01)
saveRDS(Mouse_p2g_tab,file="Mouse_p2g_tab_Oct23")

##### -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### we will write a function to output the Plot and return all the Peak to Gene links !!!!!! ########
#####
Ori_ArchR = All_project_cl
PtoG_ArchR = All_arrows_project_cl_PtoG_sub
Gene_Name = "Onecut1"
UPregion = 100000
DOWNregion = 40000
Corr_cutoff = 0.15
FDR = 0.005
PeakAvg = 2
Peak_Avg_Mat <- readRDS(file="LGS13_Peak_Avg_Oct14.rds")
FN = "13LGS_Onecut1_Oct15.pdf"

Get_p2g_fun <- function(x,corCutOff,FDRCutOff){
	corCutOff = corCutOff
	FDRCutOff = FDRCutOff
	varCutOffATAC = 0.3
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



get_TSS <- function(project,gene = "Onecut1"){
    #####
    Genes = project@geneAnnotation$genes
    #####
    Genes_sub = Genes[which(Genes$symbol == gene)]
    #####
    if(as.character(strand(Genes_sub)) == "+"){
        end(Genes_sub) = start(Genes_sub)
    }
    if(as.character(strand(Genes_sub)) == "-"){
        start(Genes_sub) = end(Genes_sub)
    }
    strand(Genes_sub) = "*"
    return(Genes_sub)
}

construction_PtoG <- function(p2g,LGS13_p2g_tab_sub,TSS,need_peaks){
    LGS13_p2g_tab_sub = LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak %in% need_peaks == T),]
    ############
    print(dim(LGS13_p2g_tab_sub))
    GR = GRanges(LGS13_p2g_tab_sub$peak)
    Mid = round(start(GR) + end(GR))/2
    ############
    TSS_loci = start(TSS)
    chr = as.character(seqnames(TSS))
    ####
    vector <- c()
    for(i in 1:length(Mid)){
        #####
        if(Mid[i] >= TSS_loci){
            v = paste0(chr,":",TSS_loci,"-",Mid[i])
        }
        if(Mid[i] < TSS_loci){
            v = paste0(chr,":",Mid[i],"-",TSS_loci)
        }
        vector <- c(vector,v)
    }
    ######
    vector <- GRanges(vector)
    vector$value = LGS13_p2g_tab_sub$Correlation
    vector$FDR = LGS13_p2g_tab_sub$FDR
    #######
    p2g$Peak2GeneLinks = vector
    ###
    return(p2g)
}

PtoG_identify_and_Plot <- function(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat){
    ############
    ##### first call PtoG ########
    ############
    print("Get_PtoG_fun!")
    LGS13_p2g_tab = Get_p2g_fun(PtoG_ArchR,corCutOff=Corr_cutoff,FDRCutOff=FDR)
    ############
    LGS13_p2g_tab_sub = LGS13_p2g_tab[which(LGS13_p2g_tab$gene == Gene_Name),]
    length(LGS13_p2g_tab_sub$peak)
    length(which(LGS13_p2g_tab_sub$peak %in% rownames(Peak_Avg_Mat) == T))
    ############
    Need_Peaks = rownames(Peak_Avg_Mat)[which(apply(Peak_Avg_Mat,1,max) > 2)] 
    LGS13_p2g_tab_sub = LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak %in% Need_Peaks == T),]
    ############
    TSS = get_TSS(PtoG_ArchR,Gene_Name)
    ############
    p2g_EXAMPLE <- getPeak2GeneLinks(
        ArchRProj = PtoG_ArchR,
        corCutOff = -10,
        resolution = 50,
        FDRCutOff = 1e-02,
        returnLoops = TRUE,
        varCutOffATAC = 0.1,
        varCutOffRNA = 0.1
    )
    #######
    p2g = construction_PtoG(p2g_EXAMPLE,LGS13_p2g_tab_sub,TSS=TSS,need_peaks=Need_Peaks)
    #######
    print("PLOT!!!")
    p <- plotBrowserTrack(
	    ArchRProj = Ori_ArchR,
        plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
        sizes = c(8, 1.5, 5, 1.5),
	    groupBy = "celltype2", 
	    geneSymbol = Gene_Name, 
	    useMatrix = "GeneIntegrationMatrix",
	    upstream = UPregion, 
	    downstream = DOWNregion,
        loops = p2g
    )
    #######
    #######
    plotPDF(plotList = p, 
        name = FN, 
        ArchRProj = Ori_ArchR, 
        addDOC = FALSE, width = 12, height = 8
    )
    #########
    return(LGS13_p2g_tab_sub)
    #########
}


Ori_ArchR = All_project_cl
PtoG_ArchR = All_arrows_project_cl_PtoG_sub
Gene_Name = "Onecut1"
UPregion = 100000
DOWNregion = 40000
Corr_cutoff = 0.15
FDR = 0.005
PeakAvg = 2
Peak_Avg_Mat <- readRDS(file="LGS13_Peak_Avg_Oct14.rds")
FN = "13LGS_Onecut1_Oct15.pdf"
PtoG_13LGS_Onecut1 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Onecut1,file="PtoG_13LGS_Onecut1_Oct16")

Gene_Name = "Onecut2"
UPregion = 30000
DOWNregion = 40000
FN = "13LGS_Onecut2_Oct15.pdf"
PtoG_13LGS_Onecut2 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Onecut2,file="PtoG_13LGS_Onecut2_Oct16")

####### ################## ######
####### ######## ####### ########

Gene_Name = "Zic3"
UPregion = 30000
DOWNregion = 40000
FN = "13LGS_Zic3_Oct15.pdf"
PtoG_13LGS_Zic3 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Zic3,file="PtoG_13LGS_Zic3_Oct16")

Gene_Name = "Sall3"
UPregion = 30000
DOWNregion = 40000
FN = "13LGS_Sall3_Oct15.pdf"
PtoG_13LGS_Sall3 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Sall3,file="PtoG_13LGS_Sall3_Oct16")

Gene_Name = "Mef2c"
UPregion = 70000
DOWNregion = 70000
FN = "13LGS_Mef2c_Oct15.pdf"
PtoG_13LGS_Mef2c <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Mef2c,file="PtoG_13LGS_Mef2c_Oct16")

Gene_Name = "Thrb"
UPregion = 300000
DOWNregion = 40000
FN = "13LGS_Thrb_Oct15.pdf"
PtoG_13LGS_Thrb <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Thrb,file="PtoG_13LGS_Thrb_Oct16")

Gene_Name = "Rxrg"
UPregion = 40000
DOWNregion = 40000
FN = "13LGS_Rxrg_Oct15.pdf"
PtoG_13LGS_Rxrg <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Rxrg,file="PtoG_13LGS_Rxrg_Oct16")

Gene_Name = "Pou2f1"
UPregion = 140000
DOWNregion = 40000
FN = "13LGS_Pou2f1_Oct15.pdf"
PtoG_13LGS_Pou2f1 <- PtoG_identify_and_Plot(Ori_ArchR,PtoG_ArchR,Gene_Name,UPregion,DOWNregion,Corr_cutoff,FDR,PeakAvg,Peak_Avg_Mat)
saveRDS(PtoG_13LGS_Pou2f1,file="PtoG_13LGS_Pou2f1_Oct16")

#########------------------------------------------------------------------------------------------------------------------------------------
#########------------------------------------------------------------------------------------------------------------------------------------

cell_index_meta = All_project_cl@cellColData
cell_index_meta_sp = split(cell_index_meta,cell_index_meta$celltype2)

cell_list <- c()
for(i in 1:length(cell_index_meta_sp)){
    tmp = cell_index_meta_sp[[i]]
    tmp_cells = rownames(tmp)
    print(length(tmp_cells))
    ####
    if(length(tmp_cells) > 2500){
        ####
        tmp_cells = tmp_cells[sample(1:length(tmp_cells),2500)]
        print(length(tmp_cells))
    }
    cell_list <- c(cell_list,tmp_cells)
}

All_arrows_project_cl_PtoG_sub = All_project_cl[cell_list]

All_arrows_project_cl_PtoG_sub <- addPeak2GeneLinks(
    ArchRProj = All_arrows_project_cl_PtoG_sub,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)


#######
####### 我们先 Plot 看一下这个东西 ############# ----------------------------------------------------------------------------------------------------------------------------------------
#######

#######
#######

setwd("/zp1/data/plyu3/LGS_13_NEW")
saveRDS(All_arrows_project_cl_PtoG_sub,file="All_arrows_project_cl_PtoG_sub_Oct14.rds")

LGS13_p2g_tab = Get_p2g_fun(All_arrows_project_cl_PtoG_sub)
####### LGS13_p2g_tab_sub = LGS13_p2g_tab[which(LGS13_p2g_tab$gene == "Onecut1"),]
LGS13_p2g_tab_sub = LGS13_p2g_tab[which(LGS13_p2g_tab$gene == "Mef2c"),]

LGS13_p2g_tab_sub$peak
####### LGS13_p2g_tab_sub = data.frame(LGS13_p2g_tab_sub)
LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak == "Itri5:96666568-96667068"),]

length(LGS13_p2g_tab_sub$peak)
length(which(LGS13_p2g_tab_sub$peak %in% rownames(Peak_Avg) == T))

#####
#####----------计算下peak在每个celltype的平均表达量 #########
#####

length(LGS13_p2g_tab$peak)
length(which(LGS13_p2g_tab$peak %in% rownames(Peak_Avg) == T))

Get_Peak_Avg_matrix <- function(All_project_cl){
    #########
    peakMatrix <- getMatrixFromProject(All_project_cl, useMatrix = "PeakMatrix")
    Peak_matrix_count = peakMatrix@assays@data[[1]]
    Peak_matrix_row = as.character(peakMatrix@rowRanges)
    rownames(Peak_matrix_count) <- Peak_matrix_row
    #colnames(Peak_matrix_count)
    ########
    cellGroups <- All_project_cl$celltype2
    averagePeakByCellType <- sapply(unique(cellGroups), function(cluster) {
        # 找到属于当前 Cluster 的细胞索引
        cellIdx <- which(cellGroups == cluster)
        # 提取属于该 Cluster 的 Peak 值并计算平均值
        Matrix::rowSums(Peak_matrix_count[, cellIdx])
    })
    #####
    factor = colSums(averagePeakByCellType) /1e6
    averagePeakByCellType = sweep(averagePeakByCellType,2,factor,FUN="/")
    #####
    return(averagePeakByCellType)
}

Peak_Avg = readRDS("LGS13_Peak_Avg_Oct14.rds")
saveRDS(Peak_Avg,file="LGS13_Peak_Avg_Oct14.rds")


[1] 395677
dim(Peak_Avg)
[1] 395677


######
Peak_Avg[which(rownames(Peak_Avg) == "Itri5:96666568-96667068"),]
###### Peak_Avg[which(rownames(Peak_Avg) == "Itri5:96667938-96668438"),]
Peak_Avg[which(rownames(Peak_Avg) == "Itri5:96665923-96666423"),]
###### Peak to bed ######

Peak_GR <- GRanges(rownames(Peak_Avg))
write.table(Peak_GR,file="LGS13_all_peaks.bed",sep="\t",quote=F,col.names=F,row.names=F)

######
###### select expressed peaks in these cell types #####
######

big = apply(Peak_Avg, 1 ,max)
k = which(big > 2)
need_peaks = rownames(Peak_Avg)[k]

#####
##### 这个PtoG一端是TSS, 另外一段是peak， valueandFDR #########
##### 构建一下 然后替换 ！！！！ #############
#####

project = All_arrows_project_cl_PtoG_sub

get_TSS <- function(project,gene = "Onecut1"){
    #####
    Genes = project@geneAnnotation$genes
    #####
    Genes_sub = Genes[which(Genes$symbol == gene)]
    #####
    if(as.character(strand(Genes_sub)) == "+"){
        end(Genes_sub) = start(Genes_sub)
    }
    if(as.character(strand(Genes_sub)) == "-"){
        start(Genes_sub) = end(Genes_sub)
    }
    strand(Genes_sub) = "*"
    return(Genes_sub)
}

###### #########
Onecut1_TSS = get_TSS(All_arrows_project_cl_PtoG_sub,"Mef2c")
TSS = Onecut1_TSS
###### #########

construction_PtoG <- function(p2g,LGS13_p2g_tab_sub,TSS,need_peaks){
    LGS13_p2g_tab_sub = LGS13_p2g_tab_sub[which(LGS13_p2g_tab_sub$peak %in% need_peaks == T),]
    ############
    print(dim(LGS13_p2g_tab_sub))
    GR = GRanges(LGS13_p2g_tab_sub$peak)
    Mid = round(start(GR) + end(GR))/2
    ############
    TSS_loci = start(TSS)
    chr = as.character(seqnames(TSS))
    ####
    vector <- c()
    for(i in 1:length(Mid)){
        #####
        if(Mid[i] >= TSS_loci){
            v = paste0(chr,":",TSS_loci,"-",Mid[i])
        }
        if(Mid[i] < TSS_loci){
            v = paste0(chr,":",Mid[i],"-",TSS_loci)
        }
        vector <- c(vector,v)
    }
    ######
    vector <- GRanges(vector)
    vector$value = LGS13_p2g_tab_sub$Correlation
    vector$FDR = LGS13_p2g_tab_sub$FDR
    #######
    p2g$Peak2GeneLinks = vector
    ###
    return(p2g)
}

p2g = construction_PtoG(p2g,LGS13_p2g_tab_sub,TSS,need_peaks)

p <- plotBrowserTrack(
	ArchRProj = All_project_cl,
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 1.5, 5, 1.5),
	groupBy = "celltype2", 
	geneSymbol = "Mef2c", 
	useMatrix = "GeneIntegrationMatrix",
	upstream = 100000, 
	downstream = 240000,
    loops = p2g
)

######
###### Body peaks
######

plotPDF(plotList = p, 
    name = "13LGS_Mef2c_Oct14.pdf", 
    ArchRProj = All_project_cl, 
    addDOC = FALSE, width = 12, height = 8
)


#####----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#####-----------13LGS samples--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#### ######
#### ######
#### ######
#### ######
#### ######






####
####
####


######
###### Old genes ##### 
######

geneSymbol = "Thrb"
upstream = 300000
downstream = 40000

geneSymbol = "Pou2f1"
upstream = 100000
downstream = 40000

geneSymbol = "Onecut1" 
upstream = 100000
downstream = 40000

geneSymbol = "Onecut2", 
upstream = 100000,
downstream = 40000,

geneSymbol = "Zic3", 
upstream = 20000,
downstream = 20000,

geneSymbol = "Sall3", 
upstream = 100000,
downstream = 50000

geneSymbol = "Mef2c"
upstream = 300000
downstream = 40000
######



load('Combine_proj_cl_202009')


p2g <- metadata(Combine_proj_cl@peakSet)$Peak2GeneLinks

Get_p2g_fun <- function(x){
	corCutOff = 0.20
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



####
################ we will find the conserved peaks ############
####

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd('/zp1/data/plyu3/LGS_13/Convert')
files = list.files()
files_need = files[grep('satsuma_summary.refined',files)]

####
####

table_list = list()

for(i in 1:length(files_need)){
	print(i)
	tmp_file = read.table(files_need[i])
	table_list = c(table_list,list(tmp_file))
}

head(table_list[[1]])


table_list_combined = do.call('rbind',table_list)
#diff1 = table_list_combined$V3 - table_list_combined$V2
#diff2 = table_list_combined$V6 - table_list_combined$V5

index = table_list_combined$V1
index_sp = strsplit(index, split=':')

chr = sapply(index_sp,function(x) x[[4]])
chr = paste('chr',chr,sep='')
start = table_list_combined$V2
end = table_list_combined$V3

library(GenomicRanges)
Mouse_alignment = GRanges(chr,IRanges(start,end),score=table_list_combined$V7)
LGS_alignment = GRanges(seqnames=table_list_combined$V4,IRanges(table_list_combined$V5,table_list_combined$V6),strand=table_list_combined$V8,score=table_list_combined$V7)

Mouse_alignment_bed <- data.frame(Mouse_alignment)
Mouse_alignment_wig <- Mouse_alignment_bed[,c(1,2,3,6)]
write.table(Mouse_alignment_wig,file="Mouse_alignment.wig",sep="\t",quote=F,row.names=F,col.names=F)


LGS13_alignment_bed <- data.frame(LGS_alignment)
LGS13_alignment_wig <- LGS13_alignment_bed[,c(1,2,3,6)]
write.table(LGS13_alignment_wig,file="LGS13_alignment.wig",sep="\t",quote=F,row.names=F,col.names=F)

setwd('/zp1/data/plyu3/LGS_13/Convert')
saveRDS(Mouse_alignment,file="Mouse_alignment_Oct17")
saveRDS(LGS_alignment,file="LGS_alignment_Oct17")


###### next we will load the mouse genome #######
######

library(BSgenome.Mmusculus.UCSC.mm10)
seq_mouse <- getSeq(BSgenome.Mmusculus.UCSC.mm10, Mouse_alignment[1])
seq_LGS13 <- getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, LGS_alignment[1])

library(BSgenome.Mmusculus.UCSC.mm10)
library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")

#start(LGS_alignment) = start(LGS_alignment)-1
#end(LGS_alignment) = end(LGS_alignment)-1
#tmp = LGS_alignment[1428]
#tmp1 = tmp
#strand(tmp1) = "+"
#seq_LGS13_p <- getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, tmp1)
#seq_LGS13_m <- getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, tmp)
#alignment1 <- pairwiseAlignment(seq_LGS13_p, seq_mouse, type = "local")
#alignment2 <- pairwiseAlignment(seq_LGS13_m, seq_mouse, type = "local")

seq2_rc <- reverseComplement(seq_LGS13_p)
alignment3 <- pairwiseAlignment(seq2_rc, seq_mouse, type = "local")
alignment <- pairwiseAlignment(seq_LGS13, seq_mouse, type = "local")

score(alignment)

start_pattern <- start(subject(alignment3))
end_pattern <- end(subject(alignment3))      
start_subject <- start(pattern(alignment3))  
end_subject <- end(pattern(alignment3))   

length_pattern <- end_pattern - start_pattern + 1
length_subject <- end_subject - start_subject + 1

####
#### how to know the overlap regions length #####
####
alignment <- pairwiseAlignment(seq_mouse, seq_mouse, type = "global")

####
#### 可以没问题 ！！！ #####
#### 搞定！！！####
####



#### Next we will compare the similarity between LGS and Mouse ########
####
#### Let us output the genome conserve sequence ######
####

#### load the mouse onecut1 peaks ####
#### load the LGS onecut1 peaks ######

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd("/zp1/data/plyu3/LGS_13_NEW/")

LGS_onecut1_peaks <- readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Onecut1_Oct16")
mouse_onecut1_peaks <- readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Onecut1_Oct16")


####### Mouse_alignment
####### LGS_alignment 




#######
#######

Compare_LGS_peaks_WITH_Mouse_peaks <- function(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment){
    ############ first we will match each mouse peak to LGS ############
    M_peak = GRanges(Mouse_peaks$peak)
    Mto_M = data.frame(findOverlaps(M_peak,Mouse_alignment))
    Mto_M_RES = data.frame(M_peak=as.character(M_peak)[Mto_M$queryHits],M_align=as.character(Mouse_alignment)[Mto_M$subjectHits])
    Mto_M_RES$LGS_align = as.character(LGS_alignment)[Mto_M$subjectHits]
    ############
    #### LGS_align to LGS peaks ######
    ############
    LGS_align_index = Mto_M_RES$LGS_align
    LGS_align_index = LGS_align_index[!duplicated(LGS_align_index)]
    L_peak = GRanges(LGS_peaks$peak)
    LGS_align_index_GR = GRanges(LGS_align_index)
    LGS_align_index_count = data.frame(findOverlaps(LGS_align_index_GR,L_peak))
    LGS_align_index_count_RES = data.frame(LGS_align=LGS_align_index[LGS_align_index_count$queryHits],L_peak=LGS_peaks$peak[LGS_align_index_count$subjectHits])
    ####
    Merge_RES = merge(Mto_M_RES,LGS_align_index_count_RES)
    ####
    library(BSgenome.Mmusculus.UCSC.mm10)
    library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")
    #####
    Total_score = c()
    for(i in 1:dim(Merge_RES)[1]){
        #########
        M_peak_tmp = GRanges(Merge_RES$M_peak[i])
        L_peak_tmp = GRanges(Merge_RES$L_peak[i])
        strand_tmp = as.character(strand(GRanges(Merge_RES$LGS_align[i])))
        #########
        M_seq = getSeq(BSgenome.Mmusculus.UCSC.mm10, M_peak_tmp)
        if(strand_tmp == "+"){
            L_seq = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, L_peak_tmp)
        }
        if(strand_tmp == "-"){
            L_seq = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, L_peak_tmp)
            L_seq <- reverseComplement(L_seq)
        }
        ########## calculate the score #########
        ##########
        alignment <- pairwiseAlignment(M_seq, L_seq, type = "local")
        start_pattern <- start(subject(alignment))
        end_pattern <- end(subject(alignment))      
        start_subject <- start(pattern(alignment))  
        end_subject <- end(pattern(alignment))   
        length_pattern <- end_pattern - start_pattern + 1
        length_subject <- end_subject - start_subject + 1
        #########
        length = min(length_pattern,length_subject)
        score = score(alignment)
        score_tmp = 0
        if(score > 50){
            #####
            score_tmp = length / 501
        }
        ####
        Total_score = c(Total_score,score_tmp)
    }
    ########
    Merge_RES$Score = Total_score
    ########
    ######## OutPut Mouse enhancer wig files ######
    ########
    M_peak_score_tab = Merge_RES[,c("M_peak","Score")]
    M_peak_score_tab_res = tapply(M_peak_score_tab$Score,M_peak_score_tab$M_peak,max)
    M_peak_score_tab_res = data.frame(M_peak = names(M_peak_score_tab_res), Score = as.numeric(M_peak_score_tab_res))
    ########
    L_peak_score_tab = Merge_RES[,c("L_peak","Score")]
    L_peak_score_tab_res = tapply(L_peak_score_tab$Score,L_peak_score_tab$L_peak,max)
    L_peak_score_tab_res = data.frame(L_peak = names(L_peak_score_tab_res), Score = as.numeric(L_peak_score_tab_res))
    #########
    file_M = paste0("Mouse_",TAG,".wig")
    M_peak_score_tab_bed = data.frame(GRanges(M_peak_score_tab_res$M_peak))[,c(1,2,3)]
    M_peak_score_tab_bed$Score = M_peak_score_tab_res$Score
    write.table(M_peak_score_tab_bed,file=file_M,sep="\t",quote=F,col.names=F,row.names=F)
    ########
    ######## OutOut LGS enhancer wig files ########
    ########
    file_M = paste0("LGS13_",TAG,".wig")
    L_peak_score_tab_bed = data.frame(GRanges(L_peak_score_tab_res$L_peak))[,c(1,2,3)]
    L_peak_score_tab_bed$Score = L_peak_score_tab_res$Score
    write.table(L_peak_score_tab_bed,file=file_M,sep="\t",quote=F,col.names=F,row.names=F)
    ########
    print("Done!!!")
}


LGS_onecut1_peaks <- readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Onecut1_Oct16")
mouse_onecut1_peaks <- readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Onecut1_Oct16")

LGS_peaks = LGS_onecut1_peaks
Mouse_peaks = mouse_onecut1_peaks
TAG = "Onecut1"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment

Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


######


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Onecut2_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Onecut2_Oct16")
TAG = "Onecut2"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Zic3_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Zic3_Oct16")
TAG = "Zic3"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Sall3_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Sall3_Oct16")
TAG = "Sall3"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Mef2c_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Mef2c_Oct16")
TAG = "Mef2c"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Thrb_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Thrb_Oct16")
TAG = "Thrb"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Rxrg_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Rxrg_Oct16")
TAG = "Rxrg"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)


LGS_peaks = readRDS("/zp1/data/plyu3/LGS_13_NEW/PtoG_13LGS_Pou2f1_Oct16")
Mouse_peaks = readRDS("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/PtoG_Mouse_Pou2f1_Oct16")
TAG = "Pou2f1"
Mouse_alignment = Mouse_alignment
LGS_alignment = LGS_alignment
Compare_LGS_peaks_WITH_Mouse_peaks(LGS_peaks,Mouse_peaks,TAG,Mouse_alignment,LGS_alignment)

#######
#######
#######
cat LGS13_Rxrg.wig LGS13_Thrb.wig LGS13_Mef2c.wig LGS13_Zic3.wig LGS13_Sall3.wig LGS13_Zic3.wig LGS13_Onecut2.wig LGS13_Onecut1.wig > LGS13_enhancer.wig

#### 1 #### 
####
####
#### 需要修改的 点：
#### 1） enhancer的多少
#### 2） alignment file 的填充 (重要)
#### 3） 增加 TSS区域 peak的比较 (重要) TSS区域需要扩大一点！！！！！ ########
#### 4） blast 没有问题 
#### 5） peak 是否需要 extend ？######
####
#### OK！！！！ #####
####


####
#### 2） alignment file 的填充 (重要) --------------------------------------------------------------------------------------------------------------------------------------------
#### 3） 不需要 填充，每一个 peak 找他对应的 alignments 就可以了 ##########
#### 直接修改一下 function #### 
#### 使用 precede() 和 follow() 函数可以分别找到下游最近和上游最近的 GRanges。弄两次就行 ########
#### 计算相似的peak只保留相似性最大的peak ########
####


ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd('/zp1/data/plyu3/LGS_13/Convert')
files = list.files()
files_need = files[grep('satsuma_summary.refined',files)]

####
####

table_list = list()

for(i in 1:length(files_need)){
	print(i)
	tmp_file = read.table(files_need[i])
	table_list = c(table_list,list(tmp_file))
}

head(table_list[[1]])


table_list_combined = do.call('rbind',table_list)

#diff1 = table_list_combined$V3 - table_list_combined$V2
#diff2 = table_list_combined$V6 - table_list_combined$V5
#

new_alignments <- data.frame(M_chr=chr,M_s=table_list_combined$V2,M_e=table_list_combined$V3,LGS_chr=table_list_combined$V4,LGS_s=table_list_combined$V5,LGS_e=table_list_combined$V6,LGS_di=table_list_combined$V8)

######
######

new_alignments$index = paste0(new_alignments$LGS_chr,"_",new_alignments$LGS_di)
new_alignments_sp = split(new_alignments,new_alignments$index)

res = sapply(new_alignments_sp,function(x) dim(x)[1])

k = which(res > 50)
names(res)[k]

new_alignments_sp_cl = new_alignments_sp[k]
sapply(new_alignments_sp_cl,function(x) dim(x)[1])

tmp_GR_tab = new_alignments_sp_cl[[1]]



######
使用 precede() 和 follow() 函数可以分别找到下游最近和上游最近的 GRanges。
######




index = table_list_combined$V1
index_sp = strsplit(index, split=':')

chr = sapply(index_sp,function(x) x[[4]])
chr = paste('chr',chr,sep='')
start = table_list_combined$V2
end = table_list_combined$V3

library(GenomicRanges)
Mouse_alignment = GRanges(chr,IRanges(start,end),score=table_list_combined$V7)
LGS_alignment = GRanges(seqnames=table_list_combined$V4,IRanges(table_list_combined$V5,table_list_combined$V6),strand=table_list_combined$V8,score=table_list_combined$V7)

Mouse_alignment_bed <- data.frame(Mouse_alignment)
Mouse_alignment_wig <- Mouse_alignment_bed[,c(1,2,3,6)]
write.table(Mouse_alignment_wig,file="Mouse_alignment.wig",sep="\t",quote=F,row.names=F,col.names=F)


LGS13_alignment_bed <- data.frame(LGS_alignment)
LGS13_alignment_wig <- LGS13_alignment_bed[,c(1,2,3,6)]
write.table(LGS13_alignment_wig,file="LGS13_alignment.wig",sep="\t",quote=F,row.names=F,col.names=F)

setwd('/zp1/data/plyu3/LGS_13/Convert')
saveRDS(Mouse_alignment,file="Mouse_alignment_Oct17")
saveRDS(LGS_alignment,file="LGS_alignment_Oct17")


####
####
####
####
####
####------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####
####
####
####
####
####
####
####
####
####
####
####
setwd("/zp1/data/plyu3/LGS_13_NEW")
Peak_Avg = readRDS("LGS13_Peak_Avg_Oct14.rds")
LGS13_total_peaks <- rownames(Peak_Avg)

setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Peak_Avg <- readRDS("Mouse_Peak_Avg_Oct14.rds")
Mouse_total_peaks <- rownames(Peak_Avg)

###########
###########
setwd('/zp1/data/plyu3/LGS_13/Convert')

Mouse_alignment <- readRDS("Mouse_alignment_Oct17")
LGS_alignment <- readRDS("LGS_alignment_Oct17")

LGS_peaks = LGS13_total_peaks
Mouse_peaks = Mouse_total_peaks

Mouse_peaks_bed = data.frame(GRanges(Mouse_peaks))
write.table(Mouse_peaks_bed[,c(1:3)],file="Mouse_peaks_total.bed",sep="\t",quote=F,col.names=F,row.names=F)

###########
########### 左右 cutoff 取 2.5 kb regions #######
###########

Refill_the_alignments <- function(Mouse_alignment,LGS_alignment){
    #########
    Mouse_alignment_s <- sort(Mouse_alignment)
    all.equal(Mouse_alignment_s,Mouse_alignment)
    #########
    Mouse_alignment_s_list = split(Mouse_alignment_s,seqnames(Mouse_alignment_s))
    #########
    diff_peaks = c()
    for(i in 1:length(Mouse_alignment_s_list)){
        ####
        print(i)
        tmp_GR = Mouse_alignment_s_list[[i]]
        tmp_s = min(start(tmp_GR))
        tmp_e = max(end(tmp_GR))
        #####
        tmp_chr_names <- as.character(seqnames(tmp_GR[1]))
        chr = sapply(strsplit(tmp_chr_names,split=":"),function(x) x[[1]])
        #####
        tmp_chr = GRanges(seqnames=chr,ranges=IRanges(start=tmp_s,end=tmp_e))
        #####
        tmp_diff = setdiff(tmp_chr,tmp_GR)
        #####
        k2 = which(width(tmp_diff) >1)
        tmp_diff = tmp_diff[k2]
        print(summary(width(tmp_diff)))
        #####
        ##### construct a table for the tmp_diff ######
        #####
        diff_peaks <- c(diff_peaks,as.character(tmp_diff))
    }
    #########
    diff_tab = data.frame(index=diff_peaks)
    ######### get left and right #####
    diff_tab$left = as.character(Mouse_alignment)[follow(GRanges(diff_tab$index), subject=Mouse_alignment)]
    diff_tab$right = as.character(Mouse_alignment)[precede(GRanges(diff_tab$index), subject=Mouse_alignment)]
    ######### get LGS left and right #######
    diff_tab$LGS_left = as.character(LGS_alignment)[match(diff_tab$left,as.character(Mouse_alignment))]
    diff_tab$LGS_right = as.character(LGS_alignment)[match(diff_tab$right,as.character(Mouse_alignment))]
    #########
    diff_tab$LGS_chr_left = sapply(strsplit(diff_tab$LGS_left,split=":"),function(x) x[[1]])
    diff_tab$LGS_chr_right = sapply(strsplit(diff_tab$LGS_right,split=":"),function(x) x[[1]])
    #### remove the chr ####
    rm1 = which(diff_tab$LGS_chr_left != diff_tab$LGS_chr_right)
    diff_tab = diff_tab[-rm1,]
    #### 
    diff_tab$LGS_strand_left = as.character(strand(GRanges(diff_tab$LGS_left)))
    diff_tab$LGS_strand_right = as.character(strand(GRanges(diff_tab$LGS_right)))
    #### remove the strand  ####
    rm2 = which(diff_tab$LGS_strand_left != diff_tab$LGS_strand_right)
    diff_tab = diff_tab[-rm2,]
    ####
    #### see the width #####
    ####
    diff_tab$right_left_1 = start(GRanges(diff_tab$LGS_right)) - start(GRanges(diff_tab$LGS_left))
    diff_tab$right_left_2 = end(GRanges(diff_tab$LGS_right)) - end(GRanges(diff_tab$LGS_left))
    ####
    rm3 = which(diff_tab$LGS_strand_left == "+" & diff_tab$right_left_1 < 0)
    rm4 = which(diff_tab$LGS_strand_left == "+" & diff_tab$right_left_2 < 0)
    rm5 = which(diff_tab$LGS_strand_left == "-" & diff_tab$right_left_1 > 0)
    rm6 = which(diff_tab$LGS_strand_left == "-" & diff_tab$right_left_1 > 0)
    ####
    rm7 = c(rm3,rm4,rm5,rm6)
    rm7 = rm7[!duplicated(rm7)]
    ####
    diff_tab_cl = diff_tab[-rm7,]
    #### Next we will merge a combined LGS regions ######
    ####
    diff_tab_cl_plus = diff_tab_cl[which(diff_tab_cl$LGS_strand_left == "+"),]
    diff_tab_cl_minus = diff_tab_cl[which(diff_tab_cl$LGS_strand_left == "-"),]
    ####
    #### don't ####
    #### next we will merge the GRanges ####
    #### head(diff_tab_cl_plus)
    LGS_s_1 = end(GRanges(diff_tab_cl_plus$LGS_left)) 
    LGS_e_1 = start(GRanges(diff_tab_cl_plus$LGS_right)) 
    #### we will remove more regions if they are no corresponding regions #####
    LGS_rm = which(LGS_s_1 > LGS_e_1)
    diff_tab_cl_plus_cl = diff_tab_cl_plus[-LGS_rm,]
    #### k = which(LGS_e_1 < LGS_s_1) head(diff_tab_cl_plus[k,])
    LGS_s_1 = end(GRanges(diff_tab_cl_plus_cl$LGS_left)) 
    LGS_e_1 = start(GRanges(diff_tab_cl_plus_cl$LGS_right)) 
    ####
    LGS_chr_1 = diff_tab_cl_plus_cl$LGS_chr
    LGS_merge = GRanges(seqnames=LGS_chr_1,IRanges(start=LGS_s_1,end=LGS_e_1),strand=diff_tab_cl_plus_cl$LGS_strand_left)
    ####
    diff_tab_cl_plus_cl$LGS_merge <- as.character(LGS_merge)
    #### Next we will see the minus strand #####
    ####
    ####
    LGS_s_2 = end(GRanges(diff_tab_cl_minus$LGS_right)) 
    LGS_e_2 = start(GRanges(diff_tab_cl_minus$LGS_left)) 
    LGS_rm2 = which(LGS_s_2 > LGS_e_2)
    diff_tab_cl_minus_cl = diff_tab_cl_minus[-LGS_rm2,]
    #####
    LGS_s_2 = end(GRanges(diff_tab_cl_minus_cl$LGS_right)) 
    LGS_e_2 = start(GRanges(diff_tab_cl_minus_cl$LGS_left)) 
    ####
    LGS_chr_2 = diff_tab_cl_minus_cl$LGS_chr
    LGS_merge2 = GRanges(seqnames=LGS_chr_2,IRanges(start=LGS_s_2,end=LGS_e_2),strand=diff_tab_cl_minus_cl$LGS_strand_left)
    ####
    diff_tab_cl_minus_cl$LGS_merge <- as.character(LGS_merge2)
    ####
    diff_tab_cl_merge_cl = rbind(diff_tab_cl_plus_cl,diff_tab_cl_minus_cl)
    #### finally, we will compare the length #####
    ####
    diff_tab_cl_merge_cl$Mouse_length = width(GRanges(diff_tab_cl_merge_cl$index))
    diff_tab_cl_merge_cl$LGS_length = width(GRanges(diff_tab_cl_merge_cl$LGS_merge))
    ####
    k_fold = which(diff_tab_cl_merge_cl$Mouse_length > 100000 | diff_tab_cl_merge_cl$LGS_length > 100000)
    diff_tab_cl_merge_cl = diff_tab_cl_merge_cl[-k_fold,]
    #########
    return(diff_tab_cl_merge_cl)
}
#######
####### merge the diff_tab_cl_merge_cl with the Ori Mapping ######
#######

Ori = data.frame(Mouse_peak = as.character(Mouse_alignment),LGS_peak=as.character(LGS_alignment))
head(Ori)
Add = data.frame(Mouse_peak = diff_tab_cl_merge_cl$index, LGS_peak = diff_tab_cl_merge_cl$LGS_merge)

Combined_alignment = rbind(Ori,Add)

setwd("/zp1/data/plyu3/LGS_13/Convert")
saveRDS(Combined_alignment,file="Combined_alignment")
#######

GR_bed = data.frame(GRanges(Combined_alignment$Mouse_peak))
write.table(GR_bed[,c(1:3)],file="M_chr.bed",sep="\t",row.names=F,col.names=F,quote=F)

GR = tmp_GR
File="M_chr1.bed"

GR_to_Bed <- function(GR,File){
    #######
    GR_bed = data.frame(GR)
    write.table(GR_bed[,c(1:3)],file=File,sep="\t",row.names=F,col.names=F,quote=F)
}

GR_to_Bed(tmp_GR,"M_chr1.bed")
GR_to_Bed(tmp_diff,"M_chr1_diff.bed")



#######---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### In this function, we will mapping all the peaks in Mouse and all the peaks in 13LGS by the mapping results #####
#######
setwd("/zp1/data/plyu3/LGS_13_NEW")
Peak_Avg = readRDS("LGS13_Peak_Avg_Oct14.rds")
LGS13_total_peaks <- rownames(Peak_Avg)

setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Peak_Avg <- readRDS("Mouse_Peak_Avg_Oct14.rds")
Mouse_total_peaks <- rownames(Peak_Avg)

######
######

M_tab <- findOverlaps(GRanges(Mouse_total_peaks),GRanges(Combined_alignment$Mouse_peak))
M_tab <- data.frame(M_tab)
M_tab$M_peaks = Mouse_total_peaks[M_tab$queryHits]
M_tab$M_Align = Combined_alignment$Mouse_peak[M_tab$subjectHits]

m = match(M_tab$M_Align,Combined_alignment$Mouse_peak)
M_tab$LGS_Align = Combined_alignment$LGS_peak[m]

M_tab = M_tab[,c(-1,-2)]

L_tab <- findOverlaps(GRanges(LGS13_total_peaks),GRanges(Combined_alignment$LGS_peak))
L_tab <- data.frame(L_tab)
L_tab$LGS_peaks = LGS13_total_peaks[L_tab$queryHits]
L_tab$LGS_Align = Combined_alignment$LGS_peak[L_tab$subjectHits]
L_tab = L_tab[,c(-1,-2)]

MergeTab = merge(M_tab,L_tab)

MergeTab = MergeTab[,c("M_Align","LGS_Align","M_peaks","LGS_peaks")]

MergeTab$Strand = as.character(strand(GRanges(MergeTab$LGS_Align)))

##### Next we will get the sequence ######
#####
library(BSgenome.Mmusculus.UCSC.mm10)
library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")

MergeTab_Plus = MergeTab[which(MergeTab$)]

M_seq = getSeq(BSgenome.Mmusculus.UCSC.mm10, GRanges(MergeTab$M_peaks))
M_seq_v = as.character(M_seq)

library(Biostrings)
MergeTab$M_seq = M_seq_v
#####
#####
MergeTab_Plus = MergeTab[which(MergeTab$Strand == "+"),]
MergeTab_Minus = MergeTab[which(MergeTab$Strand == "-"),]

L_seq_Plus = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(MergeTab_Plus$LGS_peaks))
L_seq_Plus_v = as.character(L_seq_Plus)

MergeTab_Plus$L_seq = L_seq_Plus_v


L_seq_Minus = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(MergeTab_Minus$LGS_peaks))
L_seq_Minus_Rev <- reverseComplement(L_seq_Minus)
L_seq_Minus_v = as.character(L_seq_Minus_Rev)
MergeTab_Minus$L_seq = L_seq_Minus_v

MergeTab_All = rbind(MergeTab_Plus,MergeTab_Minus)

MergeTab_All_List <- apply(MergeTab_All, 1, function(row) as.list(row))


setwd("/zp1/data/plyu3/LGS_13/Convert")
saveRDS(MergeTab_All,file="Mouse_LGS13_all_peaks_alignment.rds")

MergeTab_All <- readRDS("Mouse_LGS13_all_peaks_alignment.rds")

x = MergeTab_All[1,]
#####
compare_2_seq <- function(x){
    ######
    library(Biostrings)
    seq1 = DNAString(x$M_seq)
    seq2 = DNAString(x$L_seq)
    ######
    alignment <- pairwiseAlignment(seq1, seq2, type = "local")
    start_pattern <- start(subject(alignment))
    end_pattern <- end(subject(alignment))      
    start_subject <- start(pattern(alignment))  
    end_subject <- end(pattern(alignment))   
    length_pattern <- end_pattern - start_pattern + 1
    length_subject <- end_subject - start_subject + 1
    #########
    length = min(length_pattern,length_subject)
    score = round(score(alignment),3)
    #########
    out = paste0(length,":",score)
    #########
    return(out)
}

#####
library(pbapply)
library(parallel)
cl <- makeCluster(35)  # 这里使用2个核心


alignment_res = parLapply(cl, MergeTab_All_List, compare_2_seq)

save(alignment_res,file="alignment_res_Oct22")
#####
#####


ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd("/zp1/data/plyu3/LGS_13/Convert")

load("alignment_res_Oct22")

#####
##### OK!!! Let us add the results to each peak pairs ######
#####
strsplit(alignment_res[[1]],split=":")[[1]][1]

overlap_counts = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][1]))
overlap_score = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][2]))


MergeTab_All <- readRDS("Mouse_LGS13_all_peaks_alignment.rds")
MergeTab_All$Overlap_region = overlap_counts
MergeTab_All$Overlap_score = overlap_score

saveRDS(MergeTab_All,file="Mouse_LGS13_all_peaks_alignment.rds")


#########
######### Next we will load the signal for each peaks ###############
#########

setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")
LGS13_peaks_signals <- readRDS(file="LGS13_peaks_signals_Oct22")
colnames(LGS13_peaks_signals)

LGS13_peaks_signals$H3K27ac_signal = (LGS13_peaks_signals$H3K27ac_Rep1 + LGS13_peaks_signals$H3K27ac_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$H3K4me3_signal = (LGS13_peaks_signals$H3K4me3_Rep1 + LGS13_peaks_signals$H3K4me3_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$Otx2_signal = (LGS13_peaks_signals$Otx2_AbA_Rep1 + LGS13_peaks_signals$Otx2_AbA_Rep2 + LGS13_peaks_signals$Otx2_AbB)/3 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$Neurod1_signal = (LGS13_peaks_signals$Neurod1_Rep1 + LGS13_peaks_signals$Neurod1_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2

setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")
saveRDS(LGS13_peaks_signals,file="LGS13_peaks_signals_Oct22")
#####
#####-----Next we will process the LGS13 signals-------#########
#####


setwd("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")
Mouse_peaks_signals <- readRDS(file="Mouse_peaks_signals_Oct22")
colnames(Mouse_peaks_signals)

Mouse_peaks_signals$H3K27ac_signal = (Mouse_peaks_signals$H3K27ac_Rep1 + Mouse_peaks_signals$H3K27ac_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$H3K4me3_signal = (Mouse_peaks_signals$H3K4me3_Rep1 + Mouse_peaks_signals$H3K4me3_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$Otx2_signal = (Mouse_peaks_signals$Otx2_AbA_Rep1 + Mouse_peaks_signals$Otx2_AbA_Rep2 + Mouse_peaks_signals$Otx2_AbB_Rep1 + Mouse_peaks_signals$Otx2_AbB_Rep2)/4 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$Neurod1_signal = (Mouse_peaks_signals$Neurod1_Rep1 + Mouse_peaks_signals$Neurod1_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2


setwd("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")
saveRDS(Mouse_peaks_signals,file="Mouse_peaks_signals_Oct22")

#####-----Next we will process the Mouse signals-------#########
#####
##### we will check the signal in Mouse #####
#####

chr14:17659855-17660354
chr14:17734383-17734882


chr14:17989727-17990226

k = which(countOverlaps(GRanges(rownames(Mouse_peaks_signals)),GRanges("chr14:17989727-17990226")) > 0)
Mouse_peaks_signals[k,]

######
###### see the signal for each genes !!!!! ######
######

###### load the PtoGs !!!! #############
######


TSS and Body Peaks 
Intergenic region for each Gene

Mouse to first find each genes related peaks:

######
###### input is the annotations ######
######

setwd("/zp1/data/plyu3/LGS_13_NEW")
LGS13_Gene_Anno <- readRDS("LGS13_Gene_Anno_Oct23")
LGS13_Peak <- readRDS("LGS13_Peak_Avg_Oct14.rds")
LGS13_PtoG <- readRDS("LGS13_p2g_tab_Oct23")

######
###### This step we will add PtoG correlations #######
######

get_TSS <- function(Genes_sub){
    #####
    if(as.character(strand(Genes_sub)) == "+"){
        end(Genes_sub) = start(Genes_sub)
    }
    if(as.character(strand(Genes_sub)) == "-"){
        start(Genes_sub) = end(Genes_sub)
    }
    strand(Genes_sub) = "*"
    return(Genes_sub)
}

Function_1_find_peaks <- function(GeneName,Gtf,PtoG,PeakAvg){
    ##########
    Gene_ranges = Gtf[which(Gtf$symbol == GeneName)]
    Gene_TSS = get_TSS(Gene_ranges)
    ########## extend the TSS and merge the regions ##########
    strand(Gene_ranges) <- "*"
    start(Gene_TSS) <- start(Gene_TSS) - 2000
    end(Gene_TSS) <- end(Gene_TSS) + 2000
    ##########
    Gene_ranges_adjust = union(Gene_TSS,Gene_ranges)
    ##########
    k = which(countOverlaps(GRanges(rownames(PeakAvg)),Gene_ranges_adjust) > 0)
    Potential_peaks_TSSBody <- rownames(PeakAvg)[k]
    ##########
    ########## Next we will identify intergenetic peaks ###########
    ##########
    All_body = Gtf
    All_TSS_Plus = All_body[which(strand(All_body) == "+")]
    All_TSS_Minus = All_body[which(strand(All_body) == "-")]
    ##########
    end(All_TSS_Plus) = start(All_TSS_Plus)
    start(All_TSS_Minus) = end(All_TSS_Minus)
    All_TSS_merge = c(All_TSS_Plus,All_TSS_Minus)
    ##########
    strand(All_TSS_merge) <- "*"
    start(All_TSS_merge) <- start(All_TSS_merge) - 2000
    end(All_TSS_merge) <- end(All_TSS_merge) + 2000
    ##########
    k2 = which(countOverlaps(GRanges(rownames(PeakAvg)),Gene_ranges_adjust) == 0 & countOverlaps(GRanges(rownames(PeakAvg)),All_TSS_merge) == 0)
    ##########
    Interpeaks_total <- rownames(PeakAvg)[k2]
    ##########
    ########## target gene inter ####
    ########## 直接从 PtoG 里面找就行 #####
    ########## PtoG 500kb ######
    PtoG_cl = PtoG[which(PtoG$gene == GeneName),]
    PtoG_clcl = PtoG_cl[which(PtoG_cl$peak %in% Interpeaks_total == T),]
    ##########
    Potential_peaks_Inter = PtoG_clcl$peak
    ##########
    Total_Potential_peaks <- c(Potential_peaks_Inter,Potential_peaks_TSSBody)
    Total_Potential_peaks = Total_Potential_peaks[!duplicated(Total_Potential_peaks)]
    ##########
    ########## Next we will add the class ####
    ##########
    Total_Potential_peaks_tab = data.frame(Gene=GeneName,Peak=Total_Potential_peaks,Class="NA",PtoG="NA",MAX_Access=0)
    ##########
    class1 = which(Total_Potential_peaks_tab$Peak %in% Interpeaks_total == T)
    Total_Potential_peaks_tab$Class[class1] = "inter"
    ##########
    Gene_ranges = Gtf[which(Gtf$symbol == GeneName)]
    class2 = which(countOverlaps(GRanges(Total_Potential_peaks_tab$Peak),Gene_ranges) > 0)
    Total_Potential_peaks_tab$Class[class2] = "body"
    ##########
    Gene_TSS = get_TSS(Gene_ranges)
    strand(Gene_ranges) <- "*"
    start(Gene_TSS) <- start(Gene_TSS) - 2000
    end(Gene_TSS) <- end(Gene_TSS) + 2000
    class3 = which(countOverlaps(GRanges(Total_Potential_peaks_tab$Peak),Gene_TSS) > 0)
    Total_Potential_peaks_tab$Class[class3] = "TSS"
    ########## Next we will add PtoG #############
    ##########
    m = match(Total_Potential_peaks_tab$Peak,PtoG_cl$peak)
    Total_Potential_peaks_tab$PtoG = PtoG_cl$Correlation[m]
    ########## Next we will add the Max access ####
    ##########
    ##########
    max_access = apply(PeakAvg,1,max)
    m2 = match(Total_Potential_peaks_tab$Peak,names(max_access))
    Total_Potential_peaks_tab$MAX_Access = as.numeric(max_access)[m2]
    ##########
    Total_Potential_peaks_tab_cl = Total_Potential_peaks_tab[which(Total_Potential_peaks_tab$MAX_Access >1),]
    ###########
    rm = which(Total_Potential_peaks_tab_cl$PtoG < 0.25 & Total_Potential_peaks_tab_cl$Class == "inter")
    if(length(rm) > 0){
        Total_Potential_peaks_tab_cl = Total_Potential_peaks_tab_cl[-rm,]
    }
    ###########
    return(Total_Potential_peaks_tab_cl)
}


GeneName = "Thrb"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Thrb_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)

GeneName = "Mef2c"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Mef2c_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)

GeneName = "Rxrg"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Rxrg_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Pou2f1"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Pou2f1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Onecut1"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Onecut1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Onecut2"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Onecut2_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Zic3"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Zic3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)



GeneName = "Zic3"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Zic3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Sall3"
Gtf = LGS13_Gene_Anno
PtoG = LGS13_PtoG
PeakAvg = LGS13_Peak
LGS13_Sall3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)





setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")


Function_2_add_scores <- function(Peaks,Peak_signals){
    ##########
    m = match(Peaks$Peak,Peak_signals$peaks)
    ##########
    Peaks$H3K27ac = Peak_signals$H3K27ac_signal[m]
    Peaks$H3K4me3 = Peak_signals$H3K4me3_signal[m]
    Peaks$Otx2 = Peak_signals$Otx2_signal[m]
    Peaks$Neurod1 = Peak_signals$Neurod1_signal[m]
    ###########
    #### cutoff Otx2 > 1 ####
    #### cutoff Neurod > 1 ##
    #### cutoff H3K27ac > 1 #
    #### cutoff H3K4me3 > 5 #
    ###########
    ###########
    Peaks$epi_tag = "potential_regulatory_elements"
    k1 = which(Peaks$H3K27ac > 1 & Peaks$H3K4me3 < 5)
    Peaks$epi_tag[k1] = "Activate_enhancer"
    k2 = which(Peaks$H3K4me3 > 5)
    Peaks$epi_tag[k2] = "Promoter"
    ######
    Peaks$cut_run = ""
    ######
    k3 = which(Peaks$Otx2 > 2)
    Peaks$cut_run[k3] = paste(Peaks$cut_run[k3],"Otx2_binding")
    ######
    k4 = which(Peaks$Neurod1 > 2)
    Peaks$cut_run[k4] = paste(Peaks$cut_run[k4],"Neurod1_binding")
    ######
    ###### OK！！！####
    ######
    return(Peaks)
}

Peaks = LGS13_Thrb_peaks
LGS13_peaks_signals = readRDS(file="LGS13_peaks_signals_Oct22")
Peak_signals = LGS13_peaks_signals


Peaks = LGS13_Thrb_peaks
LGS13_Thrb_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Mef2c_peaks
LGS13_Mef2c_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Onecut1_peaks
LGS13_Onecut1_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Onecut2_peaks
LGS13_Onecut2_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Zic3_peaks
LGS13_Zic3_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Sall3_peaks
LGS13_Sall3_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Rxrg_peaks
LGS13_Rxrg_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

Peaks = LGS13_Pou2f1_peaks
LGS13_Pou2f1_peaks_res <- Function_2_add_scores(Peaks,Peak_signals)

setwd("/zp1/data/plyu3/LGS_13/Convert/")
saveRDS(LGS13_Thrb_peaks_res,file="LGS13_Thrb_peaks_res.rds")
saveRDS(LGS13_Mef2c_peaks_res,file="LGS13_Mef2c_peaks_res.rds")
saveRDS(LGS13_Onecut1_peaks_res,file="LGS13_Onecut1_peaks_res.rds")
saveRDS(LGS13_Onecut2_peaks_res,file="LGS13_Onecut2_peaks_res.rds")
saveRDS(LGS13_Zic3_peaks_res,file="LGS13_Zic3_peaks_res.rds")
saveRDS(LGS13_Sall3_peaks_res,file="LGS13_Sall3_peaks_res.rds")
saveRDS(LGS13_Rxrg_peaks_res,file="LGS13_Rxrg_peaks_res.rds")
saveRDS(LGS13_Pou2f1_peaks_res,file="LGS13_Pou2f1_peaks_res.rds")

##### #####
##### #####
##### Next for Mouse genes #####
##### #####
##### #####

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Gene_Anno <- readRDS("Mouse_Gene_Anno_Oct23")
Mouse_Peak <- readRDS("Mouse_Peak_Avg_Oct14.rds")
Mouse_PtoG <- readRDS("Mouse_p2g_tab_Oct23")



GeneName = "Thrb"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Thrb_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)

GeneName = "Mef2c"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Mef2c_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)

GeneName = "Rxrg"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Rxrg_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Pou2f1"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Pou2f1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Onecut1"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Onecut1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Onecut2"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Onecut2_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Zic3"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Zic3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)


GeneName = "Sall3"
Gtf = Mouse_Gene_Anno
PtoG = Mouse_PtoG
PeakAvg = Mouse_Peak
Mouse_Sall3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg)



setwd("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")
Mouse_peaks_signals = readRDS(file="Mouse_peaks_signals_Oct22")


Mouse_Thrb_peaks_res <- Function_2_add_scores(Mouse_Thrb_peaks,Mouse_peaks_signals)
Mouse_Mef2c_peaks_res <- Function_2_add_scores(Mouse_Mef2c_peaks,Mouse_peaks_signals)
Mouse_Onecut1_peaks_res <- Function_2_add_scores(Mouse_Onecut1_peaks,Mouse_peaks_signals)
Mouse_Onecut2_peaks_res <- Function_2_add_scores(Mouse_Onecut2_peaks,Mouse_peaks_signals)
Mouse_Zic3_peaks_res <- Function_2_add_scores(Mouse_Zic3_peaks,Mouse_peaks_signals)
Mouse_Sall3_peaks_res <- Function_2_add_scores(Mouse_Sall3_peaks,Mouse_peaks_signals)
Mouse_Rxrg_peaks_res <- Function_2_add_scores(Mouse_Rxrg_peaks,Mouse_peaks_signals)
Mouse_Pou2f1_peaks_res <- Function_2_add_scores(Mouse_Pou2f1_peaks,Mouse_peaks_signals)

setwd("/zp1/data/plyu3/LGS_13/Convert/")
saveRDS(Mouse_Thrb_peaks_res,file="Mouse_Thrb_peaks_res.rds")
saveRDS(Mouse_Mef2c_peaks_res,file="Mouse_Mef2c_peaks_res.rds")
saveRDS(Mouse_Onecut1_peaks_res,file="Mouse_Onecut1_peaks_res.rds")
saveRDS(Mouse_Onecut2_peaks_res,file="Mouse_Onecut2_peaks_res.rds")
saveRDS(Mouse_Zic3_peaks_res,file="Mouse_Zic3_peaks_res.rds")
saveRDS(Mouse_Sall3_peaks_res,file="Mouse_Sall3_peaks_res.rds")
saveRDS(Mouse_Rxrg_peaks_res,file="Mouse_Rxrg_peaks_res.rds")
saveRDS(Mouse_Pou2f1_peaks_res,file="Mouse_Pou2f1_peaks_res.rds")




#######
#######
setwd("/zp1/data/plyu3/LGS_13/Convert")
Mouse_Res = readRDS("Mouse_Thrb_peaks_res.rds")
LGS_Res = readRDS("LGS13_Thrb_peaks_res.rds")
Alignment = readRDS("Mouse_LGS13_all_peaks_alignment.rds")

#######
#######

Function_3_add_Alignments <- function(Mouse_Res,LGS_Res,Alignment){
    ##########
    ### Next we will used for mouse datasets #####
    ##########
    ##########-----for-----Mouse-------#######
    ##########
    Mouse_Res$LGS_peak_conserved = "Notfind"
    Mouse_Res$LGS_peak_overlap = 0
    Mouse_Res$LGS_peak_overlapScore = 0
    for(i in 1:length(Mouse_Res$Peak)){
        #####
        tmp_peak = Mouse_Res$Peak[i]
        k1 = which(Alignment$M_peaks == tmp_peak)
        if(length(k1) == 0){
            next
        }
        LGS_Alignment_k1 = Alignment$LGS_peaks[k1]
        ######
        k2 = which(LGS_Alignment_k1 %in% LGS_Res$Peak == T)
        if(length(k2) == 0){
            next
        }
        #######
        k3 = which(LGS_Res$Peak %in% LGS_Alignment_k1 == T)
        if(length(k3) == 0){
            next
        }
        Mouse_Res$LGS_peak_conserved[i] = paste(LGS_Res$Peak[k3], collapse = "__")
        ####### Next find the score #####
        k4 = which(Alignment$M_peaks == tmp_peak & Alignment$LGS_peaks %in% LGS_Res$Peak[k3] == T)
        RES = Alignment[k4,]
        RES = RES[,c("M_peaks","LGS_peaks","Overlap_region","Overlap_score")]
        RES_cl = RES[!duplicated(paste(RES$M_peaks,RES$LGS_peaks)),]
        #######
        if(dim(RES_cl)[1] == 1){
            Mouse_Res$LGS_peak_overlap[i] = RES_cl$Overlap_region
            Mouse_Res$LGS_peak_overlapScore[i] = RES_cl$Overlap_score
        }
        if(dim(RES_cl)[1] > 1){
            print("Mutiple!")
            rownames(RES_cl) = RES_cl$LGS_peaks
            RES_cl = RES_cl[LGS_Res$Peak[k3],]
            print(i)
            Mouse_Res$LGS_peak_overlap[i] = paste(RES_cl$Overlap_region,collapse = "__")
            Mouse_Res$LGS_peak_overlapScore[i] = paste(RES_cl$Overlap_score,collapse = "__")
        }
    }
    #########
    ####--------Next for 13LGS ######
    #########
    LGS_Res$M_peak_conserved = "Notfind"
    LGS_Res$M_peak_overlap = 0
    LGS_Res$M_peak_overlapScore = 0
    for(i in 1:length(LGS_Res$Peak)){
        #####
        tmp_peak = LGS_Res$Peak[i]
        k1 = which(Alignment$LGS_peaks == tmp_peak)
        if(length(k1) == 0){
            next
        }
        M_Alignment_k1 = Alignment$M_peaks[k1]
        ######
        k2 = which(M_Alignment_k1 %in% Mouse_Res$Peak == T)
        if(length(k2) == 0){
            next
        }
        #######
        k3 = which(Mouse_Res$Peak %in% M_Alignment_k1 == T)
        if(length(k3) == 0){
            next
        }
        LGS_Res$M_peak_conserved[i] = paste(Mouse_Res$Peak[k3], collapse = "__")
        ####### Next find the score #####
        k4 = which(Alignment$LGS_peaks == tmp_peak & Alignment$M_peaks %in% Mouse_Res$Peak[k3] == T)
        RES = Alignment[k4,]
        RES = RES[,c("M_peaks","LGS_peaks","Overlap_region","Overlap_score")]
        RES_cl = RES[!duplicated(paste(RES$M_peaks,RES$LGS_peaks)),]
        #######
        if(dim(RES_cl)[1] == 1){
            LGS_Res$M_peak_overlap[i] = RES_cl$Overlap_region
            LGS_Res$M_peak_overlapScore[i] = RES_cl$Overlap_score
        }
        if(dim(RES_cl)[1] > 1){
            print("Mutiple!")
            rownames(RES_cl) = RES_cl$M_peaks
            RES_cl = RES_cl[Mouse_Res$Peak[k3],]
            print(i)
            LGS_Res$M_peak_overlap[i] = paste(RES_cl$Overlap_region,collapse = "__")
            LGS_Res$M_peak_overlapScore[i] = paste(RES_cl$Overlap_score,collapse = "__")
        }
    }
    #######
    return(list(M=Mouse_Res,L=LGS_Res))
    #######
}

######
###### Next we will output to the excel files #######
######
setwd("/zp1/data/plyu3/LGS_13/Convert")
Mouse_Res = readRDS("Mouse_Thrb_peaks_res.rds")
LGS_Res = readRDS("LGS13_Thrb_peaks_res.rds")
Thrb_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

#######
Mouse_Res = readRDS("Mouse_Mef2c_peaks_res.rds")
LGS_Res = readRDS("LGS13_Mef2c_peaks_res.rds")
Mef2c_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

#######
Mouse_Res = readRDS("Mouse_Onecut1_peaks_res.rds")
LGS_Res = readRDS("LGS13_Onecut1_peaks_res.rds")
Onecut1_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

Mouse_Res = readRDS("Mouse_Onecut2_peaks_res.rds")
LGS_Res = readRDS("LGS13_Onecut2_peaks_res.rds")
Onecut2_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

Mouse_Res = readRDS("Mouse_Zic3_peaks_res.rds")
LGS_Res = readRDS("LGS13_Zic3_peaks_res.rds")
Zic3_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

Mouse_Res = readRDS("Mouse_Sall3_peaks_res.rds")
LGS_Res = readRDS("LGS13_Sall3_peaks_res.rds")
Sall3_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

Mouse_Res = readRDS("Mouse_Rxrg_peaks_res.rds")
LGS_Res = readRDS("LGS13_Rxrg_peaks_res.rds")
Rxrg_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)

Mouse_Res = readRDS("Mouse_Pou2f1_peaks_res.rds")
LGS_Res = readRDS("LGS13_Pou2f1_peaks_res.rds")
Pou2f1_compare = Function_3_add_Alignments(Mouse_Res,LGS_Res,Alignment)


Mouse_res_list <- list(M_Thrb = Thrb_compare$M,M_Mef2c= Mef2c_compare$M,M_Onecut1= Onecut1_compare$M,M_Onecut2= Onecut2_compare$M,M_Zic3= Zic3_compare$M,M_Sall3= Sall3_compare$M,M_Rxrg= Rxrg_compare$M,M_Pou2f1= Pou2f1_compare$M)
LGS13_res_list <- list(LGS_Thrb = Thrb_compare$L,LGS_Mef2c= Mef2c_compare$L,LGS_Onecut1= Onecut1_compare$L,LGS_Onecut2= Onecut2_compare$L,LGS_Zic3= Zic3_compare$L,LGS_Sall3= Sall3_compare$L,LGS_Rxrg= Rxrg_compare$L,LGS_Pou2f1= Pou2f1_compare$L)

library(openxlsx)
# Create a new Excel workbook
wb <- createWorkbook()

# Loop through the list and add each table to a new sheet
for (sheet_name in names(Mouse_res_list)) {
  addWorksheet(wb, sheet_name)  # Add a new sheet with the name
  writeData(wb, sheet = sheet_name, Mouse_res_list[[sheet_name]])  # Write data to the sheet
}

# Save the workbook to an Excel file
saveWorkbook(wb, file = "Mouse_8genes_peaksWithLGS13.xlsx", overwrite = TRUE)



####
wb <- createWorkbook()

# Loop through the list and add each table to a new sheet
for (sheet_name in names(LGS13_res_list)) {
  addWorksheet(wb, sheet_name)  # Add a new sheet with the name
  writeData(wb, sheet = sheet_name, LGS13_res_list[[sheet_name]])  # Write data to the sheet
}

# Save the workbook to an Excel file
saveWorkbook(wb, file = "LGS13_8genes_peaksWithLGS13.xlsx", overwrite = TRUE)



#######
#######
#######

saveRDS(Mouse_Thrb_peaks_res,file="Mouse_Thrb_peaks_res.rds")
saveRDS(Mouse_Mef2c_peaks_res,file="Mouse_Mef2c_peaks_res.rds")
saveRDS(Mouse_Onecut1_peaks_res,file="Mouse_Onecut1_peaks_res.rds")
saveRDS(Mouse_Onecut2_peaks_res,file="Mouse_Onecut2_peaks_res.rds")
saveRDS(Mouse_Zic3_peaks_res,file="Mouse_Zic3_peaks_res.rds")
saveRDS(Mouse_Sall3_peaks_res,file="Mouse_Sall3_peaks_res.rds")
saveRDS(Mouse_Rxrg_peaks_res,file="Mouse_Rxrg_peaks_res.rds")
saveRDS(Mouse_Pou2f1_peaks_res,file="Mouse_Pou2f1_peaks_res.rds")

10.181.57.113 /zp1/data/share/Isabella/Xenium5k.

#######---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Hi everyone. It looks like we have everything that we need at this point to finish up the 13-LGS paper. 
Jared and Sangeetha have recently completed Cut&Tag analysis of P2 mouse and P4/5 13-LGS retina, 
profiling multiple histone modifications, as well as binding sites for Otx2 and Neurod1. 
This should allow us to identify enhancer sites associated with cone-promoting genes in late-stage primary and neurogenic RPCs, 
as well as photoreceptor precursors, that are found in 13-LGS but not mouse.  
Specifically, we've profiled the following:

H3K4me1
H3K4me3
H3K9me9
H3K27me3
Otx2 (2 antibodies)
Neurod1
IgG control

We'd like to integrate these with the scATAC-Seq data to identify differences in regulatory sites controlling expression of cone-promoting genes 
(e.g. Onecut1/2, Zic3, Pou2f1, Sall3, Rxrg, Thrb, Mef2c) that show selectively increased expression in late-stage primary/neurogenic RPCs.  
We would like to determine whether these enhancer regions are selectively targeted by the pan-neurogenic RPC/photoreceptor factors Otx2 and Neurod1.  
We would also like to determine what other TFs are predicted to selectively target these enhancer sequences.

####
#### 我们先读取一下 chip-seq #####
####

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R

setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")
file = "Otx2_AbB.RPGC.bw"
savefile = "Otx2_AbB.norm.rds"

Read_and_Process_Wig <- function(file,savefile){
    ################
    library(rtracklayer)
    bw_GR = import(file,format = "BigWig")
    ################
    ###### gr_1bp <- unlist(tile(bw_GR, width = 1)) #######
    ###### gr_1bp$score <- rep(bw_data$score, times = width(bw_data)) #######
    saveRDS(bw_GR,file=savefile)
    print("Done!")
}

setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")

Read_and_Process_Wig(file = "Otx2_AbB.RPGC.bw",savefile = "Otx2_AbB.norm.rds")

Read_and_Process_Wig(file = "H3K27ac_Rep1.RPGC.bw",savefile = "H3K27ac_Rep1.norm.rds")
Read_and_Process_Wig(file = "H3K27ac_Rep2.RPGC.bw",savefile = "H3K27ac_Rep2.norm.rds")

Read_and_Process_Wig(file = "H3K4me3_Rep1.RPGC.bw",savefile = "H3K4me3_Rep1.norm.rds")
Read_and_Process_Wig(file = "H3K4me3_Rep2.RPGC.bw",savefile = "H3K4me3_Rep2.norm.rds")

Read_and_Process_Wig(file = "IgG_Rep1.RPGC.bw",savefile = "IgG_Rep1.norm.rds")
Read_and_Process_Wig(file = "IgG_Rep2.RPGC.bw",savefile = "IgG_Rep2.norm.rds")

Read_and_Process_Wig(file = "Otx2_AbA_Rep1.RPGC.bw",savefile = "Otx2_AbA_Rep1.norm.rds")
Read_and_Process_Wig(file = "Otx2_AbA_Rep2.RPGC.bw",savefile = "Otx2_AbA_Rep2.norm.rds")

Read_and_Process_Wig(file = "Neurod1_Rep1.RPGC.bw",savefile = "Neurod1_Rep1.norm.rds")
Read_and_Process_Wig(file = "Neurod1_Rep2.RPGC.bw",savefile = "Neurod1_Rep2.norm.rds")

setwd("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")

Read_and_Process_Wig(file = "H3K27ac_Rep1.RPGC.bw",savefile = "H3K27ac_Rep1.norm.rds")
Read_and_Process_Wig(file = "H3K27ac_Rep2.RPGC.bw",savefile = "H3K27ac_Rep2.norm.rds")

Read_and_Process_Wig(file = "H3K4me3_Rep1.RPGC.bw",savefile = "H3K4me3_Rep1.norm.rds")
Read_and_Process_Wig(file = "H3K4me3_Rep2.RPGC.bw",savefile = "H3K4me3_Rep2.norm.rds")

Read_and_Process_Wig(file = "IgG_Rep1.RPGC.bw",savefile = "IgG_Rep1.norm.rds")
Read_and_Process_Wig(file = "IgG_Rep2.RPGC.bw",savefile = "IgG_Rep2.norm.rds")

Read_and_Process_Wig(file = "Neurod1_Rep1.RPGC.bw",savefile = "Neurod1_Rep1.norm.rds")
Read_and_Process_Wig(file = "Neurod1_Rep2.RPGC.bw",savefile = "Neurod1_Rep2.norm.rds")

Read_and_Process_Wig(file = "Otx2_AbA_Rep1.RPGC.bw",savefile = "Otx2_AbA_Rep1.norm.rds")
Read_and_Process_Wig(file = "Otx2_AbA_Rep2.RPGC.bw",savefile = "Otx2_AbA_Rep2.norm.rds")

Read_and_Process_Wig(file = "Otx2_AbB_Rep1.RPGC.bw",savefile = "Otx2_AbB_Rep1.norm.rds")
Read_and_Process_Wig(file = "Otx2_AbB_Rep2.RPGC.bw",savefile = "Otx2_AbB_Rep2.norm.rds")

###############
############### ############# add the signal to the peaks ################
###############

load the mouse peaks:
load the 13LGS peaks:

ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

conda activate ArchR2
R


setwd("/zp1/data/plyu3/LGS_13_NEW")
Peak_Avg = readRDS("LGS13_Peak_Avg_Oct14.rds")
LGS13_total_peaks <- rownames(Peak_Avg)

#######
setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Peak_Avg <- readRDS("Mouse_Peak_Avg_Oct14.rds")
Mouse_total_peaks <- rownames(Peak_Avg)

###### add the file list to the functions ######
######

total_peaks = LGS13_total_peaks
signal_files = list.files("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")
signal_files = signal_files[grep("rds",signal_files)]
signal_files_N = gsub(".norm.rds","",signal_files)

setwd("/zp1/data/plyu3/LGS_13/Convert/LGSCutRun")

add_signals_by_folder <- function(total_peaks,signal_files,signal_files_N){
    library(GenomicRanges)
    #########
    total_GR = GRanges(total_peaks)
    total_GR_tab = data.frame(peaks=total_peaks)
    tmp_mat = matrix(0,nrow=length(total_peaks),ncol=length(signal_files))
    colnames(tmp_mat) <- signal_files_N
    rownames(tmp_mat) <- total_peaks
    total_GR_tab = cbind(total_GR_tab,tmp_mat)
    #########
    for(i in 1:length(signal_files)){
        #######
        print(signal_files[i])
        tmp_signal <- readRDS(signal_files[i])
        #######
        k = which(countOverlaps(tmp_signal,total_GR) > 0)
        tmp_signal_cl = tmp_signal[k]
        #######
        tmp_signal_cl_1bp <- unlist(tile(tmp_signal_cl, width = 1))
        tmp_signal_cl_1bp$score <- rep(tmp_signal_cl$score, times = width(tmp_signal_cl))
        ####### add signals ######
        f = findOverlaps(total_GR,tmp_signal_cl_1bp)
        #######
        f = data.frame(f)
        f$score = tmp_signal_cl_1bp$score[f$subjectHits]
        #######
        f_res = tapply(f$score,f$queryHits,mean)
        f_res_tab = data.frame(index=names(f_res),score = as.numeric(f_res))
        #######
        total_GR_tab[as.numeric(f_res_tab$index),i+1] = f_res_tab$score
        ########
    }
    return(total_GR_tab)
}

LGS13_peaks_signals <- add_signals_by_folder(total_peaks,signal_files,signal_files_N)
saveRDS(LGS13_peaks_signals,file="LGS13_peaks_signals_Oct22")

########

total_peaks = Mouse_total_peaks
signal_files = list.files("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")
setwd("/zp1/data/plyu3/LGS_13/Convert/MouseCutRun")
signal_files = signal_files[grep("rds",signal_files)]
signal_files_N = gsub(".norm.rds","",signal_files)
Mouse_peaks_signals <- add_signals_by_folder(total_peaks,signal_files,signal_files_N)

saveRDS(Mouse_peaks_signals,file="Mouse_peaks_signals_Oct22")

#######
#######
#######


#######--------------------------- This time we will move all the files to the new server ！！！！！ --------------------------- #####
#######--------------------------- We will move all the LGS files ！！！！ ----------------------------------------------------- #####
#######
#####---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ssh plyu3@omb1.onc.jhmi.edu
njd$rft1

#####---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#####

##### direction：/zp1/data/plyu3/Old_Server_Data/plyu3 ######
#####

rm -r /zp1/data/plyu3/PHx
rm -r

#####
##### Next we will merge the files to the omb2 server !!!! #######
#####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

#####
##### cp the folder to the omb2 server ######
#####

scp -r /zp1/data/plyu3 plyu3@omb2.onc.jhmi.edu:/zp1/data/plyu3/Old_Server_Data/

##### OK !! we need to wait the results ######

是的，使用 scp（Secure Copy Protocol）时，如果目标位置已经存在同名文件，scp 默认会直接覆盖该文件，而不会提示确认。

#####
##### 我们再开一个开始工作 #####
#####
##### 看一下 Seth 的邮件 #####
#####
##### ##### 晚上洗澡 #####
#####

##### 
##### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#####

In the previous ATAC-seq analysis (PtoG correlation, peak accessibility), I included all retinal cells from both early and late stages in Mouse and 13LGS.

I’m now revising the results with the following changes:
1. For both Mouse and 13LGS, we will only use RPC, NG, and photoreceptor cells from late-stage time points to calculate the PtoG correlation and peak accessibility.
2. The specificity scores of peaks in RPC, NG, and photoreceptor cells (compared to other cell types) will be added to the Excel file.
3. We will check enriched motifs in 13LGS-specific regulatory elements (promoters and enhancers) compared to Mouse.

###########
########### first we will revise the PtoG correlations ######
########### we will perform all the analysis on the new server #######
###########




###########
########### This time we will calculate the fold changes of each peaks compared to other cell types ####
###########


###########
###########------------- load the data on the new server -------------###########
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

###########-------------
conda activate ArchR
R
library(ArchR)

###########
###########------------- we will set the folder on the new server !!!! -----
###########

library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")

##### install the LGS genome !!!! #####
##### install.packages("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise_1.0.tar.gz")


###########
########### for Mouse ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########
########### load the Mouse ArchR project !!! ######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########
#### setwd("/zp1/data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
#### /zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
load("All_arrows_project_cl_2024")

library(ArchR)
getAvailableMatrices(All_arrows_project_cl)
All_arrows_project_cl$celltype <- factor(All_arrows_project_cl$celltype,levels=c("RPCs_S1","RPCs_S2","RPCs_S3","MG","Early_NG","Late_NG","RGC","AC_HC","Early_Cone","Cone","Early_Rod","Rod","BC"))
All_arrows_project_cl$celltype2 = All_arrows_project_cl$celltype
m = match(All_arrows_project_cl$celltype2,c("RPCs_S1","RPCs_S2","RPCs_S3","MG","Early_NG","Late_NG","RGC","AC_HC","Early_Cone","Cone","Early_Rod","Rod","BC"))
All_arrows_project_cl$celltype2 = paste0(m,"_",All_arrows_project_cl$celltype2)

####
table(All_arrows_project_cl$celltype2)
####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Meta <- All_arrows_project_cl@cellColData
saveRDS(Meta,file="Mouse_Meta_scATAC_2024Nov13")

Meta <- readRDS("Mouse_Meta_scATAC_2024Nov13")

####------------------------------------------------------------------------------------------------------------------------------------
#### add new project with arrow files: #####
####-------------------------------------------------------------------------------------------------------------------------------------

ArrowFiles <- list.files()
ArrowFiles <- ArrowFiles[grep(".arrow$",ArrowFiles)]
names = ArrowFiles
names = gsub(".arrow","",names)
names(ArrowFiles) <- names

addArchRGenome("mm10")

Mouse_Project_new_server <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Mouse_Project_new_server",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveRDS(Mouse_Project_new_server,file="Mouse_Project_new_server_Nov2024")

###### filter the project by meta ###
######
index = which(rownames(Mouse_Project_new_server@cellColData) %in% rownames(Meta) == T)

Mouse_Project_new_server_cl = Mouse_Project_new_server[rownames(Mouse_Project_new_server@cellColData)[index]]

m = match(rownames(Mouse_Project_new_server_cl@cellColData),rownames(Meta))

Mouse_Project_new_server_cl$celltype = Meta$celltype[m]
Mouse_Project_new_server_cl$celltype2 = Meta$celltype2[m]

Mouse_Project_new_server_cl <- addGroupCoverages(ArchRProj = Mouse_Project_new_server_cl, groupBy = "celltype2",threads = 20)

saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_Nov2024")

####
####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_Nov2024")

pathToMacs2 <- findMacs2()
######## BiocManager::install("BSgenome.Mmusculus.UCSC.mm10") ####

library(BSgenome.Mmusculus.UCSC.mm10)

Mouse_Project_new_server_cl <- addReproduciblePeakSet(
    ArchRProj = Mouse_Project_new_server_cl, 
    groupBy = "celltype2", 
    pathToMacs2 = pathToMacs2,
    genomeSize = 1.87e9
)

#####
#####
#####

Mouse_Project_new_server_cl <- addPeakMatrix(Mouse_Project_new_server_cl)
saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_Nov2024")

head(Mouse_Project_new_server_cl@cellColData)
table(Mouse_Project_new_server_cl$Sample)
table(Mouse_Project_new_server_cl$celltype)

Peak_Avg = Get_Peak_Avg_matrix(Mouse_Project_new_server_cl)
#### Peak_Avg = Get_Peak_Avg_matrix(All_project_cl)
dim(Peak_Avg)
saveRDS(Peak_Avg,file="Mouse_Peak_Avg_Nov18.rds")

Get_Gene_Avg_matrix <- function(All_project_cl){
    #########
    peakMatrix <- getMatrixFromProject(All_project_cl, useMatrix = "GeneIntegrationMatrix")
    Peak_matrix_count = peakMatrix@assays@data[[1]]
    Peak_matrix_row = as.character(peakMatrix@elementMetadata$name)
    rownames(Peak_matrix_count) <- Peak_matrix_row
    #colnames(Peak_matrix_count)
    ########
    cellGroups <- All_project_cl$celltype2
    averagePeakByCellType <- sapply(unique(cellGroups), function(cluster) {
        # 找到属于当前 Cluster 的细胞索引
        cellIdx <- which(cellGroups == cluster)
        # 提取属于该 Cluster 的 Peak 值并计算平均值
        Matrix::rowSums(Peak_matrix_count[, cellIdx])
    })
    #####
    factor = colSums(averagePeakByCellType) /1e6
    averagePeakByCellType = sweep(averagePeakByCellType,2,factor,FUN="/")
    #####
    return(averagePeakByCellType)
}

Gene_Avg = Get_Gene_Avg_matrix(Mouse_Project_new_server_cl)
#### Peak_Avg = Get_Peak_Avg_matrix(All_project_cl)
dim(Gene_Avg)
Gene_Avg[which(rownames(Gene_Avg) == "Trnp1"),]
saveRDS(Gene_Avg,file="Mouse_Gene_Avg_Nov18.rds")

#####
##### Gene-peak integration ###
##### 有问题 ！！！！############
##### peak numbers ######
##### Finished Creating Union Peak Set (301739)!, 1.802 mins elapsed #####
#####

table(Mouse_Project_new_server_cl$celltype2)

#####
##### reintegrate the cells ！！！！！ #########
#####

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
load("Mouse_RNA_seurat_merge_cl")

##### 这个应该没问题 ！！！ #######
#####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure2")
load("Mouse_before_seurat")

k = which(rownames(Mouse_before_seurat) %in% rownames(Mouse_RNA_seurat_merge_cl) == F)
rownames(Mouse_before_seurat)[k]

######
###### Next we will merge the Mouse_before_seurat and Mouse_RNA_seurat_merge_cl and read the raw matrix ########
######

dim(Mouse_before_seurat)







UMAPs = Mouse_RNA_seurat_merge_cl@reductions$umap@cell.embeddings
UMAPs = data.frame(UMAPs)
UMAPs$celltype = Mouse_RNA_seurat_merge_cl$celltype

write.table(UMAPs,file="Mouse_UMAPs_2025.txt",sep="\t",row.names=F,quote=F)

#######
#######








Mouse_RNA_seurat_merge_cl <- UpdateSeuratObject(RNA_seurat_merge_cl)
head(Mouse_RNA_seurat_merge_cl@meta.data)
load("RNA_seurat_meta_cl")
table(RNA_seurat_meta_cl$new_celltypes)

######
######

Mouse_RNA_seurat_merge_cl$cell_id = colnames(Mouse_RNA_seurat_merge_cl)
m = match(Mouse_RNA_seurat_merge_cl$cell_id,RNA_seurat_meta_cl$cell_id)
Mouse_RNA_seurat_merge_cl$celltype = RNA_seurat_meta_cl$new_celltypes[m]
Mouse_RNA_seurat_merge_cl$cell_id = gsub('#','.',Mouse_RNA_seurat_merge_cl$cell_id,fixed=T)
Mouse_RNA_seurat_merge_cl$time = sapply(strsplit(Mouse_RNA_seurat_merge_cl$cell_id,split=".",fixed=T),function(x) x[[1]])

######
###### OK!!! Let integrate each time points #######
######

table(Mouse_Project_new_server_cl$Sample)
table(Mouse_RNA_seurat_merge_cl$time)


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

###########-------------
conda activate ArchR
R

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2")
#### saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_cl_Nov18.rds")
#### saveRDS(Mouse_RNA_seurat_merge_cl,file="Mouse_RNA_seurat_merge_cl_Nov18.rds")

Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_cl_Nov18.rds")
Mouse_RNA_seurat_merge_cl <- readRDS("Mouse_RNA_seurat_merge_cl_Nov18.rds")

library(ArchR)
library(Seurat)
######--------------- E11 ------------------------------------------------------------------------------------------------------------------------------------------------------

E11_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E11")]
E11_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E11") == T)]
DefaultAssay(E11_rna) <- "RNA"
E11_rna <- NormalizeData(E11_rna)

######
######

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
######--------------- E11 ------------------------------------------------------------------------------------------------------------------------------------------------------

######--------------- E12 F ------------------------------------------------------------------------------------------------------------------------------------------------------
table(Mouse_Project_new_server_cl$Sample)
E12_project = Mouse_Project_new_server_cl[which(Mouse_Project_new_server_cl$Sample == "E12")]
table(Mouse_RNA_seurat_merge_cl$time)
E12_rna = Mouse_RNA_seurat_merge_cl[,which(Mouse_RNA_seurat_merge_cl$time %in% c("E12_rep1") == T)]
DefaultAssay(E12_rna) <- "RNA"
E12_rna <- NormalizeData(E12_rna)
######
######

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

######--------------- E12 ------------------------------------------------------------------------------------------------------------------------------------------------------


######--------------- E14 ------------------------------------------------------------------------------------------------------------------------------------------------------


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



##### E16 F ----------------------------------------------------------------------------
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

##### E18 F ----------------------------------------------------------------------------
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


##### P0 F ----------------------------------------------------------------------------
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


##### P2 F ----------------------------------------------------------------------------
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


##### P5 F ----------------------------------------------------------------------------
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


##### P8 F ----------------------------------------------------------------------------

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


####### we will use what ???? ####
####### 1. For both Mouse and 13LGS, we will only use RPC, NG, and photoreceptor cells from late-stage time points to calculate the PtoG correlation and peak accessibility.
########

######## switch to ArchR ####
conda activate ArchR
R
library(ArchR)

saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_cl_Nov2024")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_Nov2024")

Gene_Avg = Get_Gene_Avg_matrix(Mouse_Project_new_server_cl)
#### Peak_Avg = Get_Peak_Avg_matrix(All_project_cl)
dim(Gene_Avg)
Gene_Avg[which(rownames(Gene_Avg) == ""),]
saveRDS(Gene_Avg,file="Mouse_Gene_Avg_Nov18.rds")

Gene_Avg <- readRDS("Mouse_Gene_Avg_Nov18.rds")
Gene_Avg[which(rownames(Gene_Avg) == "Trnp1"),]

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
    dimsToUse = 1:30
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
plotPDF(p1, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Mouse_Project_new_server_cl, addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(
    ArchRProj = Mouse_Project_new_server_cl, 
    colorBy = "GeneIntegrationMatrix", 
    name = c("Opn1sw"), 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
plotPDF(p, name = "Plot-UMAP-Sample-Markers2.pdf", ArchRProj = Mouse_Project_new_server_cl, addDOC = FALSE, width = 5, height = 5)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
saveRDS(Mouse_Project_new_server_cl,file="Mouse_Project_new_server_Nov18_2024")


cell_index_meta = Mouse_Project_new_server_cl@cellColData
table(cell_index_meta$Sample)
cell_index_meta_cl = cell_index_meta[which(cell_index_meta$Sample %in% c("E18","P0","P2","P5","P8") == T),]
table(cell_index_meta_cl$Sample)
table(cell_index_meta_cl$celltype2)
cell_index_meta_cl = cell_index_meta_cl[which(cell_index_meta_cl$celltype2 %in% c("10_Cone","11_Early_Rod","12_Rod","2_RPCs_S2","3_RPCs_S3","4_MG","5_Early_NG","6_Late_NG","9_Early_Cone") == T),]
table(cell_index_meta_cl$celltype2)
getAvailableMatrices(Mouse_Project_new_server_cl)
cell_list = rownames(cell_index_meta_cl)

Mouse_Project_new_server_cl_PtoG_Sub = Mouse_Project_new_server_cl[cell_list]

#########
#########
#########

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



p1 <- plotEmbedding(ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1, name = "sub_CT.pdf", ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub, addDOC = FALSE, width = 5, height = 5)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl_PtoG_Sub <- readRDS("Mouse_Project_new_server_cl_PtoG_Sub_Nov18_2024")

Mouse_Project_new_server_cl_PtoG_Sub <- addPeak2GeneLinks(
    ArchRProj = Mouse_Project_new_server_cl_PtoG_Sub,
    reducedDims = "Harmony_Late",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)

###
saveRDS(Mouse_Project_new_server_cl_PtoG_Sub,file="Mouse_Project_new_server_cl_PtoG_Sub_Nov18_2024")


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


Mouse_Project_new_server_cl_PtoG_Sub_RES = Get_p2g_fun(Mouse_Project_new_server_cl_PtoG_Sub)
saveRDS(Mouse_Project_new_server_cl_PtoG_Sub_RES,file="Mouse_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl_PtoG_Sub_RES <- readRDS("Mouse_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")

Mouse_Project_new_server_cl_PtoG_Sub_RES[which(Mouse_Project_new_server_cl_PtoG_Sub_RES$gene=="Onecut1"),]
Mouse_Project_new_server_cl_PtoG_Sub_RES[which(Mouse_Project_new_server_cl_PtoG_Sub_RES$gene=="Onecut2"),]
Mouse_Project_new_server_cl_PtoG_Sub_RES[which(Mouse_Project_new_server_cl_PtoG_Sub_RES$gene=="Sall3"),]



#################
##################

#####----------------------------------------------------------------------------------------------------------------------------------------#####


#####
#####
#####

setwd("/zp1/data/plyu3/Old_Server_Data/LGS_13/Convert")


##### 
##### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### set the environment #####
#####

conda activate ArchR2



###########
########### for 13LGS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###########
########### 我们看一下 LGS的数据 #########
###########
########### OK!!! Done!!!! ######
###########

########### Next we will reprocess the 13LGS datasets !!!!! ###########
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

###########-------------首先需要 call peaks ######---------------------------------------------------------------
conda activate ArchR2 
R

library(ArchR)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
All_project_cl <- readRDS("All_project_cl_2024_LGS")
head(All_project_cl@cellColData)
table(All_project_cl$new_celltypes)
table(All_project_cl$timepoint)

###########-----------------------------------------------------------------------------------------------------

Meta <- All_project_cl@cellColData
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(Meta,file="LGS_Meta_scATAC_2024Nov13")
Meta <- readRDS("LGS_Meta_scATAC_2024Nov13")

####------------------------------------------------------------------------------------------------------------------------------------
#### add new project with arrow files: #####
####-------------------------------------------------------------------------------------------------------------------------------------
ArrowFiles <- list.files()
ArrowFiles <- ArrowFiles[grep(".arrow$",ArrowFiles)]
names = ArrowFiles
names = gsub(".arrow","",names)
names(ArrowFiles) <- names

#### install the LGS Genome ######

library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/ArchR_datasets')
load('genomeAnnotation_13LGS_cl')
setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/ArchR_datasets')
load('geneAnnotation_13LGS_cl')
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")

LGS13_Project_new_server <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "LGS13_Project_new_server",
  copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  geneAnnotation = geneAnnotation_13LGS,
  genomeAnnotation = genomeAnnotation_13LGS
)

saveRDS(LGS13_Project_new_server,file="LGS13_Project_new_server_Nov2024")

###### filter the project by meta ###
######
index = which(rownames(LGS13_Project_new_server@cellColData) %in% rownames(Meta) == T)

LGS13_Project_new_server_cl = LGS13_Project_new_server[rownames(LGS13_Project_new_server@cellColData)[index]]

m = match(rownames(LGS13_Project_new_server_cl@cellColData),rownames(Meta))

LGS13_Project_new_server_cl$celltype = Meta$new_celltypes[m]
LGS13_Project_new_server_cl$celltype2 = Meta$celltype2[m]

LGS13_Project_new_server_cl <- addGroupCoverages(ArchRProj = LGS13_Project_new_server_cl, groupBy = "celltype2",threads = 20)

saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_Nov2024")

####
####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

###########-------------首先需要 call peaks ######---------------------------------------------------------------
conda activate ArchR2 
R

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")

library(ArchR)
pathToMacs2 <- findMacs2()
######## BiocManager::install("BSgenome.Mmusculus.UCSC.mm10") ####

######## library(BSgenome.Mmusculus.UCSC.mm10)

sum(width(genomeAnnotation_13LGS$chromSizes))

LGS13_Project_new_server_cl <- addReproduciblePeakSet(
    ArchRProj = LGS13_Project_new_server_cl, 
    groupBy = "celltype2", 
    pathToMacs2 = pathToMacs2,
    genomeSize = 2.45e9
)

#####
##### Finished Creating Union Peak Set (428909)!, 4.561 mins elapsed. ######
#####

LGS13_Project_new_server_cl <- addPeakMatrix(LGS13_Project_new_server_cl)
saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_Nov2024")


head(LGS13_Project_new_server_cl@cellColData)
table(LGS13_Project_new_server_cl$Sample)
table(LGS13_Project_new_server_cl$celltype)

Peak_Avg = Get_Peak_Avg_matrix(LGS13_Project_new_server_cl)
#### Peak_Avg = Get_Peak_Avg_matrix(All_project_cl)
dim(Peak_Avg)
saveRDS(Peak_Avg,file="LGS13_Peak_Avg_Nov18.rds")

Peak_Avg <- readRDS("LGS13_Peak_Avg_Nov18.rds")

Peak_Avg_chr = sapply(strsplit(rownames(Peak_Avg),split=":"),function(x) x[[1]])
Peak_Avg_chr = Peak_Avg_chr[!duplicated(Peak_Avg_chr)]

LGS13_Project_new_server_cl
#####
##### ##### ##### ---------------------------------------------------------------------------------------------------------------------
#####

##### Next we will add integration RNA matrix ########
##### redo !!! ######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge")
all_seurat_RNA_big_cl_detail_new <- readRDS("all_seurat_RNA_big_cl_detail_new_202204")

head(all_seurat_RNA_big_cl_detail_new@meta.data)
table(all_seurat_RNA_big_cl_detail_new$new_celltypes)
table(all_seurat_RNA_big_cl_detail_new$timepoint)

#####
##### OK!!!! Let us integrate !!!!! #############
#####

##### UNiform ####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
library(Seurat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")
all_seurat_RNA_big_cl_detail_new <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/13LGS_RNAseq_seurat/RNA_merge/all_seurat_RNA_big_cl_detail_new_202204")

########

table(LGS13_Project_new_server_cl$Sample)
E18_1 E18_2A E18_2B  E21_1  E21_2 E21_3B  E24_1  E24_2  E26_1  E26_2  P1_1 
2306  3091   2216    747    2317  3176    7354   5605   6019   3136   3568 
P1_2  P12_1  P12_2  P17_1  P17_2  P21_1  P21_2   P4_1   P4_2   P8_1   P8_2 
6619  8195   12941   8270   6001  10999  5946    3229   6606   8921   9667 

table(all_seurat_RNA_big_cl_detail_new$timepoint)
E18  E21  E24  E26  P1  P12  P17  P21  P4  P8 
5522 3407 4659 6042 6920 8221 5277 5057 3908 4381 

######## test: ####################
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E18_1","E18_2A","E18_2B")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E18")

######## function: ##########
Integrate_Fun <- function(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples){
    ########
    ########
    tmp_project = ArchR_project[which(ArchR_project$Sample %in% ATAC_samples == T)]
    tmp_rna = RNA_seurat[,which(RNA_seurat$timepoint %in% RNA_samples == T)]
    DefaultAssay(tmp_rna) <- "RNA"
    tmp_rna <- NormalizeData(tmp_rna)
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

########
######## E18 Done !!!!--------- ######
########

ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E18_1","E18_2A","E18_2B")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E18")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

######## E21 Done !!!!! --------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E21_1","E21_2","E21_3B")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E21")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

######## E24 Done --------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E24_1","E24_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E24")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

######## E26 Done ！！！--------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("E26_1","E26_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("E26")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

######## P1 Done --------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P1_1","P1_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P1")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P4 Done !!!! ----------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P4_1","P4_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P4")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P8 Done !!!! ----------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P8_1","P8_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P8")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P12 Done ！！！----------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P12_1","P12_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P12")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

####### P17 Done！！！！----------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P17_1","P17_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P17")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)


####### P21 Done！！！----------- #####
ArchR_project = LGS13_Project_new_server_cl
ATAC_samples = c("P21_1","P21_2")
RNA_seurat = all_seurat_RNA_big_cl_detail_new
RNA_samples = c("P21")
Integrate_Fun(ArchR_project,ATAC_samples,RNA_seurat,RNA_samples)

#######
#######
####### OK!!! Let us see the total things !!! ####
####### 
#######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR 
R
library(ArchR)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")

#####
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")
#####

LGS13_Project_new_server_cl <- addIterativeLSI(
    ArchRProj = LGS13_Project_new_server_cl,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
#####
LGS13_Project_new_server_cl <- addHarmony(
    ArchRProj = LGS13_Project_new_server_cl,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

#####
LGS13_Project_new_server_cl <- addUMAP(
    ArchRProj = LGS13_Project_new_server_cl, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)


p1 <- plotEmbedding(ArchRProj = LGS13_Project_new_server_cl, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1, name = "LGS13-UMAP-Sample-Clusters.pdf", ArchRProj = LGS13_Project_new_server_cl, addDOC = FALSE, width = 5, height = 5)


p <- plotEmbedding(
    ArchRProj = LGS13_Project_new_server_cl, 
    colorBy = "GeneIntegrationMatrix", 
    name = c("Aqp4"), 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
plotPDF(p, name = "Plot-UMAP-Sample-Markers3.pdf", ArchRProj = LGS13_Project_new_server_cl, addDOC = FALSE, width = 5, height = 5)

####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Project_new_server_cl,file="LGS13_Project_new_server_Nov2024")
####
####
####

LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")

######## BiocManager::install("BSgenome.Mmusculus.UCSC.mm10") ####
########--------------- add new peaks to the LGS13 objects ------------------------------------------------------------------------------------------------------------------------------------------------------------

cell_index_meta = LGS13_Project_new_server_cl@cellColData
table(cell_index_meta$Sample)
cell_index_meta_cl = cell_index_meta[which(cell_index_meta$Sample %in% c("P1_1","P1_2","P12_1","P12_2","P17_1","P17_2","P21_1","P21_2","P4_1","P4_2","P8_1","P8_2") == T),]
table(cell_index_meta_cl$Sample)
table(cell_index_meta_cl$celltype2)
cell_index_meta_cl = cell_index_meta_cl[which(cell_index_meta_cl$celltype2 %in% c("1_Early_RPC","10_Cone","11_Rod","2_Late_RPC","3_MG","4_Early_NG","5_Late_NG","9_Photoprecursor") == T),]
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
plotPDF(p1, name = "sub_CT2.pdf", ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub, addDOC = FALSE, width = 5, height = 5)

saveRDS(LGS13_Project_new_server_cl_PtoG_Sub,file="LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")


#########----------------------------------------------------------------------------------------------------------------------------------------------------------

######### go to ArchR2 to try ####
#########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")

LGS13_Project_new_server_cl_PtoG_Sub <- readRDS("LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")

#### Mouse_Project_new_server_cl_PtoG_Sub <- readRDS("Mouse_Project_new_server_cl_PtoG_Sub_Nov2024")

LGS13_Project_new_server_cl_PtoG_Sub <- addPeak2GeneLinks(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub,
    reducedDims = "Harmony_Late",
    useMatrix = "GeneIntegrationMatrix",
    dimsToUse = 2:20,
    knnIteration = 500,
    maxDist = 500000
    #### ######
)

saveRDS(LGS13_Project_new_server_cl_PtoG_Sub,file="LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")

#
#Mouse_Project_new_server_cl_PtoG_Sub@geneAnnotation[[1]]
#LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation[[1]]
#Mouse_Project_new_server_cl_PtoG_Sub@geneAnnotation[[2]]
#LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation[[2]]
#Mouse_Project_new_server_cl_PtoG_Sub@geneAnnotation[[3]]
#LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation[[3]]
#Mouse_Project_new_server_cl_PtoG_Sub@geneAnnotation[[3]]
#LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation[[3]]
#Mouse_Project_new_server_cl_PtoG_Sub@peakSet
#LGS13_Project_new_server_cl_PtoG_Sub@peakSet
#

x = LGS13_Project_new_server_cl_PtoG_Sub

Get_p2g_fun <- function(x){
	corCutOff = 0.25
	FDRCutOff = 1e-6
	varCutOffATAC = 0.7
	varCutOffRNA = 0.3
	p2g <- metadata(x@peakSet)$Peak2GeneLinks
    ########
    metadata(p2g)
    ########
	p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
	if(!is.null(varCutOffATAC)){
    	p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
	}
	if(!is.null(varCutOffRNA)){
    	p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
	}
    ####  mRNA <- readRDS(metadata(p2g)$seRNA)
    ####  gene = mRNA@rowRanges$name
    ####  gene[grep("_",gene)]
    ####
	mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
	mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
	p2g$peak <- paste0(rowRanges(mATAC))
	p2g$gene <- rowData(mRNA)$name
	return(p2g)
}



Get_p2g_fun2 <- function(x){
	corCutOff = 0.25
	FDRCutOff = 1e-6
	varCutOffATAC = 0.7
	varCutOffRNA = 0.3
	p2g <- metadata(x@peakSet)$Peak2GeneLinks
    ########
    Peaks = metadata(p2g)$peakSet
    Genes = metadata(p2g)$geneSet
    ########
	p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
	if(!is.null(varCutOffATAC)){
    	p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
	}
	if(!is.null(varCutOffRNA)){
    	p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
	}
    ####  mRNA <- readRDS(metadata(p2g)$seRNA)
    ####  gene = mRNA@rowRanges$name
    ####  gene[grep("_",gene)]
    ####
	p2g$peak <- as.character(Peaks)[p2g$idxATAC]
	p2g$gene <- Genes$name[p2g$idxRNA]
	return(p2g)
    #### genes = p2g$gene; genes = genes[!duplicated(genes)]
}


###
###
###

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl_PtoG_Sub <- readRDS("LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")

###
### gene = LGS13_Project_new_server_cl_PtoG_Sub@geneAnnotation$genes$symbol #####
### 变成下划线了 ！！！ ##############
### "TCAF1_like" ##########
### gene = LGS13_Project_new_server_cl_PtoG_Sub_RES$gene
### gene = gene[!duplicated(gene)]
### grep("like",gene)

LGS13_Project_new_server_cl_PtoG_Sub_RES = Get_p2g_fun(LGS13_Project_new_server_cl_PtoG_Sub)
LGS13_Project_new_server_cl_PtoG_Sub_RES2 = Get_p2g_fun2(LGS13_Project_new_server_cl_PtoG_Sub)

all.equal(LGS13_Project_new_server_cl_PtoG_Sub_RES,LGS13_Project_new_server_cl_PtoG_Sub_RES2)

p2g <- getPeak2GeneLinks(
    ArchRProj = LGS13_Project_new_server_cl_PtoG_Sub,
    corCutOff = 0.1,
    resolution = 1,
    returnLoops = FALSE
)

###
saveRDS(LGS13_Project_new_server_cl_PtoG_Sub_RES,file="LGS13_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")

LGS13_Project_new_server_cl_PtoG_Sub_RES[which(LGS13_Project_new_server_cl_PtoG_Sub_RES$gene=="Onecut1"),]
LGS13_Project_new_server_cl_PtoG_Sub_RES[which(LGS13_Project_new_server_cl_PtoG_Sub_RES$gene=="Onecut2"),]
LGS13_Project_new_server_cl_PtoG_Sub_RES[which(LGS13_Project_new_server_cl_PtoG_Sub_RES$gene=="Sall3"),]

########----------OK！！！！-------- Next we will update the conserved peaks between Mouse and 13LGS ######
########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate ArchR
R
library(ArchR)
#### line:1549 ####
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")

Combined_alignment <- readRDS("Combined_alignment")

####
#### line:1549 ####
#### get the peak set from both Mouse and 13LGS #####
####


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl_PtoG_Sub <- readRDS("LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")

LGS13_PtoG_Peaks_SUB <- LGS13_Project_new_server_cl_PtoG_Sub@peakSet
LGS13_PtoG_Peaks_ALL <- LGS13_Project_new_server_cl@peakSet

LGS13_PtoG_Peaks_SUB <- as.character(LGS13_PtoG_Peaks_SUB)
LGS13_PtoG_Peaks_ALL <- as.character(LGS13_PtoG_Peaks_ALL)

all.equal(LGS13_PtoG_Peaks_SUB,LGS13_PtoG_Peaks_ALL)
saveRDS(LGS13_PtoG_Peaks_ALL,"LGS13_Project_FINALPEAKS_Nov2024")

LGS13_PtoG_Peaks_ALL <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_FINALPEAKS_Nov2024")
######
###### 429178 ####
###### for Mouse ######
######


ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate ArchR
R
library(ArchR)

#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl_PtoG_Sub <- readRDS("Mouse_Project_new_server_cl_PtoG_Sub_Nov18_2024")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_Nov2024")

Mouse_PtoG_Peaks_SUB <- Mouse_Project_new_server_cl_PtoG_Sub@peakSet
Mouse_PtoG_Peaks_ALL <- Mouse_Project_new_server_cl@peakSet

Mouse_PtoG_Peaks_SUB <- as.character(Mouse_PtoG_Peaks_SUB)
Mouse_PtoG_Peaks_ALL <- as.character(Mouse_PtoG_Peaks_ALL)

all.equal(Mouse_PtoG_Peaks_SUB,Mouse_PtoG_Peaks_ALL)
saveRDS(Mouse_PtoG_Peaks_ALL,"Mouse_Project_FINALPEAKS_Nov2024")

#########
#########
#########
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
Combined_alignment <- readRDS("Combined_alignment")

###### Map Mouse peaks with Mouse alignment #####
M_tab <- findOverlaps(GRanges(Mouse_PtoG_Peaks_ALL),GRanges(Combined_alignment$Mouse_peak))
M_tab <- data.frame(M_tab)
M_tab$M_peaks = Mouse_PtoG_Peaks_ALL[M_tab$queryHits]
M_tab$M_Align = Combined_alignment$Mouse_peak[M_tab$subjectHits]

###### Map Mouse align with LGS13 align #####
m = match(M_tab$M_Align,Combined_alignment$Mouse_peak)
M_tab$LGS_Align = Combined_alignment$LGS_peak[m]

###### remove the number indexs #####
M_tab = M_tab[,c(-1,-2)]

##### Map LGS peaks with LGS alignment #####
L_tab <- findOverlaps(GRanges(LGS13_PtoG_Peaks_ALL),GRanges(Combined_alignment$LGS_peak))
L_tab <- data.frame(L_tab)
L_tab$LGS_peaks = LGS13_PtoG_Peaks_ALL[L_tab$queryHits]
L_tab$LGS_Align = Combined_alignment$LGS_peak[L_tab$subjectHits]
L_tab = L_tab[,c(-1,-2)]

###### Remove the number indexs #####
MergeTab = merge(M_tab,L_tab)

######
MergeTab = MergeTab[,c("M_Align","LGS_Align","M_peaks","LGS_peaks")]
MergeTab$Strand = as.character(strand(GRanges(MergeTab$LGS_Align)))

#####
##### Next we will get the sequence ######
#####

library(BSgenome.Mmusculus.UCSC.mm10)
library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")

###### MergeTab_Plus = MergeTab[which(MergeTab$)]

M_seq = getSeq(BSgenome.Mmusculus.UCSC.mm10, GRanges(MergeTab$M_peaks))
M_seq_v = as.character(M_seq)

library(Biostrings)
MergeTab$M_seq = M_seq_v
#####
#####
MergeTab_Plus = MergeTab[which(MergeTab$Strand == "+"),]
MergeTab_Minus = MergeTab[which(MergeTab$Strand == "-"),]

L_seq_Plus = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(MergeTab_Plus$LGS_peaks))
L_seq_Plus_v = as.character(L_seq_Plus)
MergeTab_Plus$L_seq = L_seq_Plus_v

L_seq_Minus = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(MergeTab_Minus$LGS_peaks))
L_seq_Minus_Rev <- reverseComplement(L_seq_Minus)
L_seq_Minus_v = as.character(L_seq_Minus_Rev)
MergeTab_Minus$L_seq = L_seq_Minus_v

MergeTab_All = rbind(MergeTab_Plus,MergeTab_Minus)

MergeTab_All_List <- apply(MergeTab_All, 1, function(row) as.list(row))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
saveRDS(MergeTab_All_List,file="Mouse_LGS13_all_peaks_alignment_Nov25.rds")


#####--------------------- reload the results ------------------
MergeTab_All_List <- readRDS("Mouse_LGS13_all_peaks_alignment_Nov25.rds")

#####
compare_2_seq <- function(x){
    ######
    library(Biostrings)
    seq1 = DNAString(x$M_seq)
    seq2 = DNAString(x$L_seq)
    ######
    alignment <- pairwiseAlignment(seq1, seq2, type = "local")
    start_pattern <- start(subject(alignment))
    end_pattern <- end(subject(alignment))      
    start_subject <- start(pattern(alignment))  
    end_subject <- end(pattern(alignment))   
    length_pattern <- end_pattern - start_pattern + 1
    length_subject <- end_subject - start_subject + 1
    #########
    length = min(length_pattern,length_subject)
    score = round(score(alignment),3)
    #########
    out = paste0(length,":",score)
    #########
    return(out)
}

#####
library(pbapply)
library(parallel)
cl <- makeCluster(35) 
alignment_res = parLapply(cl, MergeTab_All_List, compare_2_seq)
save(alignment_res,file="ML_alignment_res_Nov25")
#####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate ArchR
R
library(ArchR)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
load("ML_alignment_res_Nov25")

strsplit(alignment_res[[1]],split=":")[[1]][1]

overlap_counts = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][1]))
overlap_score = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][2]))

MergeTab_All_List <- readRDS("Mouse_LGS13_all_peaks_alignment_Nov25.rds")
names(MergeTab_All_List[[1]])

M_Align <- sapply(MergeTab_All_List,function(x) x$M_Align)
LGS_Align <- sapply(MergeTab_All_List,function(x) x$LGS_Align)
M_peaks <- sapply(MergeTab_All_List,function(x) x$M_peaks)
LGS_peaks <- sapply(MergeTab_All_List,function(x) x$LGS_peaks)
Strand <- sapply(MergeTab_All_List,function(x) x$Strand)
M_seq <- sapply(MergeTab_All_List,function(x) x$M_seq)
L_seq <- sapply(MergeTab_All_List,function(x) x$L_seq)

MergeTab_All <- data.frame(M_Align=M_Align,LGS_Align=LGS_Align,M_peaks=M_peaks,LGS_peaks=LGS_peaks,Strand=Strand,M_seq=M_seq,L_seq=L_seq)

MergeTab_All$Overlap_region = overlap_counts
MergeTab_All$Overlap_score = overlap_score


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
saveRDS(MergeTab_All,file="Mouse_LGS13_all_peaks_alignment_Nov25.rds")




#########----------------------------------------------------------------------------------------------------------------------------------------------------------
######### Next we will see the Histone modifications ######
#########----------------------------------------------------------------------------------------------------------------------------------------------------------
###########

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
###########-------------首先需要 call peaks ######---------------------------------------------------------------
conda activate ArchR2 
R

#####
##### add the Histones to the new peaks !!!! #######
#####

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
conda activate ArchR
R
library(ArchR)


#####----------------------------------------------------------------------------------------------------------------------------------------------------------
##### first find the total peaks for 13LGS !!!!! 
#####-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

LGS13_Peaks_ALL <- readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_FINALPEAKS_Nov2024")

signal_files = list.files("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert/LGSCutRun")
signal_files = signal_files[grep("rds",signal_files)]
signal_files_N = gsub(".norm.rds","",signal_files)

total_peaks = LGS13_Peaks_ALL

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert/LGSCutRun")

add_signals_by_folder <- function(total_peaks,signal_files,signal_files_N){
    library(GenomicRanges)
    #########
    total_GR = GRanges(total_peaks)
    total_GR_tab = data.frame(peaks=total_peaks)
    tmp_mat = matrix(0,nrow=length(total_peaks),ncol=length(signal_files))
    colnames(tmp_mat) <- signal_files_N
    rownames(tmp_mat) <- total_peaks
    total_GR_tab = cbind(total_GR_tab,tmp_mat)
    #########
    for(i in 1:length(signal_files)){
        #######
        print(signal_files[i])
        tmp_signal <- readRDS(signal_files[i])
        #######
        k = which(countOverlaps(tmp_signal,total_GR) > 0)
        tmp_signal_cl = tmp_signal[k]
        #######
        tmp_signal_cl_1bp <- unlist(tile(tmp_signal_cl, width = 1))
        tmp_signal_cl_1bp$score <- rep(tmp_signal_cl$score, times = width(tmp_signal_cl))
        ####### add signals ######
        f = findOverlaps(total_GR,tmp_signal_cl_1bp)
        #######
        f = data.frame(f)
        f$score = tmp_signal_cl_1bp$score[f$subjectHits]
        #######
        f_res = tapply(f$score,f$queryHits,mean)
        f_res_tab = data.frame(index=names(f_res),score = as.numeric(f_res))
        #######
        total_GR_tab[as.numeric(f_res_tab$index),i+1] = f_res_tab$score
        ########
    }
    return(total_GR_tab)
}

LGS13_peaks_signals <- add_signals_by_folder(total_peaks,signal_files,signal_files_N)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
saveRDS(LGS13_peaks_signals,file="LGS13_peaks_signals_Nov2024")

######---------#######
######---------#######
######---------#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Peaks_ALL <- readRDS("Mouse_Project_FINALPEAKS_Nov2024")

total_peaks = Mouse_Peaks_ALL

signal_files = list.files("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert/MouseCutRun")
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert/MouseCutRun")

signal_files = signal_files[grep("rds",signal_files)]
signal_files_N = gsub(".norm.rds","",signal_files)
Mouse_peaks_signals <- add_signals_by_folder(total_peaks,signal_files,signal_files_N)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
saveRDS(Mouse_peaks_signals,file="Mouse_peaks_signals_Nov2024")


######---------#######
######---------Next for the average peaks accessaible score---------#######
######---------#######
######---------for the LGS samples #######------------------------------------------------------------------------------------------






setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Gene_Anno <- readRDS("LGS13_Gene_Anno_Oct23")

###### we will recall the peak average ！！！！ ######
######


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl <- readRDS("LGS13_Project_new_server_Nov2024")

All_project_cl = LGS13_Project_new_server_cl

Get_Peak_Avg_matrix <- function(All_project_cl){
    #########
    peakMatrix <- getMatrixFromProject(All_project_cl, useMatrix = "PeakMatrix")
    Peak_matrix_count = peakMatrix@assays@data[[1]]
    Peak_matrix_row = as.character(peakMatrix@rowRanges)
    rownames(Peak_matrix_count) <- Peak_matrix_row
    #colnames(Peak_matrix_count)
    ########
    cellGroups <- All_project_cl$celltype2
    averagePeakByCellType <- sapply(unique(cellGroups), function(cluster) {
        # 找到属于当前 Cluster 的细胞索引
        cellIdx <- which(cellGroups == cluster)
        # 提取属于该 Cluster 的 Peak 值并计算总的值
        Matrix::rowSums(Peak_matrix_count[, cellIdx])
    })
    #####
    factor = colSums(averagePeakByCellType) /1e6
    averagePeakByCellType = sweep(averagePeakByCellType,2,factor,FUN="/")
    #####
    return(averagePeakByCellType)
}

LGS13_Avg_PeakMat = Get_Peak_Avg_matrix(LGS13_Project_new_server_cl)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Avg_PeakMat,file="LGS13_Avg_PeakMat_Nov2024")

################### For Mouse .......... ########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_Nov2024")

Mouse_Avg_PeakMat = Get_Peak_Avg_matrix(Mouse_Project_new_server_cl)

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
saveRDS(Mouse_Avg_PeakMat,file="Mouse_Avg_PeakMat_Nov2024")

###################
################### PtoG ?????? ######################
################### Let us see the PtoG !!!! #########
###################


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Project_new_server_cl_PtoG_Sub <- readRDS("LGS13_Project_new_server_cl_PtoG_Sub_Nov2024")

LGS13_Project_new_server_cl_PtoG_Sub_RES = Get_p2g_fun(LGS13_Project_new_server_cl_PtoG_Sub)
saveRDS(LGS13_Project_new_server_cl_PtoG_Sub_RES,file="LGS13_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")

######
LGS13_PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")

######
###### This step we will add PtoG correlations #######
######


get_TSS <- function(Genes_sub){
    #####
    if(as.character(strand(Genes_sub)) == "+"){
        end(Genes_sub) = start(Genes_sub)
    }
    if(as.character(strand(Genes_sub)) == "-"){
        start(Genes_sub) = end(Genes_sub)
    }
    strand(Genes_sub) = "*"
    return(Genes_sub)
}

Function_1_find_peaks <- function(GeneName,Gtf,PtoG,PeakAvg,CTag){
    ##########
    Gene_ranges = Gtf[which(Gtf$symbol == GeneName)]
    Gene_TSS = get_TSS(Gene_ranges)
    strand(Gene_ranges) <- "*"
    ########## extend the TSS and merge the regions 500bp ##########
    start(Gene_TSS) <- start(Gene_TSS) - 500
    end(Gene_TSS) <- end(Gene_TSS) + 500
    ########## TSS 首先扩大一些 扩大 500bp ######
    k_TSS = which(countOverlaps(GRanges(rownames(PeakAvg)),Gene_TSS) > 0)
    if(length(k_TSS) > 0){
        TSS_peaks = rownames(PeakAvg)[k_TSS]
    }
    ##########
    Gene_Body = Gtf[which(Gtf$symbol == GeneName)]
    k_Body = which(countOverlaps(GRanges(rownames(PeakAvg)),Gene_Body) > 0)
    if(length(k_Body) > 0){
        Body_peaks = rownames(PeakAvg)[k_Body]
    }
    ##########
    ########## Next we will identify intergenetic peaks ###########
    ##########
    All_body = Gtf
    ########## identify all the TSS regions ####
    All_TSS_Plus = All_body[which(strand(All_body) == "+")]
    All_TSS_Minus = All_body[which(strand(All_body) == "-")]
    ##########
    end(All_TSS_Plus) = start(All_TSS_Plus)
    start(All_TSS_Minus) = end(All_TSS_Minus)
    All_TSS_merge = c(All_TSS_Plus,All_TSS_Minus)
    ##########
    strand(All_TSS_merge) <- "*"
    start(All_TSS_merge) <- start(All_TSS_merge) - 500
    end(All_TSS_merge) <- end(All_TSS_merge) + 500
    ########## All_TSS_merge[which(All_TSS_merge$symbol == "Thrb")] ####
    #All_TSS_merge_cl = All_TSS_merge[which(All_TSS_merge$symbol != GeneName)]
    #print(length(All_TSS_merge))
    #print(length(All_TSS_merge_cl))
    ##########
    Gene_ranges_all = union(All_TSS_merge,All_body)
    k_inter = which(countOverlaps(GRanges(rownames(PeakAvg)),Gene_ranges_all) == 0)
    ########## 
    Interpeaks_total <- rownames(PeakAvg)[k_inter]
    ##########
    ##########
    ########## target gene inter ####
    ########## 直接从 PtoG 里面找就行 #####
    ########## PtoG 500kb ######
    PtoG_cl = PtoG[which(PtoG$gene == GeneName),]
    PtoG_clcl = PtoG_cl[which(PtoG_cl$peak %in% Interpeaks_total == T),]
    ##########
    Potential_peaks_Inter = PtoG_clcl$peak
    ##########
    if(length(k_TSS) > 0){
        Total_Potential_peaks <- c(Potential_peaks_Inter,TSS_peaks)
    }
     if(length(k_Body) > 0){
        Total_Potential_peaks <- c(Total_Potential_peaks,Body_peaks)
    }
    Total_Potential_peaks = Total_Potential_peaks[!duplicated(Total_Potential_peaks)]
    ##########
    ########## Next we will add the class ####
    ##########
    Total_Potential_peaks_tab = data.frame(Gene=GeneName,Peak=Total_Potential_peaks,Class="NA",PtoG="NA",MAX_Access=0)
    ##########
    class1 = which(Total_Potential_peaks_tab$Peak %in% Interpeaks_total == T)
    Total_Potential_peaks_tab$Class[class1] = "inter"
    class2 = which(Total_Potential_peaks_tab$Peak %in% Body_peaks == T)
    Total_Potential_peaks_tab$Class[class2] = "body"
    class3 = which(Total_Potential_peaks_tab$Peak %in% TSS_peaks == T)
    Total_Potential_peaks_tab$Class[class3] = "tss"
    ##########
    ##########
    ########## Next we will add PtoG #############
    ##########
    m = match(Total_Potential_peaks_tab$Peak,PtoG_cl$peak)
    Total_Potential_peaks_tab$PtoG = PtoG_cl$Correlation[m]
    ########## Next we will add the Max access ####
    ##########
    ##########
    PeakAvg_cl = PeakAvg[,CTag]
    max_access = apply(PeakAvg_cl,1,max)
    ##########
    m2 = match(Total_Potential_peaks_tab$Peak,names(max_access))
    Total_Potential_peaks_tab$MAX_Access = as.numeric(max_access)[m2]
    ##########
    Total_Potential_peaks_tab_cl = Total_Potential_peaks_tab[which(Total_Potential_peaks_tab$MAX_Access > 1),]
    ###########
    keep1 = which(Total_Potential_peaks_tab_cl$Class == "tss")
    keep2 = which(Total_Potential_peaks_tab_cl$PtoG > 0.25)
    keep = c(keep1,keep2)
    keep = keep[!duplicated(keep)]
    Total_Potential_peaks_tab_cl = Total_Potential_peaks_tab_cl[keep,]
    ###########
    print(dim(Total_Potential_peaks_tab_cl))
    return(Total_Potential_peaks_tab_cl)
}

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
Gtf <- readRDS("LGS13_Gene_Anno_Oct23")
PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")
PeakAvg = readRDS("LGS13_Avg_PeakMat_Nov2024")
CTag = c("10_Cone","9_Photoprecursor","2_Late_RPC","5_Late_NG","3_MG","11_Rod")

GeneName = "Thrb"
LGS13_Thrb_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Mef2c"
LGS13_Mef2c_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Rxrg"
LGS13_Rxrg_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Pou2f1"
LGS13_Pou2f1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Onecut1"
LGS13_Onecut1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Onecut2"
LGS13_Onecut2_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Zic3"
LGS13_Zic3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Sall3"
LGS13_Sall3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)


########### Next the same gene for mouse !!!!! #######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
Mouse_Project_new_server_cl <- readRDS("Mouse_Project_new_server_Nov2024")
Mouse_Gene_Anno = Mouse_Project_new_server_cl@geneAnnotation$genes
saveRDS(Mouse_Gene_Anno,file="Mouse_Gene_Anno_Oct23")


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/Mouse_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")
PeakAvg = readRDS("Mouse_Avg_PeakMat_Nov2024")
CTag = c("3_RPCs_S3","12_Rod","10_Cone","4_MG","9_Early_Cone","11_Early_Rod")
Gtf <- readRDS("Mouse_Gene_Anno_Oct23")

GeneName = "Thrb"
Mouse_Thrb_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Mef2c"
Mouse_Mef2c_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Rxrg"
Mouse_Rxrg_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Pou2f1"
Mouse_Pou2f1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Onecut1"
Mouse_Onecut1_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Onecut2"
Mouse_Onecut2_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Zic3"
Mouse_Zic3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

GeneName = "Sall3"
Mouse_Sall3_peaks <- Function_1_find_peaks(GeneName,Gtf,PtoG,PeakAvg,CTag)

########
######## 下一步是 add epi markers ######
########


Function_2_add_scores <- function(Peaks,Peak_signals){
    ##########
    if(dim(Peaks)[1] == 0){
        return(Peaks)
    }
    m = match(Peaks$Peak,Peak_signals$peaks)
    ##########
    Peaks$H3K27ac = Peak_signals$H3K27ac_signal[m]
    Peaks$H3K4me3 = Peak_signals$H3K4me3_signal[m]
    Peaks$Otx2 = Peak_signals$Otx2_signal[m]
    Peaks$Neurod1 = Peak_signals$Neurod1_signal[m]
    ###########
    #### cutoff Otx2 > 1 ####
    #### cutoff Neurod > 1 ##
    #### cutoff H3K27ac > 1 #
    #### cutoff H3K4me3 > 5 #
    ###########
    ###########
    Peaks$epi_tag = "potential_regulatory_elements"
    k1 = which(Peaks$H3K27ac > 1 & Peaks$H3K4me3 < 5)
    Peaks$epi_tag[k1] = "Activate_enhancer"
    k2 = which(Peaks$H3K4me3 > 5)
    Peaks$epi_tag[k2] = "Promoter"
    ######
    Peaks$cut_run = ""
    ######
    k3 = which(Peaks$Otx2 > 2)
    Peaks$cut_run[k3] = paste(Peaks$cut_run[k3],"Otx2_binding")
    ######
    k4 = which(Peaks$Neurod1 > 2)
    Peaks$cut_run[k4] = paste(Peaks$cut_run[k4],"Neurod1_binding")
    ######
    ###### OK！！！####
    ######
    return(Peaks)
}


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
Mouse_peaks_signals = readRDS(file="Mouse_peaks_signals_Nov2024")

Mouse_peaks_signals$H3K27ac_signal = (Mouse_peaks_signals$H3K27ac_Rep1 + Mouse_peaks_signals$H3K27ac_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$H3K4me3_signal = (Mouse_peaks_signals$H3K4me3_Rep1 + Mouse_peaks_signals$H3K4me3_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$Otx2_signal = (Mouse_peaks_signals$Otx2_AbA_Rep1 + Mouse_peaks_signals$Otx2_AbA_Rep2 + Mouse_peaks_signals$Otx2_AbB_Rep1 + Mouse_peaks_signals$Otx2_AbB_Rep2)/4 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
Mouse_peaks_signals$Neurod1_signal = (Mouse_peaks_signals$Neurod1_Rep1 + Mouse_peaks_signals$Neurod1_Rep2)/2 - (Mouse_peaks_signals$IgG_Rep1 + Mouse_peaks_signals$IgG_Rep2)/2
saveRDS(Mouse_peaks_signals,file="Mouse_peaks_signals_Nov2024")


Peaks = Mouse_Thrb_peaks
Peak_signals = Mouse_peaks_signals

Mouse_Thrb_peaks_res <- Function_2_add_scores(Mouse_Thrb_peaks,Mouse_peaks_signals)
Mouse_Mef2c_peaks_res <- Function_2_add_scores(Mouse_Mef2c_peaks,Mouse_peaks_signals)
Mouse_Onecut1_peaks_res <- Function_2_add_scores(Mouse_Onecut1_peaks,Mouse_peaks_signals)
Mouse_Onecut2_peaks_res <- Function_2_add_scores(Mouse_Onecut2_peaks,Mouse_peaks_signals)
Mouse_Zic3_peaks_res <- Function_2_add_scores(Mouse_Zic3_peaks,Mouse_peaks_signals)
Mouse_Sall3_peaks_res <- Function_2_add_scores(Mouse_Sall3_peaks,Mouse_peaks_signals)
Mouse_Rxrg_peaks_res <- Function_2_add_scores(Mouse_Rxrg_peaks,Mouse_peaks_signals)
Mouse_Pou2f1_peaks_res <- Function_2_add_scores(Mouse_Pou2f1_peaks,Mouse_peaks_signals)

########
########
########

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
LGS13_peaks_signals <- readRDS("LGS13_peaks_signals_Nov2024")

LGS13_peaks_signals$H3K27ac_signal = (LGS13_peaks_signals$H3K27ac_Rep1 + LGS13_peaks_signals$H3K27ac_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$H3K4me3_signal = (LGS13_peaks_signals$H3K4me3_Rep1 + LGS13_peaks_signals$H3K4me3_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$Otx2_signal = (LGS13_peaks_signals$Otx2_AbA_Rep1 + LGS13_peaks_signals$Otx2_AbA_Rep2 + LGS13_peaks_signals$Otx2_AbB)/3 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2
LGS13_peaks_signals$Neurod1_signal = (LGS13_peaks_signals$Neurod1_Rep1 + LGS13_peaks_signals$Neurod1_Rep2)/2 - (LGS13_peaks_signals$IgG_Rep1 + LGS13_peaks_signals$IgG_Rep2)/2

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert/")
saveRDS(LGS13_peaks_signals,file="LGS13_peaks_signals_Nov2024")

LGS13_Thrb_peaks_res <- Function_2_add_scores(LGS13_Thrb_peaks,LGS13_peaks_signals)
LGS13_Mef2c_peaks_res <- Function_2_add_scores(LGS13_Mef2c_peaks,LGS13_peaks_signals)
LGS13_Onecut1_peaks_res <- Function_2_add_scores(LGS13_Onecut1_peaks,LGS13_peaks_signals)
LGS13_Onecut2_peaks_res <- Function_2_add_scores(LGS13_Onecut2_peaks,LGS13_peaks_signals)
LGS13_Zic3_peaks_res <- Function_2_add_scores(LGS13_Zic3_peaks,LGS13_peaks_signals)
LGS13_Sall3_peaks_res <- Function_2_add_scores(LGS13_Sall3_peaks,LGS13_peaks_signals)
LGS13_Rxrg_peaks_res <- Function_2_add_scores(LGS13_Rxrg_peaks,LGS13_peaks_signals)
LGS13_Pou2f1_peaks_res <- Function_2_add_scores(LGS13_Pou2f1_peaks,LGS13_peaks_signals)


#####
#####-----Next we will process the LGS13 signals-------#########
#####
##### 

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
Alignment <- readRDS(file="Mouse_LGS13_all_peaks_alignment_Nov25.rds")

Thrb_peaks_res_Compare = Function_3_add_Alignments(Mouse_Thrb_peaks_res,LGS13_Thrb_peaks_res,Alignment)
Onecut1_peaks_res_Compare = Function_3_add_Alignments(Mouse_Onecut1_peaks_res,LGS13_Onecut1_peaks_res,Alignment)
Onecut2_peaks_res_Compare = Function_3_add_Alignments(Mouse_Onecut2_peaks_res,LGS13_Onecut2_peaks_res,Alignment)
Zic3_peaks_res_Compare = Function_3_add_Alignments(Mouse_Zic3_peaks_res,LGS13_Zic3_peaks_res,Alignment)
Sall3_peaks_res_Compare = Function_3_add_Alignments(Mouse_Sall3_peaks_res,LGS13_Sall3_peaks_res,Alignment)
Rxrg_peaks_res_Compare = Function_3_add_Alignments(Mouse_Rxrg_peaks_res,LGS13_Rxrg_peaks_res,Alignment)
Pou2f1_peaks_res_Compare = Function_3_add_Alignments(Mouse_Pou2f1_peaks_res,LGS13_Pou2f1_peaks_res,Alignment)

LGS13_Mef2c_peaks_res
Mouse_Mef2c_peaks_res

###### This time we will output an excel file contains all these regions ########
######
Total_RES_List1 = list(LGS_Thrb = Thrb_peaks_res_Compare$L,Mouse_Thrb=Thrb_peaks_res_Compare$M)
Total_RES_List2 = list(LGS_Onecut1 = Onecut1_peaks_res_Compare$L,Mouse_Onecut1=Onecut1_peaks_res_Compare$M)
Total_RES_List3 = list(LGS_Onecut2 = Onecut2_peaks_res_Compare$L,Mouse_Onecut2=Onecut2_peaks_res_Compare$M)
Total_RES_List4 = list(LGS_Zic3 = Zic3_peaks_res_Compare$L,Mouse_Zic3=Zic3_peaks_res_Compare$M)
Total_RES_List5 = list(LGS_Sall3 = Sall3_peaks_res_Compare$L,Mouse_Sall3=Sall3_peaks_res_Compare$M)
Total_RES_List6 = list(LGS_Rxrg = Rxrg_peaks_res_Compare$L,Mouse_Rxrg=Rxrg_peaks_res_Compare$M)
Total_RES_List7 = list(LGS_Pou2f1 = Pou2f1_peaks_res_Compare$L,Mouse_Pou2f1=Pou2f1_peaks_res_Compare$M)
Total_RES_List8 = list(LGS_Mef2c = LGS13_Mef2c_peaks,Mouse_Mef2c=Mouse_Mef2c_peaks$M)

Total_RES_List = c(Total_RES_List1,Total_RES_List2,Total_RES_List3,Total_RES_List4,Total_RES_List5,Total_RES_List6,Total_RES_List7,Total_RES_List8)

library(openxlsx)
# Create a new Excel workbook
wb <- createWorkbook()
# Loop through the list and add each table to a new sheet
for (sheet_name in names(Total_RES_List)) {
  addWorksheet(wb, sheet_name)  # Add a new sheet with the name
  writeData(wb, sheet = sheet_name, Total_RES_List[[sheet_name]])  # Write data to the sheet
}

# Save the workbook to an Excel file
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Convert")
saveWorkbook(wb, file = "LGS13_and_Mouse_8genes_peaks_26Nov2024.xlsx", overwrite = TRUE)





######
######

Function_3_add_Alignments <- function(Mouse_Res,LGS_Res,Alignment){
    ##########
    ### Next we will used for mouse datasets #####
    ##########
    ##########-----for-----Mouse-------#######
    ##########
    Mouse_Res$LGS_peak_conserved = "Notfind"
    Mouse_Res$LGS_peak_overlap = 0
    Mouse_Res$LGS_peak_overlapScore = 0
    for(i in 1:length(Mouse_Res$Peak)){
        #####
        tmp_peak = Mouse_Res$Peak[i]
        k1 = which(Alignment$M_peaks == tmp_peak)
        if(length(k1) == 0){
            next
        }
        LGS_Alignment_k1 = Alignment$LGS_peaks[k1]
        ######
        k2 = which(LGS_Alignment_k1 %in% LGS_Res$Peak == T)
        if(length(k2) == 0){
            next
        }
        #######
        k3 = which(LGS_Res$Peak %in% LGS_Alignment_k1 == T)
        if(length(k3) == 0){
            next
        }
        Mouse_Res$LGS_peak_conserved[i] = paste(LGS_Res$Peak[k3], collapse = "__")
        ####### Next find the score #####
        k4 = which(Alignment$M_peaks == tmp_peak & Alignment$LGS_peaks %in% LGS_Res$Peak[k3] == T)
        RES = Alignment[k4,]
        RES = RES[,c("M_peaks","LGS_peaks","Overlap_region","Overlap_score")]
        RES_cl = RES[!duplicated(paste(RES$M_peaks,RES$LGS_peaks)),]
        #######
        if(dim(RES_cl)[1] == 1){
            Mouse_Res$LGS_peak_overlap[i] = RES_cl$Overlap_region
            Mouse_Res$LGS_peak_overlapScore[i] = RES_cl$Overlap_score
        }
        if(dim(RES_cl)[1] > 1){
            print("Mutiple!")
            rownames(RES_cl) = RES_cl$LGS_peaks
            RES_cl = RES_cl[LGS_Res$Peak[k3],]
            print(i)
            Mouse_Res$LGS_peak_overlap[i] = paste(RES_cl$Overlap_region,collapse = "__")
            Mouse_Res$LGS_peak_overlapScore[i] = paste(RES_cl$Overlap_score,collapse = "__")
        }
    }
    #########
    ####--------Next for 13LGS ######
    #########
    LGS_Res$M_peak_conserved = "Notfind"
    LGS_Res$M_peak_overlap = 0
    LGS_Res$M_peak_overlapScore = 0
    for(i in 1:length(LGS_Res$Peak)){
        #####
        tmp_peak = LGS_Res$Peak[i]
        k1 = which(Alignment$LGS_peaks == tmp_peak)
        if(length(k1) == 0){
            next
        }
        M_Alignment_k1 = Alignment$M_peaks[k1]
        ######
        k2 = which(M_Alignment_k1 %in% Mouse_Res$Peak == T)
        if(length(k2) == 0){
            next
        }
        #######
        k3 = which(Mouse_Res$Peak %in% M_Alignment_k1 == T)
        if(length(k3) == 0){
            next
        }
        LGS_Res$M_peak_conserved[i] = paste(Mouse_Res$Peak[k3], collapse = "__")
        ####### Next find the score #####
        k4 = which(Alignment$LGS_peaks == tmp_peak & Alignment$M_peaks %in% Mouse_Res$Peak[k3] == T)
        RES = Alignment[k4,]
        RES = RES[,c("M_peaks","LGS_peaks","Overlap_region","Overlap_score")]
        RES_cl = RES[!duplicated(paste(RES$M_peaks,RES$LGS_peaks)),]
        #######
        if(dim(RES_cl)[1] == 1){
            LGS_Res$M_peak_overlap[i] = RES_cl$Overlap_region
            LGS_Res$M_peak_overlapScore[i] = RES_cl$Overlap_score
        }
        if(dim(RES_cl)[1] > 1){
            print("Mutiple!")
            rownames(RES_cl) = RES_cl$M_peaks
            RES_cl = RES_cl[Mouse_Res$Peak[k3],]
            print(i)
            LGS_Res$M_peak_overlap[i] = paste(RES_cl$Overlap_region,collapse = "__")
            LGS_Res$M_peak_overlapScore[i] = paste(RES_cl$Overlap_score,collapse = "__")
        }
    }
    #######
    return(list(M=Mouse_Res,L=LGS_Res))
    #######
}














#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########----------------------------------------------------------------------------------------------------------------
#########------------###########-----------#########---------------------------------------------------------------------------------------------------------------

#########这里我们重新see prepare一下Figure2 ####
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&
###########------------- 首先需要 call peaks ##########---------------------------------------------------------------
conda activate ArchR2 
R
library(Seurat)

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')

head(LGS13_Mouse.combined@meta.data)
table(LGS13_Mouse.combined$sp)
  LGS Mouse 
24908 57892
table(LGS13_Mouse.combined$tag)
LGS@P1        LGS@P12        LGS@P17         LGS@P4         LGS@P8 
          5780           7520           4846           3501           3261 
Mouse@E18_rep2 Mouse@E18_rep3       Mouse@P0  Mouse@P2_rep2  Mouse@P2_rep3 
          9638           9732           8530           6643           8244 
      Mouse@P5  Mouse@P8_rep1  Mouse@P8_rep2 
          4359           7870           2876

#########
#########
dim(LGS13_Mouse.combined[["RNA"]]@counts)
table(LGS13_Mouse.combined$celltype)
table(LGS13_Mouse.combined$celltype2)

######### This time we will check the UMAPs for these ############
#########

scRNAseq_Output_meta_data <- function(Seurat_obj){
    Meta_tab = data.frame(Seurat_obj@meta.data)
    UMAP = Embeddings(Seurat_obj[["umap"]])
    Meta_tab = cbind(Meta_tab,UMAP)
    Meta_tab$cell_id = rownames(Meta_tab)
    return(Meta_tab)
}


All_meta_data = scRNAseq_Output_meta_data(LGS13_Mouse.combined)

##########
########## 我们需要输出一下这个UMAP到电脑上，然后再plotly里面再看一下这个图 #################
##########

colorList = c("#EF9000","#804537","#9EA220","#A74997","#61BFB9","#D11536","pink","#AAA9A9","#026AB1")

Add the colors to the meta data
table(LGS13_Mouse.combined$celltype)

##########
color_List = c(AC = "#EF9000",BC="#804537",BC_Photo_pre="pink",Cone="#9EA220",Rod="#026AB1",RPC="#61BFB9",NG="#A74997",RGC="#AAA9A9",MG="#D11536")
m = match(All_meta_data$celltype,names(color_List))
All_meta_data$color = color_List[m]

k = which(is.na(All_meta_data$color) == F)
All_meta_data = All_meta_data[k,]

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
saveRDS(All_meta_data,file="LGS_Mouse_scRNA_UMAP_all.rds")

########## On my own MAC #####
##########
library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("LGS_Mouse_scRNA_UMAP_all.rds")
####All_meta_data$color = gsub("#","FF",All_meta_data$color)
#class(All_meta_data$color)
#All_meta_data$color <- as.character(All_meta_data$color)
#table(All_meta_data$color)
k = which(All_meta_data$color == "pink")
All_meta_data$color[k] = "#FFC0CB"
table(All_meta_data$color)
#All_meta_data$color <- gsub("#","",All_meta_data$color)
color_List = c(AC = "#EF9000",BC="#804537",BC_Photo_pre="pink",Cone="#9EA220",Rod="#026AB1",RPC="#61BFB9",NG="#A74997",RGC="#AAA9A9",MG="#D11536")
All_meta_data$celltype <- factor(All_meta_data$celltype,levels=names(color_List))

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~celltype,colors=color_List,type = "scatter3d", mode = "markers",opacity = 0.75,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = 1.350966, y = 0.4084711, z = -0.2805013))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig



###### 查找视角 ！！！！！ ########
###### 查找 camera 参数 #########
library(shiny)
library(plotly)

# 示例数据

ui <- fluidPage(
  plotlyOutput("plot3d"),
  verbatimTextOutput("camera_info") # 显示当前 camera 参数
)

server <- function(input, output, session) {
  output$plot3d <- renderPlotly({
    plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
            color = ~celltype, type = "scatter3d", mode = "markers",
            opacity = 1,size=0.01)
  })
  
  # 捕获 camera 视角数据
  output$camera_info <- renderPrint({
    event_data("plotly_relayout") # 捕获 relayout 事件
  })
}

shinyApp(ui, server)

######
######
######
$scene.camera$eye$x
[1] 1.350966
$scene.camera$eye$y
[1] 0.4084711
$scene.camera$eye$z
[1] -0.2805013


#######
####### 可以，这个single cell RNAseq 结束了，下一步查找 scATACseq 的UMAP ###########
####### 这个 UMAP 小鼠和 LGS 是分开来的 ##########
####### 需要分别画一下 ！！！！ ###################
#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
LGS_late_ArchR_project_meta <- read.table("LGS_late_ArchR_project_meta.txt",sep="\t",header=T)
table(LGS_late_ArchR_project_meta$celltype)

######## 我们先画一下 LGS 的 celltype ！！！ ########
########
saveRDS(LGS_late_ArchR_project_meta,file="LGS_late_ArchR_project_meta_2024.rds")


######## ------- local MAC computer ------- ########

library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("LGS_late_ArchR_project_meta_2024.rds")
table(All_meta_data$celltype)

color_List = c(ACHC = "#EF9000",BC="#804537",Photo_BC="pink",Cone="#9EA220",Rod="#026AB1",Late_RPCs="#61BFB9",Late_NG="#A74997",RGC="#AAA9A9",MG="#D11536")
All_meta_data$celltype <- factor(All_meta_data$celltype,levels=names(color_List))

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~celltype,colors=color_List,type = "scatter3d", mode = "markers",opacity = 0.75,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = -0.8579592, y = -0.9666645, z = 1.081717))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig


$scene.camera$eye$x
[1] -0.8579592

$scene.camera$eye$y
[1] -0.9666645

$scene.camera$eye$z
[1] 1.081717


########--------------------------------------------------------------------------------------------------------------------------------------------

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
Mouse_late_ArchR_project_meta <- read.table("Mouse_late_ArchR_project_meta.txt",sep="\t",header=T)
table(Mouse_late_ArchR_project_meta$celltype)

########
######## 我们先画一下 Mouse 的 celltype ！！！ ########
########
saveRDS(Mouse_late_ArchR_project_meta,file="Mouse_late_ArchR_project_meta_2024.rds")


######## Next for the Mouse !!!! ##########
########
########


library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("Mouse_late_ArchR_project_meta_2024.rds")
table(All_meta_data$celltype)

color_List = c(ACHC = "#EF9000",BC="#804537",Photo_BC="pink",Cone="#9EA220",Rod="#026AB1",Late_RPCs="#61BFB9",Late_NG="#A74997",RGC="#AAA9A9",MG="#D11536")
All_meta_data$celltype <- factor(All_meta_data$celltype,levels=names(color_List))

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~celltype,colors=color_List,type = "scatter3d", mode = "markers",opacity = 0.75,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = 0.4808195, y = 0.976542, z = -1.112431))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig

$scene.camera$eye$x
[1] 0.3808195

$scene.camera$eye$y
[1] 0.876542

$scene.camera$eye$z
[1] -1.012431


############
############ ------------------------------- 接下来我们 plot Zic3 的基因表达情况 和 Motif binding 的情况 -------------- ###########
############
############ 首先是 RNA 表达水平的 情况 ！！！！ ##############
############

############ 需要加一列基因表达，然后比较 分开Mouse和LGS的meta data ###############
############ 不做smooth ！！！ ##########

############ load the RNA meta data !!!! #########
############

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
All_meta_data <- readRDS(file="LGS_Mouse_scRNA_UMAP_all.rds")
table(All_meta_data$sp)

All_meta_data_Mouse = All_meta_data[which(All_meta_data$sp == "Mouse"),]
All_meta_data_LGS = All_meta_data[which(All_meta_data$sp == "LGS"),]

############ load the expression matrix !!!! #######
############

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
load('LGS13_Mouse.combined')

############
DefaultAssay(LGS13_Mouse.combined) <- "RNA"
RNA_mat = LGS13_Mouse.combined[['RNA']]@data

############ add Genes #####
G_index = which(rownames(RNA_mat) == "Zic3")
C_index = match(rownames(All_meta_data_LGS),colnames(RNA_mat))
All_meta_data_LGS$Exp = RNA_mat[G_index,C_index]

########
G_index = which(rownames(RNA_mat) == "Zic3")
C_index = match(rownames(All_meta_data_Mouse),colnames(RNA_mat))
All_meta_data_Mouse$Exp = RNA_mat[G_index,C_index]

########
######## saveRDS #########
########

saveRDS(All_meta_data_LGS,file="All_meta_data_LGS_Zic3_Exp.rds")
saveRDS(All_meta_data_Mouse,file="All_meta_data_Mouse_Zic3_Exp.rds")

########
######## on the local MAC to plot the results !!!! #########
########

library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("All_meta_data_LGS_Zic3_Exp.rds")

#k = which(All_meta_data$color == "pink")
#All_meta_data$color[k] = "#FFC0CB"
#table(All_meta_data$color)
#All_meta_data$color <- gsub("#","",All_meta_data$color)
#color_List = c(AC = "#EF9000",BC="#804537",BC_Photo_pre="pink",Cone="#9EA220",Rod="#026AB1",RPC="#61BFB9",NG="#A74997",RGC="#AAA9A9",MG="#D11536")
#All_meta_data$celltype <- factor(All_meta_data$celltype,levels=names(color_List))
k = which(All_meta_data$Exp > 1.5)
All_meta_data$Exp[k] = 1.5

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~Exp,colors = c("#E8E8E8","red","purple"),type = "scatter3d", mode = "markers",opacity = 1,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = 1.350966, y = 0.4084711, z = -0.2805013))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig



#####
##### for Mouse ###
#####

library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("All_meta_data_Mouse_Zic3_Exp.rds")

###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------
###### ------------- Next we will add the Motifs scores to the Plot !!!!!! ----------------------------------------------------------------------------------------------------------------------------------


###### load the Motifs !!!! ####
######

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
FN = paste('LGS13_late','chromVAR_dev',sep='_')
load(FN)
LGS13_Z_matrix_late = chromVAR_dev@assays@data$z


setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/Human_scATACseq_new/scATAC/TRANSFAC_2018')
FN = paste('Mouse_late','chromVAR_dev',sep='_')
load(FN)
Mouse_Z_matrix_late = chromVAR_dev@assays@data$z

#######
####### load the meta datasets !!!! ##############
#######
####### Zic3 Motif: M10122
#######
#######

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load("LGS_late_ArchR_project")
head(LGS_late_ArchR_project@cellColData)
cell_id = rownames(LGS_late_ArchR_project@cellColData)

All_meta_data <- readRDS("LGS_late_ArchR_project_meta_2024.rds")
rownames(All_meta_data) <- cell_id

G_index = which(rownames(LGS13_Z_matrix_late) == "M10122")
C_index = match(rownames(All_meta_data),colnames(LGS13_Z_matrix_late))
All_meta_data$chromVAR = LGS13_Z_matrix_late[G_index,C_index]
saveRDS(All_meta_data,file="LGS_Zic3_Motif.rds")


#####

setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
All_meta_data <- readRDS("Mouse_late_ArchR_project_meta_2024.rds")
setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load("Mouse_late_ArchR_project")
head(Mouse_late_ArchR_project@cellColData)
cell_id = rownames(Mouse_late_ArchR_project@cellColData)
rownames(All_meta_data) <- cell_id


G_index = which(rownames(Mouse_Z_matrix_late) == "M10122")
C_index = match(rownames(All_meta_data),colnames(Mouse_Z_matrix_late))
All_meta_data$chromVAR = Mouse_Z_matrix_late[G_index,C_index]
saveRDS(All_meta_data,file="Mouse_Zic3_Motif.rds")


###########
###########


library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("LGS_Zic3_Motif.rds")

k1 = which(All_meta_data$chromVAR > 2)
All_meta_data$chromVAR[k1] = 2
k2 = which(All_meta_data$chromVAR < -2)
All_meta_data$chromVAR[k2] = -2

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~chromVAR,colors = c("blue","grey","red"),type = "scatter3d", mode = "markers",opacity = 1,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = -0.8579592, y = -0.9666645, z = 1.081717))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig



library(plotly)
setwd("/Users/pin/Desktop/13LGS_plotly_figures/")
All_meta_data <- readRDS("Mouse_Zic3_Motif.rds")

k1 = which(All_meta_data$chromVAR > 2)
All_meta_data$chromVAR[k1] = 2
k2 = which(All_meta_data$chromVAR < -2)
All_meta_data$chromVAR[k2] = -2

fig <- plot_ly(All_meta_data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color=~chromVAR,colors = c("blue","grey","red"),type = "scatter3d", mode = "markers",opacity = 1,size=0.01)
fig <- fig %>% layout(scene = list(xaxis = list(title = '',showticklabels = FALSE),
                                   yaxis = list(title = '',showticklabels = FALSE),
                                   zaxis = list(title = '',showticklabels = FALSE),
                                   cliponaxis = TRUE,
                                   camera = list(eye = list(x = 0.4808195, y = 0.976542, z = -1.112431))
                                  ))

fig <- fig %>% layout(plot_bgcolor='rgb(254, 247, 234)')
fig


########### Next we will Plot the Motifs by R ###########
########### for the Motif: 
########### M10122 #####.... 如何在 R 里面 Plot motif logo ############
###########

install.packages('ggseqlogo')
library('ggseqlogo')
load('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13/Raw_data/13LGS/Cone_res/PWM_list_combine_cl')


PWM_list = PWM_list_combine_cl
Name = "M10122"
Get_Motif_matrix <- function(PWM_list,Name){
    Names_list = names(PWM_list)
    K = which(Names_list == Name)
    if(length(K) == 0){print('No data !!!')}
        if(length(K) > 0){
            print('Find !!!')
            temp_data = PWM_list[[K]]
            temp_matrix = temp_data@profileMatrix
            Out_list = temp_matrix
    }
    return(Out_list)
}


Convert_to_PWM <- function(list){
    for(i in 1:length(list)){
        print(i)
        new_mat = list[[i]]
        new_mat_log = exp(new_mat)
        new_mat_log_add = new_mat_log*0.25
        list[[i]] = new_mat_log_add
    }
    return(list)
}

Zic3_Motif = Get_Motif_matrix(PWM_list,Name='M10122')
Zic3_Motif_mat = Convert_to_PWM(list(Zic3_Motif))

png('Zic3_Motif.png',height=1500,width=4000,res=72*12)
ggseqlogo(Zic3_Motif_mat,method = 'bits')
dev.off()


#################
#################
#################
################# OK!!!! we finished the Plot tasks!!!! ########
#################
#################
#################
################# Next we will call DEGs !!! ######
#################
#################
#################
################# This time we will try the DEGs analysis for this project !!! ####
#################

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

################## load the combined seurat files ########
##################


setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load('LGS13_Mouse.combined')

head(LGS13_Mouse.combined@meta.data)
table(LGS13_Mouse.combined$sp)
LGS Mouse 
24908 57892
table(LGS13_Mouse.combined$tag)

DefaultAssay(LGS13_Mouse.combined) = "RNA"
#######
####### We will call DEGs ######
#######

library(Seurat)
DefaultAssay(LGS13_Mouse.combined) <- "RNA"

#######
table(LGS13_Mouse.combined$celltype)

#######
#######
#######
LGS13_Mouse.combined$index <- paste0(LGS13_Mouse.combined$sp,":",LGS13_Mouse.combined$celltype)
table(LGS13_Mouse.combined$index)
Idents(LGS13_Mouse.combined) <- "index"

#######
#######
#######
Markers_All[[1]]$class = "LGS:RPC vs Mouse:RPC"
Markers_All[[2]]$class = "LGS:NG vs Mouse:NG"
Markers_All[[3]]$class = "LGS:BC_Photo_pre vs Mouse:BC_Photo_pre"
Markers_All[[4]]$class = "LGS:Cone vs Mouse:Cone"
Markers_All[[5]]$class = "LGS:Cone vs Mouse:Rod"



Markers1 = FindMarkers(LGS13_Mouse.combined,ident.1='LGS:RPC',ident.2='Mouse:RPC',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers1[which(rownames(Markers1) == "Otx2"),]
Markers1[which(rownames(Markers1) == "Neurod1"),]
Markers1[which(rownames(Markers1) == "Pou2f1"),]


Markers2 = FindMarkers(LGS13_Mouse.combined,ident.1='LGS:NG',ident.2='Mouse:NG',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers2[which(rownames(Markers2) == "Otx2"),]
Markers2[which(rownames(Markers2) == "Neurod1"),]
Markers2[which(rownames(Markers2) == "Pou2f1"),]


Markers3 = FindMarkers(LGS13_Mouse.combined,ident.1='LGS:BC_Photo_pre',ident.2='Mouse:BC_Photo_pre',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers3[which(rownames(Markers3) == "Onecut1"),]
Markers3[which(rownames(Markers3) == "Onecut2"),]
Markers3[which(rownames(Markers3) == "Zic3"),]

Markers4 = FindMarkers(LGS13_Mouse.combined,ident.1='LGS:Cone',ident.2='Mouse:Cone',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
Markers4[which(rownames(Markers4) == "Onecut1"),]
Markers4[which(rownames(Markers4) == "Onecut2"),]
Markers4[which(rownames(Markers4) == "Zic3"),]
Markers4[which(rownames(Markers4) == "Thrb"),]
Markers4[which(rownames(Markers4) == "Rxrg"),]
Markers4[which(rownames(Markers4) == "Sall3"),]
Markers4[which(rownames(Markers4) == "Mef2c"),]

Markers5 = FindMarkers(LGS13_Mouse.combined,ident.1='LGS:Cone',ident.2='Mouse:Rod',min.pct = 0.05,logfc.threshold = 0.25,test.use = "MAST")
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
saveRDS(Markers_All,file="Markers_All_2024_Dec.rds")

###########
########### Next we will Get average expression files for each cell types #######
###########

Celltypes = names(table(LGS13_Mouse.combined$index))
Celltypes_Need = Celltypes[c(10,7,3,4,20,17,13,14,19)]

############
RNA_Counts_Mat = LGS13_Mouse.combined[['RNA']]@counts
RNA_Counts_Mat_Meta = LGS13_Mouse.combined@meta.data

Mat = RNA_Counts_Mat
Meta = RNA_Counts_Mat_Meta

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



all_avg_mat_log[which(rownames(all_avg_mat_log) == "Thrb"),]
all_avg_mat_log[which(rownames(all_avg_mat_log) == "Zic3"),]
all_avg_mat_log[which(rownames(all_avg_mat_log) == "Rxrg"),]


######## OK!!! #####
Average_Exp_Mat = Get_Avg_Exp(Mat,Meta,Celltypes_Need)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
saveRDS(Average_Exp_Mat,file="Average_Exp_Mat_2024_Dec.rds")

Average_Exp_Mat <- readRDS("Average_Exp_Mat_2024_Dec.rds")

######
###### filter the Exp Mat by DEGs ######
######

Markers_All <- readRDS("Markers_All_2024_Dec.rds")
Markers_All_Merge = do.call("rbind",Markers_All)

Markers_All_Merge[which(Markers_All_Merge$Gene == "Onecut2"),]
Markers_All_Merge_Cl = Markers_All_Merge[which(Markers_All_Merge$p_val_adj < 1e-6 & abs(Markers_All_Merge$avg_log2FC) > 0.35),]
Average_Exp_Mat_Cl = Average_Exp_Mat[which(rownames(Average_Exp_Mat) %in% Markers_All_Merge_Cl$Gene == T),]


Markers_All_Merge_Cl = Markers_All_Merge[which(Markers_All_Merge$p_val_adj < 1e-6 & abs(Markers_All_Merge$avg_log2FC) > 0.35),]
saveRDS(kc_dat2_tab,file="kc_dat2_tab_2024.rds")
kc_dat2_tab <- readRDS("kc_dat2_tab_2024.rds")

####
Markers_All_Merge_Cl$cluster = kc_dat2_tab$cluster[match(Markers_All_Merge_Cl$Gene,kc_dat2_tab$genes)]
Markers_All_Merge_Cl = Markers_All_Merge_Cl[,c("Gene","class","cluster","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")]
####
library(writexl)
write_xlsx(Markers_All_Merge_Cl, path = "Figure2_heatmap.xlsx")




Average_Exp_Mat_Cl_Scale = t(apply(Average_Exp_Mat_Cl,1,scale))
colnames(Average_Exp_Mat_Cl_Scale) = colnames(Average_Exp_Mat_Cl)
rownames(Average_Exp_Mat_Cl_Scale) = rownames(Average_Exp_Mat_Cl) 


######

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

labels = c("Zic3","Otx2","Neurod1","Thrb","Rxrg","Onecut1","Onecut2","Sall3","Mef2c","Pou2f1")
at = match(labels,rownames(Average_Exp_Mat_Cl_Scale))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('RNA_test1.png',height=4000,width=3000,res=72*12)
Heatmap(Average_Exp_Mat_Cl_Scale, name = "XX", border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

#####
##### we will add the Genes annotations ######
#####

order = c(2,3,1,7,4,6,5,8)

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

labels = c("Zic3","Otx2","Neurod1","Thrb","Rxrg","Onecut1","Onecut2","Sall3","Mef2c","Pou2f1")
at = match(labels,rownames(Average_Exp_Mat_Cl_Scale))

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
png('RNA_test2.png',height=4000,width=3000,res=72*12)
Heatmap(Average_Exp_Mat_Cl_Scale, name = "XX", border = T,use_raster=FALSE,show_row_names=F,show_column_names=T,rect_gp = gpar(col = 'white', lwd = 0),cluster_rows = F,cluster_columns = F,col = col_fun,heatmap_legend_param = list(col_fun = col_fun,title = "",border ='black',at=c(-2,-1,0,1,2)),row_split=row_sp,column_split=col_sp) +
rowAnnotation(link = anno_mark(at = at,labels = labels),gp = gpar(fontsize = 20))
dev.off()

kc_dat2_tab = kc_dat2[order(kc_dat2$cluster,decreasing=F),]
saveRDS(kc_dat2_tab,file="kc_dat2_tab_2024.rds")

m = match(kc_dat2_tab$genes,rownames(Average_Exp_Mat_Cl))
Average_Exp_Mat_Cl_re = Average_Exp_Mat_Cl[m,]
all.equal(rownames(Average_Exp_Mat_Cl_re),kc_dat2_tab$genes)
kc_dat2_tab = cbind(kc_dat2_tab,Average_Exp_Mat_Cl_re)
saveRDS(kc_dat2_tab,file="kc_dat2_tab_2024.rds")


install.packages("writexl")
library(writexl)
write_xlsx(kc_dat2_tab, "Figure2c_MouseVS13LGS_DEGs.xlsx")

table(kc_dat2_tab$cluster)

##############
##############
##############

##############
############## Next we will prepare the PeakMatrix for the Late time point!!!! ###########
##############

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate ArchR
R

#######

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles")
PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/Finally_retinal_dev_202009/Figure2/Combine_202009/ArrowFiles/Mouse_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")
Mouse_PeakAvg = readRDS("Mouse_Avg_PeakMat_Nov2024")
Gtf <- readRDS("Mouse_Gene_Anno_Oct23")

############### see the Mouse PeakAvg ###
colSums(PeakAvg)
###############
dim(PeakAvg)


############### see the 13LGS Peak mat !!! ################
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
Gtf <- readRDS("LGS13_Gene_Anno_Oct23")
PtoG = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Project_new_server_cl_PtoG_Sub_RES_Nov18_2024")
LGS13_PeakAvg = readRDS("LGS13_Avg_PeakMat_Nov2024")

################
################
################

colnames(LGS13_PeakAvg)




















########## 比较一下 paired peaks number between the Mouse and 13LGS ########
##########

library(ggplot2)

# 创建示例数据
df <- data.frame(
  ID = rep(1:5, each = 2),         # 配对的唯一标识
  Category = rep(c("A", "B"), 5),  # 类别
  Value = c(2, 3, 4, 6, 5, 7, 8, 9, 10, 11) # 值
)

# 查看数据：
print(df)



########## 
ggplot(df, aes(x = Category, y = Value, group = ID)) +  # 按配对分组
  geom_point(aes(color = Category), size = 3) +        # 绘制点
  geom_line(aes(group = ID), color = "black", linetype = "dashed") + # 添加连接线
  theme_minimal() +
  labs(x = "Category", y = "Value", title = "Paired Points with Connecting Lines")


##########