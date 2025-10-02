#######
####### load gene-gene correlation results ########### 
#######

library(data.table)

####### Mouse_total_network_Jan2025.tsv is the arboreto output ####
####### LGS13_total_network_Jan2025.tsv is the arboreto output ####

Mouse_total_GRNs = read.table("Mouse_total_network_Jan2025.tsv",sep="\t")
LGS13_total_GRNs = read.table("LGS13_total_network_Jan2025.tsv",sep="\t")

LGS13_total_GRNs[which(LGS13_total_GRNs$V1 == "Zic3" & LGS13_total_GRNs$V2 == "Thrb"),]

#### convert .like to -like #####
#### convert .containing to -containing #########
#### LGS13_total_GRNs$V1 = gsub(".like","-like",LGS13_total_GRNs$V1)
#### LGS13_total_GRNs$V1 = gsub(".containing","-containing",LGS13_total_GRNs$V1)
#### LGS13_total_GRNs$V2 = gsub(".like","-like",LGS13_total_GRNs$V2)
#### LGS13_total_GRNs$V2 = gsub(".containing","-containing",LGS13_total_GRNs$V2)

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



#######
####### example code for cell type specific footprint signals ########
####### Add_score_to_the_motifs: calculate single-cell bw signal for each motif region
#######
########
######## peaks peak GRanges #######
######## PWM_list ######
########


Match_Motif_to_Peaks <- function(peaks,PMW_list,Genome="mm10"){
    ##########
    library(motifmatchr)
    library(parallel)
    ##########
    if(Genome=="mm10"){
        library('BSgenome.Mmusculus.UCSC.mm10.masked')
        res_list <- mclapply(PMW_list, function(x) {
            matchMotifs(
                pwms = x,
                subject = peaks,
                genome = BSgenome.Mmusculus.UCSC.mm10.masked,
                out = "positions",
                p.cutoff = 5e-05
            )
        }, mc.cores = 30)
    }
    ##########
    if(Genome=="13LGS"){
        library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")
        res_list <- mclapply(PMW_list, function(x) {
            matchMotifs(
                pwms = x,
                subject = peaks,
                genome = BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise,
                out = "positions",
                p.cutoff = 5e-05
            )
        }, mc.cores = 30)
    }
    ##########
    return(res_list)
}





########
######## Match_Motif_to_Peaks outputs #########
########

Check_normalized_Signal <- function(x,temp_bw){
    library(GenomicRanges)
	library(rtracklayer)
	temp_bw = temp_bw
	### print(summary(width(temp_bw)))
    ### for each Hints, see the footprint score ####
    ### 将 ATAC GRanges 转为覆盖度 (RleList) #######
    ### not work ####
    x = x[[1]]
    x = x[which(x$score > 2)]
    ###
    flank_size = width(x)*3
    c_start <- start(x)
    c_end   <- end(x)
    lf_start <- c_start - flank_size
    lf_end   <- c_start - 1
    rf_start <- c_end + 1
    rf_end   <- c_end + flank_size
    ####
    GR_mid = GRanges(seqnames=seqnames(x),IRanges(c_start,c_end))
    GR_left = GRanges(seqnames=seqnames(x),IRanges(lf_start,lf_end))
    GR_right = GRanges(seqnames=seqnames(x),IRanges(rf_start,rf_end))
    ####
    Get_score <- function(x,y){
        library(GenomicRanges)
        res = findOverlaps(x,y)
        res = data.frame(res)
        res$score = y$score[res$subjectHits]
        #### add sums ####
        res_sum = tapply(res$score,res$queryHits,sum)
        ####
        x$score = 0
        x$score[as.numeric(names(res_sum))] = as.numeric(res_sum)
        return(x)
    }
    #####
    GR_mid = Get_score(GR_mid,temp_bw)
    GR_left = Get_score(GR_left,temp_bw)
    GR_right = Get_score(GR_right,temp_bw)
    #####
    GR_mid$score = GR_mid$score / c(width(x))
    GR_left$score = GR_left$score / c(width(x)*3)
    GR_right$score = GR_right$score / c(width(x)*3)
    #####
    GR_mid$left = GR_left$score
    GR_mid$right = GR_right$score
    #####
    return(GR_mid)
}


########
######## file is the path of the signal file #########
######## matchtopeaks_list is the output of the Check_normalized_Signal ##########
########


Add_score_to_the_motifs <- function(file,matchtopeaks_list){
    ##########
    library(parallel)
    ##########
    library(rtracklayer)
    temp_bw = import.bw(file)
    res_list <- mclapply(matchtopeaks_list,Check_normalized_Signal,temp_bw=temp_bw,mc.cores = 30)
    ########## length(res_list) #######
    names(res_list) <- names(matchtopeaks_list)
    ####
    res_list2 = GRangesList(res_list)
    res_list3 = unlist(res_list2)
    return(res_list3)
}

##########
##########

source("/zp1/data/plyu3/All_Functions_2025/R/source_all.R")
source_all_r("/zp1/data/plyu3/All_Functions_2025/R/")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS13_Motif_matchtopeaks <- readRDS("LGS13_Motif_matchtopeaks_Jan13.rds")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Cone_res/LGS_Cone_fragments_cl_bamGR_pe_s_corrected.bw"
LGS_Cone_footprint_res = Add_score_to_the_motifs(file,LGS13_Motif_matchtopeaks)
saveRDS(LGS_Cone_footprint_res,file="LGS_Cone_footprint_res_2024")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_NG_res/LGS_Late_NG_fragments_cl_bamGR_pe_s_corrected.bw"
LGS_Late_NG_footprint_res = Add_score_to_the_motifs(file,LGS13_Motif_matchtopeaks)
saveRDS(LGS_Late_NG_footprint_res,file="LGS_Late_NG_footprint_res_2024")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_RPCs_res/LGS_Late_RPCs_fragments_cl_bamGR_pe_s_corrected.bw"
LGS_Late_RPCs_footprint_res = Add_score_to_the_motifs(file,LGS13_Motif_matchtopeaks)
saveRDS(LGS_Late_RPCs_footprint_res,file="LGS_Late_RPCs_footprint_res_2024")

file = "/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Photo_BC_res/LGS_Photo_BC_fragments_cl_bamGR_pe_s_corrected.bw"
LGS_Photo_BC_footprint_res = Add_score_to_the_motifs(file,LGS13_Motif_matchtopeaks)
saveRDS(LGS_Photo_BC_footprint_res,file="LGS_Photo_BC_footprint_res_2024")


#######
####### integrate footprint PeaktoGene and GG correlations ######
#######

Add_foot_print_to_peak <- function(footprint,TAG,Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat){
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
    LGS_motif_tf_table_cl = LGS_motif_tf_table[,c("motifs","tfs_LGS")]
    colnames(LGS_motif_tf_table_cl) = c("Motif","TF")
    LGS_motif_tf_table_cl = data.table(LGS_motif_tf_table_cl)
    ##########
    merge_table2 <- merge(
        x = merge_table,
        y = LGS_motif_tf_table_cl,
        by = c("Motif"),
        allow.cartesian = TRUE
    )
    ########## Next we will merge the gene-gene correaltion and importance ####
    total_GRNs_add = data.table(total_GRNs_add)
    ##########
    colnames(total_GRNs_add) = c("TF","gene","GG_importance","GG_Corr")
    ##########
    merge_table3 <- merge(
        x = merge_table2,
        y = total_GRNs_add,
        by = c("TF","gene")
    )
    ###########
    ########### Next we will add the average expression of the TFs ################
    ###########
    Avg_mat_sub = data.frame(TF=rownames(Avg_Mat),Exp=Avg_Mat[,TAG])
    ###
    m = match(merge_table3$TF,Avg_mat_sub$TF)
    merge_table3$TF_Exp = Avg_mat_sub$Exp[m]
    ####
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

########



setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
Peak_to_Gene_add = readRDS("Mouse_Peak_to_Gene_merge_CL_2025")
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_motif_tf_table = readRDS("LGS_motif_tf_table_2025")
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
total_GRNs_add = readRDS("Mouse_total_GRNs_add_2025Mar")
Avg_Mat = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Avg_GeneMat_Late_2025")

####
#### Avg_Mat[which(rownames(Avg_Mat) == "Mef2c"),] ####
####
footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/Mouse_Cone_footprint_res_2024")
TAG = "Cone"
Mouse_Cone_GRNs = Add_foot_print_to_peak(footprint,"Cone",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)

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

saveRDS(Mouse_Cone_GRNs,"Mouse_Cone_GRNs_2025")
saveRDS(Mouse_Rod_GRNs,"Mouse_Rod_GRNs_2025")
saveRDS(Mouse_Photopre_GRNs,"Mouse_Photopre_GRNs_2025")
saveRDS(Mouse_RPC_GRNs,"Mouse_RPC_GRNs_2025")
saveRDS(Mouse_NG_GRNs,"Mouse_NG_GRNs_2025")

#########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######### Next for the 13LGS ############
#########------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Cone_footprint_res_2024")
TAG = "Cone"
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
Peak_to_Gene_add = readRDS("LGS_13_Peak_to_Gene_merge_CL_2025")

#####
#####
head(Peak_to_Gene_add[grep("-containing",Peak_to_Gene_add$gene),])
head(Peak_to_Gene_add[grep("-like",Peak_to_Gene_add$gene),])


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
LGS_motif_tf_table = readRDS("LGS_motif_tf_table_2025")
head(LGS_motif_tf_table[grep("-containing",LGS_motif_tf_table$LGS13),])
head(LGS_motif_tf_table[grep("-like",LGS_motif_tf_table$LGS13),])


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
total_GRNs_add = readRDS("LGS13_total_GRNs_add_2025")
head(total_GRNs_add[grep("-containing",total_GRNs_add$V1),])
head(total_GRNs_add[grep("-like",total_GRNs_add$V1),])
head(total_GRNs_add[grep("-containing",total_GRNs_add$V2),])
head(total_GRNs_add[grep("-like",total_GRNs_add$V2),])


Avg_Mat = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS13_Avg_GeneMat_Late_2025")
rownames(Avg_Mat)[grep("-containing",rownames(Avg_Mat))]
rownames(Avg_Mat)[grep("-like",rownames(Avg_Mat))]


LGS13_Cone_GRNs = Add_foot_print_to_peak(footprint,"Cone",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Cone_GRNs,"LGS13_Cone_GRNs_2025")

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_NG_footprint_res_2024")
TAG = "NG"
LGS13_NG_GRNs = Add_foot_print_to_peak(footprint,"NG",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_NG_GRNs,"LGS13_NG_GRNs_2025")

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Photo_BC_footprint_res_2024")
TAG = "BC_Photo_pre"
LGS13_Photopre_GRNs = Add_foot_print_to_peak(footprint,"BC_Photo_pre",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_Photopre_GRNs,"LGS13_Photopre_GRNs_2025")

footprint = readRDS("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW/LGS_Late_RPCs_footprint_res_2024")
TAG = "RPC"
LGS13_RPC_GRNs = Add_foot_print_to_peak(footprint,"RPC",Peak_to_Gene_add,LGS_motif_tf_table,total_GRNs_add,Avg_Mat)
setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS_13_NEW")
saveRDS(LGS13_RPC_GRNs,"LGS13_RPC_GRNs_2025")

########
########----------------------------------------------------------------------------------------------------------------------------------
########



