#########
######### First load all the peaks in C2 and C3 ######
#########


########################################################################
########################################################################
########################################################################
########################################################################


squirrel_to_mouse <- read.table('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure2/squirrel_to_mouse_gene_names.csv',sep=',',header=T)
setwd('/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3')
load(file='mouse.combined')
Mouse_Gene_to_13LGS <- data.frame(Mouse_Genes = rownames(mouse.combined))
m = match(rownames(mouse.combined),squirrel_to_mouse$Gene.name)
Mouse_Gene_to_13LGS$LGS_Genes = squirrel_to_mouse$Squirrel.name[m]
Mouse_Gene_to_13LGS = Mouse_Gene_to_13LGS[-which(is.na(Mouse_Gene_to_13LGS$LGS_Genes) == T),]


kc_dat2_tab <- readRDS("kc_dat2_tab_2024.rds")
kc_dat2_tab2 = kc_dat2_tab[,c("genes","cluster")]
colnames(kc_dat2_tab2) = c("LGS_genes","cluster")
m = match(kc_dat2_tab2$LGS_genes,Mouse_Gene_to_13LGS$LGS_Genes)
kc_dat2_tab2$Mouse_genes = Mouse_Gene_to_13LGS$Mouse_Genes[m]


setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
LGS_Peak_HeatMap_RES <- readRDS(file="LGS_Peak_HeatMap_RES_2024Dec.rds")

setwd("/zp1/data/plyu3/Old_Server_Data/plyu3/LGS13_Figures/Main_Figure3")
Mouse_Peak_HeatMap_RES <- readRDS(file="Mouse_Peak_HeatMap_RES_2024Dec.rds")


LGS_Gnames = sapply(LGS_Peak_HeatMap_RES,function(x) x$Gene[1])
Mouse_Gnames = sapply(Mouse_Peak_HeatMap_RES,function(x) x$Gene[1])

#######
names(LGS_Peak_HeatMap_RES) <- LGS_Gnames
names(Mouse_Peak_HeatMap_RES) <- Mouse_Gnames

tab = data.frame(LGS=names(LGS_Peak_HeatMap_RES),Mouse=names(Mouse_Peak_HeatMap_RES) )
k = which(tab[,1] != tab[,2])

all.equal(names(LGS_Peak_HeatMap_RES),names(Mouse_Peak_HeatMap_RES))


Make_the_Table3 <- function(LGS_Peak_HeatMap_RES,Mouse_Peak_HeatMap_RES,kc_dat2_tab2){
    LGS_peaks = c()
    Mouse_peaks = c()
    ###############
    ### add the number of peaks on the kc_dat2_tab2 ########
    ###############
    head(kc_dat2_tab2)
    ###############
    for(i in 1:length(kc_dat2_tab2$LGS_genes)){
        LGS_G = kc_dat2_tab2$LGS_genes[i]
        Mouse_G = kc_dat2_tab2$Mouse_genes[i]
        #########
        k1 = which(names(LGS_Peak_HeatMap_RES) == LGS_G)
        k2 = which(names(Mouse_Peak_HeatMap_RES) == Mouse_G)
        #########
        if(length(k1) > 0){
            ######
            peaks = LGS_Peak_HeatMap_RES[[k1]]$Peak
            peaks = paste(peaks,collapse="__")
            LGS_peaks <- c(LGS_peaks,peaks)
            
        }else{
            LGS_peaks <- c(LGS_peaks,"NA")
        }
        if(length(k2) > 0){
            ######
            peaks = Mouse_Peak_HeatMap_RES[[k2]]$Peak
            peaks = paste(peaks,collapse="__")
            Mouse_peaks <- c(Mouse_peaks,peaks)
            
        }else{
            Mouse_peaks <- c(Mouse_peaks,"NA")
        }
    }
    ######
    ######
    kc_dat2_tab2$LGS_peaks = LGS_peaks
    kc_dat2_tab2$Mouse_peaks = Mouse_peaks
    ######
    return(kc_dat2_tab2)
}


LGS_Mouse_peaks_List <- Make_the_Table3(LGS_Peak_HeatMap_RES,Mouse_Peak_HeatMap_RES,kc_dat2_tab2)
LGS_Mouse_peaks_List[which(LGS_Mouse_peaks_List$LGS_genes == "Thrb"),]


#########
######### next compare these peaks #########
#########
LGS_Mouse_peaks_List_C2C3 = LGS_Mouse_peaks_List[which(LGS_Mouse_peaks_List$cluster %in% c(2,3) == T),]

k1 = which(LGS_Mouse_peaks_List_C2C3$LGS_peaks == "NA")
k2 = which(LGS_Mouse_peaks_List_C2C3$Mouse_peaks == "NA")

LGS_Mouse_peaks_List_Mouse_SP = LGS_Mouse_peaks_List_C2C3[k1,]
LGS_Mouse_peaks_List_LGS_SP = LGS_Mouse_peaks_List_C2C3[k2,]

save(LGS_Mouse_peaks_List_Mouse_SP,file="LGS_Mouse_peaks_List_Mouse_SP_2024")
save(LGS_Mouse_peaks_List_LGS_SP,file="LGS_Mouse_peaks_List_LGS_SP_2024")


#########
#########
k3 = c(k1,k2)
k3 = k3[!duplicated(k3)]
LGS_Mouse_peaks_List_Compare = LGS_Mouse_peaks_List_C2C3[-k3,]

table(LGS_Mouse_peaks_List_Compare$cluster)


Get_pairs_wise_compare_from_the_list <- function(LGS_Mouse_peaks_List){
    #########
    Out_list <- list()
    print(length(LGS_Mouse_peaks_List[,1]))
    #########
    for(i in 1:length(LGS_Mouse_peaks_List[,1])){
        if(i %% 100 == 0){
            print(i)
        }
        ###########
        peak_LGS = LGS_Mouse_peaks_List$LGS_peaks[i]
        peak_Mouse = LGS_Mouse_peaks_List$Mouse_peaks[i]
        ###########
        peak_LGS = strsplit(peak_LGS,split='__')[[1]]
        peak_Mouse = strsplit(peak_Mouse,split='__')[[1]]
        ###########
        forward = merge(data.frame(LGS_peaks = peak_LGS),data.frame(Mouse_peaks = peak_Mouse))
        reverse = merge(data.frame(LGS_peaks = peak_LGS),data.frame(Mouse_peaks = peak_Mouse))
        ########### Next we will get the seq !!! #####
        library(Biostrings)
        library(BSgenome.Mmusculus.UCSC.mm10)
        library("BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise")
        forward_LGS = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(forward$LGS_peaks))
        reverse_LGS = getSeq(BSgenome.Itridecemlineatus.Bustamante.13LGSHiRise, GRanges(reverse$LGS_peaks))
        forward_Mouse = getSeq(BSgenome.Mmusculus.UCSC.mm10, GRanges(forward$Mouse_peaks))
        reverse_Mouse = getSeq(BSgenome.Mmusculus.UCSC.mm10, GRanges(reverse$Mouse_peaks))
        ######
        reverse_Mouse <- reverseComplement(reverse_Mouse)
        ######
        forward$LGS_Seq = as.character(forward_LGS)
        reverse$LGS_Seq = as.character(reverse_LGS)
        forward$Mouse_Seq = as.character(forward_Mouse)
        reverse$Mouse_Seq = as.character(reverse_Mouse)
        #######
        ####### add tags ####
        #######
        forward$class = "forward"
        reverse$class = "reverse"
        #######
        merge = rbind(forward,reverse)
        #######
        merge$LGS_genes = LGS_Mouse_peaks_List$LGS_genes[i]
        merge$Mouse_genes = LGS_Mouse_peaks_List$Mouse_genes[i]
        merge$cluster = LGS_Mouse_peaks_List$cluster[i]
        ########
        Out_list <- c(Out_list,list(merge))
    }
    Out_list_out = do.call("rbind",Out_list)
    return(Out_list_out)
}

LGS_Mouse_peaks_List_compare = Get_pairs_wise_compare_from_the_list(LGS_Mouse_peaks_List_Compare)

saveRDS(LGS_Mouse_peaks_List_compare,file="LGS_Mouse_peaks_List_compare_2024")

MergeTab_All_List <- apply(LGS_Mouse_peaks_List_compare, 1, function(row) as.list(row))

########### ############
########### ############
########### write a funciton to blast the Peak pairs #########
########### ############
########### ############

compare_2_seq <- function(x){
    ######
    library(Biostrings)
    seq1 = DNAString(x$Mouse_Seq)
    seq2 = DNAString(x$LGS_Seq)
    ######
    alignment <- pairwiseAlignment(seq1, seq2, type = "local",gapOpening = -1, gapExtension = -0.5)
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


x = MergeTab_All_List[[10000]]

compare_2_seq2 <- function(x){
    ######
    library(Biostrings)
    seq1 = DNAString(x$Mouse_Seq)
    seq2 = DNAString(x$LGS_Seq)
    ######
    alignment <- pairwiseAlignment(seq1, seq2, type = "local",gapOpening = -2, gapExtension = -1)
    score = round(score(alignment),3)
    #########
    return(score)
}

#####
library(pbapply)
library(parallel)
cl <- makeCluster(40)  # 这里使用2个核心
alignment_res = parLapply(cl, MergeTab_All_List, compare_2_seq2)
stopCluster(cl)
save(alignment_res,file="alignment_res_2024Dec")


strsplit(alignment_res[[1]],split=":")[[1]][1]

overlap_counts = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][1]))
overlap_score = sapply(alignment_res,function(x) return(strsplit(x,split=":")[[1]][2]))

LGS_Mouse_peaks_List_compare$Overlap_region = NA
LGS_Mouse_peaks_List_compare$Overlap_score = as.numeric(alignment_res)

summary(as.numeric(LGS_Mouse_peaks_List_compare$Overlap_score))


#######
####### get the cutoff #####
#######

LGS_Mouse_peaks_List_compare = LGS_Mouse_peaks_List_compare[order(LGS_Mouse_peaks_List_compare$Overlap_score,decreasing =T),]
LGS_Mouse_peaks_List_compare = LGS_Mouse_peaks_List_compare[order(LGS_Mouse_peaks_List_compare$Overlap_score,decreasing =T),]

summary(LGS_Mouse_peaks_List_compare$Overlap_score)

library(MASS)
shapiro.test(LGS_Mouse_peaks_List_compare$Overlap_score)
fit <- fitdistr(LGS_Mouse_peaks_List_compare$Overlap_score, densfun = "normal")
pnorm(75, mean =32.313 , sd=23.878)

saveRDS(LGS_Mouse_peaks_List_compare,file="LGS_Mouse_peaks_List_compare_2024Dec16")
######### 
#########

