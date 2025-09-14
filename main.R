library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)

if(FALSE){
  remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE, force = TRUE)  
  #skip updates!
}
library(Seurat)

options(bitmapType="cairo")  #this makes DimPlot work



###################################################################################################
########################### Load and preprocess data ##############################################
###################################################################################################

### Only run if needed
if(FALSE){
  
  mat <- ReadParseBio("/home/m/mahogny/mystore/rift/P32705_1002/DGE_filtered/")  
  adata <- CreateSeuratObject(mat)
  
  ## Add metadata about infection
  adata$is_inf <- as.integer(str_split_fixed(colnames(adata),"_",2)[,1])>=7
  
  ### donor1: 1,2, 7,8
  ### donor2: 3,4, 9,10
  ### donor3: 5,6, 11,12
  
  
  adata$donor <- NA#as.integer(str_split_fixed(colnames(adata),"_",2)[,1])>=7
  adata$donor[as.integer(str_split_fixed(colnames(adata),"_",2)[,1]) %in% c(1,2,7,8)] <- "d1"
  adata$donor[as.integer(str_split_fixed(colnames(adata),"_",2)[,1]) %in% c(3,4,9,10)] <- "d2"
  adata$donor[as.integer(str_split_fixed(colnames(adata),"_",2)[,1]) %in% c(5,6,11,12)] <- "d3"
  
  
  #### Add virus counts
  get_virus_reads <- function(fvirus){
    virus1 <- str_split_fixed(str_split_fixed(readLines(fvirus),fixed("CB:Z:"),2)[,2],"\t",2)[,1]
    virus1 <- table(virus1)
    virus1 <- data.frame(bc=names(virus1), virus1=as.integer(virus1))
    
    #Get strand of read
    r1_df_raw <- read.delim(fvirus, header=FALSE, row.names=NULL)
    r1_df <- r1_df_raw %>% mutate(strand = case_when(V2 == 0 ~ "forward",
                                                   V2 == 16 ~ "reverse"))
    r1_df$V20 <- gsub('CB:Z:', '', r1_df$V20)
    r1_df <- r1_df[c("V20", "strand")]
    
    #Create columns of counts for forward and reverse strands
    r1_df <- r1_df %>%  group_by(V20) %>%summarise("forward" = sum(strand == "forward"),
                                                    "reverse" = sum(strand == "reverse"))
    names(r1_df)[1] <- "bc"
    
    #Combine strands and counts
    full_join(virus1, r1_df, by = "bc")
  }
  
  virus1 <- get_virus_reads("/home/k/kwyo0001/mystore/rift/P32705_1002/process/r1.txt")
  virus2 <- get_virus_reads("/home/k/kwyo0001/mystore/rift/P32705_1002/process/r2.txt")
  virus3 <- get_virus_reads("/home/k/kwyo0001/mystore/rift/P32705_1002/process/r3.txt")
  colnames(virus1) <- c("bc","virus_L_tot", "virus_L_fwd", "virus_L_rev") #DQ375406  L 6kb 
  colnames(virus2) <- c("bc","virus_S_tot", "virus_S_fwd", "virus_S_rev") #DQ380149  S 1.6kb   nucleocapsid
  colnames(virus3) <- c("bc","virus_M_tot", "virus_M_fwd", "virus_M_rev") #DQ380200  M 3.8kb   glyco/surface
  
  
  allvirus <- data.frame(bc=colnames(adata))
  allvirus <- merge(allvirus, virus1, all.x=TRUE)
  allvirus <- merge(allvirus, virus2, all.x=TRUE)
  allvirus <- merge(allvirus, virus3, all.x=TRUE)
  allvirus[is.na(allvirus)] <- 0
  
  #Also compute log virus counts  
  log_allvirus <- log10(1+allvirus[,-1])
  colnames(log_allvirus) <- paste0("log_",colnames(log_allvirus))
  allvirus <- cbind(allvirus, log_allvirus)
  
  all(colnames(adata) == allvirus$bc ) #yes
  
  #Add to seurat object  
  adata@meta.data <- cbind(adata@meta.data, allvirus[,-1])
  
  colSums(allvirus[,-1])
  
  if(FALSE){
    #For testing once umap etc done
    FeaturePlot(adata, features = c("log_virus_S","log_virus_M","log_virus_L"))
    ggplot(sample(adata@meta.data), aes(log_virus_S, log_virus_M, color=donor)) + geom_point()
    ggplot(adata@meta.data, aes(log_virus_L, log_virus_M, color=donor)) + geom_point()
    ggplot(adata@meta.data, aes(log_virus_S-log_virus_M, donor)) + geom_violin()
  }
  
  
  adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
  
  # Perform the filtering
  adata <- subset(adata, subset = nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 20)
  adata <- NormalizeData(adata, normalization.method = "LogNormalize")
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 3000) #### increase if not finding scaled gene
  VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Identify the most highly variable genes
  top10 <- head(VariableFeatures(adata), 30)
  # plot variable features with and without labels
  if(FALSE){
    top10
    plot1 <- VariableFeaturePlot(adata)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot2
  }
  
  
  adata <- ScaleData(adata)
  adata <- RunPCA(adata)
  adata <- FindNeighbors(adata, dims = 1:30)
  adata <- CellCycleScoring(adata, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  adata <- RunUMAP(adata, dims = 1:30)
  
  saveRDS(adata, "/home/k/kwyo0001/mystore/rift/P32705_1002/processed_adata.RDS")
  
} else {
  adata <- readRDS("/home/k/kwyo0001/mystore/rift/P32705_1002/processed_adata.RDS")
}



##########################
########################## Rename clusters
##########################


adata <- FindClusters(adata, resolution = 0.4)  #### resolution => number of clusters. higher => more clusters
DimPlot(adata, reduction = "umap")

new.cluster.ids <- c(
  "CyT",
  "CyT",
  "CyT",
  "CyT",
  "CyT_column",
  "CyT_column",
  "CyT",
  "CyT",
  "EVT & SynT",
  "EVT & SynT"
)
names(new.cluster.ids) <- levels(adata)
adata <- RenameIdents(adata, new.cluster.ids)
adata$ct <- levels(adata)
DimPlot(adata, reduction = "umap", label=TRUE, group.by = "ct")
DimPlot(adata, reduction = "umap", group.by = "Phase")
DotPlot(adata, features = c("log_virus1","log_virus2","log_virus3"), group.by = "Phase") + coord_flip()


###################################################################################################
########################### Plots #################################################################
###################################################################################################

##########################
########################## #### Response to virus
##########################

adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", group.by = "is_inf")

###2025-08-12###
#######Cell Markers - Van Riel et all + Vento#######
list_genes_cyt <- c("GATA2","GATA3","EGFR","PAGE4","TEAD4","TP63","CDH1","MKI67","ITGB6","LPCAT1","HLA-G","CSH1","ADAMTS20","ACAN","ERBB2","SERPINE2","PLAC8","PAPPA","PAPPA2","CGA")
DotPlot(adata, features = list_genes_cyt,  dot.scale = 5, cols = c("white","red")) + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem


#######IFN family#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_inf <- c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNB1","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNL4")
DotPlot(adata, features = list_genes_inf, group.by = "is_inf", dot.scale = 7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem


#######IFN related genes#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_isg <- c("IFIT1","IFIT2","IFIT3","EPSTI1","MX1","OAS1","OAS2","OAS3","APOL1","APOL16","ISG15","IFI27","PARP12","OASL","EIL2AK2","HELZ2","BST2","USP18","IRF7","IFITM1","IFITM3","SLC15A3","IL15RA","BATF2","ISG20","FCGR1A","MAP3K14","ETV7","IFNLR1","RSAD2","IFI35","TRIM25","TNFAIP3")
DotPlot(adata, features = list_genes_isg, group.by = "is_inf", dot.scale =7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem


#######Apoptosis#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_app <- c("BID","BAD","BBC3","PMAIP1","BCL2L11","BIK","BMF","BNIP3","HRK","BNIP3L","BAX","BAK1","BOK","CASP2","CASP8","CASP9","CASP10","APAF1","CASP3","CASP6","CASP7","CFLAR","BIRC2","BIRC3","XIAP","BCL2","BCL2L1","BCL2L2","MCL1","BCL2A1")
DotPlot(adata, features = list_genes_app, group.by = "is_inf", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Necrosis and Pyraptosis#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_nec <- c("RIPK1","RIPK3","MLKL","ZBP1","CASP1","CASP4","CASP5","PYCARD","GSDMD","GSDME","MEFV","IL1B","IL18","NEK7")
DotPlot(adata, features = list_genes_nec, group.by = "is_inf", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem

#######Preeclamsia#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_pree <- c("FGF14","ZNF295-AS1","MCM8","MUC21","MUC22","GPR45","TGFBRAP1","LZTS1","TMEM97P2","RUNX1","BAMBI","ADRA1D","SCN2B","SCN4B","MYCB2","INVS","SFR1","ERP44","WWTR1")
DotPlot(adata, features = list_genes_pree, group.by = "is_inf", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem


#Look at virus infection versus cell cycle
table(adata$is_inf, adata$Phase)

table(adata$is_inf, adata$Phase)/dim(adata)[2]

phaseVsInfection = table(adata$is_inf, adata$Phase) %>% as.data.frame() 
colnames(phaseVsInfection) = c("infectionStatus", "Phase", "Freq")

ggplot(phaseVsInfection, aes(x = Phase, y = Freq, fill = infectionStatus))+
  geom_bar(stat = "identity", position = "dodge")
