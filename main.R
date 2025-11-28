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
  
  mat <- ReadParseBio("~/mystore/rift/P32705_1002/DGE_filtered/")  
  adata <- CreateSeuratObject(mat)
  
  ## Add metadata about infection
  adata$is_inf <- as.integer(str_split_fixed(colnames(adata),"_",2)[,1])>=7
  
  ### donor1: 1,2, 7,8
  ### donor2: 3,4, 9,10
  ### donor3: 5,6, 11,12
  
  
  adata$donor <- NA
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
  
  virus1 <- get_virus_reads("~/mystore/rift/P32705_1002/process/r1.txt")
  virus2 <- get_virus_reads("~/mystore/rift/P32705_1002/process/r2.txt")
  virus3 <- get_virus_reads("~/mystore/rift/P32705_1002/process/r3.txt")
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
  
  saveRDS(adata, "~/mystore/rift/P32705_1002/processed_adata.RDS")
  
} else {
  adata <- readRDS("~/mystore/rift/P32705_1002/processed_adata.RDS")
}


DimPlot(adata, reduction = "umap", group.by = "Phase")
DotPlot(adata, features = c("log_virus1","log_virus2","log_virus3"), group.by = "Phase") + coord_flip()

##########################
########################## Rename clusters
##########################

adata <- FindClusters(adata, resolution = 0.4)  #### resolution => number of clusters. higher => more clusters
DimPlot(adata, reduction = "umap")

adata$ct = case_when(adata$seurat_clusters %in% c(0, 1, 2, 3, 6, 7) ~ 'CyT1',
                     adata$seurat_clusters %in% c(4, 5) ~ 'CyT2', 
                     adata$seurat_clusters %in% c(8, 9) ~ 'Juvenile STB',)
table(adata$ct)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "ct")

saveRDS(adata, "~/mystore/rift/P32705_1002/annotated_adata.RDS")



adata$inf_by_newCluster = paste(adata$is_inf, adata$ct, sep = '_')
DimPlot(adata, group.by = 'inf_by_newCluster')
table(adata$inf_by_newCluster)

####### RVFV receptor #######
FeaturePlot(adata, features = c(
  "LRP1"
))




############## Plot virus distribution
FeaturePlot(adata, features = c("log_virus1","log_virus2","log_virus3"))
DotPlot(adata, features = c("log_virus1","log_virus2","log_virus3"), group.by = "donor_is_inf") + coord_flip()
VlnPlot(adata, features = c("log_virus1","log_virus2","log_virus3"), group.by = "donor_is_inf")


##########################
########################## ridge plots
##########################


RidgePlot(adata, features = c("MKI67","TP63"), ncol = 2, group.by = "ct")
RidgePlot(adata, features = c("MKI67","TP63"), ncol = 2, group.by = "Phase")

VlnPlot(adata, features = c("MKI67","TP63"), ncol = 2, group.by = "Phase")
VlnPlot(adata, features = c("MKI67","TP63"), ncol = 2, group.by = "ct")


VlnPlot(adata, features = c("log_virus1"), ncol = 2, group.by = "Phase")
VlnPlot(adata, features = c("log_virus1"), ncol = 2, group.by = "ct")



##########################
########################## #### Custom plots for virus abundance
##########################

ggplot(adata@meta.data, aes(ct, log_virus1)) + geom_point()

ggplot(adata@meta.data, aes(ct, log_virus1)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean") + 
  theme_bw()

ggplot(adata@meta.data, aes(ct, virus1)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  theme_bw()


ggplot(adata@meta.data, aes(Phase, log_virus1)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  theme_bw()

ggplot(adata@meta.data, aes(Phase, virus1)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  theme_bw()


##########################
########################## #### Response to virus
##########################

adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", group.by = "is_inf")

list_genes_inf <- c("EGFR","TFAP2C","YAP1","TP63","ITGA6","CDX2","ELF5","BMP4","PAGE4","TEAD4","MKI67","CCNA2")
DoHeatmap(adata, features = list_genes_inf, group.by = "donor_is_inf")     #, slot = "data")# + RotatedAxis()  use data if scale problem
DotPlot(adata, features = list_genes_inf, group.by = "donor_is_inf") + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem

###2025-08-12###
#######Cell Markers - Van Riel et all + Vento#######
list_genes_cyt <- c("GATA2","GATA3","EGFR","PAGE4","TEAD4","TP63","CDH1","MKI67","ITGB6","LPCAT1","HLA-G","CSH1","ADAMTS20","ACAN","ERBB2","SERPINE2","PLAC8","PAPPA","PAPPA2","CGA")
DotPlot(adata, features = list_genes_cyt, group.by = "is_inf", dot.scale = 5, cols = c("blue","red")) + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem


###2025-08-12###
#######Cell Markers - Van Riel et all#######
list_genes_cyt <- c("TOP2A","MKI67","EGFL6","DLK1","COL1A1","PECAM1","HLA-DRA","CD14","F13A1","PTPRC","CYP19A1","HLA-G","PAPPA2","PAGE4","KRT7","CGA","CGB3","CGB5","CGB7","CGB8","SDC1")
DotPlot(adata, features = list_genes_cyt,  dot.scale = 5, cols = c("white","red")) + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem


K###2025-03-05###
#######Cell Markers#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_cyt <- c("GATA2","GATA3","KRT7","KRT8","KRT18","EGFR","ERVFRD-1","ERVV-1","ERVW-1","PAGE4","BCAM","TP63","CDH1","MKI67","ITGB6","LPCAT1","HLA-G","CSH1","ADAMTS20","ACAN","ERBB2","SERPINE2","PLAC8","PAPPA","TFAP2C","YAP1","GCM1","OVOL1","CDK1","CCNB1","RRM2","CGA","CGB3","CGB5","CGB7","CGB8","SDC1")
DotPlot(adata, features = list_genes_cyt,  dot.scale = 5, cols = c("white","red")) + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem


#######IFN family#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_inf <- c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNB1","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNL4")
DotPlot(adata, features = list_genes_inf, group.by = "is_inf", dot.scale = 7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem

list_genes_inf <- c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNB1","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNL4")
DotPlot(adata, features = list_genes_inf, group.by = "inf_by_newCluster", dot.scale = 7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######IFN related genes#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_isg <- c("IFIT1","IFIT2","IFIT3","EPSTI1","MX1","OAS1","OAS2","OAS3","APOL1","APOL16","ISG15","IFI27","PARP12","OASL","EIL2AK2","HELZ2","BST2","USP18","IRF7","IFITM1","IFITM3","SLC15A3","IL15RA","BATF2","ISG20","FCGR1A","MAP3K14","ETV7","IFNLR1","RSAD2","IFI35","TRIM25","TNFAIP3")
DotPlot(adata, features = list_genes_isg, group.by = "inf_by_newCluster", dot.scale =7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Apoptosis#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_app <- c("BID","BAD","BBC3","PMAIP1","BCL2L11","BIK","BMF","BNIP3","HRK","BNIP3L","BAX","BAK1","BOK","CASP2","CASP8","CASP9","CASP10","APAF1","CASP3","CASP6","CASP7","CFLAR","BIRC2","BIRC3","XIAP","BCL2","BCL2L1","BCL2L2","MCL1","BCL2A1")
DotPlot(adata, features = list_genes_app, group.by = "inf_by_newCluster", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Necrosis and Pyraptosis#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_nec <- c("RIPK1","RIPK3","MLKL","ZBP1","CASP1","CASP4","CASP5","PYCARD","GSDMD","GSDME","MEFV","IL1B","IL18","NEK7")
DotPlot(adata, features = list_genes_nec, group.by = "is_inf", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Wnt#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_wnt <- c("WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10A","WNT10B","WNT11","WNT16","WIF1","WLS","CASP6","FGFR3","GSN","UBQLN1","DKK1","DKK2","DKK3","DKK4","SFRP1","WIF1","SOST","IGFBP4","TPBG","APCDD1","RSPO1","RSPO2","RSPO3","NDP")
DotPlot(adata, features = list_genes_wnt, group.by = "is_inf", dot.scale = 9, cols = c("grey","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Preeclamsia#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_pree <- c("FGF14","ZNF295-AS1","MCM8","MUC21","MUC22","GPR45","TGFBRAP1","LZTS1","TMEM97P2","RUNX1","BAMBI","ADRA1D","SCN2B","SCN4B","MYCB2","INVS","SFR1","ERP44","WWTR1")
DotPlot(adata, features = list_genes_pree, group.by = "inf_by_newCluster", dot.scale = 7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



#######Inflammation#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_infl <- c("A2M","HLA-DRB1","HMGB1","HMGB2","IGHG1","IL6","IL10","IL1F10","IL1RN","IL20RB","IL25","IL2RA","IL31RA","IL36A","IL36B","IL36G","IL37","IL5RA","KDM6B","NLRP6","NOTCH1","NOTCH2","PNMA1","RAB44","RASGRP1","RBPJ","TNF","TREX1")
DotPlot(adata, features = list_genes_infl, group.by = "inf_by_newCluster", dot.scale = 7, cols = c("blue","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem


#######Rab family#######
adata$donor_is_inf <- paste(adata$is_inf, adata$donor)

DimPlot(adata, reduction = "umap", label=TRUE, group.by = "is_inf")

list_genes_rab <- c("RAB1","RAB2","RAB3","RAB4","RAB5","RAB6","RAB7","RAB8","RAB9","RAB10","RAB11","RAB12","RAB13","RAB14","RAB15","RAB16","RAB17","RAB18","RAB19","RAB20","RAB21","RAB22","RAB23","RAB24","RAB25","RAB26","RAB27","RAB28","RAB29","RAB30","RAB31","RAB32","RAB33","RAB34","RAB35","RAB36","RAB37","RAB38","RAB39","RAB40","RAB41","RAB42","RAB43","RAB44","RAB45","RAB46","RAB47","RAB48","RAB49","RAB50","RAB51","RAB52","RAB53","RAB54","RAB55","RAB56","RAB57","RAB58","RAB59","RAB60","RAB61","RAB62","RAB63","RAB64","RAB65","RAB66")
DotPlot(adata, features = list_genes_rab, group.by = "donor_is_inf", dot.scale = 9, cols = c("white","red")) + RotatedAxis()   #, slot = "data")# + RotatedAxis()  use data if scale problem



##########################
########################## #### Response to virus IF RNA DETECTED
##########################


adata$total_virus <- adata$virus1+adata$virus2+adata$virus3
median(adata$total_virus)
mean(adata$total_virus)
max(adata$total_virus)
adata$has_virus <- adata$total_virus>5


table(adata$has_virus, adata$is_inf)

DoHeatmap(adata, features = list_genes_inf, group.by = "has_virus")     #, slot = "data")# + RotatedAxis()  use data if scale problem
DotPlot(adata, features = list_genes_inf, group.by = "has_virus")     #, slot = "data")# + RotatedAxis()  use data if scale problem


plot(sort(adata$virus1+adata$virus2+adata$virus3))



##########################
########################## #### Response to virus RNA amount
##########################

#df <- as.data.frame(adata@assays$RNA@layers$scale.data)
#rownames(df) <- rownames(adata@assays$RNA@layers$scale.data)

df <- as.data.frame(t(LayerData(adata, "scale.data", features="IFIT1")))
df$virus <- adata$total_virus
ggplot(df, aes(log10(1+virus), log10(1+IFIT1))) + geom_point()



##########################
########################## dot plots
##########################


DotPlot(adata, features = c("GATA3","KRT7"), group.by = "ct") + RotatedAxis()




#DoHeatmap(adata, features = c("GATA3"), group.by = "ct", slot = "data")# + RotatedAxis()
DoHeatmap(adata, features = c("GATA3","KRT7"), group.by = "ct", slot = "data")# + RotatedAxis()



############## All cells
Idents(adata) <- adata$seurat_clusters
adata_markers <- FindAllMarkers(adata, min.pct = 0.25, logfc.threshold = 0.25)
top_gene <- adata_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)  #get top genes
View(top_gene)


############## between donors
Idents(adata) <- adata$donor
donor_markers <- FindAllMarkers(adata, min.pct = 0.25, logfc.threshold = 0.25)
donor_top_gene <- donor_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)  #get top genes per donor
View(donor_top_gene)


############## just from uninf
adata_uninf <- FindClusters(adata[,!adata$is_inf], resolution = 0.4)
DimPlot(adata_uninf, reduction = "umap", label=TRUE)

adata_markers_uninf <- FindAllMarkers(adata_uninf, min.pct = 0.25, logfc.threshold = 0.25)
top_gene_uninf <- adata_markers_uninf %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)


FeaturePlot(adata, features = c(
  "TP63",
  "RUNX3",
  "RNY1", #subset
  "SOX9" #def a subset #Hif-1α and SOX9 proteins can be used as a marker to show severity of preeclampsia and regulation of cell proliferation and angiogenesis during placental ...
))

FeaturePlot(adata, features = c(
  "KCNK13",
  "CGA",
  "AREG",
  "PTPRD",
  "H19", #lncRNA;
  "PINCR", #lncRNA;
  "HYDIN" #cilia!!!
))

####### Just infected
adata_inf <- FindClusters(adata[,adata$is_inf], resolution = 0.30)
DimPlot(adata_inf, reduction = "umap", label=TRUE)
adata_markers_inf <- FindAllMarkers(adata_inf, min.pct = 0.25, logfc.threshold = 0.25)
top_gene_inf <- adata_markers_inf %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

FeaturePlot(adata_inf, features = c(
  "PER3", #circadian
  "SYNPO2",  #defines a state. .. of what?? "Synpo2 promotes Arp2/3-dependent lamellipodia formation."
  "HYDIN", #cilia!!!
  "TMPRSS9"
))
### looks like the virus induces cell motility somehow?? SYNPO2 lamellipodia  => HYDIN  cilia

####### Just uninfected
adata_uninf <- FindClusters(adata[,!adata$is_inf], resolution = 0.2)
DimPlot(adata_uninf, reduction = "umap", label=TRUE)




FeaturePlot(adata_uninf, features = c(
  
  "LYVE1",
  
  "IGF2", #EVT
  "DSC3", #even more VCT; contrary to roser atlas :(
  "SDK1",
  "KALRN",
  "MAPK4", #VCT!
  "IGSF5" #very few SCT?
))



FeaturePlot(adata_uninf, features = c(
  "TNFRSF11A",
  "IL2RA", "IL12RB2", "IL1RL1"
  #"AIF1", "CSF1", "CD33", "CD163"
))

#
# HB, Hofbauer cells ?? fetal-derived placental macrophages


FeaturePlot(adata, features = c(
  "SH3RF3", #P/F
  "TNNI1",  #fibroblast! !
  "RTN1", #P/S
  "CACNA2D4", #HB/M
  "CACNA1C", #P/S
  "NRG2", #FB/S
  "LINC01121", #not in roser but specific
  "LINC02613" #not in roser but specific
))


FeaturePlot(adata, features = c(
  "CCND2",
  
  #two other clusters
  "PINCR",  #two clusters!! PINCR (P53-Induced Noncoding RNA) is an RNA Gene, and is affiliated with the lncRNA class.  ---- not in 10x!
  #  "REEP",
  #"ENTREP2",  P/HB?
  
  #### 5
  "NUPR1", #p/S/epi 
  "ANK1", #HB
  "H1-0", #linker histone
  
  
  #### 4
  "AREG" #specific for group 4
))


FeaturePlot(adata, features = c(
  "APOBR",#very specific; EVT subset in roser
  "ATXN1",
  "PEG10",  #super specific for VCT!!!!!
  "CDK14", #all over
  "PLAC8",#very EVT in roser
  "TRPC4",#"endothelial barrier",  -- S/P in roser
  "OAS1","HOPX",
  "MUC15",#roser, very SCT
  "MYT1",#useless
  "BST2", #p/hb/endo etc
  "MX1", #p/endo/etc; up in multiple here. several infected cell types?
  "IRF7",
  "SP110"
))
#S=stromal; #P=perivascular


#placenta infections Hoo, Ruiz-Morales et al., Cell Systems 2024  




VlnPlot(adata, features = c(
  "TP63",
  "NR2F2",
  "MCTP2",
  "EVL",
  "NEURL1",
  "IL34",
  "CACNA1C",
  "TOX3",
  "SEMA3D",
  "BATF2",
  "MX2",
  "HLA-DPB1",
  "IGSF5",
  "CGB5"
))




VlnPlot(adata, features = c(
  "TP63",
  "NR2F2",
  "MCTP2",
  "EVL",
  "NEURL1",
  "IL34",
  "CACNA1C",
  "TOX3",
  "SEMA3D",
  "BATF2",
  "MX2",
  "HLA-DPB1",
  "IGSF5",
  "CGB5"
))


#, group.by = "tree.ident")

FeaturePlot(adata, features = c(
  "TP63",
  "NR2F2",
  "MCTP2",
  "EVL",
  "NEURL1",
  "IL34",
  "CACNA1C",
  "TOX3",
  "SEMA3D",
  "BATF2",
  "MX2",
  "HLA-DPB1",
  "IGSF5",
  "CGB5"
))

#p63 maintains the proliferative CTB state, at least partially through regulation of epithelial-to-mesenchymal transition, cell adhesion, and matrix
#CTB=cytotrophoblast
#compare https://maternal-fetal-interface.cellgeni.sanger.ac.uk

#TP63 => VCT marker
#NR2F2 => multiple cell types; VCT is one
#MCTP2 => NK marker????

#NEURL1 => EVT/SCT
#IL34 many things
#CACNA1C => dS/dP =  ...
#TOX3 => epi / dS?
#SEMA3D => endo L?
#batf2 => quite mixed
#MX2 =>  all over; dM?
#CGB5 => tip of SCT
#IGSF5 => SCT



#top markers from roser!
#https://www.nature.com/articles/s41586-018-0698-6   Fig 2
FeaturePlot(adata, features = c(
  "EGFR", #mainly SCT
  "ACKR2", #EVT and SCT
  "CCL5", #t cell, dNK, 
  "PGF", 
  "MET", #SCT
  
  "HGF", #F1/F2, barely present
  "TGFB1", #EVT .. all over here
  "TGFBR2", #endo / EVT  -- can we tell EVT and SCT?
  "TGFBR1", 
  "TGFB1",
  "CXCR6"
))
#EVT => VCT => SCT


#https://www.nature.com/articles/s41586-018-0698-6


#gonadotropin: syncytiotrophoblastic cells. they seem to be a small island
FeaturePlot(adata, features = c(
  #gonadotropin
  "CGA",
  "CGB1",
  "CGB2",
  "CGB3",
  "CGB5",
  "CGB7",
  "CGB8"
))

DotPlot(adata, features = c(
  #gonadotropin
  "CGA",
  "CGB1",
  "CGB2",
  "CGB3",
  "CGB5",
  "CGB7",
  "CGB8"
)) + RotatedAxis()


FeaturePlot(adata, features = c(
  "IFNL1",
  "IFNL2",
  "IFNL3",
  "IFNL4"
))


FeaturePlot(adata, features = c(
  "TEAD4",
  "HLA-G"
))


FeaturePlot(adata, features = c(
  #  "MKI67",
  "NUPR1",
  "PINCR"
))

VlnPlot(adata, group.by =  "donor", features = c(
  "7SK"
))
FeaturePlot(adata, features = c(
  "KCNK13"
))

FeaturePlot(adata, features = c(
  "HLA-G","NOTCH1","NOTCH2","PRG2","ERBB2","TP53","TRPC4","MX1","HLA-DPB1"
))


#20241019 - CyT marker
DotPlot(adata, features = c(
  "EGFR","TFAP2C","YAP1","TP63","ITGA6","CDX2","ELF5","BMP4","PAGE4","TEAD4","MKI67","CCNA2"
))
#SynT marker
DotPlot(adata, features = c(
  "ERVFRD-1","ERVV-1","SDC1","TFAP2A","GATA2","GATA3","SLC22A11","PHLDA2","FXYD3","PNP","GCM1","ENDOU"
))
#Evt Marker
DotPlot(adata, features = c(
  "KRt77","HLA-G","ITGA1","ITGA5","NOTCH1","NOTCH2",
))



#Wnt list
FeaturePlot(adata, features = c(
  "WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10A","WNT10B","WNT11","WNT16"
))
#sorted Wnt
FeaturePlot(adata, features = c(
  "WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT9A","WNT10A","WNT11",
))
#Wnt7A related genes
FeaturePlot(adata, features = c(
  "WIF1","WLS","CASP6","FGFR3","GSN","UBQLN1"
))  
#Wnt antagonist and agonist
FeaturePlot(adata, features = c(
  "DKK1","DKK2","DKK3","DKK4","SFRP1","WIF1","SOST","WISE","IGFBP4","TPBG","APCDD1","RSPO1","RSPO2","RSPO3","NDP"
))  
#SORTED Wnt antagonist and agonist
FeaturePlot(adata, features = c(
  "DKK2","DKK3","SFRP1","WIF1","SOST","IGFBP4","TPBG","RSPO2"
))


#Immune response, Interferon
FeaturePlot(adata, features = c(
  "IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNB1","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNL4"
))
#SORT IFN
FeaturePlot(adata, features = c(
  "IFNA4","IFNA5","IFNA7","IFNA8","IFNA10","IFNA16","IFNB1","IFNK","IFNL1","IFNL2","IFNL3","IFNL4"
))
#Interferon related gene
FeaturePlot(adata, features = c(
  "MX1","OAS1","IFIT1","IFIT2","IFIT3","HLA-B","HLA-C","HLA-DPB1","HLA-DQA1","HLA-E","HLA-F","HLA-G"
))

#APOPTOSIS GENE 1
FeaturePlot(adata, features = c(
  "BID",
  "BAD",
  "BBC3",
  "PMAIP1",
  "BCL2L11",
  "BIK",
  "BMF",
  "BNIP3",
  "HRK",
  "BNIP3L",
  "BAX",
  "BAK1",
  "BOK"
))

#APOPTOSIS GENE 2
FeaturePlot(adata, features = c(
  "CASP2",
  "CASP8",
  "CASP9",
  "CASP10",
  "APAF1",
  "CASP3",
  "CASP6",
  "CASP7",
  "CFLAR",
  "BIRC2",
  "BIRC3",
  "XIAP"
))

#APOPTOSIS GENE 3
FeaturePlot(adata, features = c(
  "BCL2",
  "BCL2L1",
  "BCL2L2",
  "MCL1",
  "BCL2A1"
))

#NECROSIS GENE
FeaturePlot(adata, features = c(
  "RIPK1",
  "RIPK3",
  "MLKL",
  "ZBP1"
))

#PYROPTOSIS GENE
FeaturePlot(adata, features = c(
  "CASP1",
  "CASP4",
  "CASP5",
  "PYCARD",
  "GSDMD",
  "GSDME",
  "MEFV",
  "IL1B",
  "IL18",
  "NEK7"
))

#Inflammatory cytokine
FeaturePlot(adata, features = c(
  "IL6",
  "IL6R",
  "CLCF1",
  
))


#VEEV receptor
FeaturePlot(adata, features = c(
  "LDLRAD3",
  "VLDLR",
  "APOER2"
  
))
#another map of placenta: https://www.nature.com/articles/s41422-018-0066-y

#Look at virus infection versus cell cycle
table(adata$is_inf, adata$Phase)

table(adata$is_inf, adata$Phase)/dim(adata)[2]

phaseVsInfection = table(adata$is_inf, adata$Phase) %>% as.data.frame() 
colnames(phaseVsInfection) = c("infectionStatus", "Phase", "Freq")

ggplot(phaseVsInfection, aes(x = Phase, y = Freq, fill = infectionStatus))+
  geom_bar(stat = "identity", position = "dodge")



##########################
########################## Rename clusters, uninf, 2025-9-29
##########################
adata_uninf <- FindClusters(adata[,!adata$is_inf], resolution = 0.5)
DimPlot(adata_uninf, reduction = "umap", label=TRUE)

new.cluster.ids <- c(
  "CyT",
  "CyT",
  "CyT",
  "CyT",
  "CyT",
  "CyT_column",
  "EVT & SynT",
  "CyT_column"
)
names(new.cluster.ids) <- levels(adata_uninf)
adata_uninf <- RenameIdents(adata_uninf, new.cluster.ids)
adata_uninf$ct <- levels(adata_uninf)

#######Cell Markers - Van Riel et all + Vento#######
list_genes_cyt <- c("GATA2","GATA3","EGFR","PAGE4","TEAD4","TP63","CDH1","MKI67","ITGB6","LPCAT1","HLA-G","CSH1","ADAMTS20","ACAN","ERBB2","SERPINE2","PLAC8","PAPPA","PAPPA2","CGA")
DotPlot(adata_uninf, features = list_genes_cyt,  dot.scale = 5, cols = c("white","red")) + RotatedAxis()    #, slot = "data")# + RotatedAxis()  use data if scale problem


