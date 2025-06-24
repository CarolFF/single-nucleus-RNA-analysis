library(Seurat) 
library(Matrix)
library(ggplot2)
library(dplyr)

#data input
path = "X101SC24081731-Z01-J001/Result-X101SC24081731-Z01-J001-B1-1/2.Summary/"
sam <- list.files(path)

sceList <- list()
for (i in 1:length(sam)) {
  sceList[[i]] <- CreateSeuratObject(counts = Read10X(paste0(path,sam[i],"/filtered_feature_bc_matrix")), 
                                     project = sam[i],
                                     assay = "RNA", 
                                     min.cells = 10,
                                     min.features = 500)
}
  
sce.merge <- merge(sceList[[1]],y = sceList[-1],project = "snRNA")

head(sce.merge@meta.data)
table(sce.merge@meta.data$orig.ident)

#QC
sce.merge[["percent.mt"]] <- PercentageFeatureSet(sce.merge, pattern = "^[Mm][Tt]-")

sce.merge[["percent.rp"]] <- PercentageFeatureSet(sce.merge, pattern = "^[Rr][Pp][LlSs]")

head(sce.merge@meta.data)

p1 <- VlnPlot(sce.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("01seurat/01QC/01QC_before.pdf",plot = p1,width = 9,height = 6)
ggsave("01seurat/01QC/01QC_before.png",plot = p1,width = 9,height = 6)

sn <- subset(sce.merge, subset = nFeature_RNA < 6000 & 
                nCount_RNA < 30000 & 
                percent.mt < 8)

p2 <- VlnPlot(sn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("01seurat/01QC/02QC_after.pdf",plot = p2,width = 9,height = 6)
ggsave("01seurat/01QC/02QC_after.png",plot = p2,width = 9,height = 6)


#添加分组信息
sn@meta.data$group[sn@meta.data$orig.ident == "CKO1"] <- "CKO"
sn@meta.data$group[sn@meta.data$orig.ident == "CKO2"] <- "CKO"
sn@meta.data$group[sn@meta.data$orig.ident == "CKO3"] <- "CKO"
sn@meta.data$group[sn@meta.data$orig.ident == "FLOX1"] <- "FLOX"
sn@meta.data$group[sn@meta.data$orig.ident == "FLOX2"] <- "FLOX"
sn@meta.data$group[sn@meta.data$orig.ident == "FLOX3"] <- "FLOX"

#标准化
sn1 <- sn %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(sn)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50)

sn1 <- sn1 %>% FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.5,0.8,1.0,1.2))

print("RNA_snn_res.0.5:")
table(sn1$RNA_snn_res.0.5)

print("RNA_snn_res.0.8:")
table(sn1$RNA_snn_res.0.8)

print("RNA_snn_res.1:")
table(sn1$RNA_snn_res.1)

print("RNA_snn_res.1.2:")
table(sn1$RNA_snn_res.1.2)

sn1$RNA_snn_res.0.5 <- factor(sn1$RNA_snn_res.0.5,levels = c(0:35))
sn1$RNA_snn_res.0.8 <- factor(sn1$RNA_snn_res.0.8,levels = c(0:37))
sn1$RNA_snn_res.1 <- factor(sn1$RNA_snn_res.1,levels = c(0:38))
sn1$RNA_snn_res.1.2 <- factor(sn1$RNA_snn_res.1.2,levels = c(0:45))


col.cluster <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B', '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
                 '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764', '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
                 '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A', '#62615A','#B82129','#66762E','#7491ca')


llc_cols <- c('#b3d3e8', '#259b39', '#e4361e',
              '#ed996d', '#7491ca', '#b85919')

p3 <- DimPlot(sn1,group.by = "orig.ident",raster=FALSE,cols = llc_cols)+ggtitle("Merge samples")
ggsave("01seurat/02merge/01Merge_samples.pdf",width = 7,height = 7,plot = p3)
ggsave("01seurat/02merge/01Merge_samples.png",width = 7,height = 7,plot = p3)

p4 <- DimPlot(sn1,group.by = "group",raster=FALSE,cols = c('#F4B2BC', '#AFDCE0'))+ggtitle("Merge groups")
ggsave("01seurat/02merge/01Merge_groups.pdf",width = 7,height = 7,plot = p4)
ggsave("01seurat/02merge/01Merge_groups.png",width = 7,height = 7,plot = p4)

p7 <- DimPlot(sn1,group.by = "RNA_snn_res.0.5",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.5")
ggsave("01seurat/02merge/03cluster05.pdf",width = 7,height = 7,plot = p7)
ggsave("01seurat/02merge/03cluster05.png",width = 7,height = 7,plot = p7)

p7 <- DimPlot(sn1,group.by = "RNA_snn_res.0.8",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.8")
ggsave("01seurat/02merge/03cluster08.pdf",width = 7,height = 7,plot = p7)
ggsave("01seurat/02merge/03cluster08.png",width = 7,height = 7,plot = p7)

p7 <- DimPlot(sn1,group.by = "RNA_snn_res.1",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.1")
ggsave("01seurat/02merge/03cluster10.pdf",width = 7,height = 7,plot = p7)
ggsave("01seurat/02merge/03cluster10.png",width = 7,height = 7,plot = p7)

p7 <- DimPlot(sn1,group.by = "RNA_snn_res.1.2",raster=FALSE,label = T,cols = c(col.cluster,llc_cols))+ggtitle("RNA_snn_res.1.2")
ggsave("01seurat/02merge/03cluster12.pdf",width = 7,height = 7,plot = p7)
ggsave("01seurat/02merge/03cluster12.png",width = 7,height = 7,plot = p7)

saveRDS(sn1,"01seurat/02merge/00merge.rds")

#find markers
Idents(sn1) <- "RNA_snn_res.1.2"

markers <- FindAllMarkers(sn1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"01seurat/02merge/04markers.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"01seurat/02merge/04markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"01seurat/02merge/04markers_top30.csv")


###########top10 markers
top10 <- split(top10,top10$cluster)

for (i in 1:length(top10)){
  DotPlot(sn1,features = top10[[i]]$gene)+RotatedAxis()
  ggsave(filename = paste0("01seurat/02merge/01marker/01dotplot/dotplot_cluster_",i-1, ".png"),bg = "white")
  FeaturePlot(sn1,features = top10[[i]]$gene,raster = F,ncol = 4,cols = c("lightgrey","#ff0000"))
  ggsave(filename = paste0("01seurat/02merge/01marker/02featureplot/featureplot_cluster_",i-1, ".png"),height = 9,width = 12,bg = "white")
  
}

##########top30 GO
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(dplyr)

top30<- split(top30,top30$cluster)

for (i in 1:length(top30)){
  genes <- top30[[i]]$gene
  DEGs_entrez_id <- mapIds(x=org.Mm.eg.db,column = "ENTREZID",keys = genes,keytype = "SYMBOL")
  
  egoBP <- enrichGO(DEGs_entrez_id, OrgDb=org.Mm.eg.db, ont='BP',
                    pAdjustMethod='BH', pvalueCutoff=0.5,
                    qvalueCutoff=0.5, keyType='ENTREZID',readable = T)
  
  write.csv(egoBP@result,file = paste0("01seurat/02merge/01marker/03GO/cluster",i-1,"_BP.csv"),row.names = F)
  
  dotplot(egoBP,title= paste0('Biological process of cluster',i-1),showCategory=10) + scale_y_discrete(labels=function(enrich_go_BP) stringr::str_wrap(enrich_go_BP, width=60))
  ggsave(filename = paste0("01seurat/02merge/01marker/03GO/cluster",i-1,"_BP.png"),bg = "white",width = 8,height = 8)
  
}

##anno markers

ASPC = c("Dcn","Col1a1","Fbn1","Mmp2")
Adipocyte = c("Adipoq","Me1","Gpam","Dgat2")
Mac=c("Cd14","Cybb","Cd68","C1qb")
Mast = c("Cpa3","Hpgd","Ms4a2","Kit")
Neutrophil = c("G0s2","S100a9","Csf3r","S100a8")
DC = c("Cd1c","Flt3","Il1b","Bcl2a1")
Tcell = c("Il7r","Cd3d","Cd2","Stat4")
B = c("Ms4a1","Fcrl1","Cd79b","Igkc")
Endothelial = c("Vwf","Flt1","Emcn","Rbp7")
SMC = c("Acta2","Steap4","Gucy1a2","Rgs6")
LEC = c("Prox1","Mmrn1","Reln","Mal")

#CLUSTER
P <- FeaturePlot(sn1,features = Neutrophil,raster=FALSE,cols = c("lightgrey","#ff0000"))

ggsave("01seurat/02merge/01marker/Neutrophil.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/Neutrophil.pdf",plot = P,width = 10,height = 10)

###CLUSTER 1,3
P <- FeaturePlot(sn1,features = ASPC,raster=FALSE,cols = c("lightgrey","#ff0000"))

ggsave("01seurat/02merge/01marker/01ASPC.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/01ASPC.pdf",plot = P,width = 10,height = 10)

###CLUSTER 0,7,9,14,26,29,33
P <- FeaturePlot(sn1,features = Adipocyte,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/02Adipocyte.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/02Adipocyte.pdf",plot = P,width = 10,height = 10)


###CLUSTER 8,16,20,23,31,38
P <- FeaturePlot(sn1,features = Mac,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/03Mac.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/03Mac.pdf",plot = P,width = 10,height = 10)

##macrophage
P <- FeaturePlot(sn1,features = DC,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/06DC.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/06DC.pdf",plot = P,width = 10,height = 10)

###CLUSTER 28
P <- FeaturePlot(sn1,features = Mast,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/04Mast.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/04Mast.pdf",plot = P,width = 10,height = 10)



###CLUSTER 5,27,30,32
P <- FeaturePlot(sn1,features = Tcell,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/07Tcell.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/07Tcell.pdf",plot = P,width = 10,height = 10)

###CLUSTER 13
P <- FeaturePlot(sn1,features = B,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/08B.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/08B.pdf",plot = P,width = 10,height = 10)


###CLUSTER 12
plasma <- c("Jchain","Mzb1","Igkc","Ighm")
P <- FeaturePlot(sn1,features = plasma,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/12Plasma.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/12Plasma.pdf",plot = P,width = 10,height = 10)


###CLUSTER 2,4,6,15,17,34
P <- FeaturePlot(sn1,features = Endothelial,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/09Endothelial.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/09Endothelial.pdf",plot = P,width = 10,height = 10)


###CLUSTER 10,18,21,24,25,37
P <- FeaturePlot(sn1,features = SMC,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/10SMC.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/10SMC.pdf",plot = P,width = 10,height = 10)

###CLUSTER 11,19,22
P <- FeaturePlot(sn1,features = LEC,raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/11LEC.png",plot = P,width = 10,height = 10)
ggsave("01seurat/02merge/01marker/11LEC.pdf",plot = P,width = 10,height = 10)


P <- FeaturePlot(sn1,features = "Pdgfra",raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/01ASPC_PDGFRA.png",plot = P,width = 7,height = 7)
ggsave("01seurat/02merge/01marker/01ASPC_PDGFRA.pdf",plot = P,width = 7,height = 7)

P <- FeaturePlot(sn1,features = c("Steap4","Myocd"),raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/pericyte.png",plot = P,width = 12,height = 7)
ggsave("01seurat/02merge/01marker/pericyt.pdf",plot = P,width = 12,height = 7)

P <- FeaturePlot(sn1,features = c("Msln","Cybb"),raster=FALSE,cols = c("lightgrey","#ff0000"))
ggsave("01seurat/02merge/01marker/Msln.png",plot = P,width = 12,height = 7)
ggsave("01seurat/02merge/01marker/Msln.pdf",plot = P,width = 12,height = 7)

#anno
#############rename
celltype=data.frame(ClusterID=0:45,
                    celltype= 'Removed') 

celltype[celltype$ClusterID %in% c(0,3,13,22,32,41),2]='ASPC' 
celltype[celltype$ClusterID %in% c(1,4,5,7,31,38),2]='Adipocyte' 
celltype[celltype$ClusterID %in% c(2,6),2]='B cell' 
celltype[celltype$ClusterID %in% c(8,26,27,37),2]='EC' 
celltype[celltype$ClusterID %in% c(9,10,11,16,18,36),2]='Macrophage' 
celltype[celltype$ClusterID %in% c(17),2]='NK cell' 
celltype[celltype$ClusterID %in% c(12,15,21),2]='T cell' 
celltype[celltype$ClusterID %in% c(14),2]='Fibroblast' 
celltype[celltype$ClusterID %in% c(19),2]='Pericyte' 
celltype[celltype$ClusterID %in% c(20),2]='Cycling cell' 
celltype[celltype$ClusterID %in% c(24,30),2]='DC' 
celltype[celltype$ClusterID %in% c(25),2]='Myoepithelial' 
celltype[celltype$ClusterID %in% c(28),2]='Neuron' 
celltype[celltype$ClusterID %in% c(29),2]='LEC' 
celltype[celltype$ClusterID %in% c(34),2]='SMC'
celltype[celltype$ClusterID %in% c(40),2]='Uknown' 
celltype[celltype$ClusterID %in% c(44),2]='Mast cell' 


sn1@meta.data$celltype1 = "NA"
for(i in 1:nrow(celltype)){
  sn1@meta.data[which(sn1@meta.data$RNA_snn_res.1.2 == celltype$ClusterID[i]),'celltype1'] <- celltype$celltype[i]}

table(sn1@meta.data$celltype1)
sn1$group <- factor(sn1$group,levels = c("FLOX","CKO"))
sn1$celltype1 <- factor(sn1$celltype1,levels = c("ASPC",
                                                 "Adipocyte",
                                                 "B cell",
                                                 "NK cell",
                                                 "T cell",
                                                 "Macrophage",
                                                 "DC",
                                                 "Mast cell",
                                                 "EC",
                                                 "LEC",
                                                 "Fibroblast",
                                                 "Pericyte",
                                                 "Myoepithelial",
                                                 "SMC",
                                                 "Neuron",
                                                 "Cycling cell",
                                                 "Uknown",
                                                 "Removed"))

col.cluster <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0',
                 '#BDBF3B','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8',
                 '#F7CCDD','#A3D7E5','#B04E4F','#7B8380')

Idents(sn1) <- "celltype1"
p <- DimPlot(sn1,group.by = "celltype1",cols = col.cluster,label = T,repel = T,raster = F)
ggsave("01seurat/02merge/02anno/00UMAP_celltype1.png",width = 8,height = 8,plot = p)
ggsave("01seurat/02merge/02anno/00UMAP_celltype1.pdf",width = 8,height = 8,plot = p)

p <- DimPlot(sn1,group.by = "celltype1",split.by = "group",cols = col.cluster,label = T,repel = T,raster = F)
ggsave("01seurat/02merge/02anno/00UMAP_celltype1_group.png",width = 12,height = 8,plot = p)
ggsave("01seurat/02merge/02anno/00UMAP_celltype1_group.pdf",width = 12,height = 8,plot = p)

saveRDS(sn1,"01seurat/02merge/00merge_anno1.rds")

sn1.sub <- subset(sn1,idents = "Removed",invert = T)

sn1.sub <- sn1.sub %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(sn1.sub)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50)

Idents(sn1.sub) <- "celltype1"
p <- DimPlot(sn1.sub,group.by = "celltype1",cols = col.cluster,label = T,repel = T,raster = F)
ggsave("01seurat/02merge/02anno/01UMAP_celltype1.png",width = 8,height = 8,plot = p)
ggsave("01seurat/02merge/02anno/01UMAP_celltype1.pdf",width = 8,height = 8,plot = p)

p <- DimPlot(sn1.sub,group.by = "celltype1",split.by = "group",cols = col.cluster,label = T,repel = T,raster = F)
ggsave("01seurat/02merge/02anno/01UMAP_celltype1_group.png",width = 12,height = 8,plot = p)
ggsave("01seurat/02merge/02anno/01UMAP_celltype1_group.pdf",width = 12,height = 8,plot = p)

####proportion
##by sample
metadata <- sn1.sub@meta.data
ggplot(data = metadata, aes(x = orig.ident, fill=celltype1)) +
  geom_bar(position = 'fill', width = 0.4) +
  facet_wrap(~group,scales = 'free_x') +
  scale_fill_manual(values =  col.cluster)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))+
  labs(x="",y="",title="")

ggsave("01seurat/02merge/02anno/02prop_by_samples.pdf",width = 12,height = 9)
ggsave("01seurat/02merge/02anno/02prop_by_samples.png",width = 12,height = 9)

##by group
ggplot(data = metadata, aes(x = group, fill=celltype1)) +
  geom_bar(position = 'fill', width = 0.4) +
  scale_fill_manual(values =  col.cluster)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))+
  labs(x="",y="",title="")

ggsave("01seurat/02merge/02anno/02prop_by_group.pdf",width = 10,height = 9)
ggsave("01seurat/02merge/02anno/02prop_by_group.png",width = 10,height = 9)

df <- table(sn1.sub$celltype1,sn1.sub$orig.ident)
write.csv(df,"01seurat/02merge/02anno/02ncells_celltype1.csv")

df2 <- prop.table(df,margin=2)
write.csv(df2,"01seurat/02merge/02anno/02proportions_celltype1.csv")

saveRDS(sn1.sub,"01seurat/02merge/00merge_anno2.rds")


Idents(sn1.sub) <- "celltype1"
markers <- FindAllMarkers(sn1.sub)
write.csv(markers ,file = "01seurat/02merge/02anno/03markers~celltype1.csv")

top50 <- markers %>% group_by(cluster) %>% top_n(n=50, wt = avg_log2FC)
write.csv(top50,file = "01seurat/02merge/02anno/03markers~celltype1_top50.csv")

top30 <- markers %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)
write.csv(top50,file = "01seurat/02merge/02anno/03markers~celltype1_top50.csv")

top10 <- markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
write.csv(top10,file = "01seurat/02merge/02anno/03markers~celltype1_top10.csv")

top5 <- markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
write.csv(top10,file = "01seurat/02merge/02anno/03markers~celltype1_top5.csv")


###markers
p_heat <- DoHeatmap(subset(sn1.sub, downsample = 200), ,features = unique(top5$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("01seurat/02merge/02anno/04Marker_heatmap.pdf",plot = p_heat,width = 10,height = 10)
ggsave("01seurat/02merge/02anno/04Marker_heatmap.png",plot = p_heat,width = 10,height = 10)

##dotplot
top4 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)
write.csv(top4,"01seurat/02merge/02anno/04Marker_top4.csv")

p_dot <-DotPlot(sn1.sub, features=unique(top4$gene),cols=c('gray95','#ff0000'))+coord_flip()+RotatedAxis()+
  xlab("")+ylab("")

ggsave("01seurat/02merge/02anno/04Marker_dot.pdf",plot = p_dot,
       width = 10,height = 12)
ggsave("01seurat/02merge/02anno/04Marker_dot.png",plot = p_dot,
       width = 10,height = 12,bg = "white")

##vlnplot
markers.plot <- c("Pdgfra",
             "Adipoq",
             "Ms4a1",
             "Klrk1",
             "Themis",
             "Mrc1",
             "Flt3",
             "Cpa3",
             "Ccdc85a",
             "Reln",
             "Mgp",
             "Steap4",
             "Fhod3",
             "Myocd",
             "Robo2",
             "Top2a",
             "Timd4")

p_vln <- VlnPlot(sn1.sub,features = markers.plot,stack = T,cols = col.cluster,flip = T)+ NoLegend()+labs(x="")
ggsave("01seurat/02merge/02anno/04Marker_vlnplot2.pdf",plot = p_vln,width = 8,height = 8)
ggsave("01seurat/02merge/02anno/04Marker_vlnplot2.png",plot = p_vln,width = 8,height = 8,bg = "white")
