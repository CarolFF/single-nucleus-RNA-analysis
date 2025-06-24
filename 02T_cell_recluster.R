library(Seurat) 
library(Matrix)
library(ggplot2)
library(dplyr)

###T cell subtype
Tcell <- subset(sn1.sub,idents = "T cell")
dim(Tcell) #3254

p <- VlnPlot(Tcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("01seurat/02merge/03Tcell/01QC_T_before.pdf",plot = p,width = 9,height = 6)
ggsave("01seurat/02merge/03Tcell/01QC_T_before.png",plot = p,width = 9,height = 6)

Tcell <- subset(Tcell, subset = nFeature_RNA < 4000 & 
                  nCount_RNA < 10000)
dim(Tcell) #3221

p <- VlnPlot(Tcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("01seurat/02merge/03Tcell/01QC_T_after.pdf",plot = p,width = 9,height = 6)
ggsave("01seurat/02merge/03Tcell/01QC_T_after.png",plot = p,width = 9,height = 6)

Tcell <- Tcell %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(Tcell)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.3,0.4,0.5))

print("RNA_snn_res.0.3:")
table(Tcell$RNA_snn_res.0.3)

print("RNA_snn_res.0.4:")
table(Tcell$RNA_snn_res.0.4)

print("RNA_snn_res.0.5:")
table(Tcell$RNA_snn_res.0.5)


p<- DimPlot(Tcell,group.by = "RNA_snn_res.0.3",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.3")
ggsave("01seurat/02merge/03Tcell/02cluster03.pdf",width = 7,height = 7,plot = p)
ggsave("01seurat/02merge/03Tcell/02cluster03.png",width = 7,height = 7,plot = p)

p<- DimPlot(Tcell,group.by = "RNA_snn_res.0.4",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.4")
ggsave("01seurat/02merge/03Tcell/02cluster04.pdf",width = 7,height = 7,plot = p)
ggsave("01seurat/02merge/03Tcell/02cluster04.png",width = 7,height = 7,plot = p)

p<- DimPlot(Tcell,group.by = "RNA_snn_res.0.5",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.5")
ggsave("01seurat/02merge/03Tcell/02cluster05.pdf",width = 7,height = 7,plot = p)
ggsave("01seurat/02merge/03Tcell/02cluster05.png",width = 7,height = 7,plot = p)

saveRDS(Tcell,"01seurat/02merge/03Tcell/00Tcell.rds")

#find markers
Idents(Tcell) <- "RNA_snn_res.0.5"

markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"01seurat/02merge/03Tcell/03markers.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,"01seurat/02merge/03Tcell/03markers_top5.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"01seurat/02merge/03Tcell/03markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"01seurat/02merge/03Tcell/03markers_top30.csv")

p_heat <- DoHeatmap(Tcell, ,features = unique(top5$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("01seurat/02merge/03Tcell/04Marker_heatmap.pdf",plot = p_heat,width = 12,height = 8)
ggsave("01seurat/02merge/03Tcell/04Marker_heatmap.png",plot = p_heat,width = 12,height = 8)


#cd4&8

cdt <- subset(Tcell,idents = c(0,1,2,7))
cdt$RNA_snn_res.0.5 <- Idents(cdt)
cdt$celltype2 <- "Cd8"
cdt@meta.data[which(cdt@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype2'] <- "Cd4"

metadata <- cdt@meta.data
ggplot(data = metadata, aes(x = group, fill=celltype2)) +
  geom_bar(position = 'fill', width = 0.4) +
  #facet_wrap(~group,scales = 'free_x') +
  scale_fill_manual(values =  col.cluster)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))+
  labs(x="",y="",title="")
ggsave("01seurat/02merge/03Tcell/pro_cd4_8_optional.pdf",width = 4,height = 6)
ggsave("01seurat/02merge/03Tcell/pro_cd4_8_optional.png",width = 4,height = 6)

#Th1:

FeaturePlot(Tcell,features = c("Cxcr3","Tbx21"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("Th1.pdf",width = 12,height = 7)
ggsave("Th1.png",width = 12,height = 7)

FeaturePlot(Tcell,features = c("Trdc","Ccr7"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("γδT.pdf",width = 12,height = 7)
ggsave("γδT.png",width = 12,height = 7)

#Tn
FeaturePlot(Tcell,features = c("Ccr7","Sell"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("Tn.pdf",width = 12,height = 7)
ggsave("Tn.png",width = 12,height = 7)

FeaturePlot(Tcell,features = c("Gata3","Bcl11b","Il1rl1"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("IL2C.pdf",width = 12,height = 10)
ggsave("IL2C.png",width = 12,height = 10)

#T anno
celltype=data.frame(ClusterID=0:7,
                    celltype= 'Tcell') 
celltype[celltype$ClusterID %in% c(0,7),2]='Cd4+Tn'
celltype[celltype$ClusterID %in% c(1),2]='Cd8+Tn'
celltype[celltype$ClusterID %in% c(2),2]='Treg' 
celltype[celltype$ClusterID %in% c(3),2]='B' 
celltype[celltype$ClusterID %in% c(4),2]='Tγ/δ' 
celltype[celltype$ClusterID %in% c(5),2]='Macrophage' 
celltype[celltype$ClusterID %in% c(6),2]='IL2C' 


Tcell@meta.data$celltype2 = "NA"
for(i in 1:nrow(celltype)){
  Tcell@meta.data[which(Tcell@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype2'] <- celltype$celltype[i]}

mycol <- c('#2DA248','#DBDD8D','#A3D7E5','#8AC997','#8D574C','#DAC9AC','#F09EA1')

p<- DimPlot(Tcell,group.by = "celltype2",raster=FALSE,label = T,cols = mycol)+ggtitle("T cell subtypes")
ggsave("05Tcell_celltype2.pdf",width = 7,height = 7,plot = p)
ggsave("05Tcell_celltype2.png",width = 7,height = 7,plot = p)

saveRDS(Tcell,"00Tcell_anno1.rds")

##rm macrophage,B 
Idents(Tcell) <- "celltype2"

mac <- subset(Tcell,idents = "Macrophage")
saveRDS(mac,"00mac_from_T.rds")

B <- subset(Tcell,idents = "B")
saveRDS(B,"00B_from_T.rds")

Tsub <- subset(Tcell,idents = c("Macrophage","B"),invert = T)

##re cluster -- final res

Tsub <- Tsub %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(Tsub)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.5))

col.cluster <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B', '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
                 '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764', '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
                 '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A', '#62615A','#B82129','#66762E')

p<- DimPlot(Tsub,group.by = "RNA_snn_res.0.5",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.5")
ggsave("02anno/02cluster05.pdf",width = 7,height = 7,plot = p)
ggsave("02anno/02cluster05.png",width = 7,height = 7,plot = p)

Tsub$celltype2 <- factor(Tsub$celltype2,levels = c("Cd8+Tn",
                                                   "Cd4+Tn",
                                                   "Treg",
                                                   "Tγ/δ",
                                                   "IL2C"))

mycol <- c('#DBDD8D','#A3D7E5','#8AC997','#DAC9AC','#F09EA1')

p<- DimPlot(Tsub,group.by = "celltype2",raster=FALSE,label = T,cols = mycol)+ggtitle("T cell subtypes")
ggsave("02anno/01Tcell_celltype2.pdf",width = 7,height = 7,plot = p)
ggsave("02anno/01Tcell_celltype2.png",width = 7,height = 7,plot = p)


p<- DimPlot(Tsub,group.by = "celltype2",raster=FALSE,label = T,cols = mycol,split.by = "group")+ggtitle("T cell subtypes")
ggsave("02anno/01Tcell_celltype2~group.pdf",width = 12,height = 7,plot = p)
ggsave("02anno/01Tcell_celltype2~group.png",width = 12,height = 7,plot = p)

saveRDS(Tsub,"00Tcell_anno_final.rds")


#find markers
Idents(Tsub) <- "celltype2"

markers <- FindAllMarkers(Tsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"02anno/03markers.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,"02anno/03markers_top5.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"02anno/03markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"02anno/03markers_top30.csv")

p_heat <- DoHeatmap(Tsub ,features = unique(top5$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("02anno/04Marker_heatmap.pdf",plot = p_heat,width = 12,height = 8)
ggsave("02anno/04Marker_heatmap.png",plot = p_heat,width = 12,height = 8)


p_heat <- DoHeatmap(Tsub ,features = unique(top10$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("02anno/04Marker_heatmap_top10.pdf",plot = p_heat,width = 10,height = 10)
ggsave("02anno/04Marker_heatmap_top10.png",plot = p_heat,width = 10,height = 10)

##############BP
sce.markers=top30
head(sce.markers)

ids <- suppressWarnings(bitr(sce.markers$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db'))
head(ids)

sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

sce.markers=sce.markers[sce.markers$cluster!="none",]
dim(sce.markers)
head(sce.markers)

gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)

table(sce.markers$cluster)

gcSample # entrez id , compareCluster 
names(gcSample)

xx.BP <- compareCluster(gcSample, fun="enrichGO",OrgDb="org.Mm.eg.db" ,
                        pvalueCutoff=1,readable=TRUE,
                        ont="BP") 

BP_res <- xx.BP@compareClusterResult
head(BP_res)
write.csv(BP_res,"06Tcells_top30Markers_BP.csv")

top3 <- BP_res %>%
  group_by(Cluster) %>%
  slice_min(n = 3, order_by = pvalue)
head(top3)

write.csv(top3,"06Tcells_top30Markers_BP_top3.csv")

top3$x <- c(1,2,3,5,6,7,9,10,11,13,14,15,17,18,19)
ggplot(top3, aes(x = x, y = -log10(rev(p.adjust)), fill = rev(Cluster))) + 
  geom_bar(stat = "identity") +
  coord_flip() +  
  theme_minimal() +  
  scale_fill_manual(values = mycol) +  
  labs(x = "", y = "-Log10(p.adjust)",fill = "Celltypes")+
  theme(
    text = element_text(size = 14, family = "serif"),
    axis.ticks.length = unit(0.2, "cm"),  
    axis.ticks = element_line(size = 1),  
    plot.title = element_text(size = 14, colour = "black", hjust = 0.5),
    axis.title.y = element_text(size = 14, color = "black", vjust = 1.9, hjust = 0.5, angle = 90),
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    axis.text.x = element_text(size = 14, color = "black", vjust = 0.5, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 14, color = "black", vjust = 0.5, hjust = 1, angle = 0),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 2),
    axis.line.x = element_line(colour = "black", size = 0),
    axis.line.y = element_line(colour = "black", size = 0)
  )+
  scale_x_continuous(breaks = top3$x,labels = rev(top3$Description),position = "top")

ggsave("06Tcells_top30Markers_BP_top3_bar.pdf",width = 8,height = 8)  
ggsave("06Tcells_top30Markers_BP_top3_bar.png",width = 8,height = 8)  


FeaturePlot(Tsub,features = c("Cxcr3","Tbx21"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("02anno/Th1.pdf",width = 12,height = 7)
ggsave("02anno/Th1.png",width = 12,height = 7)

#Cd4+ recluster
cd4 <- subset(Tsub,idents = c("Treg","Cd4+Tn"))
cd4$celltype2 <- Idents(cd4)

cd4 <- cd4 %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(cd4)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.3,0.5,0.6,0.7,0.8))

cd4 <- FindClusters(cd4,resolution = c(0.6,0.7))
dir.create("03CD4")


p<- DimPlot(cd4,group.by = "RNA_snn_res.0.5",raster=FALSE,label = T,cols = col.cluster)+ggtitle("Cd4+T RNA_snn_res.0.5")
ggsave("03CD4/01umap_res05.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/01umap_res05.png",width = 7,height = 7,plot = p)

p<- DimPlot(cd4,group.by = "RNA_snn_res.0.6",raster=FALSE,label = T,cols = col.cluster)+ggtitle("Cd4+T RNA_snn_res.0.6")
ggsave("03CD4/01umap_res06.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/01umap_res06.png",width = 7,height = 7,plot = p)

p<- DimPlot(cd4,group.by = "RNA_snn_res.0.6",split.by = "group",raster=FALSE,label = T,cols = col.cluster)+ggtitle("Cd4+T RNA_snn_res.0.6")
ggsave("03CD4/01umap_res06~group.pdf",width = 12,height = 7,plot = p)
ggsave("03CD4/01umap_res06~group.png",width = 12,height = 7,plot = p)


p<- DimPlot(cd4,group.by = "RNA_snn_res.0.7",raster=FALSE,label = T,cols = col.cluster)+ggtitle("Cd4+T RNA_snn_res.0.7")
ggsave("03CD4/01umap_res07.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/01umap_res07.png",width = 7,height = 7,plot = p)

p<- DimPlot(cd4,group.by = "RNA_snn_res.0.8",raster=FALSE,label = T,cols = col.cluster)+ggtitle("Cd4+T RNA_snn_res.0.8")
ggsave("03CD4/01umap_res08.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/01umap_res08.png",width = 7,height = 7,plot = p)


FeaturePlot(cd4,features = c("Cxcr3","Tbx21"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("03CD4/Th1.pdf",width = 12,height = 7)
ggsave("03CD4/Th1.png",width = 12,height = 7)

FeaturePlot(cd4,features = c("Cd4"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("03CD4/CD4.pdf",width = 7,height = 7)
ggsave("03CD4/CD4.png",width = 7,height = 7)

saveRDS(cd4,"03CD4/00cd4_01.rds")

Idents(cd4) <- "RNA_snn_res.0.7"
markers <- FindAllMarkers(cd4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"03CD4/02markers.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,"03CD4/02markers_top5.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"03CD4/02markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"03CD4/02markers_top30.csv")

p_heat <- DoHeatmap(cd4, ,features = unique(top5$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("03CD4/03Marker_heatmap.pdf",plot = p_heat,width = 12,height = 8)
ggsave("03CD4/03Marker_heatmap.png",plot = p_heat,width = 12,height = 8)

###########top10 markers
top10 <- split(top10,top10$cluster)

for (i in 1:length(top10)){
  DotPlot(cd4,features = top10[[i]]$gene)+RotatedAxis()
  ggsave(filename = paste0("03CD4/01marker/dotplot/dotplot_cluster_",i-1, ".png"),bg = "white")
  FeaturePlot(cd4,features = top10[[i]]$gene,raster = F,ncol = 4,cols = c("lightgrey","#ff0000"))
  ggsave(filename = paste0("03CD4/01marker/featureplot/featureplot_cluster_",i-1, ".png"),height = 9,width = 12,bg = "white")
  
}


library(clusterProfiler)
library(org.Mm.eg.db) 
library(ggplot2)

top30<- split(top30,top30$cluster)

for (i in 1:length(top30)){
  genes <- top30[[i]]$gene
  DEGs_entrez_id <- mapIds(x=org.Mm.eg.db,column = "ENTREZID",keys = genes,keytype = "SYMBOL")
  
  egoBP <- enrichGO(DEGs_entrez_id, OrgDb=org.Mm.eg.db, ont='BP',
                    pAdjustMethod='BH', pvalueCutoff=0.5,
                    qvalueCutoff=0.5, keyType='ENTREZID',readable = T)
  
  write.csv(egoBP@result,file = paste0("03CD4/cluster",i-1,"_BP.csv"),row.names = F)
  
  dotplot(egoBP,title= paste0('Biological process of cluster',i-1),showCategory=10) + scale_y_discrete(labels=function(enrich_go_BP) stringr::str_wrap(enrich_go_BP, width=60))
  ggsave(filename = paste0("03CD4/cluster",i-1,"_BP.png"),bg = "white",width = 8,height = 8)
  
}


FeaturePlot(cd4,features = c("Cxcr5","Bcl6"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("03CD4/Tfh.pdf",width = 12,height = 7)
ggsave("03CD4/Tfh.png",width = 12,height = 7)

##anno
celltype=data.frame(ClusterID=0:5,
                    celltype= 'CD4') 
celltype[celltype$ClusterID %in% c(0,1),2]='Cd4+Tn'
celltype[celltype$ClusterID %in% c(2),2]='Treg'
celltype[celltype$ClusterID %in% c(3),2]='Th1' 
celltype[celltype$ClusterID %in% c(4),2]='Tfh' 
celltype[celltype$ClusterID %in% c(5),2]='CD4-(rm)' 

cd4@meta.data$celltype3 = "cd4"
for(i in 1:nrow(celltype)){
  cd4@meta.data[which(cd4@meta.data$RNA_snn_res.0.7 == celltype$ClusterID[i]),'celltype3'] <- celltype$celltype[i]}

saveRDS(cd4,"03CD4/00cd4_anno2.rds")

p<- DimPlot(cd4,group.by = "celltype3",raster=FALSE,label = T,cols = mycol)+ggtitle("CD4+T subtypes")
ggsave("03CD4/04CD4_celltype3.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/04CD4_celltype3.png",width = 7,height = 7,plot = p)


p<- DimPlot(cd4,group.by = "celltype3",raster=FALSE,label = T,cols = mycol,split.by = "group")+ggtitle("CD4+T subtypes")
ggsave("03CD4/04Tcell_celltype3~group.pdf",width = 12,height = 7,plot = p)
ggsave("03CD4/04Tcell_celltype3~group.png",width = 12,height = 7,plot = p)

#rm CD4-, cd4 final recluster
Idents(cd4) <- "celltype3"
cd4sub <- subset(cd4,idents = "CD4-(rm)",invert = T)
cd4sub$celltype3 <- Idents(cd4sub)


cd4sub <- cd4sub %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(cd4sub)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.5))

saveRDS(cd4sub,"03CD4/00cd4_final.rds")

dir.create("03CD4/CD4_final")
p<- DimPlot(cd4sub,group.by = "celltype3",raster=FALSE,label = T,cols = mycol)+ggtitle("CD4+T subtypes")
ggsave("03CD4/CD4_final/04CD4_celltype3.pdf",width = 7,height = 7,plot = p)
ggsave("03CD4/CD4_final/04CD4_celltype3.png",width = 7,height = 7,plot = p)


p<- DimPlot(cd4sub,group.by = "celltype3",raster=FALSE,label = T,cols = mycol,split.by = "group")+ggtitle("CD4+T subtypes")
ggsave("03CD4/CD4_final/04Tcell_celltype3~group.pdf",width = 12,height = 7,plot = p)
ggsave("03CD4/CD4_final/04Tcell_celltype3~group.png",width = 12,height = 7,plot = p)

Idents(cd4sub) <- "celltype3"
markers <- FindAllMarkers(cd4sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"03CD4/CD4_final/02markers.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,"03CD4/CD4_final/02markers_top5.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"03CD4/CD4_final/02markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"03CD4/CD4_final/02markers_top30.csv")

p_heat <- DoHeatmap(cd4sub,features = unique(top5$gene),group.colors = mycol)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("03CD4/CD4_final/03Marker_heatmap.pdf",plot = p_heat,width = 12,height = 8)
ggsave("03CD4/CD4_final/03Marker_heatmap.png",plot = p_heat,width = 12,height = 8)

features = c("Cd4","Lef1","Foxp3","Tbx21","Cxcr5")

for (i in 1:length(features)) {
  FeaturePlot(cd4sub,features = features[i],raster = F,cols = c("lightgrey","#ff0000"))
  ggsave(paste0("03CD4/CD4_final/",features[i],".png"),width = 7,height = 7)
  ggsave(paste0("03CD4/CD4_final/",features[i],".pdf"),width = 7,height = 7)
}

#Treg subtype
Idents(cd4) <- "celltype3"
Treg <- subset(cd4,idents = "Treg")

Treg <- Treg %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(features = rownames(Treg)) %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% 
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = c(0.2,0.3,0.4,0.5))

Treg <- FindClusters(Treg,resolution = c(0.6,0.7,0.8))

print("RNA_snn_res.0.2:")
table(Treg$RNA_snn_res.0.2)

print("RNA_snn_res.0.3:")
table(Treg$RNA_snn_res.0.3)

print("RNA_snn_res.0.4:")
table(Treg$RNA_snn_res.0.4)

print("RNA_snn_res.0.5:")
table(Treg$RNA_snn_res.0.5)

dir.create("04Treg")

p<- DimPlot(Treg,group.by = "RNA_snn_res.0.5",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.5")
ggsave("04Treg/02cluster05.pdf",width = 7,height = 7,plot = p)
ggsave("04Treg/02cluster05.png",width = 7,height = 7,plot = p)

p<- DimPlot(Treg,group.by = "RNA_snn_res.0.5",split.by = "group",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.5")
ggsave("04Treg/02cluster05~group.pdf",width = 12,height = 7,plot = p)
ggsave("04Treg/02cluster05~group.png",width = 12,height = 7,plot = p)

saveRDS(Treg,"04Treg/00Treg02.rds")

p<- DimPlot(Treg,group.by = "RNA_snn_res.0.7",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.7")
ggsave("04Treg/02cluster07.pdf",width = 7,height = 7,plot = p)
ggsave("04Treg/02cluster07.png",width = 7,height = 7,plot = p)

p<- DimPlot(Treg,group.by = "RNA_snn_res.0.8",raster=FALSE,label = T,cols = col.cluster)+ggtitle("RNA_snn_res.0.8")
ggsave("04Treg/02cluster08.pdf",width = 7,height = 7,plot = p)
ggsave("04Treg/02cluster08.png",width = 7,height = 7,plot = p)

VlnPlot(Treg,features = "Il1rl1",group.by = "RNA_snn_res.0.8")
ggsave("04Treg/Il1rl1_vln.pdf",width = 7,height = 7)
ggsave("04Treg/Il1rl1_vln.png",width = 7,height = 7)

#find markers
Idents(Treg) <- "RNA_snn_res.0.8"

markers <- FindAllMarkers(Treg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"04Treg/03markers.csv")

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5,"04Treg/03markers_top5.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,"04Treg/03markers_top10.csv")

top30 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
write.csv(top30,"04Treg/03markers_top30.csv")

p_heat <- DoHeatmap(Treg,features = unique(top10$gene),group.colors = col.cluster)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave("04Treg/04Marker_heatmap.pdf",plot = p_heat,width = 6,height = 6)
ggsave("04Treg/04Marker_heatmap.png",plot = p_heat,width = 6,height = 6)

FeaturePlot(Treg,features = c("Il1rl1","Nt5e"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("04Treg/05St2.pdf",width = 12,height = 7)
ggsave("04Treg/05St2.png",width = 12,height = 7)

VlnPlot(Treg,features = c("Il1rl1","Nt5e","Cxcr3"),split.by = "group",ncol = 2)

ggsave("04Treg/05St2_vlnplot.pdf",width = 10,height = 10)
ggsave("04Treg/05St2_vlnplot.png",width = 10,height = 10)

FeaturePlot(Treg,features = c("Foxp3","Il1rl1","Nt5e","Il2ra"),raster = F,cols = c("lightgrey","#ff0000"))
ggsave("04Treg/Foxp3_Il2ra.pdf",width = 10,height = 10)
ggsave("04Treg/Foxp3_Il2ra.png",width = 10,height = 10)

FeaturePlot(Treg,features = c("Il1rl1"),split.by = "group",raster = F,cols = c("lightgrey","#ff0000"))
ggsave("04Treg/Il1rl1_group.pdf",width = 10,height = 10)
ggsave("04Treg/Il1rl1_group.png",width = 10,height = 10)

FeaturePlot(Treg,features = c("Nt5e"),split.by = "group",raster = F,cols = c("lightgrey","#ff0000"))
ggsave("04Treg/Nt5e_group.pdf",width = 10,height = 10)
ggsave("04Treg/Nt5e_group.png",width = 10,height = 10)

VlnPlot(Treg,features = c("Foxp3","Nt5e","Cxcr3"),split.by = "group",ncol = 2)

ggsave("04Treg/05St2_vlnplot.pdf",width = 10,height = 10)
ggsave("04Treg/05St2_vlnplot.png",width = 10,height = 10)


top10 <- split(top10,top10$cluster)

for (i in 1:length(top10)){
  DotPlot(Treg,features = top10[[i]]$gene)+RotatedAxis()
  ggsave(filename = paste0("04Treg/01marker/dotplot/dotplot_cluster_",i-1, ".png"),bg = "white")
  FeaturePlot(Treg,features = top10[[i]]$gene,raster = F,ncol = 4,cols = c("lightgrey","#ff0000"))
  ggsave(filename = paste0("04Treg/01marker/featureplot/featureplot_cluster_",i-1, ".png"),height = 9,width = 12,bg = "white")
  
}
