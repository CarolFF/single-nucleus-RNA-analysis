######cell-cell communication
#CellChat-- cell-cell communication
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)

sn_merge <- readRDS("00merge_anno2.rds")

Idents(sn_merge) <- "group"

CKO <- subset(sn_merge,idents="GPSM1 ΔTreg")
FLOX <- subset(sn_merge,idents="GPSM1 f/f")

###############FLOX
FLOX.input  <- FLOX@assays$RNA@data
identity  <-  data.frame(group =FLOX$celltype2, row.names = names(FLOX$celltype2)) # create a dataframe consisting of the cell labels

FLOX <- createCellChat(object = FLOX.input)

FLOX <- addMeta(FLOX, meta = identity, meta.name = "labels")
FLOX <- setIdent(FLOX, ident.use = "labels") # set "labels" as default cell identity


groupSize <- as.numeric(table(FLOX@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse

# set the used database in the object
FLOX@DB <- CellChatDB

FLOX <- subsetData(FLOX) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel  
FLOX <- identifyOverExpressedGenes(FLOX)
FLOX <- identifyOverExpressedInteractions(FLOX)
FLOX <- projectData(FLOX, PPI.mouse)  

FLOX <- computeCommunProb(FLOX, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
FLOX <- filterCommunication(FLOX, min.cells = 10)

FLOX <- computeCommunProbPathway(FLOX)
FLOX <- aggregateNet(FLOX)

saveRDS(FLOX,"10cellchat/FLOX_cc.rds")

groupSize1 <- as.numeric(table(FLOX@idents))
mat <- FLOX@net$weight
pdf("10cellchat/01FLOX/001FLOX.pdf",width = 12,height = 12)
par(mfrow = c(6,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize1, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show.all <- FLOX@netP$pathways

dir.create("10cellchat/01FLOX/FLOX_all_pathways_com_circle") 
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual 
  netVisual(FLOX, signaling = pathways.show.all[i], out.format = c("pdf"),layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall sig
  gg <- netAnalysis_contribution(FLOX, signaling = pathways.show.all[i])
  ggsave(filename=paste0("10cellchat/01FLOX/FLOX_all_pathways_com_circle/",pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg,width = 8, height = 6)
}

###############CKO
CKO.input  <- CKO@assays$RNA@data
identity  <-  data.frame(group =CKO$celltype2, row.names = names(CKO$celltype2)) # create a dataframe consisting of the cell labels

CKO <- createCellChat(object = CKO.input)

CKO <- addMeta(CKO, meta = identity, meta.name = "labels")
CKO <- setIdent(CKO, ident.use = "labels") # set "labels" as default cell identity


groupSize <- as.numeric(table(CKO@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse

# set the used database in the object
CKO@DB <- CellChatDB

CKO <- subsetData(CKO) # subset the expression data of signaling genes for saving computation cost
CKO <- identifyOverExpressedGenes(CKO)
CKO <- identifyOverExpressedInteractions(CKO)
CKO <- projectData(CKO, PPI.mouse)  

CKO <- computeCommunProb(CKO, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
CKO <- filterCommunication(CKO, min.cells = 10)

CKO <- computeCommunProbPathway(CKO)
CKO <- aggregateNet(CKO)

saveRDS(CKO,"10cellchat/CKO_cc.rds")

groupSize1 <- as.numeric(table(FLOX@idents))
mat <- CKO@net$weight
pdf("10cellchat/02CKO/01CKO.pdf",width = 12,height = 12)
par(mfrow = c(6,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize1, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show.all <- CKO@netP$pathways

dir.create("10cellchat/02CKO/CKO_all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual 
  netVisual(CKO, signaling = pathways.show.all[i], out.format = c("pdf"),layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall sig
  gg <- netAnalysis_contribution(CKO, signaling = pathways.show.all[i])
  ggsave(filename=paste0("10cellchat/02CKO/CKO_all_pathways_com_circle/",pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg,width = 8, height = 6)
}

##############FLOX&CKO_cc
object.list <- list(FLOX = FLOX, CKO = CKO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat,"10cellchat/FLOX&CKO_cc.rds")



#####################visualize
FLOX <- readRDS("10cellchat/FLOX_cc.rds")
CKO <- readRDS("10cellchat/CKO_cc.rds")
cellchat <- readRDS("10cellchat/FLOX&CKO_cc.rds")
cellchat <- updateCellChat(cellchat)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave("10cellchat/03compare/01interaction_bar.png",width = 8,height = 6)
ggsave("10cellchat/03compare/01interaction_bar.pdf",width = 8,height = 6)

pdf("10cellchat/03compare/02interaction_circle.pdf",width = 8,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf("10cellchat/03compare/02interaction_circle_number.pdf",width = 8,height = 6)
netVisual_diffInteraction(cellchat, weight.scale = T,title.name = "Differential number of interactions")
dev.off()

pdf("10cellchat/03compare/02interaction_circle_weight.pdf",width = 8,height = 6)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",title.name = "Differential interaction strength")
dev.off()

pdf("10cellchat/03compare/03interaction_heatmap.pdf",width = 8,height = 6)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

pdf("10cellchat/03compare/04interaction_circle.pdf",width = 8,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

png("10cellchat/03compare/04interaction_circle.png")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


rankSimilarity(CKO, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#068085','#E5481A'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c('#068085','#E5481A'))
gg1 + gg2
ggsave("10cellchat/03compare/05ranket.pdf",width = 8,height = 10)
ggsave("10cellchat/03compare/05ranket.png",width = 8,height = 10)

levels(FLOX@idents)

#Identify dysfunctional signaling by comparing the communication probabities
###Treg to others
gg1 <- netVisual_bubble(cellchat, sources.use = 27, targets.use = c(2,5,11,21,23),  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in CKO", 
                        angle.x = 45, remove.isolate = T,color.text.use = T,color.text = c('#068085','#E5481A'))
gg1
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 27, targets.use = c(2,5,11,21,23),  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in CKO", 
                        angle.x = 45, remove.isolate = T,color.text.use = T,color.text = c('#068085','#E5481A'))
#> Comparing communications on a merged object
gg1 + gg2
ggsave("10cellchat/03compare/06Treg~others.png",width = 8,height = 6)
ggsave("10cellchat/03compare/06Treg~others.pdf",width = 8,height = 6)

write.csv(gg1$data,"10cellchat/03compare/06Treg~others_compare_byprop_up.csv")
write.csv(gg2$data,"10cellchat/03compare/06Treg~others_compare_byprop_down.csv")


###others to Treg
gg3 <- netVisual_bubble(cellchat, sources.use = c(2,5,11,21,23), targets.use = 27,  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in CKO", 
                        angle.x = 45, remove.isolate = T,color.text.use = T,color.text = c('#068085','#E5481A'))
gg3
#> Comparing communications on a merged object
gg4 <- netVisual_bubble(cellchat, sources.use = c(2,5,11,21,23), targets.use = 27,  comparison = c(1, 2), 
                        max.dataset = 1, title.name = "Decreased signaling in CKO", 
                        angle.x = 45, remove.isolate = T,color.text.use = T,color.text = c('#068085','#E5481A'))
#> Comparing communications on a merged object
gg3 + gg4
ggsave("10cellchat/03compare/07others~Treg.png",width = 8,height = 8)
ggsave("10cellchat/03compare/07others~Treg.pdf",width = 8,height = 8)

write.csv(gg3$data,"10cellchat/03compare/07others_to_Treg_compare_byprop_up.csv")
write.csv(gg4$data,"10cellchat/03compare/07others_to_Treg_compare_byprop_down.csv")



























