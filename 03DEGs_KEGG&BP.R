#############GO&KEGG
rm(list = ls())
gc()

library(clusterProfiler)
library(org.Mm.eg.db) 
library(ggplot2)

##all celltypes

setwd("d:/03DEGs/")

files <- list.files("./")

for (i in 1:length(files)) {
  DEGs <- read.csv(files[i])
  celltypes <- stringr::str_split(files[i], "05DEG_|\\.", simplify = TRUE)[,2]
  DEG_sig <- DEGs[abs(DEGs$avg_log2FC) >= 0.5 & DEGs$p_val_adj < 0.05,]
  DEG_sig$trend <- "up"
  DEG_sig[DEG_sig$avg_log2FC < 0,]$trend <- "down"
  group <- data.frame(gene=DEG_sig$X,change=DEG_sig$trend)
  Gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Mm.eg.db")
  #构建文件并分析
  df  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
  gcSample=split(df$ENTREZID, df$change)
  
  data_GO <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Mm.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable =T
  )
  
  dotplot(data_GO,showCategory=10)+
    xlab(paste0(celltypes))+
    scale_y_discrete(labels=function(compareCluster) stringr::str_wrap(compareCluster, width=50))
  ggsave(paste0("BP/05DEGs_",celltypes,"_BP.png"),height = 8,width = 8)
  ggsave(paste0("BP/05DEGs_",celltypes,"_BP.pdf"),height = 8,width = 8)
  write.csv(data_GO@compareClusterResult,paste0("BP/05DEGs_",celltypes,"_BP.csv"))
  
  data_kegg<-clusterProfiler::compareCluster(gcSample,fun = "enrichKEGG",  
                                      keyType = 'kegg',  
                                      organism='mmu',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1 )
  data_kegg <- setReadable(data_kegg,OrgDb = org.Mm.eg.db, keyType="ENTREZID" )
  dotplot(data_kegg,showCategory=10)+xlab(paste0(celltypes))+
    scale_y_discrete(labels=function(compareCluster) stringr::str_wrap(compareCluster, width=50))
  ggsave(paste0("05DEGs_",celltypes,"_KEGG.png"),height = 8,width = 8)
  ggsave(paste0("05DEGs_",celltypes,"_KEGG.pdf"),height = 8,width = 8)
  write.csv(data_kegg@compareClusterResult,paste0("05DEGs_",celltypes,"_KEGG.csv"))
  
}

