library(Seurat)
setwd

dir.create('results/04.SCENIC')
library(SCENIC)
library(RcisTarget)
# set.seed(123)
# subcell <- sample(colnames(sce),1000)
# colnames(sce)
# sce_tf <- sce[,subcell]

cellInfo <- data.frame(sce_tf@meta.data)
head(cellInfo)
saveRDS(cellInfo, file="results/04.SCENIC/cellInfo.Rds")

# dir.create("cisTarget_databases")
scenicOptions <- initializeScenic(org="hgnc",
                                  dbDir=,
                                  nCores=1)
saveRDS(scenicOptions, file="results/04.SCENIC/scenicOptions.Rds")

exprMat <- as.matrix(sce_tf@assays$RNA@counts)
genesKept <- geneFiltering(exprMat, scenicOptions,
                           minCountsPerGene = 10 * 0.01 * ncol(exprMat),
                           minSamples = ncol(exprMat) * 0.1)
dim(exprMat)
exprMat[1:4,1:4]
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
range(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered_log, scenicOptions, nParts=20)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
dim(exprMat_log)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)

AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
AUCmatrix[1:5,1:5]
selected_columns <- !grepl("extended", colnames(AUCmatrix))
AUCmatrix <- AUCmatrix[,selected_columns]
dim(AUCmatrix)
AUCmatrix[1:5,1:5]
# scRNAauc <- AddMetaData(sce_tf, AUCmatrix)
# scRNAauc@assays$integrated <- NULL
# sce_sub <- AddMetaData(sce_sub, AUCmatrix)

cellInfo=readRDS('results/04.SCENIC/cellInfo.Rds')
head(cellInfo)
table(cellInfo$Type)




diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=names(table(group))
  dat1=dat[dat$cluster==gr[1],-1]
  dat2=dat[dat$cluster==gr[2],-1]
  pathway=unique(colnames(dat)[-1])
  p_vale=data.frame()
  for (i in pathway){
    dd1=wilcox.test(dat1[,i],dat2[,i])$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
aucell.diff<-diff_pathway(dat=t(AUCmatrix[rownames(cell.anno),]),group=cell.anno$Type)
head(aucell.diff)
table(aucell.diff$p.value<0.05)
aucell.diff

my_mutiboxplot(dat =AUCmatrix[rownames(cell.anno),],notch = F,
               group = cell.anno$Type,group_cols = pal_jama()(7)[c(3,2)],
               ylab = 'AUCell')+theme_classic()+
  theme(text = element_text(family = 'Times'),axis.text.x = element_text(angle = 30,hjust = 1))
ggsave('results/04.SCENIC/TF_bocplot.pdf',height = 5,width = 7)

tf_regulon=readRDS("int/3.1_regulons_forAUCell.Rds")
library(clusterProfiler)
library(org.Hs.eg.db)
SYMBOL_id=tf_regulon$`XBP1 (135g)`
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
XBP1.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",ont = "BP",
                            pvalueCutoff = 0.05)
XBP1.erich.go.BP.res=XBP1.erich.go.BP@result
XBP1.erich.go.BP.res=XBP1.erich.go.BP.res[XBP1.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
XBP1.erich.go.BP.res=XBP1.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
p2=ggplot(data=XBP1.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('Enrichment analysis of XBP1 Target genes')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
dev.off()
ggsave('results/04.SCENIC/TF_XBP1_enrichment.pdf',p2,height = 4.5,width = 6)

SYMBOL_id=tf_regulon$`SPIB (344g)`
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
SPIB.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",ont = "BP",
                             pvalueCutoff = 0.05)
SPIB.erich.go.BP.res=SPIB.erich.go.BP@result
SPIB.erich.go.BP.res=SPIB.erich.go.BP.res[SPIB.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
SPIB.erich.go.BP.res=SPIB.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
p3=ggplot(data=SPIB.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('Enrichment analysis of SPIB Target genes')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
dev.off()
ggsave('results/04.SCENIC/TF_SPIB_enrichment.pdf',p3,height = 4.5,width = 6)


SYMBOL_id=tf_regulon$`RELB (281g)`
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
RELB.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",ont = "BP",
                            pvalueCutoff = 0.05)
RELB.erich.go.BP.res=RELB.erich.go.BP@result
RELB.erich.go.BP.res=RELB.erich.go.BP.res[RELB.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
RELB.erich.go.BP.res=RELB.erich.go.BP.res %>% slice_min(n =8, order_by = p.adjust)
p4=ggplot(data=RELB.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('Enrichment analysis of RELB Target genes')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
p4
ggsave('results/04.SCENIC/TF_RELB_enrichment.pdf',p4,height = 4.5,width = 6)



# SYMBOL_id=tf_regulon$`JUN (41g)`
# length(SYMBOL_id)
# entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
#                    keytype = "SYMBOL",column = "ENTREZID")
# entrez_id = na.omit(entrez_id)
# JUN.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
#                            keyType = "ENTREZID",ont = "BP",
#                            pvalueCutoff = 0.05)
# JUN.erich.go.BP.res=JUN.erich.go.BP@result
# JUN.erich.go.BP.res=JUN.erich.go.BP.res[JUN.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
# JUN.erich.go.BP.res=JUN.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
# # pdf('02.B_cell_inferCNV/chr6_enrichment.pdf',height = 6,width = 8,onefile = F)
# ggplot(data=JUN.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
#   geom_point(aes(size=Count),show.legend = T) +
#   scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
#   geom_segment(aes(yend=Description,xend=0),size=.5) +
#   labs(x='Count',y='')+theme_bw()+
#   theme(text = element_text(family = 'Times',size = 14))+
#   scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
# dev.off()


##interferon regulatory transcription factor (IRF) family.
SYMBOL_id=unique(c(tf_regulon$`IRF1 (32g)`,tf_regulon$`IRF7 (30g)`,tf_regulon$`IRF8 (12g)`))
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
IRF.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",ont = "BP",
                           pvalueCutoff = 0.05)
IRF.erich.go.BP.res=IRF.erich.go.BP@result
IRF.erich.go.BP.res=IRF.erich.go.BP.res[IRF.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
IRF.erich.go.BP.res=IRF.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
p5=ggplot(data=IRF.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+ggtitle('Enrichment analysis of IRF Target genes')+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
dev.off()
ggsave('results/04.SCENIC/TF_IRF_enrichment.pdf',p5,height = 4.5,width = 6)

pdf('results/04.SCENIC/Fig4.pdf',height = 15,width = 12,onefile = F)
mg_merge_plot(p1,mg_merge_plot(p2,p3,p4,p5,nrow=2,ncol=2,labels = LETTERS[2:5]),
              nrow=2,labels = c('A',''),heights = c(1,2))
dev.off()
