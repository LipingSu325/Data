library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source

library(Seurat)
# library(SeuratObject)
library(Matrix)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)


dir.create('results/01.cell_annotation')
data.sample=read.delim('origin_datas/GEO/GSE182434_cell_annotation.txt.gz')
head(data.sample)
colnames(data.sample)[3]='Samples'
table(data.sample$Tissue)
table(data.sample$TumorNormal)
data.sample=data.sample[data.sample$Tissue!='FL',]
table(data.sample$TumorNormal,data.sample$Tissue)

data.count=fread('origin_datas/GEO/GSE182434_raw_count_matrix.txt.gz')
dim(data.count)
data.count=as.data.frame(data.count)
rownames(data.count)=data.count[,1]
length(rownames(data.count))
length(unique(rownames(data.count)))
data.count=data.count[,-1]
data.count[1:5,1:5]


data.count=data.count[,data.sample$ID]
dim(data.count)

library(Seurat)
datalist <- list()
for (i in unique(data.sample$Samples)) {
  ce <- data.sample$ID[data.sample$Samples == i]
  my.count <- data.count[, ce]
  datalist[[i]] <- CreateSeuratObject(counts=my.count, project = i,min.cells = 3,min.features = 200)
  datalist[[i]]$Samples=data.sample$Samples[data.sample$Samples == i]
  datalist[[i]]$Patient=data.sample$Patient[data.sample$Samples == i]
  datalist[[i]]$CellType=data.sample$CellType[data.sample$Samples == i]
  datalist[[i]]$TumorNormal=data.sample$TumorNormal[data.sample$Samples == i]
  datalist[[i]]$COO =data.sample$COO [data.sample$Samples == i]
  rm(my.count)
}


for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  rm(sce)
}

sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  14285
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<6000 & percent.mt<10)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))

sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples", lambda=0.5, max.iter.harmony=30)
ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
dev.off()
sce <- RunUMAP(sce, dims=1:20, reduction="harmony")
after_batch=DimPlot(sce,group.by='Samples',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
after_batch


library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")
sce <- FindClusters(object = sce,resolution = .1)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
sce=subset(sce, subset=seurat_clusters%in% c(0:7,9))
my.cols=brewer.pal(11,"Set3")[-9]

seurat_clusters_umap=DimPlot(sce,group.by='seurat_clusters',reduction="umap",label = T,pt.size = 0.2,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
seurat_clusters_umap+
  
  DimPlot(sce,group.by='CellType',reduction="umap",label = T,pt.size = 0.2,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())

table(sce$seurat_clusters,sce$CellType)

Logfc = 0.25
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
Idents(sce)<-'seurat_clusters'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.cell_annotation/seurat_clusters_degs.csv')
# 
# 
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top5$gene)
length(unique(Top5$gene))
DotPlot(object = sce, features = unique(Top5$gene),
        cols=c("snow", "blue"),scale = T,col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()
table(sce$seurat_clusters,sce$CellType)
#


marker <- data.frame(cluster = 0:9,cell = 0:9)
marker[marker$cluster %in% c(0),2] <- 'Naive T cells'
marker[marker$cluster %in% c(1),2] <- 'CD8+ T cells'
marker[marker$cluster %in% c(2,4,7,9),2] <- 'B cells'
marker[marker$cluster %in% c(3),2] <- 'Tregs'
marker[marker$cluster %in% c(5),2] <- 'Myeloid cells'
marker[marker$cluster %in% c(6),2] <- 'NK cells'
marker[marker$cluster %in% c(8),2] <- ' '
marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols=brewer.pal(11,"Set3")[-9]
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="umap",label = F,pt.size = 0.2,cols =my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap=LabelClusters(cell_type_umap,id = 'cell_type',family='Times')
cell_type_umap
ggsave('results/01.cell_annotation/cell_type_umap.pdf',cell_type_umap,height = 5,width = 5.5)
table(sce$cell_type)

# 0：CD4+T    	IL7R,LTB
# 1：CD8+T		GZMK，CST7，CD8A，IFNG
# 3: Tregs		FOXP3,CTLA4,TIGIT,ICOS
# 2,7: B cell		SPIB,CD79A,MS4A1,JCHAIN
# 4:Myeloid cell	CST3,LYZ,C1QB
# 5: NK cell		XCL1,GNLY,XCL2,KLRD1

marker_gene=c('CD79A','MS4A1',
              'CST3','LYZ','C1QB',
              'GZMK','CST7','CD8A','IFNG',
              'FOXP3','TNFRSF4','ICOS','CTLA4',
              'NKG7','XCL1','GNLY','KLRC1' ,
              'LTB','TCF7','CCR7'
              
)
Idents(sce)='cell_type'
marker.dot=DotPlot(object = sce, features = marker_gene,cols=c("#80B1D3", "#D95F02"),scale = T,col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')#+coord_flip()
marker.dot
ggsave('results/01.cell_annotation/marker.dot.pdf',marker.dot,height = 5,width = 8)


saveRDS(sce,file = 'results/01.cell_annotation/sce.rds')

table(sce$Samples,sce$TumorNormal)
cell_freq1=data.frame(t(prop.table(table(sce$cell_type,sce$Patient),margin=2)),
                      tissue=rep(c('Tumor','Tonsil'),c(4,1)))
cell_freq1
colnames(cell_freq1)[1:3]<-c('Patient','cell_type','Freq')
cell_freq1$patinet=substr(cell_freq1$Patient,1,8)
cell_prop_fig1=ggplot(cell_freq1,aes(x=patinet,y=Freq,fill=cell_type))+
  scale_fill_manual(values = my.cols)+
  facet_grid(~tissue,scales = 'free',space='free')+
  geom_bar(position = "fill",stat="identity")+
  xlab('')+ylab('Proportion')+theme_bw()+
  theme(text = element_text(family = 'Times',size=12),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig1
ggsave('results/01.cell_annotation/cell_prop_fig1.pdf',cell_prop_fig1,height = 4,width = 6)


meta.data=sce@meta.data
head(meta.data)
# meta.data=meta.data[meta.data$TumorNormal=='Tumor',]
bar = meta.data %>% group_by(TumorNormal, cell_type) %>% count()
Type_label = c('Tumor','Normal')
bar$TumorNormal = factor(bar$TumorNormal, levels=Type_label)
bar = bar %>% group_by(TumorNormal) %>% mutate(percent=100*n/sum(n))

cell_prop_fig2=ggplot(data=bar, aes(x=cell_type, y=percent, fill=TumorNormal,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FB8072","#BEBADA"))+theme_classic()+
  ggtitle("Percent(%)")+xlab('')+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"),
        legend.text=element_text(family = 'Times', size=12), 
        legend.title=element_blank(), text = element_text(family = 'Times'))
ggsave('results/01.cell_annotation/cell_prop_fig2.pdf',cell_prop_fig2,height = 4,width = 7)
