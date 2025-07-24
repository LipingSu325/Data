dir.create('results/03.B_cell_subcluster')
table(sce$cell_type)
my.data=subset(sce,cell_type=='B cells')
table(my.data$TumorNormal)
# Normal  Tumor 
# 1505   3630 
my.data <- NormalizeData(my.data)
my.data <- FindVariableFeatures(my.data, selection.method = "vst", nfeatures = 2000)
my.data <- ScaleData(my.data, features = rownames(my.data))
my.data <- RunPCA(my.data, features = VariableFeatures(my.data))
my.data = RunHarmony(my.data, group.by.vars="Samples", lambda=0.5, max.iter.harmony=30)
ElbowPlot(my.data,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
dev.off()
my.data <- RunUMAP(my.data, dims=1:15, reduction="harmony")

my.data <- FindNeighbors(my.data, dims = 1:15, reduction="harmony")
my.data <- FindClusters(object = my.data,resolution = .2)

DefaultAssay(my.data) <- "RNA"
colnames(my.data@meta.data)
table(my.data@meta.data$seurat_clusters)
my.data=subset(my.data,seurat_clusters%in%c(0,1,2,3,4,6))
p=DimPlot(my.data,group.by='seurat_clusters',reduction="umap",label = T,pt.size = 0.4)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
ggsave('results/03.B_cell_subcluster/subcluster_umap.pdf',p,height = 5,width = 5.5)

marker <- data.frame(cluster = 0:6,cell = 0:6)
marker[marker$cluster %in% c(0),2] <- 'subcluster 1'
marker[marker$cluster %in% c(1),2] <- 'subcluster 2'
marker[marker$cluster %in% c(2),2] <- 'subcluster 3'
marker[marker$cluster %in% c(3),2] <- 'subcluster 4'
marker[marker$cluster %in% c(4),2] <- 'subcluster 5'
marker[marker$cluster %in% c(5),2] <- ''
marker[marker$cluster %in% c(6),2] <- 'subcluster 6'
marker
my.data@meta.data$subcluster <- sapply(my.data@meta.data$seurat_clusters,function(x){marker[x,2]})
my.cols2=brewer.pal(9,"Set2")
cell_type_umap2=DimPlot(my.data,group.by='subcluster',reduction="umap",label = F,pt.size = 0.2,cols =my.cols2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap2=LabelClusters(cell_type_umap2,id = 'subcluster',family='Times')
cell_type_umap2
pdf('results/03.B_cell_subcluster/subcluster_umap2.pdf',height = 7,width = 7.5)
cell_type_umap2
dev.off()

bar = my.data@meta.data %>% group_by(TumorNormal, subcluster) %>% count()
Type_label = c('Tumor','Normal')
bar$TumorNormal = factor(bar$TumorNormal, levels=Type_label)
bar = bar %>% group_by(TumorNormal) %>% mutate(percent=100*n/sum(n))

pdf('results/03.B_cell_subcluster/subcluster_barplot2.pdf',height = 7,width = 8.5)
ggplot(data=bar, aes(x=subcluster, y=percent, fill=TumorNormal,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FB8072","#BEBADA"))+theme_classic()+
  ggtitle("Percent(%)")+xlab('')+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"),
        legend.text=element_text(family = 'Times', size=12), 
        legend.title=element_blank(), text = element_text(family = 'Times'))
dev.off()

Idents(my.data)<-'subcluster'
Bcell_DEGs <- FindAllMarkers(object = my.data,logfc.threshold = 0.25, min.pct = 0.25,only.pos = T)
head(Bcell_DEGs)
Bcell_DEGs["pct.diff"]=Bcell_DEGs$pct.1-Bcell_DEGs$pct.2
Bcell_DEGs <- Bcell_DEGs[Bcell_DEGs$p_val_adj<0.05,]
table(Bcell_DEGs$cluster)
length(unique(Bcell_DEGs$gene))#762
head(Bcell_DEGs)
table(Bcell_DEGs$cluster)
write.xlsx(Bcell_DEGs,'results/03.B_cell_subcluster/Bcell_subcluster_DEGs.xlsx',overwrite = T)
Top5 <- Bcell_DEGs %>% group_by(cluster) %>% slice_min(n =5, order_by = pct.2)
length(Top5$gene)
length(unique(Top5$gene))
# DotPlot(object = my.data, features = unique(Top5$gene),cols=c("snow", "blue"),scale = T)+
#   RotatedAxis()+ ggtitle("Marker Genes")+
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlab('')+ylab('')+coord_flip()

DotPlot(object = my.data, features = unique(Top5$gene),cols=c("snow", "blue"),scale = T,group.by = 'subcluster')+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()

VlnPlot(my.data, features = c('CR2'),pt.size = 0,group.by = 'subcluster') #B cell 0
VlnPlot(my.data, features = c("TOP2A",'MKI67'),pt.size = 0,group.by = 'subcluster')#Proliferative B cells 2
VlnPlot(my.data, features = c('LMO2'),pt.size = 0,group.by = 'subcluster')#Germinal center B cell	3
VlnPlot(my.data, features = c('SDC1'),pt.size = 0,group.by = 'subcluster')#Plasma cell	6


p=VlnPlot(my.data, features = c('CR2',"TOP2A",'MKI67','LMO2','SDC1'),pt.size = 0,group.by = 'subcluster',ncol = 5)#Plasma cell	6
ggsave('results/03.B_cell_subcluster/marker_vlnplot2.pdf',height = 5,width = 18)

Bcell_DEGs.list=split(x=Bcell_DEGs,f=Bcell_DEGs$cluster)
Bcell_DEGs.list=sapply(Bcell_DEGs.list, function(x){subset(x,select='gene',drop=TRUE)})
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
subcluster.name=names(Bcell_DEGs.list)
subcluster.name=subcluster.name[order(subcluster.name)]

p=list()
for (s in subcluster.name) {
  entrez_id = mapIds(x = org.Hs.eg.db, keys = Bcell_DEGs.list[[s]],keytype = "SYMBOL",column = "ENTREZID")
  entrez_id = na.omit(entrez_id)
  c0.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05)
  c0.erich.go.BP.res=c0.erich.go.BP@result
  write.xlsx(c0.erich.go.BP.res,paste0('results/03.B_cell_subcluster/',s,'_enrichment.xlsx'))
  c0.erich.go.BP.res=c0.erich.go.BP.res[c0.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
  c0.erich.go.BP.res=c0.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust) 
  p[[which(subcluster.name==s)]]=ggplot(data = c0.erich.go.BP.res, 
                                        mapping = aes(x=Count,y=reorder(Description,Count), fill = -log10(p.adjust))) +
    geom_bar(stat="identity")+ theme_bw()+
    scale_fill_gradient(low="grey80",high =c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")[which(subcluster.name==s)])+
    labs(y=NULL,title = "subcluster ")+ggtitle(s)+
    theme(text = element_text(family = 'Times',size = 14))+
    scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
  
}


subcluster.enrichment.fig=mg_merge_plot(p,ncol=3,nrow = 2)
pdf('results/03.B_cell_subcluster/subcluster_enrichment.pdf',height = 10,width = 18,onefile = F)
subcluster.enrichment.fig
dev.off()


table(my.data$subcluster,my.data$TumorNormal)
table(my.data$subcluster,my.data$COO)
table(my.data$subcluster,my.data$Patient)
