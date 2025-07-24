# sce=readRDS('results/01.cell_annotation/sce.rds')
# my.data=readRDS('results/03.B_cell_subcluster/my.data.rds')
# table(sce$cell_type)
# sce$group=sce$cell_type
# sce$group[colnames(my.data)]=my.data$subcluster
# table(sce$group)
# colnames(sce@meta.data)
# my.data1=subset(sce,TumorNormal=='Tumor' & group!='B cells')
# table(my.data1$group)
# saveRDS(my.data1,file = 'results/my.data.cellchat.rds')

dir.create('results/05.cellchat')
my.data1=readRDS('results/05.cellchat/my.data.cellchat.rds')
library(CellChat)
table(my.data1$group)
DefaultAssay(my.data1) <- 'RNA'
cellchat <- createCellChat(object = my.data1, meta = my.data1@meta.data, group.by = "group")
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat@data.signaling
#
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)

library(NMF)
# selectK(cellchat, pattern = "outgoing")
# selectK(cellchat, pattern = "incoming")

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
dev.off()
saveRDS(cellchat,file='results/05.cellchat/cellchat.RData')

cellchat=readRDS('results/05.cellchat/cellchat.RData')
groupSize <- as.numeric(table(cellchat@idents))

cellchat@netP$pathways
length(cellchat@netP$pathways)
#19

par(mfrow=c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")
dev.off()

table(my.data1$group)
pdf('results/05.cellchat/cellchat_all_pathway.pdf',height = 9,width = 7)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,sources.use = c('subcluster 3'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'))

netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,
                 sources.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'),
                 targets.use = c('subcluster 3'))
dev.off()
cellchat@netP$pathways

##IL16-CD4
pdf('results/05.cellchat/cellchat_IL16_pathway.pdf',height = 3,width = 5,onefile = F)
netVisual_bubble(cellchat, signaling = c('IL16','IL2','IL1'),sources.use = c('subcluster 3'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'))
dev.off()

pdf('results/05.cellchat/cellchat_LT_pathway.pdf',height = 3,width = 5,onefile = F)
netVisual_bubble(cellchat, signaling = c('LT'),sources.use = c('subcluster 3'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'))
dev.off()

pdf('results/05.cellchat/cellchat_GALECTIN_pathway.pdf',height = 3,width = 5,onefile = F)
netVisual_bubble(cellchat, signaling = c("GALECTIN" ),sources.use = c('subcluster 3'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'))
dev.off()

pdf('results/05.cellchat/cellchat_MIF_pathway.pdf',height = 3,width = 5,onefile = F)
netVisual_bubble(cellchat, signaling = c("MIF" ),sources.use = c('subcluster 3'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells'))
dev.off()


#1,11
pdf('results/05.cellchat/MIF_CD8+T.pdf',height = 3,width = 6,onefile = F)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathway[1],
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells')[1])
dev.off()
pdf('results/05.cellchat/LT_CD8+T.pdf',height = 3,width = 6,onefile = F)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathway[11],
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells')[1])
dev.off()

pdf('results/05.cellchat/MIF_NK.pdf',height = 3,width = 6,onefile = F)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathway[1],
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells')[3])
dev.off()

pdf('results/05.cellchat/LT_NK.pdf',height = 3,width = 6,onefile = F)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathway[11],
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells')[3])
dev.off()

#8
pdf('results/05.cellchat/fig5d.pdf',height = 6,width = 4,onefile = F)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathway,
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6')[2],
                 targets.use = c('CD8+ T cells','NK cells'))
dev.off()


netVisual_bubble(cellchat, signaling = cellchat@netP$pathway[11],
                 sources.use = c('subcluster 1','subcluster 3','subcluster 4','subcluster 6'),
                 targets.use = c('CD8+ T cells','Naive T cells','NK cells','Tregs','Myeloid cells')[4])
