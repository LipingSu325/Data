dir.create('results/02.B_cell_inferCNV')
#infercnv#######
setwd
library(Seurat)
my.data=readRDS
counts_matrix<-my.data@assays$RNA@counts
counts_matrix=data.frame(as.matrix(counts_matrix))
counts_matrix[1:5,1:5]

gene_pos = data.table::fread("infercnv/infercnv_genes.txt", data.table = F)

cell_anno=my.data@meta.data
cell_anno=data.frame(cell=rownames(cell_anno),type=cell_anno$TumorNormal)
head(cell_anno)


write.table(cell_anno,'infercnv/cell_anno.txt',quote = F,sep = '\t',row.names = F,col.names = F)
write.table(gene_pos,'infercnv/gene_order_file.txt',quote = F,sep = '\t',row.names = F,col.names = F)

library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file='infercnv/cell_anno.txt',
                                    ref_group_names=c("Normal"),
                                    gene_order_file='infercnv/gene_order_file.txt', 
                                    delim="\t")

options(scipen = 100)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv/output_dir", 
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=T,
                             HMM_type="i3",
                             analysis_mode="subclusters",
                             HMM_report_by="subcluster",
                             output_format="pdf",
                             no_prelim_plot=T,
                             num_threads = 8)
# #59800000
pred_cnv_genes=read.delim('results/02.B_cell_inferCNV/HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat')
head(pred_cnv_genes)
pred_cnv_genes$cell_type=str_split_fixed(pred_cnv_genes$cell_group_name,'[.]',2)[,2]
pred_cnv_genes$chr=str_split_fixed(pred_cnv_genes$gene_region_name,'-',2)[,1]
head(pred_cnv_genes)
table(pred_cnv_genes$state)
table(pred_cnv_genes$cell_type)
pred_cnv_genes$tissue=str_split_fixed(pred_cnv_genes$cell_type,'_',2)[,1]
# chr6p_del_genes=subset(pred_cnv_genes,tissue=='Tumor' & chr=='chr6' & end<59800000 & state==1)
# chr6q_del_genes=subset(pred_cnv_genes,tissue=='Tumor' & chr=='chr6' & start>59800000 & end<170805979 & state==1)
chr6_del_genes=subset(pred_cnv_genes,tissue=='Tumor' & chr=='chr6' & end<170805979 & state==1)
length(unique(chr6_del_genes$gene))
# chr19q_amp_genes=subset(pred_cnv_genes,tissue=='Tumor' & chr=='chr19' & start>26200000 & end<58617616 & state==3)
chr1q_amp_genes=subset(pred_cnv_genes,tissue=='Tumor' & chr=='chr1' & start>123400000 & end<248956422 & state==3)
length(unique(chr1q_amp_genes$gene))
library(clusterProfiler)
library(org.Hs.eg.db)
SYMBOL_id=unique(chr6_del_genes$gene)
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
chr6.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",ont = "BP",
                        pvalueCutoff = 0.05)
chr6.erich.go.BP.res=chr6.erich.go.BP@result
write.xlsx(chr6.erich.go.BP.res,'results/02.B_cell_inferCNV/chr6_enrichment.xlsx',overwrite = T)
chr6.erich.go.BP.res=chr6.erich.go.BP.res[chr6.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
chr6.erich.go.BP.res=chr6.erich.go.BP.res %>% slice_min(n =10, order_by = p.adjust)
pdf('02.B_cell_inferCNV/chr6_enrichment.pdf',height = 6,width = 8,onefile = F)
ggplot(data=chr6.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#80B1D3", high = "#FDB462")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
dev.off()
library(clusterProfiler)
library(org.Hs.eg.db)
SYMBOL_id=unique(chr1q_amp_genes$gene)
length(SYMBOL_id)
entrez_id = mapIds(x = org.Hs.eg.db,keys = SYMBOL_id,
                   keytype = "SYMBOL",column = "ENTREZID")
entrez_id = na.omit(entrez_id)
chr1.erich.go.BP = enrichGO(gene = entrez_id,OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",ont = "BP",
                       pvalueCutoff = 0.05)
chr1.erich.go.BP.res=chr1.erich.go.BP@result
write.xlsx(chr1.erich.go.BP.res,'results/02.B_cell_inferCNV/chr1q_enrichment.xlsx',overwrite = T)
chr1.erich.go.BP.res=chr1.erich.go.BP.res[chr1.erich.go.BP.res$p.adjust<0.05,c('Description','GeneRatio','Count','p.adjust')]
chr1.erich.go.BP.res=chr1.erich.go.BP.res %>% slice_min(n =5, order_by = p.adjust)
pdf('02.B_cell_inferCNV/chr1q_enrichment.pdf',height = 6,width = 8,onefile = F)
ggplot(data=chr1.erich.go.BP.res,aes(x=Count,y=reorder(Description,Count), color = -log10(p.adjust))) +
  geom_point(aes(size=Count),show.legend = T) +
  scale_color_gradient(low = "#BEBADA", high = "#FB8072")+
  geom_segment(aes(yend=Description,xend=0),size=.5) +
  labs(x='Count',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 35))
dev.off()
