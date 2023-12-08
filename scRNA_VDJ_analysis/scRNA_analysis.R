library(readxl)
library(ggsignif)
library(magrittr)
library(ggrepel)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(dplyr)
library(patchwork)
### Set working dir ###
pro.dir=getwd() #用当前目录作为工作目录
sam<-c("NP-gB","NP-gH")
merge.list=c()
for(i in sam){
  #### Read data ####
  cym1<-Read10X(paste0(pro.dir,"/",i,"_filtered_feature_bc_matrix/")) #读取CellCosmo的output矩阵
  obj<-CreateSeuratObject(cym1,project = i)
  HB.genes <- c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.HB"]]<-PercentageFeatureSet(obj, features=HB.genes) 
  obj <- subset(obj, subset = nFeature_RNA > 200  & 
                  percent.mt < 20 & percent.HB < 5)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  merge.list=c(merge.list,obj)
}
#### Merge for integration ####
mouse <- merge(merge.list[[1]],merge.list[2:length(merge.list)])
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mouse)
mouse <- ScaleData(mouse, features = all.genes)
mouse <- RunPCA(mouse, features = VariableFeatures(mouse))
mouse <- RunUMAP(mouse, dims = 1: 30)
mouse <- RunTSNE(mouse, dims = 1: 30)
#### Clustering & markers ####
mouse <- FindNeighbors(mouse, dims = 1:30)
mouse <- FindClusters(mouse, resolution = 0.5)
all.markers <- FindAllMarkers(mouse, only.pos = TRUE,logfc.threshold = 0.5,min.pct = 0.25)
write.csv(all.markers15,"ALL_Markers.csv")
## Assign Celltypes Manually ##
########Visualisation#######
dir.create("../Annotation",recursive = T)
setwd("../Annotation")

col1=c("#FFADAD","#0F7B6D","#FDFFB6","#6940A6","#DFAB00","#9BF6FF","#A0C4FF","#FFC6FF","#DF3E3E")
pdf("Annotation_umap.pdf",width = 5,height = 5)
DimPlot(mouse, reduction = "umap", label = T,cols = c(col1),
        pt.size = 0.1,label.size = 4,repel = TRUE,group.by = 'annotation')
dev.off()
pdf("Annotation_umap_split.pdf",width = 12)
show(DimPlot(mouse, reduction = "umap", label = T,split.by = 'orig.ident',cols = c(col1),
             pt.size = 0.1,label.size = 4,repel = TRUE,group.by = 'annotation',))
dev.off()
pdf("Annotation_umap_no_label.pdf",width = 5,height = 5)
DimPlot(mouse, reduction = "umap", label = F,cols = c(col1),
        pt.size = 0.1,label.size = 4,repel = TRUE,group.by = 'annotation')
dev.off()
pdf("Annotation_umap_split_no_label.pdf",width = 12)
show(DimPlot(mouse, reduction = "umap", label = F,split.by = 'orig.ident',cols = c(col1),
             pt.size = 0.1,label.size = 4,repel = TRUE,group.by = 'annotation',))
dev.off()
pdf("Annotation_tsne.pdf")
DimPlot(mouse, reduction = "tsne", label = T,cols = c(col1),
        pt.size = 0.6,label.size = 4,repel = TRUE,group.by = 'annotation')
dev.off()
### Markers ###
Idents(mouse)=mouse$annotation
cell.markers=FindAllMarkers(mouse,only.pos = T)
write.csv(cell.markers,"CellMarkers.csv")
ma_p<-cell.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
col2=brewer.pal(9,"Reds")
pdf("Dot_marker_annotation_mouse.pdf",height = 10)
DotPlot(mouse,features = unique(ma_p$gene),
        group.by = 'annotation',cols = col2[c(1,7)])+
  xlab("Markers")+coord_flip()+
  ylab("Celltypes")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

####### VDJ analysis ########
#read data gB
tcrgB<-read.csv("NP-gB_vdj_b/filtered_contig_annotations_NP-gB.csv")
tcrgH<-read.csv("NP-gH_vdj_b/filtered_contig_annotations_NP-gH.csv")
seq_stat_gb=as.data.frame(table(tcrgB$cdr3_nt)) 
seq_stat_gb=seq_stat_gb[order(seq_stat_gb$Freq,decreasing = T),]
seq_stat_gb=rbind(seq_stat_gb,data.frame(Var1="Singleton",Freq=sum(seq_stat_gb$Freq==1)))
seq_stat_gb=filter(seq_stat_gb,seq_stat_gb$Freq!=1)
write.csv(seq_stat_gb,"Sequence_stat_gb.csv")
ig_gb=tcrgB[,c("barcode","c_gene")]
ig_gb$igtype="IgK"
ig_gb$igtype[grep(pattern = "^IGHG",x=ig_gb$c_gene)]="IgG"
ig_gb$igtype[grep(pattern = "^IGHA",x=ig_gb$c_gene)]="IgA"
ig_gb$igtype[grep(pattern = "^IGHM",x=ig_gb$c_gene)]="IgM"
ig_gb$igtype[grep(pattern = "^IGHD",x=ig_gb$c_gene)]="IgD"
ig_gb$igtype[grep(pattern = "^IGLC",x=ig_gb$c_gene)]="Igλ"
ig_stat_gb=as.data.frame(table(ig_gb$igtype)) 
write.csv(ig_stat_gb,"IG_stat_gb.csv")

seq_stat_gh=as.data.frame(table(tcrgH$cdr3_nt)) 
seq_stat_gh=seq_stat_gh[order(seq_stat_gh$Freq,decreasing = T),]
seq_stat_gh=rbind(seq_stat_gh,data.frame(Var1="Singleton",Freq=sum(seq_stat_gh$Freq==1)))
seq_stat_gh=filter(seq_stat_gh,seq_stat_gh$Freq!=1)
write.csv(seq_stat_gh,"Sequence_stat_gh.csv")

ig_gh=tcrgH[,c("barcode","c_gene")]
ig_gh$igtype="IgK"
ig_gh$igtype[grep(pattern = "^IGHG",x=ig_gh$c_gene)]="IgG"
ig_gh$igtype[grep(pattern = "^IGHA",x=ig_gh$c_gene)]="IgA"
ig_gh$igtype[grep(pattern = "^IGHM",x=ig_gh$c_gene)]="IgM"
ig_gh$igtype[grep(pattern = "^IGHD",x=ig_gh$c_gene)]="IgD"
ig_gh$igtype[grep(pattern = "^IGLC",x=ig_gh$c_gene)]="Igλ"
ig_stat_gh=as.data.frame(table(ig_gh$igtype)) 
write.csv(ig_stat_gh,"IG_stat_gh.csv")

##############
df1=data.frame(sample="NP-gB",j_gene=tcrgB$j_gene)
df2=data.frame(sample="NP-gH",j_gene=tcrgH$j_gene)
df=rbind(df1,df2) 
df=df[df$j_gene!="",]
pr_t <- prop.table(table(df$j_gene,df$sample))
pr_t<-as.data.frame(pr_t)
pr_t$Freq <-pr_t$Freq * 100
colnames(pr_t)<-c("Jgene","Group","Freq")
write.csv(pr_t,"Proportion_Jgene_usage.csv")
pdf("Proportion_Jgene_usage.pdf",width = 10)
ggplot(pr_t, aes(x = Jgene , y = Freq, fill =Group)) + 
  geom_bar(stat = 'identity',width = 0.8,position = position_dodge()) +#
  xlab("Celltypes")+
  ylab("Percentage(%)")+
  scale_fill_manual(values=c(color4))+theme_bw()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
dev.off()
#######
tcrgB<-tcrgB[!duplicated(tcrgB$barcode),]
tcrgB$barcode=paste0(tcrgB$barcode,"_1")
tcrgB <- tcrgB[,c("barcode","raw_clonotype_id","chain","c_gene")]
names(tcrgB)[names(tcrgB)=="raw_clonotype_id"]<- "clonotype_id"
clonogB<-read.csv("../../NP-gB_vdj_b/clonotypes.csv")
tcrgB <- merge(tcrgB,clonogB[,c("clonotype_id","cdr3s_aa")])
names(tcrgB)[names(tcrgB)=="cdr3s_aa"]<- "cdr3s_aa"

tcrgH<-tcrgH[!duplicated(tcrgH$barcode),]
tcrgH$barcode=paste0(tcrgH$barcode,"_2")
tcrgH <- tcrgH[,c("barcode","raw_clonotype_id","chain","c_gene")]
names(tcrgH)[names(tcrgH)=="raw_clonotype_id"]<- "clonotype_id"
clonogH<-read.csv("NP-gH_vdj_b/clonotypes.csv")
tcrgH <- merge(tcrgH,clonogH[,c("clonotype_id","cdr3s_aa")])
names(tcrgH)[names(tcrgH)=="cdr3s_aa"]<- "cdr3s_aa"

tcr=rbind(tcrgB,tcrgH)
tcr <-tcr[,c(2,1,3,5,4)]
rownames(tcr)<-tcr[,1]
tcr[,1]<-NULL
colnames(tcr)<- paste0("b",colnames(tcr))
mouse<-AddMetaData(mouse,tcr)
mouse$annotation=droplevels(mouse$annotation)
#####
pr_t <- prop.table(table(mouse$bchain, mouse$annotation), margin = 2)
pr_t<-as.data.frame(pr_t) %>% na.omit() %>% droplevels()
pr_t$Freq <-pr_t$Freq * 100
colnames(pr_t)<-c("Chain","CellType","Freq")
write.csv(pr_t,paste0("proportion_Bchain.csv"))
pdf("Proportion_Bchain.pdf")
ggplot(pr_t, aes(x = CellType , y = Freq, fill =Chain)) + 
  geom_bar(stat = 'identity',width = 0.8) +#
  xlab("Celltypes")+
  ylab("Percentage(%)")+
  scale_fill_manual(values=c(color4))+theme_bw()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
dev.off()
#######
pr_t <- prop.table(table(mouse$bc_gene, mouse$annotation), margin = 2)
pr_t<-as.data.frame(pr_t) %>% na.omit() %>% droplevels()
levels(pr_t$Var1)[1]="No type"
pr_t$Var1[which(pr_t$Var1=="")]="No type"
pr_t$Freq <-pr_t$Freq * 100
colnames(pr_t)<-c("C_gene","CellType","Freq")
write.csv(pr_t,paste0("proportion_c_gene.csv"))
pdf("Proportion_c_gene.pdf")
ggplot(pr_t, aes(x = CellType , y = Freq, fill =C_gene)) + 
  geom_bar(stat = 'identity',width = 0.8) +#
  xlab("Celltypes")+
  ylab("Percentage(%)")+
  scale_fill_manual(values=c(color1,color2))+theme_bw()+
  theme(text=element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank())
dev.off()
#### monocle3 ####
library(monocle3)
cds <- as.cell_data_set(mouse)
cds <- cluster_cells(cds,cluster_method = 'louvain')
cds <- learn_graph(cds)
cds@clusters$UMAP$clusters=Idents(obj)
pdf(paste("trajectory.pdf",sep = ""))
show(plot_cells(cds, label_groups_by_cluster = F, label_leaves = T, 
                label_branch_points = F,trajectory_graph_segment_size = 0.7,
                group_label_size = 3,cell_size = 0.6)+scale_color_manual(values=c(color1,color2,color3))
)
dev.off()
max.avp <- which.max(unlist(FetchData(obj, root)))
max.avp <- colnames(obj)[max.avp]
cds <- order_cells(cds,root_cells = max.avp)
cds@clusters$UMAP$clusters=mouse$annotation
cds <- order_cells(cds)
pdf(paste0("trajectory_mouse_pseudotime.pdf"))
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T, 
           label_branch_points = F,trajectory_graph_segment_size = 0.5,group_label_size = 3,cell_size = 1)
dev.off()
pdf("trajectory_mouse.pdf")
plot_cells(cds, label_groups_by_cluster = F, label_leaves = T,
           label_branch_points = F,trajectory_graph_segment_size = 0.7,
           group_label_size = 4,cell_size = 0.6,graph_label_size = 1)+scale_color_manual(values=col1)
dev.off()
##########RNA Velocity############
dir.create("../Velocity",recursive = T)
setwd("../Velocity")
library(velocyto.R)
library(umap)
########## Read data 
#NP-gB
ldata <- read.loom.matrices("velocyto/NP-gB.loom")

head(colnames(ldata$spliced))
mouse$orig.ident[1]
mouse$orig.ident[3000]
colnames(ldata$spliced)<-gsub("sample_alignments_6UPJ7:","",colnames(ldata$spliced)) %>% gsub("x","-1_1",.)
colnames(ldata$unspliced)<-colnames(ldata$spliced)
colnames(ldata$ambiguous)<-colnames(ldata$spliced)
#Filter loom data with seurat cells
validcells <- intersect(colnames(ldata$spliced),colnames(mouse))
ldata$spliced<-ldata$spliced[,validcells]
ldata$unspliced<-ldata$unspliced[,validcells]
ldata$ambiguous<-ldata$ambiguous[,validcells]
sp_7 <- ldata$spliced #60
unsp_7 <- ldata$unspliced #

#NP-gH
ldata <- read.loom.matrices("velocyto/NP-gH.loom")

head(colnames(ldata$spliced))
mouse$orig.ident[1]
mouse$orig.ident[3000]
colnames(ldata$spliced)<-gsub("sample_alignments_7BFEK:","",colnames(ldata$spliced)) %>% gsub("x","-1_2",.)
colnames(ldata$unspliced)<-colnames(ldata$spliced)
colnames(ldata$ambiguous)<-colnames(ldata$spliced)
#Filter loom data with seurat cells
validcells <- intersect(colnames(ldata$spliced),colnames(mouse))
ldata$spliced<-ldata$spliced[,validcells]
ldata$unspliced<-ldata$unspliced[,validcells]
ldata$ambiguous<-ldata$ambiguous[,validcells]

sp_14 <- ldata$spliced
unsp_14 <- ldata$unspliced
##MERGE
sp <-cbind(sp_7,sp_14)
unsp<-cbind(unsp_7,unsp_14)
umap <- mouse@reductions$umap@cell.embeddings
#Calculate cell distance & evaluate velocity
cell.dist <- as.dist(1-armaCor(t(umap)))

rvel <-gene.relative.velocity.estimates(sp,unsp,deltaT2 = 2,kCells = 10,
                                        
                                        cell.dist = cell.dist,fit.quantile = 0.02,n.cores = 32)

saveRDS(rvel,"Relative_velo.RDS")

#Draw the plot
gg <- DimPlot(mouse, reduction = "umap", label = F,cols = c(col1),
              pt.size = 0.1,label.size = 4,repel = TRUE,group.by = 'annotation')
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors)<-rownames(umap)

pdf("RNA_Velocity_NP.pdf")

p4 <- show.velocity.on.embedding.cor(umap,rvel,n = 50,scale = "sqrt",cell.colors = ac(colors,alpha = 0.5),cex=0.8,
                                     
                                     arrow.scale = 5,min.grid.cell.mass=1.0,show.grid.flow = T,n.cores = 12,
                                     
                                     grid.n=40,arrow.lwd=0.5,do.par=F,cell.border.alpha=0.05,min.arrow.size = 0.05,
                                     
                                     main="Cell Velocity NP",return.details =F)

dev.off( )
