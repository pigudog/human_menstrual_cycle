rm(list = ls())
library(readr)
library(rhdf5)
library(tidyverse)
library(patchwork)
library(stringr)
library(Matrix)
library(hdf5r)
suppressMessages(library(plotly))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(gplots))
suppressMessages(library(genefilter))
suppressMessages(library(scran))
library(DoubletFinder)

#####################################################################3
# create 10x
mentral_10x_day_donor_ctype <- read_csv("GSE111976_summary_10x_day_donor_ctype.csv")
mentral_10x_donor_phase <- read_csv("GSE111976_summary_10x_donor_phase.csv")
GSE111976_ct_endo_10x <- readRDS("~/mentral/GSE111976/scRNA/GSE111976_ct_endo_10x.rds")
cellinfo <- subset(mentral_10x_day_donor_ctype, select = c("day", "donor", "cell_type"))
phase <- cellinfo$donor
secretory_early  <- grep("14", phase, fixed = TRUE) 
phase[secretory_early] = "secretory_early"
proliferative  <- grep("19", phase, fixed = TRUE) 
phase[proliferative] = "proliferative"
secretory_mid  <- grep("20", phase, fixed = TRUE) 
phase[secretory_mid] = "secretory_mid"
secretory_earlymid  <- grep("29", phase, fixed = TRUE) 
phase[secretory_earlymid] = "secretory_early-mid"
secretory_earlymid  <- grep("39", phase, fixed = TRUE) 
phase[secretory_earlymid] = "secretory_early-mid"
secretory_earlymid  <- grep("41", phase, fixed = TRUE) 
phase[secretory_earlymid] = "secretory_early-mid"
secretory_late  <- grep("57", phase, fixed = TRUE) 
phase[secretory_late] = "secretory_late"
secretory_early  <- grep("58", phase, fixed = TRUE) 
phase[secretory_early] = "secretory_early"
secretory_early  <- grep("60", phase, fixed = TRUE) 
phase[secretory_early] = "secretory_early"
proliferative  <- grep("63", phase, fixed = TRUE) 
phase[proliferative] = "proliferative"
table(phase)
metadata_10x <- cbind(cellinfo,phase)
dim(metadata_10x)
dim(GSE111976_ct_endo_10x)
# > dim(metadata_10x)
# [1] 71032     4
# > dim(GSE111976_ct_endo_10x)
# [1] 33538 71032

# addmeatadata
mentral_10x <- CreateSeuratObject(counts = GSE111976_ct_endo_10x,project = "10x")
mentral_10x <- AddMetaData(object = mentral_10x,     #seurat object
                    metadata = metadata_10x$day,    #add needed metadata
                    col.name = "day")  #ranme
mentral_10x <- AddMetaData(object = mentral_10x,    
                           metadata = metadata_10x$donor,    
                           col.name = "donor")  
mentral_10x <- AddMetaData(object = mentral_10x,    
                           metadata = metadata_10x$cell_type,    
                           col.name = "cell_type")  
mentral_10x <- AddMetaData(object = mentral_10x,    
                           metadata = metadata_10x$phase,    
                           col.name = "phase")  

save(mentral_10x,file = "mentral_10x_raw.Rdata")


# create C1
mentral_C1_day_donor_ctype <- read_csv("GSE111976_summary_C1_day_donor_ctype.csv")
mentral_C1_donor_phase <- read_csv("GSE111976_summary_C1_donor_phase.csv")
GSE111976_ct_endo_C1 <- read_csv("~/mentral/GSE111976/scRNA/GSE111976_ct.csv")
dim(GSE111976_ct_endo_C1)
dim(mentral_C1_day_donor_ctype)
# Create the Seurat object with all the data (unfiltered)
raw_data <- as.data.frame(GSE111976_ct_endo_C1)
rownames(raw_data)<-raw_data[,1]
raw_data <-raw_data[,2:2149]
mentral_C1 <- CreateSeuratObject(counts = raw_data)

# Creat metadata
meta_data <- as.data.frame(mentral_C1_day_donor_ctype)
# add rownames to metadta 
rownames(meta_data) <- meta_data$cell_name
meta_data<-meta_data[,-1]
meta_data<-meta_data[,1:3]
phase <- meta_data$donor
for (i in 1:length(phase)){
  for (j in 1:length(mentral_C1_donor_phase$donor)){
    if (phase[i] ==mentral_C1_donor_phase$donor[j] ){
      phase[i]<-mentral_C1_donor_phase$phase_canonical[j]
      break
    }
  }
}
table(phase)
metadata_C1 <- cbind(meta_data,phase)
# add metadata to Seurat object 
mentral_C1 <- AddMetaData(object = mentral_C1, metadata = metadata_C1)
save(mentral_C1,file = "mentral_C1_raw.Rdata")

#############################################################################################
rm(list = ls())
# For article repeatability, we set up random seeds
set.seed(1314)
# check the different between C1 and 10x DATA
load("~/mentral/GSE111976/scRNA/mentral_10x_raw.Rdata")
load("~/mentral/GSE111976/scRNA/mentral_C1_raw.Rdata")

# check the cell_type before
table(mentral_10x$cell_type)
table(mentral_C1$cell_type)

# check the 
table(mentral_10x$phase)
table(mentral_C1$phase)

# let's start it with 10x!!!
#####################################################################
### Standard pre-processing workflow
# 1. SoupX -- without 10x raw data
# 2. Screen QC
##mitochondrion, ribosome、erythrocyte
mentral_10x[["percent.mt"]] <- PercentageFeatureSet(mentral_10x, pattern = "^MT-")
mentral_10x[["percent.rb"]] <- PercentageFeatureSet(mentral_10x, pattern = "^RP[SL]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(mentral_10x))
mentral_10x[["percent.HB"]]<-PercentageFeatureSet(mentral_10x, features=HB.genes)
# cannot find ercc(External RNA Control Consortium): one of spike in
table(grepl("^ERCC-",rownames(mentral_10x)))
# FALSE 
# 33538

### QC violin picture
# possible theme
theme.set2 = theme(axis.title.x=element_blank())
# the elements of picture
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "donor"
# vlnplot before qc
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(mentral_10x, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("./preanalysis10x/vlnplot_before_qc.png", plot = violin, width = 12, height = 8)
### set QC strandards
minGene=100
maxGene=10000
maxUMI=100000
pctMT=30
pctHB=1
### vlnplot after QC
mentral_10x <- subset(mentral_10x, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB <pctHB)
dim(mentral_10x)
# [1] 33538 66730
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(mentral_10x, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("./preanalysis10x/vlnplot_after_qc.png", plot = violin, width = 12, height = 8) 


save(mentral_10x,file = 'mentral_10x_after_screen.Rdata')

# 3. harmony for integrate
library(harmony)
library(tidyverse)
library(patchwork)
### standardize data with SCT v2!
mentral_10x <- SCTransform(mentral_10x, vst.flavor = "v2")
DefaultAssay(mentral_10x)
# [1] "SCT"

### PCA
mentral_10x <- RunPCA(mentral_10x, npcs=50, verbose=FALSE)
table(mentral_10x$donor)
mentral_10x <- RunHarmony(mentral_10x, 
                    group.by.vars="orig.ident",
                    reduction = "pca", 
                    assay.use="SCT", 
                    max.iter.harmony = 20) 
# group.by.vars: The parameter is to set which group to integrate by
# max.iter.harmony: Set the number of iterations, default is 10. When running RunHarmony, results will indicate how many iterations have elapsed before convergence is complete.

elowplot <- ElbowPlot(mentral_10x, ndims = 50)
ggsave("./preanalysis10x/elbow_plot.png", plot = elowplot, width = 12, height = 8) 

pc.num=1:30
mentral_10x <- RunUMAP(mentral_10x, reduction="harmony", dims=pc.num)
mentral_10x <- FindNeighbors(mentral_10x, reduction = "harmony",
                       dims = pc.num)
mentral_10x<-FindClusters(mentral_10x,resolution=0.8)
p <- DimPlot(mentral_10x, group.by = "donor")
ggsave("./preanalysis10x/UMAP_donor_primality_0.8.png", p, width = 8, height = 6)
p <- DimPlot(mentral_10x, group.by = "cell_type")
ggsave("./preanalysis10x/UMAP_cell_type_primality_0.8.png", p, width = 8, height = 6)
p <- DimPlot(mentral_10x, group.by = "phase")
ggsave("./preanalysis10x/UMAP_phase_primality_0.8.png", p, width = 8, height = 6)
p <- DimPlot(mentral_10x, group.by = "seurat_clusters")
ggsave("./preanalysis10x/UMAP_seurat_clusters_primality_0.8.png", p, width = 8, height = 6)
p <- DimPlot(mentral_10x, group.by = "donor", split.by = "donor", ncol = 4)
ggsave("./preanalysis10x/UMAP_donor_primality_split_0.8.png", p, width = 18, height = 12)
save(mentral_10x,file = 'mentral_10x_after_harmony.Rdata')

###########################################################################################
# Find HVG and cluster
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
rm(list=ls())
options(stringsAsFactors = F)
mentral_10x.markers <- FindAllMarkers(object = mentral_10x,assay = "SCT", only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
# wriet to decidual/HVGmarkers_res_0.5
write.csv(mentral_10x.markers,file=paste0("./preanalysis10x",'/HVGmarkers_res_0.8.csv'))
save(mentral_10x,mentral_10x.markers,file = 'mentral_10x_HVG.Rdata')

library(dplyr) 
p<-DimHeatmap(object = mentral_10x,dims = 1:15)
ggsave(p,filename=paste0("preanalysis10x",'/mentral_10x_dimheatmap.pdf'),width = 24,height = 18)
top10 <- mentral_10x.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p <- DoHeatmap(mentral_10x,top10$gene,size=3)
ggsave(p,filename=paste0("preanalysis10x",'/top10_heatmap.png'),width = 24,height = 18)
# Specify genes  - Mac
genes_to_check = c( "FCN1","MS4A7","CD14")
# featureplot
p <- FeaturePlot(mentral_10x, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis10x",'/fearureplot_sepcify_mac.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(mentral_10x, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis10x",'/dotplot_sepcify_specific_MAC.pdf'),width = 16,height = 12)
# 13 macrophage


#################################################################################################
# cell annotation
# Specify genes  - Mac stromal
# Stromal fibroblast - ADGRL4 COL5A1 LUM COL6A3 CRISPLD2 COL6A1
# epithelium TACSTD2 UCA1 SDC4 KLF5 WFDC2 EPCAM
## Unciliated epithelium(glandular)
## Unciliated epithelium(luminal)
## Ciliated epithelium FOXJ1 DYNLRB2 SNTN C9orf24 CDHR3 DYDC2 
# Lymphocyte PTPRC CCL5 STK17B
# Macropage LYZ IL1B CD14 HLA-DQA1 MSA4A6A AIF1
# Smooth muscle cell ACTA2 MCAM BGN NOTCH3 GUCY1A2 RGS5
genes_to_check = c("ADGRL4","COL5A1","LUM","COL6A3","CRISPLD2","COL6A1",  # Stromal fibroblast
                   "TACSTD2","UCA1","SDC4","KLF5","WFDC2","EPCAM",        # Unciliated epithelium
                   "FOXJ1","","DYNLRB2","SNTN","C9orf24","CDHR3","DYDC2", # Ciliated epithelium
                   "PTPRC","CCL5","STK17B",                               # Lymphocyte
                   "LYZ","IL1B","CD14","HLA-DQA1","MSA4A6A","AIF1",       # Macrophage
                   "ACTA2","MCAM","BGN","NOTCH3","GUCY1A2","RGS5"         # Smooth muscle cell
                   )
# featureplot
p <- FeaturePlot(mentral_10x, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis10x",'/fearureplot_sepcify.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(mentral_10x, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis10x",'/dotplot_sepcify_specific.png'),width = 16,height = 12)


#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:25,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,4,6,13,14,15,19,20,22),2]='Stromal fibroblast' 
celltype[celltype$ClusterID %in% c(1,2,3,5,7,8,10,11,12,16,18,21,24,25),2]='Unciliated epithelium'
celltype[celltype$ClusterID %in% c(9,17),2]='Ciliated epithelium'
celltype[celltype$ClusterID %in% c(8,10),2]='Lymphocyte'
celltype[celltype$ClusterID %in% c(23),2]='Macrophage'
celltype[celltype$ClusterID %in% c(15),2]='Smooth muscle cell'


head(celltype)
celltype 
table(celltype$celltype)

new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(mentral_10x)
mentral_10x <- RenameIdents(mentral_10x, new.cluster.ids)
save(mentral_10x,file = 'mentral_10x_after_cluster.Rdata')
# legend
p<-DimPlot(mentral_10x, reduction = "umap", label = TRUE, pt.size = 0.5) 
ggsave(p,filename=paste0("preanalysis10x",'/umap_firstcluster.png'),width = 16,height = 12)

# percentage - donor
sample_table <- as.data.frame(table(mentral_10x@meta.data$donor,mentral_10x@active.ident))
names(sample_table) <- c("Samples","celltype","CellNumber")
# color
library(ggsci)
library(Cairo)
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
ggsave(plot_sample,filename=paste0("preanalysis10x",'/percentage_donor_ctype.png'),width = 6,height = 8)

# percentage -phase
sample_table <- as.data.frame(table(mentral_10x@meta.data$phase,mentral_10x@active.ident))
names(sample_table) <- c("Phase","celltype","CellNumber")

plot_sample<-ggplot(sample_table,aes(x=Phase,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
ggsave(plot_sample,filename=paste0("preanalysis10x",'/percentage_phase_ctype.png'),width = 6,height = 8)

