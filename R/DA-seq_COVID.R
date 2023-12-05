
library(Seurat)
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

source("convenience.R")

python2use <- "C:/Program Files/python.exe"

GPU <- 3

##为FIt-SNE R包装器设置路径
fitsneR <- "D:/download/FIt-SNE-1.1.0/fast_tsne.R"
#fitsneR <- "D:/mingw64/bin/FIt-SNE/fast_tsne.R"

## Load data
main_S <- readRDS(url("https://ndownloader.figshare.com/files/22927382"))
##main_S <- readRDS("covid_nbt_main.rds")

DefaultAssay(main_S) <- "integrated"

head(main_S@meta.data)

table(main_S@meta.data$severity)

##检查
names(main_S@commands)

## analysis
#对数据进行表达特征的缩放，是的具有相似的尺度，消除尺度差异。
main_S <- ScaleData(main_S)

#进行主成分分析，主成分为90.verbose=F表示不显示详细信息
main_S <- RunPCA(main_S, npcs = 90, verbose = F)

#使用tSNE进行降维处理
main_S <- runFItSNE(
  main_S, dims.use = 1:90, seed.use = 3, fast.R.path = fitsneR,
  ann_not_vptree = FALSE, nthreads = 12
)

TSNEPlot(main_S, group.by = "severity")
TSNEPlot(main_S, group.by = "celltype", label = T)


## Separate immune cells from patient samples
# get immune cells
sort(unique(main_S@meta.data$celltype))

epi_type <- c("Basal","Ciliated","Ciliated-diff","FOXN4","Ionocyte","IRC",
              "Secretory","Secretory-diff","Squamous","outliers_epithelial","unknown_epithelial")

immune_type <- setdiff(sort(unique(main_S@meta.data$celltype)), epi_type)

immune_S <- subset(x = main_S, cells = which(main_S@meta.data$celltype %in% immune_type))


# remove control cells
immune_S <- subset(immune_S,
                   cells = which(immune_S@meta.data$severity != "control"))

immune_S@meta.data$severity <- factor(as.character(immune_S@meta.data$severity))

#可视化
TSNEPlot(immune_S, group.by = "severity")

TSNEPlot(immune_S, group.by = "celltype", label = T)

#为每个细胞类型分配一个数字标签

immune_S@meta.data$celltype_num <- as.numeric(factor(immune_S@meta.data$celltype))

#### DA-seq on immune cells
##获取样本标签
table(main_S@meta.data$patient)

table(main_S@meta.data$sample)

moderate_labels <- unique(main_S@meta.data$sample[main_S@meta.data$severity == "moderate"])

critical_labels <- unique(main_S@meta.data$sample[main_S@meta.data$severity == "critical"])

##找差异丰度细胞
da_cells_immune <- getDAcells(
  X = immune_S@reductions$pca@cell.embeddings,
  cell.labels = immune_S@meta.data$sample,
  labels.1 = moderate_labels, labels.2 = critical_labels,
  k.vector = seq(200,3200,200), 
  plot.embedding = immune_S@reductions$tsne@cell.embeddings,
)

da_cells_immune <- updateDAcells(
  da_cells_immune, pred.thres = c(-0.8,0.8), 
  plot.embedding = immune_S@reductions$tsne@cell.embeddings, size = 0.01
)

da_cells_immune$pred.plot
da_cells_immune$da.cells.plot

## DA regions

da_regions_immune <- getDAregion(
  X = immune_S@reductions$pca@cell.embeddings, da.cells = da_cells_immune,
  cell.labels = immune_S@meta.data$sample,
  labels.1 = moderate_labels, labels.2 = critical_labels,
  resolution = 0.01,
  plot.embedding = immune_S@reductions$tsne@cell.embeddings, size = 0.1
)
da_regions_immune$da.region.plot
da_regions_immune$DA.stat

n_da_immune <- length(unique(da_regions_immune$da.region.label)) - 1

##本地的marker
# set cluster overlap for each DA
da_clusters_immune <- c(
  "1" = "Neu", "2" = "nrMa", "4" = "CTL", "5" = "Neu"
)

# Seurat negbinom to find local markers
immune_S <- addDAslot(immune_S, da.regions = da_regions_immune)
local_markers_immune <- list()
for(i in 1:n_da_immune){
  if(!as.character(i) %in% names(da_clusters_immune)){next()}
  local_markers_immune[[as.character(i)]] <- SeuratLocalMarkers(
    object = immune_S, da.region.to.run = i, cell.label.slot = "celltype", 
    cell.label.to.run = da_clusters_immune[as.character(i)],
    assay = "RNA", test.use = "negbinom", min.diff.pct = 0.09, only.pos = T
  )
  local_markers_immune[[as.character(i)]]$pct.diff <- local_markers_immune[[as.character(i)]]$pct.1 -
    local_markers_immune[[as.character(i)]]$pct.2
}

lapply(c(1,2,4,5), function(x){
  write.table(
    local_markers_immune[[as.character(x)]], paste0("markers/covidChua_local_DA",x,".txt"), 
    sep = "\t", quote = F, row.names = T, col.names = T
  )
})



# liao et al.数据分析 ---------------------------------------------------------
#由网页下载转为绝对路径
liao_S <- readRDS("DA-seq/nCoV.rds")

##创建聚类
epi_clusters <- c(13,16,25,28,31)

macro_clusters <- as.character(c(0,1,2,3,4,5,7,8,
                                 10,11,12,18,21,22,23,26))

t_clusters <- as.character(c(6,9,14))

cd8t_clusters <- as.character(c(9,14))

#计算模块分数
liao_immune_S <- subset(liao_S, 
                        cells = which(!liao_S@meta.data$cluster %in% epi_clusters & liao_S$group != "HC"))

liao_immune_S@meta.data$severity <- gsub("O", "moderate", 
                                         liao_immune_S@meta.data$group1)

liao_immune_S@meta.data$severity <- gsub("S/C", "critical", 
                                         liao_immune_S@meta.data$severity)

da_gene_modules <- lapply(local_markers_immune, FUN = function(x){
  rownames(x)[1:min(100,nrow(x))]
})

names(da_gene_modules) <- paste0("DA", names(da_gene_modules))

liao_immune_S <- AddModuleScore(
  object = liao_immune_S, features = da_gene_modules, 
  assay = "RNA", name = names(da_gene_modules)
)

for(i in 1:4){
  colnames(liao_immune_S@meta.data)[grep(names(da_gene_modules)[i],
                                         colnames(liao_immune_S@meta.data))] <-
    names(da_gene_modules)[i]
}


##产生图表
library(scales)

da_cols <- hue_pal()(n_da_immune)

da_order <- order(da_regions_immune$da.region.label)

tsne_embedding <- immune_S@reductions$tsne@cell.embeddings

## TSNE plots
gg1 <- plotCellLabel(
  tsne_embedding, label = immune_S@meta.data$severity, size = 0.01, do.label = F
) + theme_tsne

ggsave(gg1, filename = "figs/covidChua_a.png", 
       width = 50, height = 50, units = "mm", dpi = 1200)

ggsave(g_legend(gg1, legend.position = "right"),
       filename = "figs/covidChua_a_legend.pdf", 
       width = 0.7, height = 0.3, dpi = 1200)


