
# package准备： --------------------------------------------------------------

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)


# 对数据进行预处理 ----------------------------------------------------------------
data("sim_trajectory", package = "miloR")

##创建一个singlecellexperiment对象
traj_sce <- sim_trajectory[['SCE']]

colData(traj_sce) <- DataFrame(traj_meta)

##
logcounts(traj_sce) <- log(counts(traj_sce) + 1)
traj_sce <- runPCA(traj_sce, ncomponents=30)
traj_sce <- runUMAP(traj_sce)

plotUMAP(traj_sce)


# 创建一个Milo对象 --------------------------------------------------------------
##从singlecellexperiment转化为milo对象
traj_milo <- Milo(traj_sce)

reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")

traj_milo


# 创建一个KNN 图像 --------------------------------------------------------------

traj_milo <- buildGraph(traj_milo, k = 10, d = 30)


# 定义一个具有代表性的邻域 ------------------------------------------------------------
##在Milo指南中有对每个函数意思进行标注
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

##在得到邻域之后，最好检查一下邻域的大小
plotNhoodSizeHist(traj_milo)


# 对领域细胞进行计数 ---------------------------------------------------------------
##对每个邻域的每个样本进行计数
traj_milo <- countCells(traj_milo, 
                        meta.data = data.frame(colData(traj_milo)),
                        samples="Sample")

head(nhoodCounts(traj_milo))


# 差异丰度检测 ------------------------------------------------------------------
##我们需要思考实验设计，
##并且设计的矩阵必须将感兴趣的条件和样本匹配，此时condition成为了协变量

traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Condition")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$Sample

##重新排序行名去匹配milo中的列

traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]

traj_design

##使用这个统计数据，需要在Milo对象中存储最近邻居之间的距离。
traj_milo <- calcNhoodDistance(traj_milo, d=30)

rownames(traj_design) <- traj_design$Sample
da_results <- testNhoods(traj_milo, design = ~ Condition, design.df = traj_design)

##计算每个邻域的Fold-change和修正的p值，这表明条件之间是否存在显著的丰度差异。
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 


# 可视化邻域及展示差异丰度 ------------------------------------------------------------

traj_milo <- buildNhoodGraph(traj_milo)

plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

#对da_results数据按照"SpatialFDR"列进行排序，并返回排序后的前几行数据。

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# 检测差异分析的结果 ---------------------------------------------------------------
#1.先检查未校正P值的分布，以验证检验是否平衡。
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

#火山图可视化结果（每个点代表一个邻域）
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)



embryo_milo <- buildNhoodGraph(embryo_milo)
umap_pl <- plotReducedDim(embryo_milo, dimred = "umap", colour_by="stage", text_by = "celltype", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results, layout="umap",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")



