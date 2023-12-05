r

# packages: ---------------------------------------------------------------

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(MouseGastrulationData)
#library(Matrix)


# 数据准备： -------------------------------------------------------------------
##为了使数据计算更快，调取部分数据样本（胚胎），
##并且数据通过主成分分析降维以及FastMNN进行批次矫正、数据整合
###（fastMNN算法由batchelor包提供，安装seurat包时已经自动安装过此包。
###要单独使用fastMNN算法整合数据，需要使用SeuratWrappers包）

select_samples <- c(2,  3,  6, 15,
                    10, 14, 20, 30
)
#EmbryoAtlasData小鼠原肠胚形成时间数据
embryo_data = EmbryoAtlasData(samples = select_samples)

embryo_data


# 数据可视化 -------------------------------------------------------------------
##重新计算该细胞子集的UMAP嵌入以可视化数据。

#去除含有缺失值的样本，只保留非缺失的样本
embryo_data <- embryo_data[,apply(reducedDim(embryo_data, "pca.corrected"), 
                                  1, function(x) !all(is.na(x)))]

embryo_data <- runUMAP(embryo_data, dimred = "pca.corrected", name = 'umap')

plotReducedDim(embryo_data, colour_by="stage", dimred = "umap") 



# 差异丰度检测---创建一个Milo对象 -----------------------------------------------------
##用milo对象可以储存邻域KNN图信息
embryo_milo <- Milo(embryo_data)

embryo_milo


# 创建KNN-graph -------------------------------------------------------------
#需要添加KNN到milo中，构建一个MNN矫正的PCA维度
##参数定义：d代表降维数；k代表邻域细胞数量

embryo_milo <- buildGraph(embryo_milo, k = 30, d = 30, reduced.dim = "pca.corrected")


# 在KNN中定义一个具有代表性的邻域 -------------------------------------------------------
#定义一个单元的邻域，即索引，作为KNN图中通过一条边连接到索引单元的单元组。
#******采用Gut2015.KNN采样方法，对代表性细胞的子集进行采样作为指标。
#参数定义：prop代表开始时随机取样的细胞比例，refined
embryo_milo <- makeNhoods(embryo_milo, prop = 0.1, k = 30, d=30, 
                          refined = TRUE, reduced_dims = "pca.corrected")

#用此检测k值是否合适
#（如果直方图的峰值较为集中且不过于偏离某个特定的邻域大小，
#说明k值选择合理，细胞的邻域大小分布较为均匀）
plotNhoodSizeHist(embryo_milo)


# 邻域细胞计数 ------------------------------------------------------------------
#需要计算每个邻域的每个样本的细胞数量
#（milo靠检测重复样本中细胞数量变化来检测差异丰度）

embryo_milo <- countCells(embryo_milo,
                          meta.data = as.data.frame(colData(embryo_milo)),
                          sample="sample")

head(nhoodCounts(embryo_milo))



# 定义实验设计 ------------------------------------------------------------------
#利用GLM模型，使用edgeR的负二项GLM
#设计实验用于进行差异表达分析
embryo_design <- data.frame(colData(embryo_milo))[,c("sample", "stage", 
                                                     "sequencing.batch")]

#转变批次信息为因子
embryo_design$sequencing.batch <- as.factor(embryo_design$sequencing.batch) 

embryo_design <- distinct(embryo_design)

rownames(embryo_design) <- embryo_design$sample

embryo_design


# 计算邻域的连通性 ------------------------------------
#在milo对象中储存最近邻域的距离
#这将是整个过程中最耗时的步骤：(需要对格式进行处理，)
embryo_milo <- calcNhoodDistance(embryo_milo,
                                 d=30, 
                                 reduced.dim = "pca.corrected")

# 检测差异丰度 ------------------------------------------------------------------
#测试实验阶段之间的差异，同时考虑到技术批次之间的可变性
da_results <- testNhoods(embryo_milo,
                         design = ~ sequencing.batch + stage, 
                         design.df = embryo_design,
                         reduced.dim = "pca.corrected")

head(da_results)

##计算每个邻域的Fold-change和修正的p值，这表明条件之间是否存在显著的丰度差异。
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

# 检测差异分析的结果 ---------------------------------------------------------------
#1.先检查未校正P值的分布，以验证检验是否平衡。
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

#火山图可视化结果（每个点代表一个邻域）
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

#可视化与单个细胞嵌入相关的DA结果，每个节点代表一个邻域，edges表示两个领域共同都有的细胞，

embryo_milo <- buildNhoodGraph(embryo_milo)

#polt单细胞的umap。
#使用umap降维，颜色由stage编码，数据点文本标签显示celltype，fill=none表示不显示填充的图例
umap_pl <- plotReducedDim(embryo_milo, dimred = "umap",
                          colour_by="stage", 
                          text_by = "celltype", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results, layout="umap",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

#观察DA在某些特定的细胞中明显。

da_results <- annotateNhoods(embryo_milo, da_results, coldata_col = "celltype")

head(da_results)

#如果邻域具有相同的趋势，可以通过设置cell-type_fraction值，将具有混淆的类型排除
#此处函数可以展示cell-type_fraction的分布趋势
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)


da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, 
                              "Mixed", da_results$celltype)

#DA折叠变化在不同细胞类型中的分布
plotDAbeeswarm(da_results, group.by = "celltype")


# 找差异丰度种群的markers ---------------------------------------------------------

#一旦发现具有显著的差异丰度，如果要进行在邻域中找到细胞特有的基因，
#1.首先需要基于生物学问题定义邻域
embryo_milo <- logNormCounts(embryo_milo)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)

da_nhood_markers <- findNhoodGroupMarkers(embryo_milo, da_results, 
                                          subset.row = rownames(embryo_milo)[1:10])

da_nhood_markers <- findNhoodGroupMarkers(embryo_milo, da_results, 
                                          subset.row = rownames(embryo_milo)[1:10], 
                                          aggregate.samples = TRUE,
                                          sample_col = "sample")

head(da_nhood_markers)


##邻域自动分组
#在邻域细胞类型非常不同的情况下，使用社区检测将邻域划分为组，基于(1)两个邻域之间共享单元的数量;
#(2) DA小区褶皱变化方向;(3)折痕变化的差异。
## Run buildNhoodGraph to store nhood adjacency matrix
embryo_milo <- buildNhoodGraph(embryo_milo)

## Find groups
da_results <- groupNhoods(embryo_milo, da_results, max.lfc.delta = 2)
head(da_results)
#检测分组情况：
plotNhoodGroups(embryo_milo, da_results, layout="umap") 

plotDAbeeswarm(da_results, "NhoodGroup")

#试图改变函数的对应参数，检测效果：
plotDAbeeswarm(groupNhoods(embryo_milo, da_results, max.lfc.delta = 0.5), 
               group.by = "NhoodGroup") + ggtitle("max LFC delta=0.5")

plotDAbeeswarm(groupNhoods(embryo_milo, da_results, max.lfc.delta = 1) , 
               group.by = "NhoodGroup") + ggtitle("max LFC delta=1")

plotDAbeeswarm(groupNhoods(embryo_milo, da_results, max.lfc.delta = 2)   , 
               group.by = "NhoodGroup") + ggtitle("max LFC delta=2")


# 寻找邻域基因特征 ----------------------------------------------------------------
#用groupnhood对社区进行了分组，
## Exclude zero counts genes
keep.rows <- rowSums(logcounts(embryo_milo)) != 0
embryo_milo <- embryo_milo[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(embryo_milo)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)

#运行findNhoodGroupMarkers来测试每个邻居群体的一对一差异基因表达
nhood_markers <- findNhoodGroupMarkers(embryo_milo, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, 
                                       sample_col = "sample")

head(nhood_markers)

gr2_markers <- nhood_markers[c("logFC_2", "adj.P.Val_2")] 
colnames(gr2_markers) <- c("logFC", "adj.P.Val")

head(gr2_markers[order(gr2_markers$adj.P.Val), ])


