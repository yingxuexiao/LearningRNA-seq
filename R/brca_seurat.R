library(Seurat)

#改变工作的环境只包含三个10X文件
getwd()
data_dir <- "C:/Users/Administrator/Documents/GitHub/LearningRNA-seq/BrCa"

#检查该路径是否只有三个文件
list.files(data_dir)


# 利用Seurat正式读取文件 ----------------------------------------------------------
##Read10X命令读取三个文件，得到一个带行名（基因名）及列名（细胞名）的count的矩阵
count_matrix <- Read10X(data.dir = data_dir)

# 创建Seurat对象
seurat_object = CreateSeuratObject(counts = count_matrix,
                                   min.cells = 3, 
                                   min.features = 200, 
                                   project = "seurat_object")

seurat_object = CreateSeuratObject(counts = count_matrix,
                                   min.cells = 0, 
                                   min.features = 0, 
                                   project = "brca")
brca <- seurat_object
brca[["percent.mt"]] <- PercentageFeatureSet(brca, pattern = "^MT-")

VlnPlot(brca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##过滤
brca <- subset(brca, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
brca

#标准化
brca <- NormalizeData(brca, normalization.method = "LogNormalize", scale.factor = 1e4)


#识别高变基因
brca <- FindVariableFeatures(brca, selection.method = "vst", nfeatures = 2000) #返回两千个高变基因

#提取前10的高变基因
top10 <- head(VariableFeatures(brca), 10)
top10

展示高变基因
plot1 <- VariableFeaturePlot(brca)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#归一化
brca <- ScaleData(brca, vars.to.regress = "percent.mt")

#降维
brca <- RunPCA(brca, features = VariableFeatures(object = brca)) ##默认会输出5个主成分
brca <- FindNeighbors(brca, dims = 1:10)
brca <- FindClusters(brca, resolution = 0.5)

#查看前5分析细胞聚类数ID
head(Idents(brca), 5)

table(brca@meta.data$seurat_clusters)

#可视化
brca <- RunUMAP(brca, dims = 1:10)

p1 <- DimPlot(brca, reduction = "umap", label = T, label.size = 5)

brca <- RunTSNE(brca, dims = 1:10)

p2 <- DimPlot(brca, reduction = "tsne", label = T, label.size = 5)
p1+p2
