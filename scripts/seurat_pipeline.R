######################################## Tutoial ########################################
# https://satijalab.org/seurat/articles/pbmc3k_tutorial


######################################## IO ########################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("\n\n\tUsage: Rscript seurat_pipeline.R <pipseeker.dir> <seurat.outdir>\n\n")
}
pipseeker.dir <- args[1]
seurat.outdir <- args[2]
# pipseeker.dir <- "/data/mengxf/Project/KML250418_scRNA/work/250421/pipseeker/HEK-3T3-Human-Mouse-Mixture-T2_S3_results/raw_matrix"
# seurat.outdir <- "/data/mengxf/Project/KML250418_scRNA/work/250421/seurat"

dir.create(seurat.outdir, showWarnings = FALSE, recursive = TRUE)
mtrcs.vlnplot.before <- file.path(seurat.outdir, "qc_vlnplot_filtered_before.png")
mtrcs.vlnplot.after <- file.path(seurat.outdir, "qc_vlnplot_filtered_after.png")
mtrcs.rltn.scatter.before <- file.path(seurat.outdir, "qc_metrics_relationships_before.png")
mtrcs.rltn.scatter.after <- file.path(seurat.outdir, "qc_metrics_relationships_after.png")
vrbl.ftr.figure <- file.path(seurat.outdir, "slct_variable_feature.png")
pca.dim.scatter <- file.path(seurat.outdir, "pca_dimensional_reduction_genes_scatter.png")
pca.dim.heatmap <- file.path(seurat.outdir, "pca_dimensional_reduction_genes_heatmap.png")
pca.dim.elbowplot <- file.path(seurat.outdir, "pca_relevant_dimensions_elbow.png")
pca.dim.rdctn.figure <- file.path(seurat.outdir, "pca_dimensional_reduction_plot.png")
umap.dim.rdctn.figure <- file.path(seurat.outdir, "umap_dimensional_reduction_plot.png")
tsne.dim.rdctn.figure <- file.path(seurat.outdir, "tsne_dimensional_reduction_plot.png")
sig.gene.vlnplot <- file.path(seurat.outdir, "mrkr_significant_genes_vlnplot.png")
mrkr.on.umap.figure <- file.path(seurat.outdir, "mrkr_on_umap_plot.png")
mrkr.clstr.heatmap <- file.path(seurat.outdir, "mrkr_cluster_heatmap.png")
seurat.object.rds <- file.path(seurat.outdir, "seurat_object.rds")
pca.detail.file <- file.path(seurat.outdir, "pca_detail.txt")
mrkr.tbl.file <- file.path(seurat.outdir, "mrkr_detail.tsv")


######################################## Library ########################################
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


######################################## Function ########################################
load.dataset <- function(data.dir) {
  # 加载 PBMC 数据集
  pbmc.data <- Read10X(data.dir = data.dir)
  # 使用原始(非标准化)数据初始化 Seurat 对象
  # * project 不同样本可以不同命名, 后期可以合并
  data <- CreateSeuratObject(counts = pbmc.data, project = "SeuratProject", min.cells = 3, min.features = 200)
  # [[操作符可以向对象 metadata 添加列, 添加线粒体百分比
  # * 可以添加核糖体 RNA 和红细胞等
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  return(data)
}

plot.qc.metrics <- function(outfile, data) {
  p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(outfile, p, width = 12, height = 12, dpi = 300)
}

plot.qc.relationship <- function(data, outfile) {
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(outfile, plot1 + plot2, width = 12, height = 12, dpi = 300)
}

normlize.select.feature <- function(data) {
  # 规范化数据, 校正至每个细胞10000X
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  # 特征选择, 选择变异系数较大的前 2000 个基因
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  return(data)
}

plot.variable.feature <- function(data, outfile) {
  # 标注差异最大的 10 个基因
  top10 <- head(VariableFeatures(data), 10)
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(data)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(outfile, plot1 + plot2, width = 12, height = 12, dpi = 300)
}

run.pca <- function(data) {
  # PCA 前线性转换(缩放)数据
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  # 运行 PCA
  data <- RunPCA(data, features = VariableFeatures(object = data))
  return(data)
}

print.pca.detail <- function(data, pcafile) {
  # 前10 个 PCA 结果明细
  sink(pcafile)
  print(data[["pca"]], dims = 1:5, nfeatures = 5)
  sink()
}

plot.pca <- function(data, scttrplt, htmplt, elbowplt) {
  # PCA 贡献度肘图
  p <- ElbowPlot(data)
  ggsave(elbowplt, p, width = 8, height = 8, dpi = 300, bg = "white")
  # PCA 变异来源基因点图
  p <- VizDimLoadings(data, dims = 1:2, reduction = "pca")
  ggsave(scttrplt, p, width = 12, height = 12, dpi = 300)
  # PCA 变异来源基因热图
  png(htmplt, width = 12, height = 12, res = 300, units = "in")
  DimHeatmap(data, dims = 1:9, cells = 500, balanced = TRUE)
  dev.off()
}

run.clstr <- function(data) {
  data <- FindNeighbors(data, dims = 1:10)
  data <- FindClusters(data, resolution = 0.5)
  data <- RunUMAP(data, dims = 1:10)
  data <- RunTSNE(data, dims = 1:10)
  return(data)
}

plot.clstr.dim <- function(data, pcaplot, umaplot, tsneplot) {
  # 聚类后 PCA 图
  p <- DimPlot(data, reduction = "pca")
  ggsave(pcaplot, p, width = 12, height = 12, dpi = 300)
  # UMAP
  p <- DimPlot(data, reduction = "umap")
  ggsave(umaplot, p, width = 12, height = 12, dpi = 300)
  # tSNE
  p <- DimPlot(data, reduction = "tsne")
  ggsave(tsneplot, p, width = 12, height = 12, dpi = 300)
}

plot.mrkr.vln.umap <- function(mrkrs, data, vlnplot, mrkr.umap.plot) {
  # 差异基因在不同细胞分类中的表达小提琴图
  # 每个 cluster 选第 1 个
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    pull(gene) -> sig.genes
  # 显著基因在各 cluster 间表达水平小提琴图
  p <- VlnPlot(data, features = sig.genes)
  ggsave(sig.gene.vlnplot, p, width = 12, height = 12, dpi = 300)
  # 显著基因在 UMAP 聚类图上的位置
  p <- FeaturePlot(data, features = sig.genes)
  ggsave(mrkr.on.umap.figure, p, width = 12, height = 12, dpi = 300)
}

plot.mrkr.heatmap <- function(mrkrs, data, htmplt) {
  # 显著基因与 cluster 热图
  mrkrs %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  p <- DoHeatmap(data, features = top10$gene)
  ggsave(htmplt, p, width = 12, height = 12, dpi = 300)
}


######################################## Main ########################################
# 加载数据
pbmc <- load.dataset(pipseeker.dir)

# 质控流程
# 过滤前, QC 参数小提琴图
plot.qc.metrics(mtrcs.vlnplot.before, pbmc)
# 过滤前, QC 参数相关性图
plot.qc.relationship(pbmc, mtrcs.rltn.scatter.before)

# 细胞过滤
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 &
  nFeature_RNA < 8000 & percent.mt < 25)
# 过滤后, QC 参数小提琴图
plot.qc.metrics(mtrcs.vlnplot.after, pbmc)
# 过滤后, QC 参数相关性图
plot.qc.relationship(pbmc, mtrcs.rltn.scatter.after)

# 标准化数据 + 特征选择
pbmc <- normlize.select.feature(pbmc)
# 变量特征图
plot.variable.feature(pbmc, vrbl.ftr.figure)

# PCA
pbmc <- run.pca(pbmc)
plot.pca(pbmc, pca.dim.scatter, pca.dim.heatmap, pca.dim.elbowplot)

# 聚类 + UMAP + tSNE
pbmc <- run.clstr(pbmc)
plot.clstr.dim(pbmc, pca.dim.rdctn.figure, umap.dim.rdctn.figure, tsne.dim.rdctn.figure)

# 差异表达分析
# 使用所有 cluster 查找 Markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 1)
write.table(
  pbmc.markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster")],
  mrkr.tbl.file,
  sep = "\t",
  row.names = F
)
plot.mrkr.vln.umap(pbmc.markers, pbmc, sig.gene.vlnplot, mrkr.on.umap.figure)
plot.mrkr.heatmap(pbmc.markers, pbmc, mrkr.clstr.heatmap)

# 保存数据
# PCA
print.pca.detail(pbmc, pca.detail.file)
# Seurat 对象
saveRDS(pbmc, file = seurat.object.rds)
