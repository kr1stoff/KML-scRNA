# SingleR 细胞类型注释
# https://lishensuo.github.io/posts/bioinfo/019%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--singler%E7%BB%86%E8%83%9E%E7%B1%BB%E5%9E%8B%E6%B3%A8%E9%87%8A/


######################################## IO ########################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("\n\n\tUsage: Rscript singleR <threads> <pbmc.rds> <singler.ref.rds> <outdir>\n\n")
}
threads <- as.numeric(args[1])
pbmc.rds <- args[2]
singler.ref.rds <- args[3]
outdir <- args[4]

# threads <- 16
# pbmc.rds <- "/data/mengxf/Project/KML250418_scRNA/work/250421/seurat/seurat_object.rds"
# # * 无需联网下载. 保存注释文件然后导入
# # ref <- HumanPrimaryCellAtlasData()
# # saveRDS(ref, file = "/data/mengxf/Database/scRNA/singleR/HumanPrimaryCellAtlasData.rds")
# singler.ref.rds <- "/data/mengxf/Database/scRNA/singleR/HumanPrimaryCellAtlasData.rds"
# outdir <- "/data/mengxf/Project/KML250418_scRNA/work/250421/singleR"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
umap.celltype.figure <- file.path(outdir, "umap_singler_celltype.png")
score.heatmap.figure <- file.path(outdir, "score_heatmap.png")
seurat.singler.object.rds <- file.path(outdir, "seurat_singler_object.rds")


######################################## Library ########################################
library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)


######################################## Main ########################################
# 加载 Seurat 对象
pbmc <- readRDS(pbmc.rds)
# 获取标准化后的数据矩阵
norm.count <- GetAssayData(pbmc, layer = "data")
# 导入 SingleR 参考数据集
ref <- readRDS(singler.ref.rds)
# 运行 SingleR 预测细胞类型
pred <- SingleR(
    test = norm.count,
    ref = ref,
    clusters = pbmc$seurat_clusters,
    labels = ref$label.main,
    num.threads = threads
)
# 将预测结果添加到 Seurat 对象中
pbmc$celltype <- pred$labels[match(pbmc$seurat_clusters, rownames(pred))]

# 输出注释细胞类型后的 UMAP 图
p <- DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave(umap.celltype.figure, p, width = 12, height = 12, dpi = 300)

# 输出注释细胞类型后的热图
png(score.heatmap.figure, width = 12, height = 12, res = 300, units = "in")
plotScoreHeatmap(pred)
dev.off()

# 保存 Seurat 对象
saveRDS(pbmc, file = seurat.singler.object.rds)
