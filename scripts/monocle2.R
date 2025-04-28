# 细胞轨迹
# https://mp.weixin.qq.com/s/wzn5Jdf5DR1LrhXecYg7Yw?mpshare=1&scene=1&srcid=0425BCnQ1vJqJAayIiEEkHbP&sharer_shareinfo=3bd9163fcb7dcd454867dc1e0fe016df&sharer_shareinfo_first=3bd9163fcb7dcd454867dc1e0fe016df&from=industrynews&color_scheme=light#rd
# https://cole-trapnell-lab.github.io/monocle-release/docs/


######################################## IO ########################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("\n\n\tUsage: Rscript monocle2.R <threads> <seurat.singler.object.rds> <outdir>\n\n")
}
threads <- as.numeric(args[1])
seurat.singler.object <- args[2]
outdir <- args[3]

# threads <- as.numeric("16")
# seurat.singler.object <- "/data/mengxf/Project/KML250418_scRNA/work/250421/singleR/seurat_singler_object.rds"
# outdir <- "/data/mengxf/Project/KML250418_scRNA/work/250421/monocle2"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
clstr.pseudo.fig <- file.path(outdir, "cluster_pseudotime.png")
genes.jitter.fig <- file.path(outdir, "genes_jitter.png")
pseudo.heatmap.fig <- file.path(outdir, "pseudotime_heatmap.png")


######################################## Library ########################################
library(monocle)
library(Seurat)
library(ggplot2)


######################################## Main ########################################
# 读取 singleR 注释获得 Seurat 对象
sce <- readRDS(seurat.singler.object)
Idents(sce) <- sce$celltype

# 数据准备
# expr.matrix
mono.tj <- sce
mono.mtrx <- as(as.matrix(GetAssayData(mono.tj, slot = "counts")), "sparseMatrix")
# featuredata
feature.ann <- data.frame(gene_id = rownames(mono.mtrx), gene_short_name = rownames(mono.mtrx))
rownames(feature.ann) <- rownames(mono.mtrx)
mono.fd <- new("AnnotatedDataFrame", data = feature.ann)
# phenodata
sample.ann <- mono.tj@meta.data
rownames(sample.ann) <- colnames(mono.mtrx)
mono.pd <- new("AnnotatedDataFrame", data = sample.ann)

# 开始分析
# 构建 CellDataSet 对象
# ! lowerDetectionLimit 参数提速, 默认 0.1, 增加数值可以加速
mono.cds <- newCellDataSet(
    mono.mtrx,
    phenoData = mono.pd,
    featureData = mono.fd,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 0.1
)
# 查看 phenodata、featuredat
# head(pData(mono.cds))
# head(fData(mono.cds))
# 使用 估计尺度因子+估计离散度 预处理数据
mono.cds <- estimateSizeFactors(mono.cds)
mono.cds <- estimateDispersions(mono.cds)
# 筛选基因
disp.tbl <- dispersionTable(mono.cds)
# ! mean_expression 提速 默认 0.1, 增加数值可以加速
unsup_clustering_genes <- subset(disp.tbl, mean_expression >= 0.1)
mono.cds <- setOrderingFilter(mono.cds, unsup_clustering_genes$gene_id)
# DDRtree 降维
mono.cds <- reduceDimension(
    mono.cds,
    max_components = 2,
    method = "DDRTree"
)

# 计算 pseudotime 值
mono.cds <- orderCells(mono.cds)
# head(pData(mono.cds))

# 高变基因
disp.genes <- subset(disp.tbl, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff.test <- differentialGeneTest(
    mono.cds[disp.genes, ],
    cores = threads,
    fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
sig.gene.names <- row.names(subset(diff.test, qval < 1e-50))

# 可视化
# Cluster/Pseudotime 轨迹分布图
p1 <- plot_cell_trajectory(mono.cds, color_by = "celltype", size = 1, show_backbone = TRUE)
p2 <- plot_cell_trajectory(mono.cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
ggsave(clstr.pseudo.fig, p1 + p2, width = 12, height = 12, dpi = 300)

# 基因抖动图
p <- plot_genes_jitter(mono.cds[sig.gene.names[1:5], ], grouping = "State", color_by = "State")
ggsave(genes.jitter.fig, p, width = 12, height = 12, dpi = 300)

# 拟时序相关基因热图
png(pseudo.heatmap.fig, width = 9, height = 18, res = 300, units = "in")

plot_pseudotime_heatmap(
    mono.cds[sig.gene.names, ],
    num_clusters = 4,
    show_rownames = T, return_heatmap = T
)
dev.off()
