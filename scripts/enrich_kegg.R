######################################## IO ########################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("\n\n\tUsage: Rscript enrich_kegg.R <markers.detail.table> <outdir>\n\n")
}
genes.tbl.file <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


######################################## library ########################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# ! 不要动环境, rstudio-server 导入会报错, 但是在命令行里面不会报错
library(pathview)


######################################## Function ########################################
# 生成 geneFC 对象
generate.geneFC <- function(gntbl) {
    data <- read.delim(genes.tbl.file)
    # 不去重会报错, 不同 cluster 会有很多重复基因
    data <- data[!duplicated(data$gene), ]
    rownames(data) <- data$gene
    cnvrt.ids <- bitr(data$gene,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db,
    )
    geneFC <- data[cnvrt.ids$SYMBOL, "avg_log2FC"]
    names(geneFC) <- cnvrt.ids$ENTREZID
    return(geneFC)
}


# 进行 KEGG 富集分析
analyza.enrich.kegg <- function(entrezids, outdir) {
    # ! rstudio-server 跑不了，命令行可以
    enrich.kegg.res <- enrichKEGG(
        gene = entrezids,
        organism = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    write.table(
        enrich.kegg.res@result,
        file = file.path(outdir, "detail.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
    )
    return(enrich.kegg.res)
}


# 绘制 KEGG 富集分析结果
kegg.plot <- function(enrich.kegg.res, outdir) {
    # 绘制条形图
    png(
        file = file.path(outdir, "barplot.png"),
        width = 12,
        height = 18,
        units = "in",
        res = 300
    )
    print(barplot(
        enrich.kegg.res,
        drop = TRUE,
        showCategory = 20,
        main = "30min"
    ))
    dev.off()
    # 绘制点图
    png(
        file = file.path(outdir, "bubble.png"),
        width = 12,
        height = 18,
        units = "in",
        res = 300
    )
    print(dotplot(enrich.kegg.res))
    dev.off()
}


# 绘制 KEGG 通路图
kegg.pathway.plot <- function(geneFC, pwid, outdir) {
    # 绘制通路图, 画富集最高的通路
    wd0 <- getwd()
    setwd(outdir)
    pathview(
        gene.data = geneFC,
        pathway.id = pwid,
        species = "hsa",
        out.suffix = pwid
    )
    setwd(wd0)
}


######################################## main ########################################
geneFC <- generate.geneFC(genes.tbl.file)
enrich.kegg.res <- analyza.enrich.kegg(names(geneFC), outdir)
# 没有显著的富集通路就跳过
if (dim(enrich.kegg.res)[1] == 0) next
kegg.plot(enrich.kegg.res, outdir)
kegg.pathway.plot(geneFC, enrich.kegg.res$ID[1], outdir)
