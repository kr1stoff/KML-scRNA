######################################## IO ########################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("\n\n\tUsage: Rscript enrich_GO.R <markers.detail.table> <outdir>\n\n")
}
genes.tbl.file <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


######################################## library ########################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


######################################## Function ########################################
# enrichGO 可视化
gp.plot <- function(enrich.go.res, outdir) {
    # 绘制条形图
    png(
        file = file.path(outdir, "barplot.png"),
        width = 12,
        height = 18,
        units = "in",
        res = 300
    )
    print(
        barplot(
            enrich.go.res,
            drop = TRUE,
            showCategory = 10,
            split = "ONTOLOGY"
        ) + facet_grid(ONTOLOGY ~ ., scale = "free")
    )
    dev.off()

    # 绘制点图
    png(
        file = file.path(outdir, "bubble.png"),
        width = 12,
        height = 18,
        units = "in",
        res = 300
    )
    print(
        dotplot(enrich.go.res, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
    )
    dev.off()
}

# enrichGO 富集分析
analyza_enrich_go <- function(gntbl, outdir) {
    data <- read.delim(gntbl)
    enrich.go.res <- enrichGO(
        gene = data$gene,
        keyType = "SYMBOL",
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = FALSE
    )
    write.table(
        enrich.go.res,
        file = file.path(outdir, "detail.tsv"),
        sep = "\t",
        quote = F,
        row.names = F
    )
    return(enrich.go.res)
}


######################################## main ########################################
enrich.go.res <- analyza_enrich_go(genes.tbl.file, outdir)
tempres <- gp.plot(enrich.go.res, outdir)
