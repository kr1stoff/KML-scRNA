WORKDIR="/data/mengxf/Project/KML250418_scRNA/work/250421"
THREADS=32

# 质控
mkdir -p ${WORKDIR}/fastqc
mamba -n basic run fastqc -t $THREADS -f fastq -o ${WORKDIR}/fastqc \
    /data/rawdata/supplier/ilmn_250418_scRNA_demo/illumina_single_cell_3RNA_prep/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R1_001.fastq.gz \
    /data/rawdata/supplier/ilmn_250418_scRNA_demo/illumina_single_cell_3RNA_prep/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R2_001.fastq.gz &
# GC, Q20, Q30
mkdir ${WORKDIR}/fastq_count
/data/mengxf/Software/GitHub/fastq_count/fastq_count \
    -output ${WORKDIR}/fastq_count/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R1_001.tsv \
    /data/rawdata/supplier/ilmn_250418_scRNA_demo/illumina_single_cell_3RNA_prep/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R1_001.fastq.gz &
/data/mengxf/Software/GitHub/fastq_count/fastq_count \
    -output ${WORKDIR}/fastq_count/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R2_001.tsv \
    /data/rawdata/supplier/ilmn_250418_scRNA_demo/illumina_single_cell_3RNA_prep/HEK-3T3-Human-Mouse-Mixture-T2_S3_L001_R2_001.fastq.gz &

# 二级分析
# L00[12] 2 个 lane 的数据，pipseeker 也可以识别
mkdir ${WORKDIR}/pipseeker
/data/mengxf/Software/pipseeker-v3.3.0-linux/pipseeker full \
    --chemistry V \
    --fastq /data/rawdata/supplier/ilmn_250418_scRNA_demo/illumina_single_cell_3RNA_prep/HEK-3T3-Human-Mouse-Mixture-T2_S3 \
    --star-index-path /data/mengxf/Database/scRNA/fluentbio/pipseeker-gex-reference-GRCh38-2022.04 \
    --annotation /data/mengxf/Database/scRNA/fluentbio/human-pbmc-references/references/human-pbmc-v4.csv \
    --output-path ${WORKDIR}/pipseeker \
    --threads $THREADS

# 三级分析
# Seurat: 细胞质控 + 过滤 + 标准化 + 特征选择 + PCA/UMAP/tSNE 降维 + 差异基因 (marker)
mamba -n scrna run Rscript /data/mengxf/Project/KML250418_scRNA/scripts/seurat_pipeline.R \
    ${WORKDIR}/pipseeker/raw_matrix \
    ${WORKDIR}/seurat

# SingleR: 细胞注释
mkdir -p ${WORKDIR}/singleR
mamba -n scrna run Rscript /data/mengxf/Project/KML250418_scRNA/scripts/singleR.R \
    $THREADS \
    /data/mengxf/Project/KML250418_scRNA/work/250421/seurat/seurat_object.rds \
    /data/mengxf/Database/scRNA/singleR/HumanPrimaryCellAtlasData.rds \
    ${WORKDIR}/singleR

# Monocle2: 细胞轨迹
mamba -n scrna run Rscript /data/mengxf/Project/KML250418_scRNA/scripts/monocle2.R \\
$THREADS \
    ${WORKDIR}/singleR/seurat_singler_object.rds \
    ${WORKDIR}/monocle2

# clusterProfiler: 富集分析
# GO
mamba -n R.rnaseq run Rscript /data/mengxf/Project/KML250418_scRNA/scripts/enrich_go.R \
    ${WORKDIR}/seurat/mrkr_detail.tsv \
    ${WORKDIR}/enrich/GO
# KEGG
mamba -n R.rnaseq run Rscript /data/mengxf/Project/KML250418_scRNA/scripts/enrich_kegg.R \
    ${WORKDIR}/seurat/mrkr_detail.tsv \
    ${WORKDIR}/enrich/KEGG
