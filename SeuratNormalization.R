# install seurat if needed
#install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("multtest")

library(Seurat)
library(ggplot2)
library(sctransform)
library("Matrix")
library(monocle)

#redo data import

## load data
Data <- Read10X(data.dir = "filtered_feature_bc_matrix");  
pbmc <- CreateSeuratObject(counts = Data$`Gene Expression`)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt") # remove the 'novel/unmapped transcript' and mitochrondia transcripts
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 20000)

# export sample ID (required)
sampleID =  pbmc@assays$RNA@data@Dimnames
sampleID = sampleID[[2]]
df = data.frame(sampleID)
write.csv(df,'sampleID.csv', row.names = FALSE)

# ---------------- processing ----------------------------------
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))


##-------------- output normalized expression ---------------------
writeMM(pbmc$RNA@data, "NormalizedExpression.txt")
writeMM(pbmc$RNA@counts, "AdjustedCount.txt")


## ----------------- write the filtered gene list --------------------------
geneList = pbmc@assays$RNA@data@Dimnames
geneList = geneList[[1]]
write.csv(geneList, "filterGeneList.csv")

#-------------------------ATAC peaks------------------------------------
#-----------------------------------------------------------------------------
pbmc <- CreateSeuratObject(counts = Data$Peaks, assay="ATAC")

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nCount_ATAC>5000)

# export sample ID (required)
sampleID =  pbmc@assays$ATAC@data@Dimnames
sampleID = sampleID[[2]]
df = data.frame(sampleID)
write.csv(df,'Peak_sampleID.csv', row.names = FALSE)

# ---------------- processing ----------------------------------
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

## normalization
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunLSI(pbmc, n=50, scale.max = NULL)


##-------------- output normalized expression ---------------------
writeMM(pbmc$ATAC@data, "PeaK_NormalizedExpression.txt")
writeMM(pbmc$ATAC@counts, "Peak_AdjustedCount.txt")


## ----------------- write the filtered gene list --------------------------
geneList = pbmc@assays$ATAC@data@Dimnames
geneList = geneList[[1]]
write.csv(geneList, "peak_filterGeneList.csv")


