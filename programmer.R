# Load Libraries
library(tidyverse)
library(Seurat)
library(tximport)
library(EnsDb.Hsapiens.v86)

# Import Alevin matrix
file_path = file.path("/projectnb/bf528/users/lava_lamp/project_4/daisy_pr4/1_data_curator/salmon_output/alevin/quants_mat.gz")
txi <- tximport(file_path, type="alevin")

# Convert GeneID to Gene Name
txi[["counts"]]@Dimnames[[1]] <- gsub("\\..*","", txi[["counts"]]@Dimnames[[1]])
geneids <- txi[["counts"]]@Dimnames[[1]]
genenames <- ensembldb::select(EnsDb.Hsapiens.v86, keys= geneids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

# Match GeneIDs and Gene Names
geneids_df <- data.frame(geneids)
genes_merged <- merge(x = geneids_df, y = genenames, by.x = "geneids", by.y = "GENEID", all.x = TRUE, na.fill = TRUE)

# Remove rows with NA gene names
ds <- as.data.frame(txi$counts)
ds <- ds[!is.na(match(geneids,genenames$GENEID)),]

# Assign Symbols + Collapse Duplicate Genes
ds$symbol <- genenames$SYMBOL
ds <- ds %>% group_by(symbol) %>% summarise_all(funs(sum))

# Change row names and remove symbol column used in previous step
ds <- as.data.frame(ds)
row.names(ds) <- ds$symbol
ds <- ds[,colnames(ds)!="symbol"]

# Create Seurat Object
cells <- CreateSeuratObject(counts = ds)

# Data QC - Create Mitochondrial Percent
cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")

# Data QC - nFeatures, nCount and PercentMT
png("dataqc_vln1.png")
VlnPlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("dataqc_scatter.png", width = 800)
plot1 <- FeatureScatter(cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

cells_subset <- subset(cells, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)

# Normalize Data
cells_subset <- NormalizeData(cells_subset)

# Identification of highly variable features (feature selection)
cells_subset <- FindVariableFeatures(cells_subset, selection.method = "vst", nfeatures = 2000)

# Identify 10 most highly variable genes
top10_variable_features <- head(VariableFeatures(cells_subset), 10)

# Plot variable features with and without labels
png("feature_selection.png", width = 1000)
plot1 <- VariableFeaturePlot(cells_subset)
plot2 <- LabelPoints(plot = plot1, points = top10_variable_features, repel = TRUE)
plot2
dev.off()

# Scale Data
all.genes <- rownames(cells_subset)
cells_subset <- ScaleData(cells_subset, features = all.genes)

# Perform Linear Dimensional Reduction
cells_subset <- RunPCA(cells_subset, features = VariableFeatures(object = cells_subset))
png("pca.png")
DimPlot(cells_subset, reduction = "pca")
dev.off()

# Determine Dimensionality of Dataset
cells_subset <- JackStraw(cells_subset, num.replicate = 100)
cells_subset <- ScoreJackStraw(cells_subset, dims = 1:20)

png("jackstraw.png")
JackStrawPlot(cells_subset, dims = 1:15)
dev.off()

png("elbow.png")
ElbowPlot(cells_subset)
dev.off()

# Cluster Cells
cells_subset <- FindNeighbors(cells_subset, dims = 1:10)
cells_subset <- FindClusters(cells_subset, resolution = 0.5)

# Run non-linear dimensionality reduction
cells_subset <- RunUMAP(cells_subset, dims = 1:10)
png("umap.png")
DimPlot(cells_subset, reduction = "umap")
dev.off()

# Save RDS
saveRDS(cells_subset, file = "programmer_output.rds")

# Cells Per Cluster
cluster_numbers <- table(cells@meta.data$seurat_clusters)
cluster_numbers <- as.data.frame(cluster_numbers)
cluster_numbers %>% ggplot(mapping = aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity") +
  ggtitle("Distribution of Cells Across Clusters") +
  xlab("Cluster ID") +
  ylab("Number of Cells")
