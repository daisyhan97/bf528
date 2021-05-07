# Load Packages
library(tidyverse)
library(Seurat)
set.seed(3423)

# Load Data
cells <- readRDS("/projectnb2/bf528/users/lava_lamp/project_4/daisy_pr4/2_programmer/programmer_output.rds")

# Identify cluster biomarkers
cells_markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- cells_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Marker genes used in original analysis
marker_genes <- read_csv("MarkerGenes.csv")
allgenes <- unlist(strsplit(marker_genes$Genes, ","))

# Visualize Clusters, by Cluster and by Biomarker
png("umap.png")
DimPlot(cells, reduction = "umap")
dev.off()
png("features.png", width = 1000, height = 1000)
FeaturePlot(cells, features = allgenes)
dev.off()

# Heatmap of Clusters
png("heatmap_clusters.png", height = 800, width = 1000)
DoHeatmap(cells, features = top_markers$gene) + NoLegend()
dev.off()
  # Notable Pairs: (0,1), (3,6), (4,8)

# Violin Plots
VlnPlot(cells, features = c("GCG")) # Alpha
VlnPlot(cells, features = c("INS")) # Beta
VlnPlot(cells, features = c("SST")) # Delta
VlnPlot(cells, features = c("KRT19")) # Ductal 
VlnPlot(cells, features = c("PPY")) # Gamma
VlnPlot(cells, features = c("CPA1")) #Acinar
VlnPlot(cells, features = c("PDGFRB")) # Stellate
VlnPlot(cells, features = c("VWF", "PECAM1","CD34")) # Vascular
VlnPlot(cells, features = c("SDS", "CD163")) # Macrophage

# Cell Assignments
# 0, 1 - Alpha (GCG)
# 2 - Gamma (PPY)
# 3, 6 - Beta (INS)
# 4, 8 - Acinar (CPA1)
# 5 - Ductal (KRT19)
# 7 - Stellate (PDGFRB)
# 9 - Vascular (VWF)
# 10 - Macrophage (SDS)

# Assign Cell Type to Clusters
new_cluster_ids <- c("Alpha", "Alpha", "Delta", "Beta", "Acinar", "Ductal", "Beta", "Stellate", "Acinar", "Vascular", "Macrophage")
names(new_cluster_ids) <- levels(cells)
cells <- RenameIdents(cells, new_cluster_ids)

# Clustered UMAP
png("umap_celltype.png")
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

# Clustered Heatmap
png("heatmap_celltype.png", height = 1000, width = 1200)
DoHeatmap(cells, features = top_markers$gene) + NoLegend()
dev.off

# Novel Marker Genes
novel_markers <- cells_markers %>% filter(!gene %in% allgenes)
top_novel_markers <- novel_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
png("heatmap_novel.png", height = 1000, width = 1200)
DoHeatmap(cells, features = top_novel_markers$gene) + NoLegend()
dev.off()

# Recreate Heatmap
heatmap_genes <- c("GCG", "INS", "PPY", "SST", "GHRL", "PRSS1", "KRT19", "SPARC", "VWF", "RGS5", "PDGFRA", "SOX10", "SDS", "TPSAB1", "TRAC")
png("heatmap_recreation.png", height = 800, width = 1000)
DoHeatmap(cells, features = heatmap_genes) + NoLegend()
dev.off()

# Save Seurat Object
saveRDS(cells, file = "analyst_output.rds")

# Save Novel Genes
write_csv(novel_markers, "marker_genes.csv")