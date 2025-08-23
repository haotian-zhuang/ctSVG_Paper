suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
f <- "human_pancreas.rds"
counts <- readRDS(paste0("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/filterexpr/", f))
coord <- readRDS(paste0("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/centroid/", f))
cc <- readRDS(paste0("cellclu/", f))
if (length(cc) > 1e4) {
  set.seed(2025)
  cc <- sample(cc, 1e4)
}

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]

npcs = 10
resolution = 0.1
d <- CreateSeuratObject(counts)
d <- NormalizeData(d, verbose = FALSE)
d <- FindVariableFeatures(d, verbose = FALSE)
d <- ScaleData(d, verbose = FALSE)
d <- RunPCA(d, verbose = FALSE, npcs = npcs)
d <- FindNeighbors(d, verbose = FALSE, reduction = "pca", dims = 1:npcs)
d <- FindClusters(d, verbose = FALSE, resolution = resolution, random.seed = 0)
d <- RunUMAP(d, verbose = FALSE, reduction = "pca", dims = 1:npcs)

df <- data.frame(x = coord[, 1], y = coord[, 2], cluster = d$seurat_clusters)
ggplot(data = df, aes(x = x, y = y)) +
  geom_point(aes(color = cluster)) +
  theme_minimal() +
  coord_fixed() +
  theme(legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave(filename = paste0("svg/simu/", sub(".rds", "/simdata/spatial_cluster.jpg", f)),
       dpi = 500, width = 5.5, height = 4)

saveRDS(counts, file = paste0("svg/simu/", sub(".rds", "/simdata/counts.rds", f)))
saveRDS(coord, file = paste0("svg/simu/", sub(".rds", "/simdata/coord.rds", f)))
saveRDS(d, file = paste0("svg/simu/", sub(".rds", "/simdata/seurat.rds", f)))
saveRDS(d$seurat_clusters, file = paste0("svg/simu/", sub(".rds", "/simdata/cellclu.rds", f)))
