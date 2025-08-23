suppressPackageStartupMessages(library(nnSVG))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    expr <- log2(t(t(counts) / colSums(counts) * 1e4) + 1)
    
    res <- nnSVG(input = expr, spatial_coords = coord,
                 n_threads = 20)
    res <- as.data.frame(res)
    res <- res[order(res[, "rank"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
    
    res_nnsvg <- do.call(rbind, sapply(levels(cc), function(k) {
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
    
    saveRDS(res_nnsvg, file = paste0("svg/simu/", f, "/", signal, "/method/nnsvg_global/res.rds"))
  }
}
