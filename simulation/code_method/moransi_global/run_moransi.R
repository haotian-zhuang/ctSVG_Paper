suppressPackageStartupMessages(library(Seurat))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    expr <- log1p(t(t(counts) / colSums(counts) * 1e4))
    
    gene_split <- split(rownames(expr), sample(1:10, nrow(expr), replace = TRUE))
    names(gene_split) <- NULL
    res <- do.call(rbind, parallel::mclapply(gene_split, function(gs) {
      RunMoransI(data = expr[gs, ], pos = coord)
    }, mc.cores = 10))
    
    # res <- RunMoransI(data = expr, pos = coord)
    res <- res[order(res[, "p.value"], -res[, "observed"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    
    res_moransi <- do.call(rbind, sapply(levels(cc), function(k) {
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
    
    saveRDS(res_moransi, file = paste0("svg/simu/", f, "/", signal, "/method/moransi_global/res.rds"))
  }
}
