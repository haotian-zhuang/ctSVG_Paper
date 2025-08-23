suppressPackageStartupMessages(library(Giotto))
suppressPackageStartupMessages(library(Matrix))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    
    obj <- createGiottoObject(raw_exprs = counts, spatial_locs = coord)
    obj <- normalizeGiotto(gobject = obj)
    obj <- createSpatialNetwork(gobject = obj)
    res <- binSpect(obj, bin_method = "kmeans", do_parallel = TRUE, cores = 2, verbose = FALSE)
    res <- as.data.frame(res)
    res <- res[order(res[, "p.value"], -res[, "score"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    
    res_giotto <- do.call(rbind, sapply(levels(cc), function(k) {
      data.frame(cluster = k, gene = res$genes, res)
    }, simplify = FALSE))
    
    saveRDS(res_giotto, file = paste0("svg/simu/", f, "/", signal, "/method/giotto_global/res.rds"))
  }
}
