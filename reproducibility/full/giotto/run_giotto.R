suppressPackageStartupMessages(library(Giotto))
suppressPackageStartupMessages(library(Matrix))

for (f in list.files("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/filterexpr/", pattern = ".rds")) {
  print(f)
  
  counts <- readRDS(paste0("svg/crop/data/counts/", f))
  coord <- readRDS(paste0("svg/crop/data/coord/", f))
  cc <- readRDS(paste0("svg/crop/data/cellclu/", f))

  cc <- cc[cc %in% names(which(table(cc) >= 1e3))]
  cc <- droplevels(cc)
  
  counts <- counts[, names(cc)]
  coord <- coord[names(cc), ]
  
  sccell <- split(names(cc), cc)
  
  res_giotto <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    # i <- names(cc)[cc == k]
    i <- sccell[[k]] #
    counts_clu <- counts[, i]
    obj <- createGiottoObject(raw_exprs = counts[rowMeans(counts_clu > 0) > 0.01, i], spatial_locs = coord[i, ])
    obj <- normalizeGiotto(gobject = obj)
    obj <- createSpatialNetwork(gobject = obj)
    res <- binSpect(obj, bin_method = "kmeans", do_parallel = FALSE, verbose = FALSE)
    res <- as.data.frame(res)
    res <- res[order(res[, "p.value"], -res[, "score"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    data.frame(cluster = k, gene = res$genes, res)
  }, simplify = FALSE))
  saveRDS(res_giotto, file = paste0("svg/crop/full/giotto/", f))
}
