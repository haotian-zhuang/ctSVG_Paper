suppressPackageStartupMessages(library(MERINGUE))

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
  set.seed(2025)
  sccell <- lapply(sccell, function(k) {
    idx <- sample(length(k), length(k) * 0.8)
    k[idx]
  })
  
  res_meringue <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    # i <- names(cc)[cc == k]
    i <- sccell[[k]] #
    counts_clu <- counts[, i]
    expr <- normalizeCounts(counts = counts[rowMeans(counts_clu > 0) > 0.01, i], log = FALSE, verbose = FALSE)
    w <- getSpatialNeighbors(pos = coord[i, ])
    res <- getSpatialPatterns(mat = expr, weight = w)
    res <- res[order(res[, "p.value"], -res[, "observed"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
  saveRDS(res_meringue, file = paste0("svg/crop/down/meringue/", f))
}
