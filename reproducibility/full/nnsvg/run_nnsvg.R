suppressPackageStartupMessages(library(nnSVG))

for (f in list.files("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/filterexpr/", pattern = ".rds")) {
  print(f)
  
  counts <- readRDS(paste0("svg/crop/data/counts/", f))
  coord <- readRDS(paste0("svg/crop/data/coord/", f))
  cc <- readRDS(paste0("svg/crop/data/cellclu/", f))
  
  cc <- cc[cc %in% names(which(table(cc) >= 1e3))]
  cc <- droplevels(cc)
  
  counts <- counts[, names(cc)]
  coord <- coord[names(cc), ]
  expr <- log2(t(t(counts) / colSums(counts) * 1e4) + 1)
  
  sccell <- split(names(cc), cc)
  
  res_nnsvg <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    # i <- names(cc)[cc == k]
    i <- sccell[[k]] #
    expr_clu <- expr[, i]
    res <- nnSVG(input = expr[rowMeans(expr_clu > 0) > 0.01, i], spatial_coords = coord[i, ],
                 order = "Sum_coords", n_threads = 20)
    res <- as.data.frame(res)
    res <- res[order(res[, "rank"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
  saveRDS(res_nnsvg, file = paste0("svg/crop/full/nnsvg/", f))
}
