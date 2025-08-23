suppressPackageStartupMessages(library(SPARK))
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
  counts <- as(counts, "sparseMatrix")
  
  sccell <- split(names(cc), cc)
  
  res_sparkg <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    # i <- names(cc)[cc == k]
    i <- sccell[[k]] #
    counts_clu <- counts[, i]
    res <- sparkg_sw(counts = counts[rowMeans(counts_clu > 0) > 0.01, i], location = coord[i, ])
    res <- res[["res_mtest"]]
    res <- res[order(res[, "combinedPval"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "combinedPval"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
  saveRDS(res_sparkg, file = paste0("svg/crop/full/sparkg/", f))
}
