suppressPackageStartupMessages(library(Seurat))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
expr <- log1p(t(t(counts) / colSums(counts) * 1e4))

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_moransi <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    i <- names(cc)[cc == k]
    expr_clu <- expr[, i]
    res <- RunMoransI(data = expr[rowMeans(expr_clu > 0) > 0.01, i], pos = coord[i, ])
    res <- res[order(res[, "p.value"], -res[, "observed"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_moransi, file = "svg/scdesign3/method/moransi/res.rds")
