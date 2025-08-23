suppressPackageStartupMessages(library(MERINGUE))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_meringue <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    i <- names(cc)[cc == k]
    counts_clu <- counts[, i]
    expr <- normalizeCounts(counts = counts[rowMeans(counts_clu > 0) > 0.01, i], log = FALSE, verbose = FALSE)
    w <- getSpatialNeighbors(pos = coord[i, ])
    res <- getSpatialPatterns(mat = expr, weight = w)
    res <- res[order(res[, "p.value"], -res[, "observed"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_meringue, file = "svg/scdesign3/method/meringue/res.rds")
