suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

counts <- readRDS("svg/sim_mouse_brain_v2/simdata/sim_counts.rds")
coord <- readRDS("svg/sim_mouse_brain_v2/simdata/coord.rds")
cc <- readRDS("svg/sim_mouse_brain_v2/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- as(counts, "sparseMatrix")

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_sparkg <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    i <- names(cc)[cc == k]
    counts_clu <- counts[, i]
    res <- sparkg_sw(counts = counts[rowMeans(counts_clu > 0) > 0.01, i], location = coord[i, ])
    res <- res[["res_mtest"]]
    res <- res[order(res[, "combinedPval"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "combinedPval"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_sparkg, file = "svg/sim_mouse_brain_v2/method/sparkg/res.rds")
