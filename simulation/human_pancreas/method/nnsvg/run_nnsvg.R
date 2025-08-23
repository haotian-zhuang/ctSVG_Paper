suppressPackageStartupMessages(library(nnSVG))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
expr <- log2(t(t(counts) / colSums(counts) * 1e4) + 1)

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_nnsvg <- do.call(rbind, sapply(levels(cc), function(k) {
    print(k)
    i <- names(cc)[cc == k]
    expr_clu <- expr[, i]
    if (length(i) > 65) {
      res <- nnSVG(input = expr[rowMeans(expr_clu > 0) > 0.01, i], spatial_coords = coord[i, ],
                   n_threads = 20)
    } else {
      res <- nnSVG(input = expr[rowMeans(expr_clu > 0) > 0.01, i], spatial_coords = coord[i, ],
                   order = "Sum_coords", n_threads = 20)
    }
    res <- as.data.frame(res)
    res <- res[order(res[, "rank"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_nnsvg, file = "svg/scdesign3/method/nnsvg/res.rds")
