source("ctSVG/function.R")

counts <- readRDS("svg/sim_mouse_brain_v2/simdata/sim_counts.rds")
coord <- readRDS("svg/sim_mouse_brain_v2/simdata/coord.rds")
cc <- readRDS("svg/sim_mouse_brain_v2/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
expr <- log1p(t(t(counts) / colSums(counts) * 1e4))

sccell <- split(names(cc), cc)
sccell <- lapply(sccell, function(k) { k[!detectOutlier(coord = coord[k, ])] })

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_pretsa <- do.call(rbind, sapply(names(sccell), function(k) {
    print(k)
    i <- sccell[[k]]
    expr_clu <- expr[, i]
    res <- spatialTest(expr = expr[rowMeans(expr_clu > 0) > 0.01, i], coord = coord[i, ], coord_permute = NULL, knot = FALSE)
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_pretsa, file = "svg/sim_mouse_brain_v2/method/pretsa/res.rds")
