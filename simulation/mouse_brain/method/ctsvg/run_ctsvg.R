source("ctSVG/function.R")

counts <- readRDS("svg/sim_mouse_brain_v1/simdata/sim_counts.rds")
coord <- readRDS("svg/sim_mouse_brain_v1/simdata/coord.rds")
cc <- readRDS("svg/sim_mouse_brain_v1/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
expr <- log1p(t(t(counts) / colSums(counts) * 1e4))

sccell <- split(names(cc), cc)
sccell <- lapply(sccell, function(k) { k[!detectOutlier(coord = coord[k, ])] })
clu.list <- readRDS("svg/sim_mouse_brain_v1/simdata/clusterlist.rds")
jac.mat <- readRDS("svg/sim_mouse_brain_v1/simdata/jaccardmat.rds")

print(start_time <- Sys.time())
gc1 <- gc(reset = TRUE)
runtime = system.time({
  res_ctsvg <- do.call(rbind, sapply(names(sccell), function(k) {
    print(k)
    top.index <- order(jac.mat[k, ], decreasing = TRUE)[1:100]
    cell.index <- clu.list[[k]][top.index]
    
    i <- sccell[[k]]
    coord_permute <- lapply(cell.index, function(j) {
      coord.tmp <- coord[j, ]
      rownames(coord.tmp) <- sample(rownames(coord.tmp), replace = FALSE)
      coord.tmp
    })
    
    expr_clu <- expr[, i]
    res <- spatialTest(expr = expr[rowMeans(expr_clu > 0) > 0.01, ], coord = coord[i, ], coord_permute = coord_permute, knot = FALSE)
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
})
gc2 <- gc()

time = runtime[["elapsed"]]
memory = sum(gc2[, 6] - gc1[, 6])
saveRDS(res_ctsvg, file = "svg/sim_mouse_brain_v1/method/ctsvg/res.rds")
