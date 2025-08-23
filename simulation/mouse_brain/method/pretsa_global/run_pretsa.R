source("ctSVG/function.R")

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]
expr <- log1p(t(t(counts) / colSums(counts) * 1e4))

res <- spatialTest(expr = expr, coord = coord, coord_permute = NULL, knot = FALSE)

res_pretsa <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = rownames(res), res)
}, simplify = FALSE))

saveRDS(res_pretsa, file = "svg/scdesign3/method/pretsa_global/res.rds")
