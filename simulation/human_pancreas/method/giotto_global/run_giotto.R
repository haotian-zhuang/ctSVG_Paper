suppressPackageStartupMessages(library(Giotto))
suppressPackageStartupMessages(library(Matrix))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]

obj <- createGiottoObject(raw_exprs = counts, spatial_locs = coord)
obj <- normalizeGiotto(gobject = obj)
obj <- createSpatialNetwork(gobject = obj)
res <- binSpect(obj, bin_method = "kmeans", do_parallel = FALSE, verbose = FALSE)
res <- as.data.frame(res)
res <- res[order(res[, "p.value"], -res[, "score"]), ]
res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")

res_giotto <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = res$genes, res)
}, simplify = FALSE))

saveRDS(res_giotto, file = "svg/scdesign3/method/giotto_global/res.rds")
