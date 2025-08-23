suppressPackageStartupMessages(library(Seurat))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]
expr <- log1p(t(t(counts) / colSums(counts) * 1e4))

res <- RunMoransI(data = expr, pos = coord)
res <- res[order(res[, "p.value"], -res[, "observed"]), ]
res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")

res_moransi <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = rownames(res), res)
}, simplify = FALSE))

saveRDS(res_moransi, file = "svg/scdesign3/method/moransi_global/res.rds")
