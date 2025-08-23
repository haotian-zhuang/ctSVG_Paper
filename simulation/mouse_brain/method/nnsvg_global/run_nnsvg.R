suppressPackageStartupMessages(library(nnSVG))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]
expr <- log2(t(t(counts) / colSums(counts) * 1e4) + 1)

res <- nnSVG(input = expr, spatial_coords = coord,
             n_threads = 20)
res <- as.data.frame(res)
res <- res[order(res[, "rank"]), ]
res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")

res_nnsvg <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = rownames(res), res)
}, simplify = FALSE))

saveRDS(res_nnsvg, file = "svg/scdesign3/method/nnsvg_global/res.rds")
