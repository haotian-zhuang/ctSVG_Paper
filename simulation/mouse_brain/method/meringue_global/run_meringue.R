suppressPackageStartupMessages(library(MERINGUE))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]

expr <- normalizeCounts(counts = counts, log = FALSE, verbose = FALSE)
w <- getSpatialNeighbors(pos = coord)
res <- getSpatialPatterns(mat = expr, weight = w)
res <- res[order(res[, "p.value"], -res[, "observed"]), ]
res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")

res_meringue <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = rownames(res), res)
}, simplify = FALSE))

saveRDS(res_meringue, file = "svg/scdesign3/method/meringue_global/res.rds")
