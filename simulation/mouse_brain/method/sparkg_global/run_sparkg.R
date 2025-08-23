suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

counts <- readRDS("svg/scdesign3/simdata/sim_counts.rds")
coord <- readRDS("svg/scdesign3/simdata/coord.rds")
cc <- readRDS("svg/scdesign3/simdata/cellclu.rds")

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]
counts <- counts[rowMeans(counts > 0) > 0.01, ]
counts <- as(counts, "sparseMatrix")

res <- sparkg_sw(counts = counts, location = coord)
res <- res[["res_mtest"]]
res <- res[order(res[, "combinedPval"]), ]
res[, "fdr"] <- stats::p.adjust(res[, "combinedPval"], method = "fdr")

res_sparkg <- do.call(rbind, sapply(levels(cc), function(k) {
  data.frame(cluster = k, gene = rownames(res), res)
}, simplify = FALSE))

saveRDS(res_sparkg, file = "svg/scdesign3/method/sparkg_global/res.rds")
