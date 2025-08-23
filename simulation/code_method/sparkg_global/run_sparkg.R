suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
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
    
    saveRDS(res_sparkg, file = paste0("svg/simu/", f, "/", signal, "/method/sparkg_global/res.rds"))
  }
}
