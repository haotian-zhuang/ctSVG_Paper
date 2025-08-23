suppressPackageStartupMessages(library(SPARK))
suppressPackageStartupMessages(library(Matrix))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1", "v2", "v3", "v4", "v5")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- as(counts, "sparseMatrix")
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_sparkg <- do.call(rbind, sapply(levels(cc), function(k) {
        print(k)
        i <- names(cc)[cc == k]
        counts_clu <- counts[, i]
        res <- sparkg_sw(counts = counts[rowMeans(counts_clu > 0) > 0.01, i], location = coord[i, ])
        res <- res[["res_mtest"]]
        res <- res[order(res[, "combinedPval"]), ]
        res[, "fdr"] <- stats::p.adjust(res[, "combinedPval"], method = "fdr")
        data.frame(cluster = k, gene = rownames(res), res)
      }, simplify = FALSE))
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_sparkg, file = paste0("svg/simu/", f, "/", signal, "/method/sparkg/res.rds"))
  }
}
