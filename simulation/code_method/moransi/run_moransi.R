suppressPackageStartupMessages(library(Seurat))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    expr <- log1p(t(t(counts) / colSums(counts) * 1e4))
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_moransi <- do.call(rbind, sapply(levels(cc), function(k) {
        print(k)
        i <- names(cc)[cc == k]
        expr_clu <- expr[, i]
        res <- RunMoransI(data = expr[rowMeans(expr_clu > 0) > 0.01, i], pos = coord[i, ])
        res <- res[order(res[, "p.value"], -res[, "observed"]), ]
        res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
        data.frame(cluster = k, gene = rownames(res), res)
      }, simplify = FALSE))
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_moransi, file = paste0("svg/simu/", f, "/", signal, "/method/moransi/res.rds"))
  }
}
