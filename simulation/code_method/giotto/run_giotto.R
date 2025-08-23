suppressPackageStartupMessages(library(Giotto))
suppressPackageStartupMessages(library(Matrix))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_giotto <- do.call(rbind, sapply(levels(cc), function(k) {
        print(k)
        i <- names(cc)[cc == k]
        counts_clu <- counts[, i]
        obj <- createGiottoObject(raw_exprs = counts[rowMeans(counts_clu > 0) > 0.01, i], spatial_locs = coord[i, ])
        obj <- normalizeGiotto(gobject = obj)
        obj <- createSpatialNetwork(gobject = obj)
        res <- binSpect(obj, bin_method = "kmeans", do_parallel = TRUE, cores = 2, verbose = FALSE)
        res <- as.data.frame(res)
        res <- res[order(res[, "p.value"], -res[, "score"]), ]
        res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
        data.frame(cluster = k, gene = res$genes, res)
      }, simplify = FALSE))
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_giotto, file = paste0("svg/simu/", f, "/", signal, "/method/giotto/res.rds"))
  }
}
