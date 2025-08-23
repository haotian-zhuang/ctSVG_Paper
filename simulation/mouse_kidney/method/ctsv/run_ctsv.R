suppressPackageStartupMessages(library(CTSV))
suppressPackageStartupMessages(library(SpatialExperiment))

for (f in "mouse_kidney") {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    spe <- SpatialExperiment(assays = list(counts = counts), spatialCoords = coord)
    W <- sapply(sort(unique(cc)), function(x) as.numeric(cc == x))
    dimnames(W) <- list(names(cc), sort(unique(cc)))
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_ctsv <- CTSV(spe = spe, W = W, num_core = 20)
      saveRDS(res_ctsv, file = paste0("svg/simu/", f, "/", signal, "/method/ctsv/raw_res.rds"))
      
      pval_ctsv <- res_ctsv[["pval"]]
      rownames(pval_ctsv) <- rownames(counts)
      ctsv1 <- pval_ctsv[, 1:(ncol(pval_ctsv) / 2)]
      ctsv2 <- pval_ctsv[, (ncol(pval_ctsv) / 2 + 1):ncol(pval_ctsv)]
      colnames(ctsv1) <- colnames(ctsv2) <- levels(cc)
      
      res_ctsv <- do.call(rbind, sapply(levels(cc), function(k) {
        res <- data.frame(cluster = k, gene = rownames(pval_ctsv),
                          p1 = ctsv1[, k], p2 = ctsv2[, k])
        res[, "pval"] <- pmin(res[, "p1"], res[, "p2"])
        res <- res[order(res[, "pval"]), ]
        res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
        res
      }, simplify = FALSE))
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_ctsv, file = paste0("svg/simu/", f, "/", signal, "/method/ctsv/res.rds"))
  }
}
