suppressPackageStartupMessages(library(spacexr))

for (f in "human_pancreas") {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    W <- sapply(sort(unique(cc)), function(x) as.numeric(cc == x))
    dimnames(W) <- list(names(cc), sort(unique(cc)))
    
    spatialRNA <- SpatialRNA(coords = as.data.frame(coord), counts = counts)
    reference <- Reference(counts = counts, cell_types = cc, min_UMI = -Inf)
    myRCTD <- create.RCTD(spatialRNA = spatialRNA, reference = reference, max_cores = 1,
                          gene_cutoff = -Inf, fc_cutoff = -Inf, gene_cutoff_reg = -Inf, fc_cutoff_reg = -Inf, UMI_min = -Inf,
                          UMI_max = Inf, counts_MIN = -Inf, UMI_min_sigma = -Inf, CELL_MIN_INSTANCE = -Inf)
    myRCTD@config[["MIN_OBS"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_BULK"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_REG"]] <- -Inf
    myRCTD@config[["CONFIDENCE_THRESHOLD"]] <- -Inf
    # myRCTD <- run.RCTD(RCTD = myRCTD)
    myRCTD@config$RCTDmode <- "full"
    myRCTD <- import_weights(myRCTD = myRCTD, weights = W)
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_cside <- run.CSIDE.nonparam(myRCTD = myRCTD, cell_type_threshold = -Inf, gene_threshold = -Inf,
                                      doublet_mode = FALSE, fdr = Inf)
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_cside, file = paste0("svg/simu/", f, "/", signal, "/method/cside/raw_res.rds"))
  }
}
