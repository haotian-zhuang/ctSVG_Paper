suppressPackageStartupMessages(library(spVC))
suppressPackageStartupMessages(library(Triangulation))

for (f in "mouse_kidney") {
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
    
    boundary <- readRDS(paste0("svg/real/spvc/local/boundary/", f, ".rds"))
    Tr.cell <- TriMesh(boundary, n = 2)
    V <- as.matrix(Tr.cell$V)
    Tr <- as.matrix(Tr.cell$Tr)
    
    print(start_time <- Sys.time())
    gc1 <- gc(reset = TRUE)
    runtime = system.time({
      res_spvc <- test.spVC(Y = counts, X = W, S = coord, V = V, Tr = Tr, twostep = FALSE,
                            para.cores = 20, filter.min.nonzero = 0, filter.spot.counts = 0)
    })
    gc2 <- gc()
    
    time = runtime[["elapsed"]]
    memory = sum(gc2[, 6] - gc1[, 6])
    saveRDS(res_spvc, file = paste0("svg/simu/", f, "/", signal, "/method/spvc/raw_res_full.rds"))
  }
}
