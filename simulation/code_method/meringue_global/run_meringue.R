suppressPackageStartupMessages(library(MERINGUE))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    
    expr <- normalizeCounts(counts = counts, log = FALSE, verbose = FALSE)
    w <- getSpatialNeighbors(pos = coord)
    
    gene_split <- split(rownames(expr), sample(1:10, nrow(expr), replace = TRUE))
    names(gene_split) <- NULL
    res <- do.call(rbind, parallel::mclapply(gene_split, function(gs) {
      getSpatialPatterns(mat = expr[gs, ], weight = w)
    }, mc.cores = 10))
    
    # res <- getSpatialPatterns(mat = expr, weight = w)
    res <- res[order(res[, "p.value"], -res[, "observed"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "p.value"], method = "fdr")
    
    res_meringue <- do.call(rbind, sapply(levels(cc), function(k) {
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
    
    saveRDS(res_meringue, file = paste0("svg/simu/", f, "/", signal, "/method/meringue_global/res.rds"))
  }
}
