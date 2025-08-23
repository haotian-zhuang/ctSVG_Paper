source("ctSVG/function.R")

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    expr <- log1p(t(t(counts) / colSums(counts) * 1e4))
    
    res <- spatialTest(expr = expr, coord = coord, coord_permute = NULL, knot = FALSE)
    
    res_pretsa <- do.call(rbind, sapply(levels(cc), function(k) {
      data.frame(cluster = k, gene = rownames(res), res)
    }, simplify = FALSE))
    
    saveRDS(res_pretsa, file = paste0("svg/simu/", f, "/", signal, "/method/pretsa_global/res.rds"))
  }
}
