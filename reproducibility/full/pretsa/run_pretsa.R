source("ctSVG/function.R")

for (f in list.files("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/filterexpr/", pattern = ".rds")) {
  print(f)
  
  counts <- readRDS(paste0("svg/crop/data/counts/", f))
  coord <- readRDS(paste0("svg/crop/data/coord/", f))
  cc <- readRDS(paste0("svg/crop/data/cellclu/", f))
  
  cc <- cc[cc %in% names(which(table(cc) >= 1e3))]
  cc <- droplevels(cc)
  
  counts <- counts[, names(cc)]
  coord <- coord[names(cc), ]
  expr <- log1p(t(t(counts) / colSums(counts) * 1e4))
  
  sccell <- split(names(cc), cc)
  
  res_pretsa <- do.call(rbind, sapply(names(sccell), function(k) {
    print(k)
    i <- sccell[[k]]
    expr_clu <- expr[, i]
    res <- spatialTest(expr = expr[rowMeans(expr_clu > 0) > 0.01, i], coord = coord[i, ], coord_permute = NULL, knot = FALSE)
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
  saveRDS(res_pretsa, file = paste0("svg/crop/full/pretsa/", f))
}
