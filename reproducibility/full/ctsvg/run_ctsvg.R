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
  
  clu.list <- readRDS(paste0("ctSVG/clusterlist/", f))
  jac.mat <- readRDS(paste0("ctSVG/jaccardmat/", f))
  
  clures <- do.call(rbind, sapply(names(sccell), function(k) {
    print(k)
    top.index <- order(jac.mat[k, ], decreasing = TRUE)[1:100]
    cell.index <- clu.list[[k]][top.index]
    
    i <- sccell[[k]]
    coord_permute <- lapply(cell.index, function(j) {
      j <- j[j %in% names(cc)] #
      coord.tmp <- coord[j, ]
      rownames(coord.tmp) <- sample(rownames(coord.tmp), replace = FALSE)
      coord.tmp
    })
    
    expr_clu <- expr[, i]
    res <- spatialTest(expr = expr[rowMeans(expr_clu > 0) > 0.01, ], coord = coord[i, ], coord_permute = coord_permute, knot = FALSE)
    data.frame(cluster = k, gene = rownames(res), res)
  }, simplify = FALSE))
  saveRDS(clures, file = paste0("svg/crop/full/ctsvg/", f))
}
