suppressPackageStartupMessages(library(Seurat))
source("ctSVG/function.R")

f <- "mouse_intestine.rds"
d <- readRDS(paste0("svg/simu/", sub(".rds", "/simdata/seurat.rds", f)))
sc <- Idents(d)
sccell <- split(names(sc), sc)

clu.list <- vector(mode = "list", length = length(sccell))
names(clu.list) <- names(sccell)

n.max <- 1e3
resolution = 0.1
jac.mat <- matrix(nrow = length(sccell), ncol = n.max)
rownames(jac.mat) <- names(sccell)

Sys.time()
for (i in seq_len(n.max)) {
  d.tmp <- FindClusters(d, verbose = FALSE, resolution = resolution, random.seed = i)
  sc.tmp <- Idents(d.tmp)
  sccell.tmp <- split(names(sc.tmp), sc.tmp)
  jac.tmp <- PairWiseJaccardSets(sc, sc.tmp)
  
  jac.mat[, i] <- apply(jac.tmp, 1, max)
  match.tmp <- apply(jac.tmp, 1, which.max)
  for (j in names(match.tmp)) {
    clu.list[[j]] <- c(clu.list[[j]], list(sccell.tmp[[match.tmp[j]]]))
  }
  if (i %% 20 == 0) {
    print(paste0("Repeat clustering: ", i))
    print(rowSums(jac.mat > 0.5, na.rm = TRUE))
  }
}

saveRDS(clu.list, file = paste0("svg/simu/", sub(".rds", "/simdata/clusterlist.rds", f)))
saveRDS(jac.mat, file = paste0("svg/simu/", sub(".rds", "/simdata/jaccardmat.rds", f)))
