suppressPackageStartupMessages(library(CTSV))
suppressPackageStartupMessages(library(SpatialExperiment))

f <- "mouse_kidney.rds"
print(f)

counts <- readRDS(paste0("svg/crop/data/counts/", f))
coord <- readRDS(paste0("svg/crop/data/coord/", f))
cc <- readRDS(paste0("svg/crop/data/cellclu/", f))

cc <- cc[cc %in% names(which(table(cc) >= 1e3))]
cc <- droplevels(cc)

sccell <- split(names(cc), cc)
set.seed(2025)
sccell <- lapply(sccell, function(k) {
  idx <- sample(length(k), length(k) * 0.8)
  k[idx]
})

i <- unlist(sccell)
cc_sub <- cc[i]
counts_sub <- counts[, names(cc_sub)]
coord_sub <- coord[names(cc_sub), ]
counts_sub <- counts_sub[rowMeans(counts_sub > 0) > 0.01, ]

print(dim(counts_sub))
print(table(cc_sub))

spe <- SpatialExperiment(assays = list(counts = counts_sub), spatialCoords = coord_sub)
W <- sapply(levels(cc), function(x) as.numeric(cc_sub == x))
dimnames(W) <- list(names(cc_sub), levels(cc))
res_ctsv <- CTSV(spe = spe, W = W, num_core = 20)

pval_ctsv <- res_ctsv[["pval"]]
rownames(pval_ctsv) <- rownames(counts_sub)
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
saveRDS(res_ctsv, file = paste0("svg/crop/down/ctsv/", f))
