suppressPackageStartupMessages(library(spacexr))

f <- "mouse_embryo.rds"
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
cc <- cc[i]

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

res_cside <- run.CSIDE.nonparam(myRCTD = myRCTD, cell_type_threshold = -Inf, gene_threshold = -Inf,
                                doublet_mode = FALSE, fdr = Inf)
saveRDS(res_cside, file = paste0("svg/crop/down/cside/raw_res/", f))
