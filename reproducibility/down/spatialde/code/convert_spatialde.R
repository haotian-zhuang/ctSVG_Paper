for (f in list.files("svg/crop/down/spatialde/raw_res/", pattern = ".csv.gz")) {
  print(f)
  
  cc <- readRDS(paste0("svg/crop/data/cellclu/", sub(".csv.gz", ".rds", f)))
  raw_res <- read.csv(paste0("svg/crop/down/spatialde/raw_res/", f), row.names = 1, check.names = FALSE)
  raw_res$cluster <- as.character(raw_res$cluster)
  
  res_spatialde <- do.call(rbind, sapply(levels(cc), function(k) {
    res <- raw_res[raw_res$cluster == k, ]
    res <- res[!duplicated(res$g), ]
    rownames(res) <- res$g
    res <- res[order(res[, "pval"], -res[, "LLR"]), ]
    res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
    data.frame(gene = res$g, res)
  }, simplify = FALSE))
  saveRDS(res_spatialde, file = paste0("svg/crop/down/spatialde/", sub(".csv.gz", ".rds", f)))
}
