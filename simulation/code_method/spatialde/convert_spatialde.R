for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    raw_res <- read.csv(paste0("svg/simu/", f, "/", signal, "/method/spatialde/raw_res.csv"), row.names = 1, check.names = FALSE)
    raw_res$cluster <- as.character(raw_res$cluster)
    
    res_spatialde <- do.call(rbind, sapply(levels(cc), function(k) {
      res <- raw_res[raw_res$cluster == k, ]
      res <- res[!duplicated(res$g), ]
      rownames(res) <- res$g
      res <- res[order(res[, "pval"], -res[, "LLR"]), ]
      res[, "fdr"] <- stats::p.adjust(res[, "pval"], method = "fdr")
      data.frame(gene = res$g, res)
    }, simplify = FALSE))
    
    saveRDS(res_spatialde, file = paste0("svg/simu/", f, "/", signal, "/method/spatialde/res.rds"))
  }
}
