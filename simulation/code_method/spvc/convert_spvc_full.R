for (f in c("human_crc", "human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    raw_res <- readRDS(paste0("svg/simu/", f, "/", signal, "/method/spvc/raw_res_full.rds"))
    res <- do.call(rbind, sapply(names(raw_res$results.full), function(g) {
      pval <- raw_res$results.full[[g]][["p.value"]]
      data.frame(cluster = names(pval), gene = g, pval = pval)
    }, simplify = FALSE))
    rownames(res) <- NULL
    
    res <- res[grepl("gamma", res$cluster), ]
    res <- res[res$cluster != "gamma_0", ]
    res$cluster <- sub("gamma_X", "", res$cluster)
    res$cluster <- as.character(as.numeric(res$cluster) - 1)
    
    res_spvc <- do.call(rbind, sapply(sort(unique(res$cluster)), function(k) {
      res_clu <- res[res$cluster == k, ]
      res_clu <- res_clu[order(res_clu[, "pval"]), ]
      res_clu[, "fdr"] <- stats::p.adjust(res_clu[, "pval"], method = "fdr")
      res_clu
    }, simplify = FALSE))
    rownames(res_spvc) <- NULL
    
    saveRDS(res_spvc, file = paste0("svg/simu/", f, "/", signal, "/method/spvc/res.rds"))
  }
}
