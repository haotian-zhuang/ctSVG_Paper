for (f in list.files("svg/crop/full/cside/raw_res/", pattern = ".rds")) {
  print(f)
  
  raw_res <- readRDS(paste0("svg/crop/full/cside/raw_res/", f))
  res_cside <- do.call(rbind, sapply(names(raw_res@de_results[["all_gene_list"]]), function(k) {
    res_clu <- raw_res@de_results[["all_gene_list"]][[k]]
    res_clu <- data.frame(cluster = k, gene = rownames(res_clu), pval = res_clu[, "p_val"])
    res_clu <- res_clu[order(res_clu[, "pval"]), ]
    res_clu[, "fdr"] <- stats::p.adjust(res_clu[, "pval"], method = "fdr")
    res_clu
  }, simplify = FALSE))
  rownames(res_cside) <- NULL
  saveRDS(res_cside, file = paste0("svg/crop/full/cside/", f))
}
