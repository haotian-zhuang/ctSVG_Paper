suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(PRROC))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
pal <- c(brewer.pal(8, "Set1")[c(1:5)], "gold", brewer.pal(8, "Set1")[c(7:8)], brewer.pal(8, "Dark2")[c(1:3)],
         brewer.pal(8, "Set1")[c(2:5)], "gold", brewer.pal(8, "Set1")[c(7:8)])
names(pal) <- c("ctSVG", paste0(c("Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE"), "-CT"),
                "CTSV", "C-SIDE", "spVC",
                "Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE")

method <- c("ctsvg", "giotto", "meringue", "moransi", "nnsvg", "pretsa", "sparkg", "spatialde",
            "ctsv", "cside", "spvc",
            "giotto_global", "meringue_global", "moransi_global", "nnsvg_global", "pretsa_global", "sparkg_global", "spatialde_global")
names(method) <- c("ctSVG", paste0(c("Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE"), "-CT"),
                   "CTSV", "C-SIDE", "spVC",
                   "Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE")

pval_name <- c("pval.parametric", "p.value", "p.value", "p.value", "pval", "pval", "combinedPval", "pval",
               "pval", "pval", "pval",
               "p.value", "p.value", "p.value", "pval", "pval", "combinedPval", "pval")
names(pval_name) <- method

ds <- sub(".rds", "", list.files("ctSVG/ctSVG", pattern = ".rds"))
names(ds) <- c("human colon cancer", "human lung cancer", "human pancreas",
               "mouse brain", "mouse embryo", "mouse intestine", "mouse kidney")

signal <- c("v1")
names(signal) <- c("v1")

df <- data.frame()
for (d in ds) {
  for (s in signal) {
    for (m in method) {
      res <- readRDS(paste0("svg/simu/", d, "/", s, "/method/", m, "/res.rds"))
      colnames(res)[colnames(res) == pval_name[m]] <- "pval"
      if (m == "ctsvg") {
        colnames(res)[3] <- "fdr"
      }
      res <- res[!is.na(res$pval), ]
      res$svg <- factor(ifelse(res$fdr < 0.05, 1, 0))
      
      for (k in unique(res$cluster)) {
        genelist <- readRDS(paste0("svg/simu/", d, "/", s, "/genelist/", k, ".rds"))
        genelist <- factor(ifelse(genelist == "svg", 1, 0))
        res.tmp <- res[res$cluster == k, ]
        rownames(res.tmp) <- res.tmp$gene
        
        if (nrow(res.tmp) < 1e3) {
          next
        }
        
        genelist <- genelist[rownames(res.tmp)]

        auprc <- pr.curve(scores.class0 = -res.tmp[genelist == 1, "pval"],
                          scores.class1 = -res.tmp[genelist == 0, "pval"])[["auc.integral"]]
        
        df <- rbind(df, data.frame(auprc = auprc,
                                   method = names(method)[method == m],
                                   cluster = paste0("cell type ", (as.numeric(k) + 1)),
                                   data = names(ds)[ds == d],
                                   signal = names(signal)[signal == s]))
      }
    }
  }
}
df$method <- factor(df$method, levels = names(method))

auprc_order <- df %>%
  group_by(method) %>%
  summarise(auprc = median(auprc)) %>%
  arrange(desc(auprc)) %>%
  pull(method)
df$method <- factor(df$method, levels = rev(auprc_order))

ggplot(df, aes(x = auprc, y = method, fill = method)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "AUPRC", y = NULL) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
ggsave("svg/simu/auprc.pdf", width = 10, height = 4)
