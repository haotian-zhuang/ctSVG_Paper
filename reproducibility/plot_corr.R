suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
pal <- c(brewer.pal(8, "Set1")[c(1:5)], "gold", brewer.pal(8, "Set1")[c(7:8)],
         brewer.pal(8, "Dark2")[c(1:3)],
         brewer.pal(8, "Set1")[c(2:5)], "gold", brewer.pal(8, "Set1")[c(7:8)])
names(pal) <- c("ctSVG", paste0(c("Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE"), "-CT"),
                "CTSV", "C-SIDE", "spVC",
                "Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE")

ds <- list.files("ctSVG/ctSVG", pattern = ".rds")
names(ds) <- c("human colon cancer", "human lung cancer", "human pancreas",
               "mouse brain", "mouse embryo", "mouse intestine", "mouse kidney")

method <- c("ctsvg", "giotto", "meringue", "moransi", "nnsvg", "pretsa", "sparkg", "spatialde",
            "ctsv", "cside")
names(method) <- c("ctSVG", paste0(c("Giotto", "MERINGUE", "Moran's I", "nnSVG", "PreTSA", "SPARK", "SpatialDE"), "-CT"),
                   "CTSV", "C-SIDE")

df <- data.frame()
for (m in method) {
  for (d in intersect(list.files(paste0("svg/crop/full/", m), pattern = ".rds"),
                      list.files(paste0("svg/crop/down/", m), pattern = ".rds"))) {
    sub1 <- readRDS(paste0("svg/crop/full/", m, "/", d))
    sub2 <- readRDS(paste0("svg/crop/down/", m, "/", d))
    
    # C-SIDE filtered out genes
    clu <- intersect(names(which(table(sub1$cluster) >= 1e3)),
                     names(which(table(sub2$cluster) >= 1e3)))
    sub1 <- sub1[sub1$cluster %in% clu, ]
    sub2 <- sub2[sub2$cluster %in% clu, ]
    
    corr <- sapply(unique(sub1$cluster), function(i) {
      s1 <- sub1[sub1$cluster == i, ]
      s2 <- sub2[sub2$cluster == i, ]
      r1 <- setNames(seq_len(nrow(s1)), nm = s1$gene)
      r2 <- setNames(seq_len(nrow(s2)), nm = s2$gene)
      gene <- intersect(names(r1), names(r2))
      r1 <- r1[gene]
      r2 <- r2[gene]
      cor(r1, r2, method = "spearman")
    })
    
    df <- rbind(df, data.frame(corr = corr,
                               cluster = unique(sub1$cluster),
                               method = names(method)[method == m],
                               data = names(ds)[ds == d]))
  }
}

corr_order <- df %>%
  group_by(method) %>%
  summarise(corr = median(corr)) %>%
  arrange(desc(corr)) %>%
  pull(method)

df$method <- factor(df$method, levels = rev(corr_order))
ggplot(df, aes(x = corr, y = method, fill = method)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Spearman correlation", y = NULL) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
ggsave("svg/crop/correlation.pdf", width = 10, height = 2.5)
