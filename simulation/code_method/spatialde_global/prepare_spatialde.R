suppressPackageStartupMessages(library(data.table))

for (f in c("human_lungcancer", "human_pancreas", "mouse_brain", "mouse_embryo", "mouse_intestine", "mouse_kidney")) {
  for (signal in c("v1")) {
    print(paste0(f, " ", signal))
    
    counts <- readRDS(paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
    coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
    cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))
    
    counts <- counts[, names(cc)]
    coord <- coord[names(cc), ]
    counts <- counts[rowMeans(counts > 0) > 0.01, ]
    
    counts <- data.frame(t(counts), check.names = FALSE)
    sample_info <- data.frame(coord, total_counts = rowSums(counts), cc, check.names = FALSE)
    fwrite(counts, paste0("svg/simu/", f, "/", signal, "/method/spatialde_global/counts.csv.gz"), compress = "gzip", row.names = TRUE)
    fwrite(sample_info, paste0("svg/simu/", f, "/", signal, "/method/spatialde_global/sample_info.csv.gz"), compress = "gzip", row.names = TRUE)
  }
}
