suppressPackageStartupMessages(library(data.table))

for (f in list.files("/hpc/group/jilab/zj/visiumHD/analysis/visiumHD/visa/filterexpr/", pattern = ".rds")) {
  print(f)
  
  counts <- readRDS(paste0("svg/crop/data/counts/", f))
  coord <- readRDS(paste0("svg/crop/data/coord/", f))
  cc <- readRDS(paste0("svg/crop/data/cellclu/", f))
  
  cc <- cc[cc %in% names(which(table(cc) >= 1e3))]
  cc <- droplevels(cc)
  
  counts <- counts[, names(cc)]
  coord <- coord[names(cc), ]
  
  counts <- data.frame(t(counts), check.names = FALSE)
  sample_info <- data.frame(coord, total_counts = rowSums(counts), cc, check.names = FALSE)
  fwrite(counts, paste0("svg/crop/full/spatialde/data/counts/", sub(".rds", ".csv.gz", f)), compress = "gzip", row.names = TRUE)
  fwrite(sample_info, paste0("svg/crop/full/spatialde/data/sample_info/", sub(".rds", ".csv.gz", f)), compress = "gzip", row.names = TRUE)
}
