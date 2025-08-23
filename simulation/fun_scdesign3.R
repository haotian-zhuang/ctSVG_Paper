suppressPackageStartupMessages(library(scDesign3))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

counts <- readRDS(paste0("svg/simu/", f, "/simdata/counts.rds"))
coord <- readRDS(paste0("svg/simu/", f, "/simdata/coord.rds"))
cc <- readRDS(paste0("svg/simu/", f, "/simdata/cellclu.rds"))

counts <- counts[, names(cc)]
coord <- coord[names(cc), ]

gene_detection <- sapply(levels(cc), function(k) {
  rowMeans(counts[, cc == k] > 0)
})
min_detection <- apply(gene_detection, 1, min)
topgene <- rownames(counts)[order(min_detection, decreasing = TRUE)][1:2000]

sim_counts <- do.call(cbind, sapply(levels(cc), function(k) {
  print(k)
  i <- names(cc)[cc == k]
  example_sce <- SingleCellExperiment(list(counts = counts[topgene, i]),
                                      colData = data.frame(row.names = rownames(coord[i, ]),
                                                           spatial1 = coord[i, "row"],
                                                           spatial2 = coord[i, "col"]))
  
  set.seed(1)
  example_data <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = NULL,
    spatial = c("spatial1", "spatial2"),
    other_covariates = NULL,
    corr_by = "ind"
  )
  
  example_marginal <- fit_marginal(
    data = example_data,
    predictor = "gene",
    mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 50)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 20,
    usebam = FALSE
  )
  
  set.seed(1)
  example_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 20,
    input_data = example_data$dat
  )
  
  example_para <- extract_para(
    sce = example_sce,
    marginal_list = example_marginal,
    n_cores = 20,
    family_use = "nb",
    new_covariate = example_data$newCovariate,
    data = example_data$dat
  )
  
  dev_explain <- sapply(example_marginal, function(x){
    sum = summary(x$fit)
    return(sum$dev.expl)
  })
  dev_ordered <- order(dev_explain, decreasing = TRUE)
  num_de <- round(length(dev_ordered)/10)
  ordered <- dev_explain[dev_ordered]
  if (signal == "v1") {
    de_idx <- names(ordered)[1:num_de] # top 10%
    non_de_idx <- names(ordered)[-(1:num_de)]
  } else if (signal == "v2") {
    de_idx <- names(ordered)[(1:num_de) + 2*num_de] # top 20% - 30%
    non_de_idx <- names(ordered)[-((1:num_de) + 2*num_de)]
  } else if (signal == "v3") {
    de_idx <- names(ordered)[(1:num_de) + 4*num_de] # top 40% - 50%
    non_de_idx <- names(ordered)[-((1:num_de) + 4*num_de)]
  } else if (signal == "v4") {
    de_idx <- names(ordered)[(1:num_de) + 6*num_de] # top 60% - 70%
    non_de_idx <- names(ordered)[-((1:num_de) + 6*num_de)]
  } else if (signal == "v5") {
    de_idx <- names(ordered)[(1:num_de) + 8*num_de] # top 80% - 90%
    non_de_idx <- names(ordered)[-((1:num_de) + 8*num_de)]
  }
  ###
  de_label <- c(rep("svg", length(de_idx)), rep("nonsvg", length(non_de_idx)))
  names(de_label) <- c(de_idx, non_de_idx)
  saveRDS(de_label, file = paste0("svg/simu/", f, "/", signal, "/genelist/", k, ".rds"))
  ###
  non_de_mat <- apply(example_para$mean_mat[,non_de_idx], 2, function(x){
    avg <- (max(x)+min(x))/2
    new_mean <- rep(avg, length(x))
    return(new_mean)
  })
  example_para$mean_mat[,non_de_idx] <- non_de_mat
  
  set.seed(1)
  example_newcount <- simu_new(
    sce = example_sce,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    n_cores = 20,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = example_data$newCovariate,
    important_feature = rep(TRUE, dim(example_sce)[1]),
    filtered_gene = example_data$filtered_gene
  )
  
  simu_sce <- SingleCellExperiment(list(counts = example_newcount), colData = example_data$newCovariate)
  logcounts(simu_sce) <- log1p(counts(simu_sce))
  
  de_genes = de_idx[1:5]
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  expre = lapply(de_genes, function(x){
    curr = as.matrix(counts(simu_sce)[x,])
    curr = log1p(curr)
    return(rescale(curr))
  })
  long = do.call(rbind, expre)
  long = as.data.frame(long)
  colnames(long) <- "Expression"
  long$gene = do.call(c, lapply(de_genes, function(x){rep(x,dim(expre[[1]])[1])}))
  long$x = rep(loc[,1],5)
  long$y = rep(loc[,2],5)
  as_tibble(long, rownames = "Cell") %>% ggplot(aes(x = x, y = y, color = Expression)) +geom_point(size = 0.1)+facet_grid(~gene)+ scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) + coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  ggsave(paste0("svg/simu/", f, "/", signal, "/plot_de/", k, ".jpg"), dpi = 500, width = 10, height = 6)
  
  non_de_genes = non_de_idx[1:5]
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  expre = lapply(non_de_genes, function(x){
    curr = as.matrix(counts(simu_sce)[x,])
    curr = log1p(curr)
    return(rescale(curr))
  })
  long = do.call(rbind, expre)
  long = as.data.frame(long)
  colnames(long) <- "Expression"
  long$gene = do.call(c, lapply(non_de_genes, function(x){rep(x,dim(expre[[1]])[1])}))
  long$x = rep(loc[,1],5)
  long$y = rep(loc[,2],5)
  as_tibble(long, rownames = "Cell") %>% ggplot(aes(x = x, y = y, color = Expression)) +geom_point(size = 0.1)+facet_grid(~gene)+ scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) + coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  ggsave(paste0("svg/simu/", f, "/", signal, "/plot_nonde/", k, ".jpg"), dpi = 500, width = 10, height = 6)
  
  return(counts(simu_sce))
}, simplify = FALSE))
sim_counts <- sim_counts[, names(cc)]
saveRDS(sim_counts, file = paste0("svg/simu/", f, "/", signal, "/sim_counts.rds"))
