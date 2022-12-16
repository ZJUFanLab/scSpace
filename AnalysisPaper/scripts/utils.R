# custom functions for downstream analysis


run_RCTD <- function(sc_meta,
                     sc_data,
                     st_meta,
                     st_data,
                     celltype_index){
  ### Create the Reference object
  cell_types <- as.factor(sc_meta[[celltype_index]]); names(cell_types) <- colnames(sc_data)
  reference <- Reference(sc_data, cell_types, require_int = FALSE)
  
  coords <- data.frame(x = st_meta$xcoord, y = st_meta$ycoord)
  rownames(coords) <- colnames(st_data)
  
  ### Create SpatialRNA object
  puck <- SpatialRNA(coords, st_data, require_int = FALSE)
  
  myRCTD <- create.RCTD(puck, reference, CELL_MIN_INSTANCE = 5)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights <- normalize_weights(results$weights) 
  norm_weights <- data.frame(as.matrix(norm_weights))
  
  return(norm_weights)
}


run_Seurat <- function(sc_meta,
                       sc_data,
                       res){
  sc_seu <- CreateSeuratObject(sc_data, meta.data = sc_meta, verbose = FALSE)
  sc_seu <- NormalizeData(sc_seu, verbose = FALSE)
  sc_seu <- FindVariableFeatures(sc_seu, nfeatures = 2000, verbose = FALSE)
  sc_seu <- ScaleData(sc_seu, features = rownames(sc_seu), verbose = FALSE)
  sc_seu <- RunPCA(sc_seu, verbose = FALSE)
  sc_seu <- RunTSNE(sc_seu, dims = 1:20, verbose = FALSE)
  sc_seu <- FindNeighbors(sc_seu)
  sc_seu <- FindClusters(sc_seu, resolution = res)
  return(sc_seu)
}


cal_dist <- function(meta,
                     coord_index,
                     group_by,
                     selected_type,
                     ignore_select_type = FALSE){
  dist <- as.matrix(dist(meta[, coord_index]))
  
  ct_type <- sort(unique(meta[[group_by]]))
  
  if(ignore_select_type){
    ct_type <- setdiff(ct_type, selected_type)
  }
  
  select_ct_index <- which(meta[[group_by]] == selected_type)
  
  pseudo_dist_table <- data.frame()
  
  message('Beginning normalized distance calculating...')
  pb <- progress::progress_bar$new(format = ' Calculating [:bar] :percent eta: :eta', 
                                   total = length(ct_type), clear = FALSE, width = 60, 
                                   complete = "+", incomplete = "-")
  
  # dist
  for (i in ct_type) {
    ct_index <- which(meta[[group_by]] == i)
    dist_list <- c()
    if(i == selected_type){
      for (j in ct_index) {
        col_index <- j
        dist_list <- c(dist_list, dist[ct_index, col_index])
        ct_index <- setdiff(ct_index, col_index)
      }
    }
    else{
      for (j in select_ct_index) {
        col_index <- j
        dist_list <- c(dist_list, dist[ct_index, col_index])
      }
    }
    
    dist_table <- data.frame(dist = dist_list,
                             group = rep(paste0(selected_type, '_', i), length(dist_list)))
    
    pseudo_dist_table <- rbind(pseudo_dist_table, dist_table)
    
    pb$tick()
  }
  print('Normalized distance calculating done.')
  
  pseudo_dist_table$dist <- pseudo_dist_table$dist / max(pseudo_dist_table$dist)
  return(pseudo_dist_table)
  
}


plotEnrichment_new <- function(pathway, stats,
                               gseaParam=1,
                               ticksSize=0.2) {
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    # geom_point(color="blue", size=0.1) +
    # geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    # geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="#059669", size = 1.5) + theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    
    labs(x="rank", y="enrichment score")
  g
}


run_DRSC <- function(data,
                     meta,
                     cluster_num){
  
  
  seu <- CreateSeuratObject(data, meta.data = meta)
  seu <- NormalizeData(seu, verbose = F)
  # choose 480 spatially variable features
  seu <- FindSVGs(seu, nfeatures = 480)
  ### Given K
  seu <- DR.SC(seu, K=cluster_num, platform = 'seqfish', verbose=F)
  
  res <- seu@meta.data
  
  return(res)
}


run_BayesSpace <- function(data,
                           meta,
                           cluster_num){
  
  
  sce <- SingleCellExperiment(assays=list(counts=as(as.matrix(data), "dgCMatrix")), colData=meta)
  
  set.seed(102)
  sce <- spatialPreprocess(sce, platform="Visium", 
                           n.PCs=15, n.HVGs=2000, log.normalize=TRUE)
  set.seed(149)
  sce <- spatialCluster(sce, q=cluster_num, platform="Visium", d=15,
                        init.method="mclust", model="t", gamma=2,
                        nrep=10000, burn.in=100,
                        save.chain=TRUE)
  
  res <- data.frame(sce@colData)
  
  return(res)
}


# evaluate clustering results in benchmarking step
eva_function <- function(clu_obj,
                         sim_id,
                         method_choose,
                         all_cluster = TRUE,
                         select_spatial_cluster = NULL){
  if(all_cluster){
    print('calculating all clusters...')
    meta <- clu_obj
  }
  else{
    print('calculating spatial clustering...')
    meta <- clu_obj
    meta <- meta[meta$Group %in% select_spatial_cluster, ]
  }
  
  meta[[method_choose]] <- as.factor(meta[[method_choose]])
  
  ari <- mclust::adjustedRandIndex(meta$Group, meta[[method_choose]])
  nmi <- NMI::NMI(data.frame(cell = 1:nrow(meta), clu = meta$Group), 
                  data.frame(cell = 1:nrow(meta), clu = meta[[method_choose]]))$value
  
  
  
  res_obj <- data.frame(
    sim = sim_id,
    method = method_choose,
    ari = ari,
    nmi = nmi)
  
  return(res_obj)
}








 # old spatial-informed clustering for scSpace
spa_cluster <- function(sc_seu, 
                        coord_index,
                        seeds = 100,
                        res,
                        Ks, 
                        Kg){
  knn_idx <- data.frame(BiocNeighbors::findKNN(sc_seu@meta.data[,coord_index], k = Ks)$index)
  N <- nrow(knn_idx)
  W <- matrix(0, N, N)
  
  diag(W) <- 1
  for (i in 1:nrow(W)) {
    W[i, as.numeric(knn_idx[i,])] <- 1
  }
  
  colnames(W) <- rownames(W) <- colnames(sc_seu)
  
  sc_pca_ori <- sc_seu@reductions$pca@cell.embeddings
  
  
  # create gene expression k-nearest neighbors
  knn_idx <- BiocNeighbors::findKNN(sc_pca_ori, k = Kg)$index
  
  # create weights from binary adjacency matrix W
  sp_graph <- igraph::graph_from_adjacency_matrix(W)
  sp_graph <- igraph::simplify(sp_graph)
  sp_weight <- matrix(ncol = ncol(knn_idx), nrow = nrow(knn_idx))
  for (i in 1:nrow(knn_idx)) {
    to_idx <- as.numeric(knn_idx[i,])
    sp_weight[i,] <- igraph::distances(sp_graph, rownames(W)[i], colnames(W)[to_idx])
  }
  
  weight <- 1/(0 + as.vector(sp_weight)) + 0
  
  # graph-based clustering with spatial weights
  df <- data.frame(from = rep(1:nrow(knn_idx), Kg),
                   to = as.vector(knn_idx),
                   weight = weight)
  
  g <- igraph::graph_from_data_frame(df, directed = FALSE)
  g <- igraph::simplify(g)
  
  set.seed(seeds)
  cluster <- leidenAlg::leiden.community(g, resolution = res, n.iterations = 2)
  sc_seu$scSpace <- cluster$membership
  
  return(sc_seu)
}



