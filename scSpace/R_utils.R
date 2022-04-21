#' This is utils R script
#' 


# data processing
data_process <- function(sc_data, 
                         sc_meta,
                         normalize = TRUE,
                         n_features = 2000,
                         npcs = 50,
                         dims = 20){
  seu_obj <- Seurat::CreateSeuratObject(sc_data, meta.data = sc_meta, verbose = FALSE)
  if(normalize){
    seu_obj <- Seurat::NormalizeData(seu_obj, verbose = FALSE)
  }
  seu_obj <- Seurat::FindVariableFeatures(seu_obj, nfeatures = n_features, verbose = FALSE)
  seu_obj <- Seurat::ScaleData(seu_obj, verbose = FALSE)
  seu_obj <- Seurat::RunPCA(seu_obj, npcs = npcs, verbose = FALSE)
  seu_obj <- Seurat::RunTSNE(seu_obj, dims = 1:dims)
  
  return(seu_obj)
}


# search leiden resolution of target cluster number
search_res <- function(spa_graph, 
                       target_num, 
                       start = 0.4, 
                       step = 0.1, 
                       n_seed = 123){
  set.seed(n_seed)
  res <- start
  message(paste0("Start at res = ", res, ", step = ", step))
  
  original_num <- length(table(leidenAlg::leiden.community(spa_graph, resolution = res)$membership))
  message(paste0("Res = ", res, ", number of clusters = ", original_num))
  
  while (original_num != target_num) {
    
    if(original_num < target_num){
      sign <- 1
    } else{
      sign <- (-1)
    }
    set.seed(n_seed)
    new_num <- length(table(leidenAlg::leiden.community(spa_graph, resolution = (res + step*sign))$membership))
    message(paste0("Res = ", (res + step*sign), ", number of clusters = ", new_num))
    
    if(new_num == target_num){
      res <- res + step*sign
      # message(paste0("Recommended res = ", res))
      break
    } else{
      if(new_num < target_num){
        new_sign <- 1
      } else{
        new_sign <- (-1)
      }
      if(new_sign == sign){
        res <- res + step*sign
        message(paste0("Res changed to ", res))
        original_num <- new_num
      } else{
        step <- step / 2
        message(paste0("Step changed to ", step))
      }
    }
    
  }
  
  message(paste0("Recommended res = ", res))
  return(res)
}


# spatial-informed clustering
spa_cluster <- function(seu_obj,
                        Ks = 10,
                        Kg = 20,
                        alpha = 0,
                        beta = 0,
                        target_num = 0,
                        res = 0.4,
                        step = 0.1,
                        n_seed = 123){
  
  # create gene expression graph
  knn_idx <- data.frame(BiocNeighbors::findKNN(seu_obj@meta.data[,c('Pseudo_space1','Pseudo_space2')], k = Ks)$index)
  N <- nrow(knn_idx)
  W <- matrix(0, N, N)
  
  diag(W) <- 1
  for (i in 1:nrow(W)) {
    W[i, as.numeric(knn_idx[i,])] <- 1
  }
  
  colnames(W) <- rownames(W) <- colnames(seu_obj)
  
  pca_matrix <- seu_obj@reductions$pca@cell.embeddings
  
  knn_idx <- BiocNeighbors::findKNN(pca_matrix, k = Kg)$index
  
  # create space graph
  sp_graph <- igraph::graph_from_adjacency_matrix(W)
  sp_graph <- igraph::simplify(sp_graph)
  sp_weight <- matrix(ncol = ncol(knn_idx), nrow = nrow(knn_idx))
  for (i in 1:nrow(knn_idx)) {
    to_idx <- as.numeric(knn_idx[i,])
    sp_weight[i,] <- igraph::distances(sp_graph, rownames(W)[i], colnames(W)[to_idx])
  }
  
  weight <- 1 / (alpha + as.vector(sp_weight)) + beta
  
  # create spatial-weighted gene expression graph
  df <- data.frame(from = rep(1:nrow(knn_idx), Kg),
                   to = as.vector(knn_idx),
                   weight = weight)
  
  g <- igraph::graph_from_data_frame(df, directed = FALSE)
  g <- igraph::simplify(g)
  
  if(target_num == 0){
    message(paste0('Unsupervised clustering with res = ', res, '...'))
    set.seed(n_seed)
    spa_clu <- leidenAlg::leiden.community(g, resolution = res)
    seu_obj$scSpace <- spa_clu$membership
  } else if(target_num <= 1){
    stop('Target cluster number must be greater than 1!')
    } else{
      message(paste0('Unsupervised clustering with target cluster num = ', target_num, ', sarching appropriate res firstly...'))
      res_new <- search_res(spa_graph = g,
                            target_num = target_num,
                            start = res,
                            step = step,
                            n_seed = n_seed)
      set.seed(n_seed)
      spa_clu <- leidenAlg::leiden.community(g, resolution = res_new)
      seu_obj$scSpace <- spa_clu$membership
    }
  
  return(seu_obj)
}


# saving
save_data <- function(seu_obj,
                      project_name){
  sc_data <- data.frame(seu_obj@assays$RNA@counts)
  sc_meta <- seu_obj@meta.data
  sc_meta$tsne1 <- seu_obj@reductions$tsne@cell.embeddings[, 1]
  sc_meta$tsne2 <- seu_obj@reductions$tsne@cell.embeddings[, 2]
  
  output_file <- paste0('./data/', project_name, '/output_data')
  if (!dir.exists(output_file)){
    dir.create(output_file)}
  
  write.csv(sc_data, file = paste0(output_file, '/output_sc_data.csv'))
  write.csv(sc_meta, file = paste0(output_file, '/output_sc_meta.csv'))
  
  message('Saving data OK!')

}




