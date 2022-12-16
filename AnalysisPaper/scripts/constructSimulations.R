# simulating paired scRNA-seq and ST data using splatter.
library(splatter)

# generate gene expression data and spatial coordinates
sim_data_construct <- function(gene_num = 5000,
                               batch_cell_num_1 = 1000,
                               batch_cell_num_2 = 1000,
                               group_prob,
                               de_prob,
                               max_range = 20,
                               seed = 123
){
  
  params <- newSplatParams(seed = seed)
  params <- setParams(
    params, 
    update = list(nGenes = gene_num,
                  batchCells = c(batch_cell_num_1, batch_cell_num_2),
                  group.prob = group_prob,
                  de.prob = de_prob)
  )
  
  message('Simulating data...')
  sim <- splatSimulate(params, method = "group", verbose = FALSE)
  counts <- data.frame(sim@assays@data@listData[["counts"]])
  meta <- data.frame(sim@colData)
  meta$Group <- as.character(meta$Group)
  
  sc_meta <- meta[meta$Batch == 'Batch1', ]
  sc_data <- counts[, rownames(sc_meta)]
  message('ScRNA-seq data simulation done...')
  
  st_meta <- meta[meta$Batch == 'Batch2', ]
  st_data <- counts[, rownames(st_meta)]
  st_meta$pseudo_x <- 0
  st_meta$pseudo_y <- 0
  
  cell_type_num <- length(group_prob)
  cell_type_name <- names(table(st_meta$Group))
  cell_num_per_ct <- as.vector(table(st_meta$Group))
  
  x_center <- sample(-max_range:max_range, cell_type_num)
  y_center <- sample(-max_range:max_range, cell_type_num)
  
  for (i in 1:cell_type_num) {
    st_meta[st_meta$Group == cell_type_name[i],]$pseudo_x <- rnorm(cell_num_per_ct[i], x_center[i])
    st_meta[st_meta$Group == cell_type_name[i],]$pseudo_y <- rnorm(cell_num_per_ct[i], y_center[i])
  }
  
  message('SrT data simulation done...')
  message('Simulating done.')
  
  return(list(sc_data = sc_data,
              sc_meta = sc_meta,
              st_data = st_data,
              st_meta = st_meta))
}

# generate proportion of each celltype randomly
generate_ct_prob <- function(ct_num, 
                             seed = 123){
  set.seed(seed)
  ct_prob <- runif(ct_num, min = 0, max = 1)
  ct_prob <- ct_prob / sum(ct_prob)
  return(ct_prob)
}

# generate group differential expression of each celltype randomly, 0.2 for main celltypes and 0.01 for spatial-subtypes
generate_de_prob <- function(ct_num,
                             sub_ct_num){
  de_prob <- c(rep(0.2, (ct_num - sub_ct_num)), rep(0.01, sub_ct_num))
  return(de_prob)
}