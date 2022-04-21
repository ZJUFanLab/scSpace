#' This R script is for scSpace cluster
#' @return This R script returns out data saved in folders "output_data"
#' @import  processed data
#' @export: output data
#' @example demo

# parameters
option_list <- list(optparse::make_option('--project', type = 'character', 
                                          default = 'demo', help = 'project name'),
                    optparse::make_option('--sc_data', type = 'character', 
                                          default = 'demo_sc_data', help = 'sc_data file'),
                    optparse::make_option('--sc_meta', type = 'character', 
                                          default = 'demo_sc_meta', help = 'sc_meta file'),
                    optparse::make_option('--sub_cluster', type = 'logical', 
                                          default = FALSE, help = 'True: sub-clustering; False: clustering'),
                    optparse::make_option('--idents', type = 'character', 
                                          default = 'Group', help = 'cell type annotation'),
                    optparse::make_option('--select_celltype', type = 'character', 
                                          default = 'Group1', help = 'cell type for sub-clustering'),
                    optparse::make_option('--normalize', type = 'logical', 
                                          default = TRUE, help = 'True: normalize the data'),
                    optparse::make_option('--n_features', type = 'integer', 
                                          default = 2000, help = 'Number of features to select as top variable features'),
                    optparse::make_option('--n_pcs', type = 'integer', 
                                          default = 50, help = 'Total Number of PCs to compute'),
                    optparse::make_option('--n_dims', type = 'integer', 
                                          default = 20, help = 'Total Number of dims to use'),
                    optparse::make_option('--Ks', type = 'integer', 
                                          default = 10, help = 'Number of nearest cells in pseudo space'),
                    optparse::make_option('--Kg', type = 'integer', 
                                          default = 20, help = 'Number of nearest cells in gene expression'),
                    optparse::make_option('--sp_alpha', type = 'integer', 
                                          default = 0, help = 'parameter for spatial weight'),
                    optparse::make_option('--sp_beta', type = 'integer', 
                                          default = 0, help = 'parameter for spatial weight'),
                    optparse::make_option('--n_seed', type = 'integer', 
                                          default = 12312, help = 'random seed'),
                    optparse::make_option('--target_num', type = 'integer', 
                                          default = 0, help = 'target cluster number for spatial-informed cluster'),
                    optparse::make_option('--res', type = 'double', 
                                          default = 0.5, help = 'start for searching res'),
                    optparse::make_option('--step', type = 'double', 
                                          default = 0.1, help = 'iteration step for searching res'))


args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

source('R_utils.R')

message("=== scSpace clustering ===")
message('Loading data...')
sc_data <- read.csv(file = paste0('./data/', args$project, '/', args$sc_data, '.csv'), row.names = 1)
sc_meta <- read.csv(file = paste0('./data/', args$project, '/', args$sc_meta, '.csv'), row.names = 1)
pseudo_space <- read.csv(file = paste0('./data/', args$project, '/processed_data/pseudo_space.csv'), row.names = 1)
message('Data have been loaded...')

rownames(sc_meta) <- colnames(sc_data)
sc_meta$Pseudo_space1 <- pseudo_space[, 1]
sc_meta$Pseudo_space2 <- pseudo_space[, 2]

if(args$sub_cluster){
  message(paste0('Spatial-informed sub-clustering of ', args$select_celltype, '...'))
  sc_meta <- sc_meta[sc_meta[args$idents] == args$select_celltype, ]
  sc_data <- sc_data[, rownames(sc_meta)]
}

message('Data processing...')
sc_seu <- data_process(sc_data = sc_data,
                       sc_meta = sc_meta,
                       normalize = args$normalize,
                       n_features = args$n_features,
                       npcs = args$n_pcs,
                       dims = args$n_dims)

message('Begin spatial-informed clustering...')
sc_seu <- spa_cluster(seu_obj = sc_seu,
                      Ks = args$Ks,
                      Kg = args$Kg,
                      alpha = args$sp_alpha,
                      beta = args$sp_beta,
                      target_num = args$target_num,
                      res = args$res,
                      step = args$step,
                      n_seed = args$n_seed)

message('Begin saving data...')
save_data(seu_obj = sc_seu,
          project_name = args$project)


