#' This R script is for pre-processing step.
#' @return This R script returns processed data saved in folders "processed_data"
#' @import original scRNA-seq and spatial transcriptomics data
#' @export: processed data
#' @example demo

# parameters
option_list <- list(optparse::make_option('--project', type = 'character', 
                                          default = 'demo', help = 'project name'),
                    optparse::make_option('--sc_data', type = 'character', 
                                          default = 'demo_sc_data', help = 'sc_data file'),
                    optparse::make_option('--sc_meta', type = 'character', 
                                          default = 'demo_sc_meta', help = 'sc_meta file'),
                    optparse::make_option('--st_data', type = 'character', 
                                          default = 'demo_st_data', help = 'st_data file'),
                    optparse::make_option('--st_meta', type = 'character', 
                                          default = 'demo_st_meta', help = 'st_meta file'),
                    optparse::make_option('--normalize', type = 'logical', 
                                          default = TRUE, help = 'True: normalize the data'),
                    optparse::make_option('--st_type', type = 'character', 
                                          default = 'st', help = 'st: ST, 10x, or slide-seq; image: MERFISH, seqFISH, or STARmap'),
                    optparse::make_option('--n_features', type = 'integer', 
                                          default = 2000, help = 'Number of features to select as top variable features'))


args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

message("=== Data Pre-processing ===")
message('Loading data...')
sc_data <- read.csv(file = paste0('./data/', args$project, '/', args$sc_data, '.csv'), row.names = 1)
sc_meta <- read.csv(file = paste0('./data/', args$project, '/', args$sc_meta, '.csv'), row.names = 1)
st_data <- read.csv(file = paste0('./data/', args$project, '/', args$st_data, '.csv'), row.names = 1)
st_meta <- read.csv(file = paste0('./data/', args$project, '/', args$st_meta, '.csv'), row.names = 1)
message('Data have been loaded...')

rownames(sc_meta) <- colnames(sc_data)
rownames(st_meta) <- colnames(st_data)

inter_gene <- intersect(rownames(sc_data), rownames(st_data))
message(paste0('There are ', length(inter_gene), 
               ' common genes in scRNA-seq and spatial transcriptomics data...'))
st_data <- st_data[inter_gene,]
sc_data <- sc_data[inter_gene,]

sc_seu <- Seurat::CreateSeuratObject(sc_data, meta.data = sc_meta, verbose = FALSE)
st_seu <- Seurat::CreateSeuratObject(st_data, meta.data = st_meta, verbose = FALSE)

if(args$normalize){
  message('Normalizing the data...')
  sc_seu <- Seurat::NormalizeData(sc_seu, verbose = FALSE)
  st_seu <- Seurat::NormalizeData(st_seu, verbose = FALSE)
  message('Data have been normalized...')
} else {
  message('Data have been normalized...')
}

if(args$st_type == 'st'){
  if(nrow(st_data) < args$n_features){
    stop('There are too few genes in ST data, please check again!')
  } else{
    message(paste0('Select ', args$n_features, ' features as top variable features...'))
    sc_seu <- Seurat::FindVariableFeatures(sc_seu, nfeatures = args$n_features, verbose = FALSE)
    st_seu <- Seurat::FindVariableFeatures(st_seu, nfeatures = args$n_features, verbose = FALSE)
    hvg_gene <- unique(c(c(sc_seu@assays$RNA@var.features), c(st_seu@assays$RNA@var.features)))
  }
} else if(args$st_type == 'image'){
  if(nrow(st_data) > args$n_features){
    message(paste0('Select ', args$n_features, ' features as top variable features...'))
    sc_seu <- Seurat::FindVariableFeatures(sc_seu, nfeatures = args$n_features, verbose = FALSE)
    st_seu <- Seurat::FindVariableFeatures(st_seu, nfeatures = args$n_features, verbose = FALSE)
    hvg_gene <- unique(c(c(sc_seu@assays$RNA@var.features), c(st_seu@assays$RNA@var.features)))
  } else{
    hvg_gene <- inter_gene
  }
} else {
  stop('Please select the correct type of spatial transcriptomic data, `st` or `image`!')
}


sc_data_output <- data.frame(sc_seu@assays$RNA@data[hvg_gene, ])
st_data_output <- data.frame(st_seu@assays$RNA@data[hvg_gene, ])
sc_meta_output <- data.frame(sc_seu@meta.data)
st_meta_output <- data.frame(st_seu@meta.data)

message('Data have been pre-processed, begin saving...')

output_file <- paste0('./data/', args$project, '/processed_data')
if (!dir.exists(output_file)){
  dir.create(output_file)}

write.csv(sc_data_output, file = paste0(output_file, '/sc_data_processed.csv'))
write.csv(st_data_output, file = paste0(output_file, '/st_data_processed.csv'))
write.csv(sc_meta_output, file = paste0(output_file, '/sc_meta_processed.csv'))
write.csv(st_meta_output, file = paste0(output_file, '/st_meta_processed.csv'))

message('Pre-processed have been done!')

