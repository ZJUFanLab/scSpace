# this R script is for scSpace spatial-informed clustering step

source('../scripts/utils.R')
# parameters
sc_data_file <- '~/workspace/scSpace/data/[9]human_covid19/original_data/sc_data_sampled.csv'
sc_meta_file <- '~/workspace/scSpace/data/[9]human_covid19/original_data/sc_meta_sampled.csv'
pseudo_coords_file <- '~/workspace/scSpace/data/[9]human_covid19/pseudo_space/pseudo_coords.csv'
n_features <- 2000
Ks <- 10
Kg <- 20

# load data
print('load data...')
sc_data <- read.csv(file = sc_data_file, row.names = 1)
sc_meta <- read.csv(file = sc_meta_file, row.names = 1)
pseudo_coords <- read.csv(file = pseudo_coords_file, row.names = 1)

rownames(sc_meta) <- colnames(sc_data)
sc_meta$pseudo_x <- pseudo_coords[,1]
sc_meta$pseudo_y <- pseudo_coords[,2]

mdm_meta <- sc_meta[sc_meta$cell_type_fine == 'Monocyte-derived macrophages' | sc_meta$cell_type_fine == 'Transitioning MDM', ]
mdm_data <- sc_data[, rownames(mdm_meta)]
sc_seu <- SeuratObject::CreateSeuratObject(mdm_data, meta.data = mdm_meta, verbose = FALSE)
sc_seu <- FindVariableFeatures(sc_seu, nfeatures = n_features, verbose = FALSE)
sc_seu <- ScaleData(sc_seu)
sc_seu <- RunPCA(sc_seu)

sc_seu <- spa_cluster(
  sc_seu = sc_seu,
  coord_index = c('pseudo_x', 'pseudo_y'),
  Ks = 10,
  Kg = 20,
  res = 0.6
)
