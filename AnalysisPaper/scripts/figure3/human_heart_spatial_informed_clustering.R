# this R script is for scSpace spatial-informed clustering step

source('../scripts/utils.R')
# parameters
sc_data_file <- '~/workspace/spatial_cluster/human_heart/sc_data.csv'
sc_meta_file <- '~/workspace/spatial_cluster/human_heart/sc_meta.csv'
pseudo_coords_file <- '~/workspace/spatial_cluster/human_heart/pseudo_coords.csv'
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

sc_seu <- SeuratObject::CreateSeuratObject(sc_data, meta.data = sc_meta, verbose = FALSE)
sc_seu <- NormalizeData(sc_seu, verbose = FALSE)
sc_seu <- FindVariableFeatures(sc_seu, nfeatures = n_features, verbose = FALSE)
sc_seu <- ScaleData(sc_seu)
sc_seu <- RunPCA(sc_seu)

sc_seu <- spa_cluster(
  sc_seu = sc_seu,
  coord_index = c('pseudo_x', 'pseudo_y'),
  Ks = 10,
  Kg = 20,
  res = 0.35
)



