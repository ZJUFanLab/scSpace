if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")

if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")

if (!requireNamespace("leidenAlg", quietly = TRUE))
  install.packages("leidenAlg")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('BiocNeighbors', force = T)

if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")
