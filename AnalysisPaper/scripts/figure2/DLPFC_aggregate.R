slice_id <- c(151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676)
data_list <- list(
  test_slice = slice_id, 
  train_slice = list(c(151508, 151509, 151510), c(151507, 151509, 151510), c(151507, 151508, 151510), c(151507, 151508, 151509),
                     c(151670, 151671, 151672), c(151669, 151671, 151672), c(151669, 151670, 151672), c(151669, 151670, 151671),
                     c(151674, 151675, 151676), c(151673, 151675, 151676), c(151673, 151674, 151676), c(151673, 151674, 151675)))


# pseudo space aggregate
for (i in 1:length(slice_id)) {
  test_slice <- data_list$test_slice[i]
  original_coord <- read.csv(paste0('~/workspace/scSpace/data/DLPFC_new/data/', test_slice, '_meta.csv'), row.names = 1)
  
  train_list_tmp <- data_list$train_slice[[i]]
  
  p1 <- read.csv(paste0('~/workspace/scSpace/data/DLPFC_new/result/', 
                        test_slice, '_pseudo_space(trained_on_', train_list_tmp[1], ').csv'), row.names = 1)
  p2 <- read.csv(paste0('~/workspace/scSpace/data/DLPFC_new/result/', 
                        test_slice, '_pseudo_space(trained_on_', train_list_tmp[2], ').csv'), row.names = 1)
  p3 <- read.csv(paste0('~/workspace/scSpace/data/DLPFC_new/result/', 
                        test_slice, '_pseudo_space(trained_on_', train_list_tmp[3], ').csv'), row.names = 1)
  
  p_new <- (p1 + p2 + p3) / 3
  
  original_coord$pseudo_space1 <- p_new$X0
  original_coord$pseudo_space2 <- p_new$X1
  
  write.csv(original_coord, file = paste0('~/workspace/scSpace/data/DLPFC_new/result/', test_slice, '_pseudo_space_aggregate.csv'))
  
}