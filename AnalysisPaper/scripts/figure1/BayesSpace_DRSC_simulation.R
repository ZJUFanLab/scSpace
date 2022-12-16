source('../scripts/utils.R')


################## DR.SC ##################
sim_id <- paste0('sim', 1:140)
time_df <- data.frame(matrix(ncol = 1, nrow = 140))
rownames(time_df) <- sim_id; colnames(time_df) <- 'run_time'
res_list <- list()
for (i in 1:length(sim_id)) {
  if(i <= 50){
    data <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/original_data/', sim_id[i], '_sc_data.csv'), row.names = 1)
    meta <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/original_data/', sim_id[i], '_sc_meta.csv'), row.names = 1)
    pseudo_space <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/pseudo_space/pseudo_coords_', sim_id[i], '.csv'), row.names = 1)
  } else if(i > 50){
    data <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/original_data/', sim_id[i], '_sc_data.csv'), row.names = 1)
    meta <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/original_data/', sim_id[i], '_sc_meta.csv'), row.names = 1)
    pseudo_space <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/result/', sim_id[i], '_pseudo_space.csv'), row.names = 1)
  }
  
  meta$row <- pseudo_space[, 1]
  meta$col <- pseudo_space[, 2]
  
  ct_num <- length(unique(meta$Group))
  start_time <- Sys.time()
  res <- run_DRSC(data = data, meta = meta, cluster_num = ct_num)
  end_time <- Sys.time()
  time_df$run_time[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  res_list[[sim_id[i]]] <- res
}

save(res_list, file = '~/workspace/scSpace/data/bayesspace/drsc_sim140.Rdata')
write.csv(time_df, file = '~/workspace/scSpace/data/bayesspace/drsc_sim140_runtime.csv')


################## BayesSpace ##################
sim_id <- paste0('sim', 1:140)
time_df <- data.frame(matrix(ncol = 1, nrow = 140))
rownames(time_df) <- sim_id; colnames(time_df) <- 'run_time'
res_list <- list()
for (i in 1:length(sim_id)) {
  if(i <= 50){
    data <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/original_data/', sim_id[i], '_sc_data.csv'), row.names = 1)
    meta <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/original_data/', sim_id[i], '_sc_meta.csv'), row.names = 1)
    pseudo_space <- read.csv(paste0('~/workspace/scSpace/data/[0]simulated_data/pseudo_space/pseudo_coords_', sim_id[i], '.csv'), row.names = 1)
  } else if(i > 50){
    data <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/original_data/', sim_id[i], '_sc_data.csv'), row.names = 1)
    meta <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/original_data/', sim_id[i], '_sc_meta.csv'), row.names = 1)
    pseudo_space <- read.csv(paste0('~/workspace/scSpace/data/add_subcluster/result/', sim_id[i], '_pseudo_space.csv'), row.names = 1)
  }
  
  meta$row <- pseudo_space[, 1]
  meta$col <- pseudo_space[, 2]
  
  ct_num <- length(unique(meta$Group))
  start_time <- Sys.time()
  res <- run_BayesSpace(data = data, meta = meta, cluster_num = ct_num)
  end_time <- Sys.time()
  time_df$run_time[i] <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  res_list[[sim_id[i]]] <- res
}

save(res_list, file = '~/workspace/scSpace/data/bayesspace/bayesspace_sim140.Rdata')
write.csv(time_df, file = '~/workspace/scSpace/data/bayesspace/bayesspace_sim140_runtime.csv')



