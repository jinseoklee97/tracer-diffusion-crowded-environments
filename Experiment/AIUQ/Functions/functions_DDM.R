
processing <- function(I, params = list()) {
  
  model = list()
  
  if (is.null(params$pxsz)) model$pxsz = 1 else model$pxsz = params$pxsz
  if (is.null(params$nframes)) model$nframes = dim(I)[3] else model$nframes = params$nframes
  if (is.null(params$nx)) model$nx = dim(I)[1] else model$nx = params$nx
  if (is.null(params$ny)) model$ny = dim(I)[2] else model$ny = params$ny
  if (!is.null(params$mindt)) model$mindt = params$mindt else stop("Please specify mininum lag time.")
  if (!is.null(params$q)) model$q = params$q else stop("Please specify q.")
  if (!is.null(params$len_q)) model$len_q = params$len_q else stop("Please specify len_q.")
  if (!is.null(params$I_q_matrix)) model$I_q_matrix = params$I_q_matrix else stop("Please specify I_q_matrix.")
  if (!is.null(params$q_ori_ring_loc_index)) model$q_ori_ring_loc_index = params$q_ori_ring_loc_index else stop("Please specify q_ori_ring_loc_index.")
  if (!is.null(params$index_dt_selected)) model$index_dt_selected = params$index_dt_selected else model$index_dt_selected = 1:2
  
  if (!is.null(params$I_o_q_2)) {
    model$I_o_q_2 = params$I_o_q_2
    model$I_o_q_2_on = 0
  } else {
    model$I_o_q_2_on = 1
  }
  
  stopifnot(is.numeric(model$pxsz), is.numeric(model$mindt), 
            is.numeric(model$nframes), is.numeric(model$nx), is.numeric(model$ny))
  
  model$dt = model$mindt * (1:(model$nframes - 1))
  model$ndt = (model$nframes - 1):1
  
  if (min(model$nx, model$ny) %% 2 == 0) {
    model$sz = min(model$nx, model$ny) - 1
  } else {
    model$sz = min(model$nx, model$ny)
  }
  
  
  model$Dqt = matrix(nrow = model$len_q, ncol = length(model$index_dt_selected))
  
  for (idx in seq_along(model$index_dt_selected)) {
    k = model$index_dt_selected[idx]
    
    avg = 0
    for (j in 1:model$ndt[k]) {
      avg = avg + abs(model$I_q_matrix[,j+k] - model$I_q_matrix[,j])^2/model$sz^2
    }
    avg = avg/model$ndt[k]
    
    avg_over_q_per_dt = numeric(length = model$len_q)
    for(i in 1:model$len_q){
      avg_over_q_per_dt[i]=mean(avg[model$q_ori_ring_loc_index[[i]]])
    }
    
    model$Dqt[,idx] = avg_over_q_per_dt
  }
  
  
  
  ## I_o_q_2 
  if (model$I_o_q_2_on == 1) {
    avg_I_2_over_time = 0
    
    for (t in 1:model$nframes) {
      avg_I_2_over_time = avg_I_2_over_time + abs(model$I_q_matrix[,t])^2/model$sz^2
    }
    avg_I_2_over_time = avg_I_2_over_time/model$nframes
    
    model$I_o_q_2 = numeric(length = model$len_q)
    for(i in 1:model$len_q){
      model$I_o_q_2[i]=mean(avg_I_2_over_time[model$q_ori_ring_loc_index[[i]]])
    }
  }
  return(model)
}





## analysis no truncation 
analysis <- function(params = list()) {
  
  model = list()
  
  if (!is.null(params$dt)) model$dt = params$dt
  
  if (!is.null(params$index_dt_selected)) {
    model$index_dt_selected = params$index_dt_selected
  } else {
    stop("Index_dt_delected cannot be NULL.")
  }
  
  # Dqt 
  ddm_ori = params$Dqt
  q_ori = params$q
  I_o_2_ori = matrix(params$I_o_q_2, nrow=1)
  
  # A_hat, B_hat 
  if (!is.null(params$sigma_2_0_hat)) {
    model$sigma_2_0_hat = params$sigma_2_0_hat
    model$B_hat = 2 * model$sigma_2_0_hat
  } else {
    model$B_hat = min(mean(ddm_ori[nrow(ddm_ori),]), min(ddm_ori[,1]))
    model$sigma_2_0_hat = model$B_hat/2
  }
  model$A_hat = 2 * (I_o_2_ori - model$sigma_2_0_hat)
  
  # q truncation
  q_max_num = min(100, length(q_ori))
  q = q_ori[1:q_max_num]
  A_hat = as.numeric(model$A_hat)[1:q_max_num]
  
  
  model$Dqt_q_max = ddm_ori[1:q_max_num, , drop=FALSE]
  
  n_q = length(q)
  n_time = ncol(model$Dqt_q_max)
  
  model$index_cut = rep(0, n_q)
  q_index_selected = c()
  
  # original selection conditions
  n_q_min = 10
  n_q_max = min(n_q, 83)
  
  
  Dqt_select_dt = model$Dqt_q_max[, model$index_dt_selected, drop=FALSE]
  
  # f(q,t)
  fqt = 1 - (Dqt_select_dt - model$B_hat) / A_hat
  fqt[fqt <= 0] = NA
  
  # MSD(q,t)
  Qmat = matrix(q, nrow=length(q), ncol=length(model$index_dt_selected))
  MSD_q_dt = 4 * log(1/fqt) / (Qmat^2)
  
  # apply filtering index_cut
  MSD = rep(NA, length(model$index_dt_selected))
  for (ti in model$index_dt_selected) {
    MSD[ti] = median(MSD_q_dt[, ti], na.rm=TRUE)
  }
  
  if (length(na.omit(MSD)) < 2) {
    stop("Fewer than two valid MSD values were obtained. Please increase the number of lag times selected.")
  }
  
  model$MSD = MSD
  
  return(model)
}
