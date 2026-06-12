library(data.table)
library(fftwtools)
library(geometry) 
library(plot3D)
library(SuperGauss)
library(lhs)
library(RobustGaSP)
library(AIUQ)
library(bioimagetools)
library(EBImage)
library(ijtiff)
library(here)

#file_name = "SP 0.72 AF_500frame"
file_name = "40%_1"
file_path = here("Videos", paste0(file_name, ".tif"))

source(here("Functions", "functions_nonparametric.R"))
source(here("Functions", "functions_DDM.R"))

run_timed_block = T

if (run_timed_block) {
  t_start_total = Sys.time()
  
  M=100 # no. of particles
  mindt = 2 # mindt
  pxsz =0.293  # pixel size 
  
  intensity_ori = bioimagetools::readTIF(file_path, as.is=TRUE)
  
  # check if intensity_ori is normalized
  rng = range(intensity_ori, na.rm = TRUE)
  
  if (!(isTRUE(all.equal(rng[1], 0)) && isTRUE(all.equal(rng[2], 1)))) {
    intensity_ori = (intensity_ori - rng[1]) / (rng[2] - rng[1])
  }
  

  ori_format = 'SST_array' ##real data 
  q_thr = 0.999
  
  if(ori_format=='T_SS_mat'|ori_format=='S_ST_mat'|ori_format=='SST_array'){
    intensity=intensity_format_transform(ori_format=ori_format,intensity_ori,edge=0,exclude=0) ##R
  }
  data_fft=FFT2D(intensity=intensity,pxsz=pxsz,mindt=mindt) 
  index_list=Get_q_ring_loc(data_fft$sz,data_fft$len_q)
  ini_est_list=Get_A_B_ini_num_q_max_est(sz=data_fft$sz,len_q=data_fft$len_q,len_t=data_fft$len_t,
                                         I_q_matrix=data_fft$I_q_matrix,
                                         q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,threshold=q_thr,beta=0) 
                                   
  
  # several subsetting
  num_selected_upper = 20
  q_index_selected_here=unique(ceiling(exp(seq(from=log(1),to=log(ini_est_list$num_q_max),by=ceiling(log(ini_est_list$num_q_max))/num_selected_upper))))
  num_iteration_max = 100
  
  subsample_t = 5
  
  # construct initial parameters using DDM approach ---- 
  param_list_proc = list()
  param_list_proc$mindt = mindt
  param_list_proc$pxsz = pxsz
  param_list_proc$len_q = data_fft$len_q
  param_list_proc$q = data_fft$q
  param_list_proc$I_q_matrix = data_fft$I_q_matrix
  param_list_proc$I_o_q_2 = ini_est_list$I_o_q_2_ori
  param_list_proc$q_ori_ring_loc_index = index_list$q_ori_ring_loc_index

  FFT_model = processing(intensity_ori, param_list_proc)
  GP_model = analysis(FFT_model)


  # tune the number of initial parameters and the scaling factor applied to the slope of initials
  num_estim = 6
  slope_factor = 0.9
  
  
  design_optimization='log_equal_space' 
  model_name='direct_nonparametric'
  initial_param=get_initial_param_nonparametric(model_name=model_name,d_input=data_fft$d_input,num_estim=num_estim,sigma_0_2_ini=min(ini_est_list$I_o_q_2_ori),
                                                design_optimization=design_optimization)
  
  param_ini_DDM = vector(length = num_estim+1)
  param_ini_DDM[1] = log(GP_model$MSD[1]) 
  param_ini_DDM[length(param_ini_DDM)] = initial_param$param_initial[num_estim+1]
  for (i in 1:num_estim){
    param_ini_DDM[i] = param_ini_DDM[1] + slope_factor * (initial_param$input_training[i]-log(mindt))*((log(GP_model$MSD[2])-log(GP_model$MSD[1]))/(log(2*mindt)-log(mindt)))
  }
  
  
  # optimization ---- 
  time_record_stage_1_DDM=system.time({
    m_param_optim_no_restart_direct_DDM = optim(param_ini_DDM, log_lik_param_nonparametric_ori, I_q_cur=data_fft$I_q_matrix,B_cur=NA,
                                                num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
                                                q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
                                                sz=data_fft$sz,len_t=data_fft$len_t,
                                                d_input=data_fft$d_input,q=data_fft$q,model_name="direct_nonparametric",
                                                input_training=initial_param$input_training,
                                                q_index_selected=q_index_selected_here,
                                                subsample_t=subsample_t,
                                                control = list(fnscale=-1,maxit=100), method = 'L-BFGS-B'
    )
    
    theta_est_optim_no_restart_direct_DDM=exp(m_param_optim_no_restart_direct_DDM$par[-length(m_param_optim_no_restart_direct_DDM$par)])
    MSD_est_optim_no_restart_direct_DDM=Get_MSD_nonparametric(theta=theta_est_optim_no_restart_direct_DDM,d_input=data_fft$d_input,input_training=initial_param$input_training,model_name="direct_nonparametric") 
  })[3]
  
  time_record_stage_2_DDM=system.time({
    m_param_optim_no_restart_direct_DDM_fine = optim(m_param_optim_no_restart_direct_DDM$par, log_lik_param_nonparametric_ori, I_q_cur=data_fft$I_q_matrix,B_cur=NA,
                                                     num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
                                                     q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
                                                     sz=data_fft$sz,len_t=data_fft$len_t,
                                                     d_input=data_fft$d_input,q=data_fft$q,model_name="direct_nonparametric",
                                                     input_training=initial_param$input_training,
                                                     q_index_selected=q_index_selected_here,
                                                     control = list(fnscale=-1,maxit=60), method = 'L-BFGS-B'
    )  
    
    theta_est_optim_no_restart_direct_DDM_fine=exp(m_param_optim_no_restart_direct_DDM_fine$par[-length(m_param_optim_no_restart_direct_DDM_fine$par)])
    MSD_est_optim_no_restart_direct_DDM_fine=Get_MSD_nonparametric(theta=theta_est_optim_no_restart_direct_DDM_fine,d_input=data_fft$d_input,input_training=initial_param$input_training,model_name="direct_nonparametric") 
  })[3]
  
 
  t_end_total = Sys.time()
  time_no_uncertainty = as.numeric(difftime(t_end_total, t_start_total, units = "secs"))
}

plot(log10(data_fft$d_input),log10(MSD_est_optim_no_restart_direct_DDM_fine),type='l',col='cyan', lwd = 2,
     xlab = expression(Delta*t), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">")) 

# Uncertainty ---- 
compute_param_uncertainty=T

time_uncertainty = system.time({
  if(compute_param_uncertainty){
    
    param_uq_range=param_uncertainty_nonparametric(param_est=m_param_optim_no_restart_direct_DDM_fine$par,I_q_cur=data_fft$I_q_matrix,B_cur=NA,
                                                   num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
                                                   q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
                                                   sz=data_fft$sz,len_t=data_fft$len_t,
                                                   d_input=data_fft$d_input,q=data_fft$q,model_name=model_name,input_training=initial_param$input_training,M=M,
                                                   q_index_selected=q_index_selected_here,
                                                   num_iteration_max=num_iteration_max, estimation_method = "asymptotics") ##M_B is the number of bootstrap
    
    for(i_p in 1:length(m_param_optim_no_restart_direct_DDM_fine$par)){
      param_uq_range[,i_p]=c(min(m_param_optim_no_restart_direct_DDM_fine$par[i_p],param_uq_range[1,i_p]),
                             max(m_param_optim_no_restart_direct_DDM_fine$par[i_p],param_uq_range[2,i_p]))
      
    }
    
    theta_lower=exp(param_uq_range[1,-length(m_param_optim_no_restart_direct_DDM_fine$par)])
    theta_upper=exp(param_uq_range[2,-length(m_param_optim_no_restart_direct_DDM_fine$par)])
    
    MSD_lower=Get_MSD_nonparametric(theta=theta_lower,d_input=data_fft$d_input,input_training=initial_param$input_training,model_name=model_name)
    MSD_upper=Get_MSD_nonparametric(theta=theta_upper,d_input=data_fft$d_input,input_training=initial_param$input_training,model_name=model_name)
    
  }
  
})[3]




plot(log10(data_fft$d_input),log10(MSD_est_optim_no_restart_direct_DDM_fine),type='l',col='cyan', lwd = 2,
     xlab = expression(Delta*t), ylab = expression("<" * Delta*x^2 * "(" * Delta*t * ")" * ">")) 
lines(log10(data_fft$d_input), log10(MSD_lower),type='l',col='gray')
lines(log10(data_fft$d_input), log10(MSD_upper),type='l',col='gray')


record_data = T

if(record_data){
  d_input_MSD_record=cbind(data_fft$d_input, MSD_est_optim_no_restart_direct_DDM_fine, MSD_lower, MSD_upper)
  colnames(d_input_MSD_record)=c('d_input', 'n_par_AIUQ', 'lower', 'upper')
  dir.create("Results", showWarnings = FALSE)
  
  write.csv(d_input_MSD_record,
            paste0("Results/", file_name, "_AIUQ_nonparametric.csv"),
            row.names = FALSE)
}


