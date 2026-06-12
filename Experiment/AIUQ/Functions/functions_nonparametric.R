###function for ab initio SAM

##
intensity_format_transform<-function(ori_format,intensity,edge=0,exclude=0){
  if(ori_format=='SST_array'){ ##two space indices and one time indice, tiff
    #sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    if(dim(intensity)[1]%%2==0){ ##is it correct if it is a rectangle
      sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = dim(intensity)[1]
    }
    len_t = dim(intensity)[3] #relevant to format
    
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_transform[,i] = as.vector(intensity[1:sz,1:sz,i])  # 1:479; 481:480*2-1,...
    }
    #mid_size = dim(intensity)[1]/2 ##is it useful? 
  }else if(ori_format=='S_ST_mat'){ ##time,  space and time as a matrix, matlab simulation 
    intensity=as.matrix(intensity)
    
    #sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    if(dim(intensity)[1]%%2==0){
      sz = dim(intensity)[1]-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = dim(intensity)[1]
    }
    
    len_t = dim(intensity)[2]/dim(intensity)[1] #relevant to format
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_transform[,i]= as.vector(intensity[1:sz,(1+(sz+1)*(i-1)):(sz+(sz+1)*(i-1))]) 
    }
    
    
  }else if(ori_format=='T_SS_mat'){ ##simulation for R, square matrix, need to think about 
    intensity=as.matrix(intensity)
    
    if(sqrt(dim(intensity)[2])%%2==0){
      sz = sqrt(dim(intensity)[2])-1 #pixel dimension; total # of pixels in each image = sz^2
    }else{
      sz = sqrt(dim(intensity)[2])
    }
    len_t = dim(intensity)[1]
    
    #I_q = array(NaN,dim = c(sz,sz,len_t)) #Fourier transformed intensity in t, not dt
    intensity_transform=matrix(NA,sz^2,len_t)
    for(i in 1:len_t){
      intensity_mat=matrix(intensity[i,],sz+1,sz+1 ) 
      intensity_transform[,i]= as.vector(intensity_mat[1:sz,1:sz] )     
    }
  }
  return(intensity_transform)
}

##having really use edge and exclude 
FFT2D<-function(intensity,pxsz,mindt){
  sz=sqrt(dim(intensity)[1])
  len_t=dim(intensity)[2]
  I_q_matrix=matrix(NA,sz^2,len_t)
  for(i in 1:len_t){
    I_q_matrix[,i] = as.vector(fftw2d(matrix(intensity[,i],sz,sz)))  #
  }
    
  
  ans_list=list()
  ans_list$sz=sz
  ans_list$len_q=length(1:((sz-1)/2))
  ans_list$len_t=len_t
  ans_list$I_q_matrix=I_q_matrix
  ans_list$q= (1:((sz-1)/2))*2*pi/(sz*pxsz)
  ans_list$input= mindt*(1:(len_t))  ##t
  ans_list$d_input = ans_list$input[1:length(ans_list$input)]-ans_list$input[1] ##delta t, including zero 
  
  return(ans_list)
}


## Define fftshift: 
# Function that used to swap the first quadrant with the third, and the second quadrant
# with the forth. So the zero-frequency component will be at the center of the array.
fftshift <- function(input_matrix, dim = -1) {
  
  rows <- dim(input_matrix)[1]    
  cols <- dim(input_matrix)[2]    
  
  swap_up_down <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
  }
  
  swap_up_down_reverse <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half):rows), (1:cols)], input_matrix[1:(rows_half-1), (1:cols)]))
  }
  
  swap_left_right_reverse <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half):cols)], input_matrix[1:rows, 1:(cols_half-1)]))
  }
  
  
  if (dim == -1) {
    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))
  }
  else if (dim == 1) {
    return(swap_up_down(input_matrix))
  }
  else if (dim == 2) {
    return(swap_left_right(input_matrix))
  }else if(dim == 3){
    input_matrix <- swap_up_down_reverse(input_matrix)
    return(swap_left_right_reverse(input_matrix))
    
  }
  else {
    stop("Invalid dimension parameter")
  }
}

Get_q_ring_loc<-function(sz,len_q){
  v = (-(sz-1)/2):((sz-1)/2)
  x = matrix(rep(v,each = sz), byrow = FALSE,nrow=sz)
  y = matrix(rep(v,each = sz), byrow = TRUE,nrow=sz)
  
  theta_q = geometry::cart2pol(x, y)
  #dim(theta_q)
  q_ring_num = theta_q[,(sz+1):dim(theta_q)[2]]
  q_ring_num = round(q_ring_num) ###after transformed q number
  #len_q = length(1:((sz-1)/2))
  nq_index = vector(mode = "list")
  for(i in 1:len_q){
    nq_index[[i]] = which(q_ring_num==i)
  }
  #q_ring_num[nq_index[[3]]] #check whether index is correct 
  #> which(q_ring_num==0)
  #[1] 4901
  
  q_ori_ring_loc=fftshift(q_ring_num, dim = 3)  ###this is the location for original 
  #> which(q_ori_ring_loc==0)
  #[1] 1
  #image2D(q_ori_ring_loc)
  ###only use 1 zero
  # if(length(which(q_ring_num==0))>1){
  #   index_0_num=floor(length(q_ring_num)/2)+1
  #   q_ring_num_0_index=which(q_ring_num==0)
  #   selected_0_index= which(q_ring_num_0_index==index_0_num)
  #   q_ring_num[q_ring_num_0_index[-selected_0_index]]=1 ###change to 1
  #   
  #   q_ori_ring_loc_index=which(q_ori_ring_loc==0)
  #   q_ori_ring_loc[q_ori_ring_loc_index[-1]]=1
  # }
    
  #max(q_ring_num)
  q_ori_ring_loc_index=as.list(1:len_q)
  total_q_ori_ring_loc_index=NULL
  for(i in 1:len_q){
    q_ori_ring_loc_index[[i]]=which(q_ori_ring_loc==i)
    total_q_ori_ring_loc_index=c(total_q_ori_ring_loc_index,  q_ori_ring_loc_index[[i]])
  }
  
  

  #image2D(q_ring_num)
  #image2D(q_ori_ring_loc)
  ans_list=list()
  ans_list$q_ring_loc=q_ring_num
  ans_list$q_ori_ring_loc=q_ori_ring_loc
  ans_list$q_ori_ring_loc_index=q_ori_ring_loc_index ##original index 
  ans_list$total_q_ori_ring_loc_index=total_q_ori_ring_loc_index ##original index 
  

  return(ans_list)
}


##q_ori_ring_loc_index already contains isotropic information 
Get_A_B_ini_num_q_max_est<-function(sz,len_q,len_t, 
                                    I_q_matrix,q_ori_ring_loc_index,threshold=0.999,beta=0.5,anisotropic=T){
  
  #avg_I_2_ori=matrix(0,nrow = sz, ncol = sz)
  avg_I_2_ori=0
  
  for(i in 1:len_t){
    #avg_I_2_ori=avg_I_2_ori+abs(I_q[,,i])^2/sz^2
    avg_I_2_ori=avg_I_2_ori+abs(I_q_matrix[,i])^2/sz^2

  }
  avg_I_2_ori = avg_I_2_ori/len_t
  
  I_o_q_2_ori=rep(NA,len_q)
  for(i in 1:len_q){
    I_o_q_2_ori[i]=mean(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
  }
  
  ##change to minimum?
  # I_o_q_2_ori_last = I_o_q_2_ori[len_q]
  # 
  # B_est_ini = 2*I_o_q_2_ori_last
  # A_est_ini = abs(2*(I_o_q_2_ori - I_o_q_2_ori_last))

  # Feb 9, 2026
  I_o_q_2_ori_min = min(I_o_q_2_ori)

  B_est_ini = max(2*abs(I_o_q_2_ori_min),10^{-20})
  A_est_ini = 2*(I_o_q_2_ori - B_est_ini/2)

  # Mar 07, 2026
  q_reach_thr = which(cumsum(A_est_ini) / sum(A_est_ini) >= threshold)[1]
  
  tail_start = max(ceiling(0.9 * len_q), q_reach_thr + 1)
  
  if (tail_start > len_q) {
    num_q_max = q_reach_thr
  } else {
    tail_idx_select = which(A_est_ini[tail_start:len_q] >= 0.001) 
    
    if (length(tail_idx_select) == 0) {
      num_q_max = q_reach_thr
    } else {
      num_q_max = tail_start - 1 + max(tail_idx_select)
    }
  }
  
  
  if(num_q_max/len_q<=beta){
    num_q_max=ceiling(beta*len_q)
  }
  
  ##get unique index, this maybe improved
  q_ori_ring_loc_unique_index=as.list(1:len_q)
  for(i in 1:len_q){
    unique_val=unique(avg_I_2_ori[q_ori_ring_loc_index[[i]]])
    unique_val=unique_val[1:(length(q_ori_ring_loc_index[[i]])/2)]
    index_selected=NULL
    for(j in 1:length(unique_val)){
       index_selected=c(index_selected,which(avg_I_2_ori==unique_val[j])[1])
    }
    q_ori_ring_loc_unique_index[[i]]=index_selected ###q_ori_ring_loc_index[[i]][index_selected]
  }
  
  total_q_ori_ring_loc_unique_index=NULL
  for(i in 1:len_q){
    total_q_ori_ring_loc_unique_index=c(total_q_ori_ring_loc_unique_index,     q_ori_ring_loc_unique_index[[i]])
  }
  
  ans_list=list()
  ans_list$A_est_ini=A_est_ini
  ans_list$B_est_ini=B_est_ini
  ans_list$num_q_max=num_q_max
  ans_list$I_o_q_2_ori=I_o_q_2_ori
  ans_list$avg_I_2_ori=avg_I_2_ori
  ans_list$q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index
  ans_list$total_q_ori_ring_loc_unique_index=total_q_ori_ring_loc_unique_index
  
  #for anisotropic
  if(anisotropic){
    #image2D(t(q_ori_ring_loc[1:len_q,1:len_q]))
    q1_unique_index=as.list(1:len_q)
    q2_unique_index=as.list(1:len_q)
    
    #sub_q_ori_ring_loc=q_ori_ring_loc[1:len_q,1:len_q]
    for(i in 1:len_q){
      #index_here=which(sub_q_ori_ring_loc==i)
      #col_here=floor(index_here/len_q)+1
      #row_here=index_here%%len_q
      index_here=q_ori_ring_loc_unique_index[[i]]
      total_num_unique_index_here=length( q_ori_ring_loc_unique_index[[i]])
      q1_unique_index[[i]]=q2_unique_index[[i]]=rep(NA,total_num_unique_index_here)
      for(j in 1:(total_num_unique_index_here)){
        q2_unique_index[[i]][j]=floor((index_here[j]-1)/sz) ##could contain zero
        
        left_here=(index_here[j]-1)%%sz
        if(left_here<=len_q){
          q1_unique_index[[i]][j]=(index_here[j]-1)%%sz ##could contain zero
        }else{
          q1_unique_index[[i]][j]=sz-(index_here[j]-1)%%sz-1 ##could contain zero
        }
        # if(j<=length(row_here)){
        #   q1_unique_index[[i]][j]=row_here[j]
        #   q2_unique_index[[i]][j]=col_here[j]
        # }else{
        #   j_here=2*length(row_here)-j ###length(row_here)-(j-)
        #   q1_unique_index[[i]][j]=row_here[j_here]
        #   q2_unique_index[[i]][j]=col_here[j_here]
        # }
      }
    }
    ans_list$q1_unique_index=q1_unique_index
    ans_list$q2_unique_index=q2_unique_index
    
  }
  
  

  return(ans_list)
}

Get_MSD_nonparametric<-function(theta,d_input,input_training,model_name){
  if(model_name=='direct_nonparametric'){
    ##range.par can be adjust
    #100*d_input[2]
    range_par=(max(input_training)-min(input_training))
    
    ##trend
    #X=cbind(rep(1,length(input_training)),as.numeric(input_training))
    #X_testing=cbind(rep(1,length(d_input[-1])),as.numeric(log(d_input[-1])))
    #
    #m_direct=rgasp(design=input_training,response=log(theta),
    #               range.par = range_par,trend=as.matrix(X)) ##input_training is log
    #m_direct_log_MSD=predict(m_direct,testing_input=as.matrix(log(d_input[-1])),testing_trend=X_testing)$mean

    ##no trend
    m_direct=rgasp(design=input_training,response=log(theta),range.par = range_par) ##input_training is log 
    
    ##small nugget
    #m_direct=rgasp(design=input_training,response=log(theta),range.par = range_par,nugget=10^{-4}) ##input_training is log
    
    m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])))$mean
    
    MSD=c(0,exp(m_direct_log_MSD))

  }
  # else if(model_name=='direct_additive_nonparametric'){
  #   ##range.par can be adjust
  #   #100*d_input[2]
  #   range_par=(max(input_training)-min(input_training))
  #   output=log(cumsum(theta))
  #   m_direct=rgasp(design=input_training,response=output,range.par = range_par) ##input_training is log
  #   m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])))$mean
  #   MSD=c(0,exp(m_direct_log_MSD))
  # }
  else if(model_name=='varying_power_nonparametric'){
    #if(input_training[1]==0){
      ###this is just for time to be 0 to 1, we can always input 0 to 1 and transform later
      range_par=(max(input_training)-min(input_training)) ##input training is in the log space
  
      alpha_training=2*theta[-1]/(1+theta[-1])
      m_varying_power=rgasp(design=input_training[-1],response=(alpha_training),range.par = range_par) ##input_training is log
      m_varying_power_pred_alpha=predict(m_varying_power,as.matrix(log(d_input[-1])))$mean
      log_msd=log(theta[1])+(m_varying_power_pred_alpha)*log(d_input[-1])
  
      MSD=c(0,exp(log_msd))
    #}else{ ##need to make sure no 0 is in training
      
   # }
    
  }
  else if(model_name=='varying_power_nonparametric_msd'){
    range_par=(max(input_training)-min(input_training)) ##input training is in the log space
    log_MSD_training=rep(NA,length(input_training))
    log_MSD_training[1]=log(theta[1])
    
    log_MSD_training[2:length(input_training)]= log_MSD_training[1]+2*theta[-1]/(1+theta[-1])*(input_training[-1])

    m_direct=rgasp(design=input_training,response=(log_MSD_training),range.par = range_par) ##input_training is log
    m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])))$mean
    
    #X=cbind(rep(1,length(input_training)),input_training)
    #X_testing=cbind(rep(1,length(d_input[-1])),log(d_input[-1]))
    
    
    #m_direct=rgasp(design=input_training,response=(log_MSD_training),
    #               range.par = range_par*0.1,trend=X) ##input_training is log
    #m_direct_log_MSD=predict(m_direct,as.matrix(log(d_input[-1])),testing_trend=X_testing)$mean
    
    MSD=c(0,exp(m_direct_log_MSD))
    
  }
  return(MSD)
}


get_circ_firstrow = function(p, first_row){
  second_half = rev(first_row[2:p])
  circ_firstrow = c(first_row, second_half)
  return(circ_firstrow)
}

###changed July 15, 2025
log_lik_param_nonparametric_ori<-function(param,I_q_cur,B_cur,num_q_cur,I_o_q_2_ori,
                                          d_input,q_ori_ring_loc_unique_index,sz,len_t,q,model_name,input_training,
                                          q_index_selected=1:num_q_cur,subsample_t=1,whitening=F
){ #recompute_num_q_cur=F
  #print(param)
  #theta=exp(param)
  p=length(param)-1
  #theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
  
  if(whitening){
    range_par=(max(input_training)-min(input_training))
    ##no trend
    m_direct=rgasp(design=input_training,response=param[-(p+1)],range.par = range_par) ##input_training is log
    log_theta=m_direct@L%*%param[-(p+1)]
    theta=exp(log_theta)
  }else{
    theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
    
  }
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }
  
  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
  
  
  
  ##the model is defined by MSD
  #MSD=Get_MSD_nonparametric(theta,d_input,input_training,model_name)
  
  subsample_d_index=seq(1,length(d_input),subsample_t)
  d_input_subsample=d_input[subsample_d_index]
  
  MSD=Get_MSD_nonparametric(theta,d_input_subsample,input_training,model_name)
  
  
  log_lik_sum=0
  len_t_subsample=length(MSD)
  
  NTz <- NormalToeplitz$new(len_t_subsample) ##maybe create it outside?
  
  eta=B_cur/4 ##nugget
  
  num_q_select=length(q_index_selected)
  
  i_log_lik_record=rep(NA,num_q_select)
  
  for(i_q_selected in q_index_selected){
    #print(q_selected)
    
    #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    
    
    q_selected=q[i_q_selected]
    #beta_q = (D*q[i_q_selected]^2)
    sigma_2=A_cur[i_q_selected]/4
    
    acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
    acf[1] = acf[1]+eta
    acf=as.numeric(acf)
    
    i_log_lik_record[i_q_selected] = sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    
    # if (is.nan(i_log_lik)) {
    #   ## try new estimator: get the first row of the new estimator
    #   circ_mat = generate_circ_mat(length(acf), acf)
    #   lambda = fft(circ_mat[1,]) ## eigen values
    #   lambda_truncated = ifelse(Re(lambda) <= 0, 10^{-4}, lambda)
    #   firstrow_reconstructed = Re(fft(lambda_truncated, inverse = TRUE) / (2*length(acf) - 1)) # the first row of reconstructed circulant matrix
    #   i_log_lik = sum(NTz$logdens(z = t(output_re), acf = firstrow_reconstructed[1:length(acf)])) + sum(NTz$logdens(z = t(output_im), acf = firstrow_reconstructed[1:length(acf)]))
    # }
    
    #log_lik_sum=log_lik_sum+sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    log_lik_sum = log_lik_sum +  i_log_lik_record[i_q_selected]
    
  }
  
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index))*log(2*pi) ##add 2pi
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index)*len_t)*log(2*pi) ##add 2pi? we may not need to add these constant, we should remove them later
  if(is.nan(log_lik_sum)){
    log_lik_sum=-10^15
    #log_lik_sum=-10^15-10^15*runif(1)
    
    
  }
  
  #print(c((param),log_lik_sum))
  #print(c(log_lik_sum))
  return(log_lik_sum)
}



## For BFGS More-Thuente line search (minimization problem)
# neg_log_lik_param_nonparametric_ori<-function(param,I_q_cur,B_cur,num_q_cur,I_o_q_2_ori,
#                                               d_input,q_ori_ring_loc_unique_index,sz,len_t,q,model_name,input_training,
#                                               q_index_selected=1:num_q_cur
# ){ 
#   p=length(param)-1
#   theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
#   if(is.na(B_cur)){ ##this fix the dimension 
#     sigma_2_0_hat=exp(param[p+1]) ##noise
#     B_cur=2*sigma_2_0_hat
#   }
#   
#   A_cur = abs(2*(I_o_q_2_ori - B_cur/2)) 
#   
#   
#   
#   ##the model is defined by MSD
#   MSD=Get_MSD_nonparametric(theta,d_input,input_training,model_name)
#   log_lik_sum=0
#   
#   NTz <- NormalToeplitz$new(len_t) ##maybe create it outside?
#   
#   eta=B_cur/4 ##nugget
#   
#   num_q_select=length(q_index_selected)
#   
#   i_log_lik_record=rep(NA,num_q_select)
#   
#   for(i_q_selected in q_index_selected){
#     #print(q_selected)
#     
#     output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
#     output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
#     
#     
#     q_selected=q[i_q_selected]
#     #beta_q = (D*q[i_q_selected]^2) 
#     sigma_2=A_cur[i_q_selected]/4
#     
#     acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
#     acf[1] = acf[1]+eta
#     acf=as.numeric(acf)
#     
#     i_log_lik_record[i_q_selected] = sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
#     
#     # if (is.nan(i_log_lik)) {
#     #   ## try new estimator: get the first row of the new estimator
#     #   circ_mat = generate_circ_mat(length(acf), acf)
#     #   lambda = fft(circ_mat[1,]) ## eigen values
#     #   lambda_truncated = ifelse(Re(lambda) <= 0, 10^{-4}, lambda)
#     #   firstrow_reconstructed = Re(fft(lambda_truncated, inverse = TRUE) / (2*length(acf) - 1)) # the first row of reconstructed circulant matrix 
#     #   i_log_lik = sum(NTz$logdens(z = t(output_re), acf = firstrow_reconstructed[1:length(acf)])) + sum(NTz$logdens(z = t(output_im), acf = firstrow_reconstructed[1:length(acf)]))
#     # }
#     
#     #log_lik_sum=log_lik_sum+sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
#     log_lik_sum = log_lik_sum +  i_log_lik_record[i_q_selected]
#     
#   }
#   
#   #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index))*log(2*pi) ##add 2pi
#   #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index)*len_t)*log(2*pi) ##add 2pi? we may not need to add these constant, we should remove them later
#   if(is.nan(log_lik_sum)){
#     log_lik_sum=-10^15
#     #log_lik_sum=-10^15-10^15*runif(1)
#     
#     
#   }
#   
#   #print(c((param),log_lik_sum))
#   print(c(-log_lik_sum))
#   return(-log_lik_sum) 
# }

neg_log_lik_param_nonparametric_ori<-function(param,I_q_cur,B_cur,num_q_cur,I_o_q_2_ori,
                                          d_input,q_ori_ring_loc_unique_index,sz,len_t,q,model_name,input_training,
                                          q_index_selected=1:num_q_cur,subsample_t=1,whitening=F
){ #recompute_num_q_cur=F
  #print(param)
  #theta=exp(param)
  p=length(param)-1
  #theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
  
  if(whitening){
    range_par=(max(input_training)-min(input_training))
    ##no trend
    m_direct=rgasp(design=input_training,response=param[-(p+1)],range.par = range_par) ##input_training is log
    log_theta=m_direct@L%*%param[-(p+1)]
    theta=exp(log_theta)
  }else{
    theta=exp(param[-(p+1)]) ##first p parameters are parameters in ISF
    
  }
  if(is.na(B_cur)){ ##this fix the dimension
    sigma_2_0_hat=exp(param[p+1]) ##noise
    B_cur=2*sigma_2_0_hat
  }
  
  A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
  
  
  
  ##the model is defined by MSD
  #MSD=Get_MSD_nonparametric(theta,d_input,input_training,model_name)
  
  subsample_d_index=seq(1,length(d_input),subsample_t)
  d_input_subsample=d_input[subsample_d_index]
  
  MSD=Get_MSD_nonparametric(theta,d_input_subsample,input_training,model_name)
  
  
  log_lik_sum=0
  len_t_subsample=length(MSD)
  
  NTz <- NormalToeplitz$new(len_t_subsample) ##maybe create it outside?
  
  eta=B_cur/4 ##nugget
  
  num_q_select=length(q_index_selected)
  
  i_log_lik_record=rep(NA,num_q_select)
  
  for(i_q_selected in q_index_selected){
    #print(q_selected)
    
    #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
    output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],subsample_d_index])/(sz)
    
    
    q_selected=q[i_q_selected]
    #beta_q = (D*q[i_q_selected]^2)
    sigma_2=A_cur[i_q_selected]/4
    
    acf = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
    acf[1] = acf[1]+eta
    acf=as.numeric(acf)
    
    i_log_lik_record[i_q_selected] = sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    
    # if (is.nan(i_log_lik)) {
    #   ## try new estimator: get the first row of the new estimator
    #   circ_mat = generate_circ_mat(length(acf), acf)
    #   lambda = fft(circ_mat[1,]) ## eigen values
    #   lambda_truncated = ifelse(Re(lambda) <= 0, 10^{-4}, lambda)
    #   firstrow_reconstructed = Re(fft(lambda_truncated, inverse = TRUE) / (2*length(acf) - 1)) # the first row of reconstructed circulant matrix
    #   i_log_lik = sum(NTz$logdens(z = t(output_re), acf = firstrow_reconstructed[1:length(acf)])) + sum(NTz$logdens(z = t(output_im), acf = firstrow_reconstructed[1:length(acf)]))
    # }
    
    #log_lik_sum=log_lik_sum+sum(NTz$logdens(z = t(output_re), acf = acf))+sum(NTz$logdens(z = t(output_im), acf = acf))
    log_lik_sum = log_lik_sum +  i_log_lik_record[i_q_selected]
    
  }
  
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index))*log(2*pi) ##add 2pi
  #log_lik_sum=log_lik_sum-0.5*sum(lengths(q_ori_ring_loc_unique_index)*len_t)*log(2*pi) ##add 2pi? we may not need to add these constant, we should remove them later
  if(is.nan(log_lik_sum)){
    log_lik_sum=-10^15
    #log_lik_sum=-10^15-10^15*runif(1)
    
    
  }
  
  #print(c((param),log_lik_sum))
  #print(c(log_lik_sum))
  return(-log_lik_sum)
}



###two initial starts?

# get_initial_param_nonparametric <- function(model_name,d_input,num_estim=6,sigma_0_2_ini=NA){
#   
#   d_int=log(data_fft$d_input[length(data_fft$d_input)])/(num_estim-1)
#   input_index=rep(NA,num_estim)
#   input_index[1]=2
#   input_seq_fit=d_int*(1:(num_estim-1))
#   for(i in 1:(num_estim-1)){
#     input_index[i+1]=which(abs(log(data_fft$d_input)-input_seq_fit[i])==min(abs(log(data_fft$d_input)-input_seq_fit[i])))
#   }
#   input_training=log(data_fft$d_input[input_index])
#   
#   if(model_name=='direct_nonparametric'){
#     #param_initial=matrix(NA,1,num_estim+1) #include B
#     #param_initial[1,]=log(c(1,sigma_0_2_ini))#method='L-BFGS-B',
#     
#     param_initial=log(c((exp(input_training)*0.1),sigma_0_2_ini)) ##
#   }else if(model_name=='varying_power_nonparametric'){
#     ##this is somewhat ad-hoc which start from a FBM with some random parameter
#     param_initial1=log(c(0.1,rep(.1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
#     #param_initial2=log(c(0.5,rep(0.5,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
#     param_initial2=log(c(0.1,rep(1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
#     
#     param_initial=cbind(param_initial1,param_initial2)
#   }else if(model_name=='varying_power_msd_nonparametric'){
#     param_initial=log(c(0.1,rep(1,length(input_training)-1),sigma_0_2_ini)) ##ad hoc
#     
#   }
#   return.list=list()
#   return.list$input_training=input_training
#   return.list$param_initial=param_initial
#   
#   return(return.list)
# }

get_initial_param_nonparametric <- function(model_name,d_input,num_estim=6,sigma_0_2_ini=NA,design_optimization="log_equal_space"){
  if(design_optimization=='equal_space'){
    input_index=rep(NA,num_estim)
    input_index[1]=2 ##perhaps start from the beginning, this is like half or full
    d_int=floor(length(data_fft$d_input)/(num_estim-1))
    input_index[2:num_estim]=floor(1:(num_estim-1))*d_int
  }else if(design_optimization=='log_equal_space'){ ##log equal space 
    d_int=(log(data_fft$d_input[length(data_fft$d_input)])-log(data_fft$d_input[2]))/(num_estim-1)
    
    input_index=rep(NA,num_estim)
    
    input_index[1]=2 ##perhaps start from the beginning, first one
    input_seq_fit=log(data_fft$d_input[2])+d_int*(1:(num_estim)) 
    for(i in 2:(length(input_index))){
      input_index[i]=which(abs(log(data_fft$d_input)-input_seq_fit[i-1])==min(abs(log(data_fft$d_input)-input_seq_fit[i-1]))) 
    }
  }else if(design_optimization=='log_equal_space_middle'){ ##middle 
    skip=8
    d_int=(log(data_fft$d_input[length(data_fft$d_input)])-log(data_fft$d_input[2]))/(num_estim-1+skip)
    
    input_index=rep(NA,num_estim)
    
    input_index[1]=2 ##perhaps start from the beginning, first one
    input_seq_fit=log(data_fft$d_input[2])+d_int*(skip:(num_estim+skip-1)) ##the subsetting is like skip some in the log space? this is almost like equally space then, perhaps try equally space
    for(i in 2:(length(input_index))){
      input_index[i]=which(abs(log(data_fft$d_input)-input_seq_fit[i-1])==min(abs(log(data_fft$d_input)-input_seq_fit[i-1])))
    }
  }
  input_training=log(data_fft$d_input[input_index])
  
  if(model_name=='direct_nonparametric'){
    #param_initial=matrix(NA,1,num_estim+1) #include B
    #param_initial[1,]=log(c(1,sigma_0_2_ini))#method='L-BFGS-B',
    
    param_initial=log(c((exp(input_training)*0.1),sigma_0_2_ini)) ##
  }else if(model_name=='varying_power_nonparametric'){
    ##this is somewhat ad-hoc which start from a FBM with some random parameter
    param_initial1=log(c(0.1,rep(.1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
    #param_initial2=log(c(0.5,rep(0.5,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
    param_initial2=log(c(0.1,rep(1,length(input_training)-1)+0.001*runif(length(input_training)-1),sigma_0_2_ini)) ##
    
    param_initial=cbind(param_initial1,param_initial2)
  }else if(model_name=='varying_power_msd_nonparametric'){
    param_initial=log(c(0.1,rep(1,length(input_training)-1),sigma_0_2_ini)) ##ad hoc
    
  }
  return.list=list()
  return.list$input_training=input_training
  return.list$param_initial=param_initial
  return.list$input_index=input_index
  
  return(return.list)
}


### 
# Get_MSD_grad_nonparameteric<-function(theta,d_input,model_name,input_training){
# 
#   #if(model_name=='varying_power_nonparametric'){
# 
#     range_par=(max(input_training)-min(input_training)) ##input training is in the log space
#     
#     alpha_training=2*theta[-1]/(1+theta[-1])
#     m_varying_power=rgasp(design=input_training[-1],response=(alpha_training),range.par = range_par) ##input_training is log
#     m_varying_power_pred_alpha=predict(m_varying_power,as.matrix(log(d_input[-1])))$mean
#     
#     
#     MSD_grad=matrix(NA,length(d_input),length(theta))
#     MSD_grad[,1]=c(0,exp(m_varying_power_pred_alpha*log(d_input[-1])))
#     beta=theta[1]
#     for(i in 2:length(theta)){
#       alpha_training_deriv=rep(0,length(alpha_training))
#       alpha_training_deriv[i-1]=1
#       m_varying_power_deriv=rgasp(design=input_training[-1],response=alpha_training_deriv,range.par = range_par) ##input_training is log
#       
#       m_varying_power_pred_alpha_deriv_i=predict(m_varying_power_deriv,as.matrix(log(d_input[-1])))$mean
#       
#       MSD_grad[,i]= beta*c(0,log(d_input[-1])*(d_input[-1]^alpha_training[i-1])*m_varying_power_pred_alpha_deriv_i)
#     }
#   #}###need to add grad of CFBM
#   
#   return(MSD_grad)
# }

# Get_grad_trans_nonparametric<-function(theta,d_input,model_name){
#  #if(model_name=='varying_power_nonparametric'){
#     beta=theta[1]
#     alpha=2*theta[-1]/(1+theta[-1])
#     grad_trans= c(beta, alpha*(1-alpha/2))
#     #grad_trans = c(beta, (2-alpha)^2/2)
#   #}
#   return(grad_trans)
# }



# param_uncertainty_nonparametric<-function(param_est, optimization_method, I_q_cur,B_cur=NA,num_q_cur,I_o_q_2_ori,
#                                           q_ori_ring_loc_unique_index,
#                                           sz,len_t,d_input,q, q_index_selected_here,
#                                           model_name,input_training,estimation_method='asymptotics',M,num_iteration_max,lmm,lower_bound=NA){
#   #q_lower=q-min(q)/2
#   q_lower=q-min(q)
#   
#   #param_ini=param_est
#   if (optimization_method == "optim") {
#   m_param_lower = try(optim(param_est,log_lik_param_nonparametric_ori,I_q_cur=data_fft$I_q_matrix,B_cur=NA,
#                             num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
#                             #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
#                             q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
#                             method='L-BFGS-B',
#                             control = list(fnscale=-1,maxit=num_iteration_max),sz=data_fft$sz,len_t=data_fft$len_t,q_index_selected=q_index_selected_here,
#                             d_input=d_input,q=q_lower,model_name=model_name,input_training=input_training),TRUE)
#   } else if (optimization_method == "BFGS") {
#     step_sizes = 10^c(-9, -6, -3, -1, -0.3, -0.05, -0.001, 0) 
#     
#     m_param_lower = BFGS(f = log_lik_param_nonparametric_ori, x0 = param_est, I_q_cur=data_fft$I_q_matrix,B_cur=NA,
#                          num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
#                          q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
#                          sz=data_fft$sz,len_t=data_fft$len_t,
#                          d_input=data_fft$d_input,q=q_lower,model_name=model_name,
#                          input_training=input_training,q_index_selected=q_index_selected_here,
#                          max_iter = num_iteration_max, step_sizes = step_sizes, lag_step = 20, lmm = 20)
#     
#   }
#   
#   
#   # if(class(m_param_lower)[1]=="try-error"){
#   #   compute_twice=T
#   #   m_param_lower = try(optim(param_est+c(rep(0.5,p),0),log_lik_param,I_q_cur=I_q_cur,B_cur=NA,
#   #                             num_q_cur=num_q_cur,I_o_q_2_ori=I_o_q_2_ori,
#   #                             #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
#   #                             q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
#   #                             method='L-BFGS-B',lower=lower_bound,
#   #                             control = list(fnscale=-1,maxit=num_iteration_max),sz=sz,len_t=len_t,
#   #                             d_input=d_input,q=q_lower,model_name=model_name),TRUE)
#   #   
#   # }
#   
#   #q_upper=q+min(q)/2
#   q_upper=q+min(q)
#   
#   if (optimization_method == "optim") {
#   m_param_upper = try(optim(param_est,log_lik_param_nonparametric_ori,I_q_cur=data_fft$I_q_matrix,B_cur=NA,
#                             num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
#                             #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
#                             q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
#                             method='L-BFGS-B',
#                             control = list(fnscale=-1,maxit=num_iteration_max),sz=data_fft$sz,len_t=data_fft$len_t,q_index_selected=q_index_selected_here,
#                             d_input=d_input,q=q_upper,model_name=model_name,input_training=input_training),TRUE)
#   } else if(optimization_method == "BFGS"){
#     step_sizes = 10^c(-9, -6, -3, -1, -0.3, -0.05, -0.001, 0) 
#     
#     m_param_upper = BFGS(f = log_lik_param_nonparametric_ori, x0 = param_est, I_q_cur=data_fft$I_q_matrix,B_cur=NA,
#                          num_q_cur=ini_est_list$num_q_max,I_o_q_2_ori=ini_est_list$I_o_q_2_ori,
#                          q_ori_ring_loc_unique_index=ini_est_list$q_ori_ring_loc_unique_index,
#                          sz=data_fft$sz,len_t=data_fft$len_t,
#                          d_input=data_fft$d_input,q=q_upper,model_name=model_name,
#                          input_training=input_training,q_index_selected=q_index_selected_here,
#                          max_iter = num_iteration_max, step_sizes = step_sizes, lag_step = 20, lmm = 20)
#     
#   }
#   
#   # if(class(m_param_upper)[1]=="try-error"){
#   #   compute_twice=T
#   #   m_param_upper = try(optim(param_est+c(rep(0.5,p),0),log_lik_param,I_q_cur=I_q_cur,B_cur=NA,
#   #                             num_q_cur=num_q_cur,I_o_q_2_ori=I_o_q_2_ori,
#   #                             #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
#   #                             q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
#   #                             method='L-BFGS-B',lower=lower_bound,
#   #                             control = list(fnscale=-1,maxit=num_iteration_max),sz=sz,len_t=len_t,
#   #                             d_input=d_input,q=q_upper,model_name=model_name),TRUE)
#   #   
#   # }
#   
#   p=length(param_est)-1
#   param_range=matrix(NA,2,p+1)
#   for(i in 1:(p+1) ){
#     if(optimization_method == "optim"){
#     param_range[1,i]=min(m_param_lower$par[i],m_param_upper$par[i])
#     param_range[2,i]=max(m_param_lower$par[i],m_param_upper$par[i])
#     } else if(optimization_method == "BFGS"){
#       param_range[1,i]=min(m_param_lower$x[i],m_param_upper$x[i])
#       param_range[2,i]=max(m_param_lower$x[i],m_param_upper$x[i])
#     }
#   }
#   half_length_param_range_fft=(param_range[2,]-param_range[1,])/2
#   
#   
#   if(estimation_method=='asymptotics'){
#     theta=exp(param_est)
#     p=length(param_est)-1
#     theta=exp(param_est[-(p+1)]) ##first p parameters are parameters in ISF
#     if(is.na(B_cur)){ ##this fix the dimension
#       sigma_2_0_hat=exp(param_est[p+1]) ##noise
#       B_cur=2*sigma_2_0_hat
#     }
#     #A_cur = 2*(I_o_q_2_ori - B_cur/2)
#     
#     A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
#     eta=B_cur/4 ##nugget
#     
#     MSD=Get_MSD_nonparametric(theta,d_input=d_input,model_name=model_name,input_training=input_training)
#     
#     MSD_grad=Get_MSD_grad_nonparameteric(theta,d_input,model_name,input_training=input_training)
#     grad_trans=Get_grad_trans_nonparametric(theta,d_input,model_name)
#     #Hessian_inv_list=as.list(1:num_q_cur)
#     #Hessian_inv_sum=0
#     Hessian_list=as.list(1:num_q_cur)
#     Hessian_sum=0
#     
#     for(i_q_selected in 1: num_q_cur){
#       #print(q_selected)
#       
#       #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
#       #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
#       
#       
#       #q_selected=q[i_q_selected]
#       q_selected=q[i_q_selected]
#       
#       sigma_2=A_cur[i_q_selected]/4
#       
#       acf0 = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
#       acf=acf0
#       acf[1] = acf[1]+eta ##for grad this is probably no adding
#       acf=as.numeric(acf)
#       
#       Tz <- Toeplitz$new(len_t,acf=acf) ##maybe create it outside?
#       Hessian=matrix(NA,p+1,p+1) ##last one is
#       acf_grad=matrix(NA,len_t,p+1)
#       for(i_p in 1:p){
#         acf_grad[,i_p]=-acf0*q_selected^2/4* MSD_grad[,i_p]*grad_trans[i_p]
#       }
#       #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
#       acf_grad[,p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
#       acf_grad[1,p+1]= acf_grad[1,p+1]+0.5
#       acf_grad[,p+1]= acf_grad[,p+1]*sigma_2_0_hat
#       
#       for(i_p in 1:(p+1) ){
#         for(j_p in 1:(p+1) ){
#           Hessian[i_p,j_p]=Tz$trace_hess(as.numeric(acf_grad[,i_p]), as.numeric(acf_grad[,j_p]) )
#         }
#       }
#       
#       #Hessian_list[[i_q_selected]]=Hessian
#       Hessian_sum=Hessian_sum+Hessian*length(q_ori_ring_loc_unique_index[[i_q_selected]])
#       ###a litte more conservation is to say they are perfectly correlated in a ring
#       #Hessian_sum=Hessian_sum+Hessian
#       
#       
#       #Hessian_inv_list[[i_q_selected]]=solve(Hessian)
#       #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]*length(q_ori_ring_loc_unique_index[[i_q_selected]])
#       ##a litte more conservation is to say they are perfectly correlated in a ring
#       #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]
#       
#       
#       #cor(Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][2],]),Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][1],]))
#     }
#     
#     Hessian_sum=Hessian_sum*M/sum(lengths(q_ori_ring_loc_unique_index[1:num_q_cur]))
#     sd_theta_B=sqrt(diag(solve(Hessian_sum) ))
#     #sd_theta_B=sqrt(diag(Hessian_inv_sum/sum(lengths(q_ori_ring_loc_unique_index))^2 ))
#     
#     param_range[1,]=param_range[1,]-sd_theta_B*qnorm(0.975) ##this is 10 times larger to account for not estimating A and B correct and other misspecification
#     param_range[2,]=param_range[2,]+sd_theta_B*qnorm(0.975)
#     half_length_param_range_est=sd_theta_B*qnorm(0.975)
#     
#   }
#   return(param_range)
# }


param_uncertainty_nonparametric<-function(param_est,I_q_cur,B_cur=NA,num_q_cur,I_o_q_2_ori,
                                          q_ori_ring_loc_unique_index,
                                          sz,len_t,d_input,q,
                                          q_index_selected, 
                                          model_name,input_training,estimation_method='asymptotics',M,num_iteration_max,lower_bound=NA){
  
  #no need
  #t_min = min(d_input[d_input > 0])
  #d_input_rescale = d_input / t_min           
  #input_training_rescale = input_training - log(t_min) 
  
  q_lower=q-min(q)
  
  #param_ini=param_est
  # m_param_lower_stage_1 = try(optim(param_est,log_lik_param_nonparametric_ori,I_q_cur=data_fft$I_q_matrix,B_cur=NA,
  #                           num_q_cur=num_q_cur,I_o_q_2_ori=I_o_q_2_ori,
  #                           #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
  #                           q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
  #                           method='L-BFGS-B',
  #                           control = list(fnscale=-1,maxit=100),sz=sz,len_t=len_t, 
  #                           d_input=d_input,q=q_lower,model_name=model_name,input_training=input_training,
  #                           q_index_selected = q_index_selected, subsample_t=subsample_t
  #                           ),TRUE)
  
  m_param_lower = optim(param_est, log_lik_param_nonparametric_ori, I_q_cur=I_q_cur,B_cur=NA,
                                               num_q_cur=num_q_max,I_o_q_2_ori=I_o_q_2_ori,
                                               q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                                               sz=sz,len_t=len_t,
                                               d_input=d_input,q=q_lower,model_name="direct_nonparametric",
                                               input_training=input_training,
                                               q_index_selected=q_index_selected_here,
                                               control = list(fnscale=-1,maxit=100), method = 'L-BFGS-B')  
  
  
  q_upper=q+min(q)
  
  # m_param_upper_stage_1 = try(optim(param_est,log_lik_param_nonparametric_ori,I_q_cur=I_q_cur,B_cur=NA,
  #                           num_q_cur=num_q_cur,I_o_q_2_ori=I_o_q_2_ori,
  #                           #q_ori_ring_loc_index=index_list$q_ori_ring_loc_index,
  #                           q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
  #                           method='L-BFGS-B',
  #                           control = list(fnscale=-1,maxit=100),sz=sz,len_t=len_t,
  #                           d_input=d_input,q=q_upper,model_name=model_name,input_training=input_training,
  #                           q_index_selected = q_index_selected, subsample_t=subsample_t
  #                           ),TRUE)
  
  m_param_upper = optim(param_est, log_lik_param_nonparametric_ori, I_q_cur=I_q_cur,B_cur=NA,
                        num_q_cur=num_q_max,I_o_q_2_ori=I_o_q_2_ori,
                        q_ori_ring_loc_unique_index=q_ori_ring_loc_unique_index,
                        sz=sz,len_t=len_t,
                        d_input=d_input,q=q_upper,model_name="direct_nonparametric",
                        input_training=input_training,
                        q_index_selected=q_index_selected_here,
                        control = list(fnscale=-1,maxit=100), method = 'L-BFGS-B')  
  
  ##here this is log theta, 
  p=length(param_est)-1
  param_range=matrix(NA,2,p+1)
  for(i in 1:(p+1) ){
    param_range[1,i]=min(m_param_lower$par[i],m_param_upper$par[i])
    param_range[2,i]=max(m_param_lower$par[i],m_param_upper$par[i])
    
  }
  half_length_param_range_fft=(param_range[2,]-param_range[1,])/2
  
  ##this asymptotics is on log theta as well
  if(estimation_method=='asymptotics'){
    
    p=length(param_est)-1
    theta=exp(param_est[-(p+1)]) ##first p parameters are parameters in ISF
    if(is.na(B_cur)){ ##this fix the dimension
      sigma_2_0_hat=exp(param_est[p+1]) ##noise
      B_cur=2*sigma_2_0_hat
    }
    #A_cur = 2*(I_o_q_2_ori - B_cur/2)
    
    A_cur = abs(2*(I_o_q_2_ori - B_cur/2))
    eta=B_cur/4 ##nugget
    
    #why rescaled? only MSD_grad matters, MSD and grad_trans do not matter
    #MSD=Get_MSD_nonparametric(theta,d_input=d_input_rescale,model_name=model_name,input_training=input_training_rescale)
    #MSD_grad_save=Get_MSD_grad_nonparameteric(theta,d_input_rescale,model_name,input_training=input_training_rescale)
    #grad_trans_save=Get_grad_trans_nonparametric(theta,d_input_rescale,model_name)
    grad_trans=rep(1,p) ##derivative on the param not theta
    MSD=Get_MSD_nonparametric(theta,d_input=d_input,model_name=model_name,input_training=input_training)
    ##numerical gradient
    delta=0.001 ## delta for numerical difference
    MSD_grad=matrix(NA,length(MSD),p)
    for(i_p in 1:p){
      param_change=param_est
      param_change[i_p]=  param_change[i_p]+delta
      theta_change=exp(param_change[-(p+1)]) ##first p parameters are parameters in ISF
      # if((B_cur==2*sigma_2_0_hat)){ ##this fix the dimension
      #   sigma_2_0_change=exp(param_change[p+1]) ##noise
      #   B_change=2*sigma_2_0_change
      # }
      MSD_change=Get_MSD_nonparametric(theta_change,d_input=d_input,model_name=model_name,input_training=input_training)
      MSD_grad[,i_p]=(MSD_change-MSD)/delta
  
    }
    #MSD=Get_MSD_nonparametric(theta,d_input=d_input,model_name=model_name,input_training=input_training)
    #MSD_grad=Get_MSD_grad_nonparameteric(theta,d_input,model_name,input_training=input_training)
    #grad_trans=Get_grad_trans_nonparametric(theta,d_input,model_name)
    
    #Hessian_inv_list=as.list(1:num_q_cur)
    #Hessian_inv_sum=0
    Hessian_list=as.list(1:num_q_cur)
    Hessian_sum=0
    
    for(i_q_selected in 1: num_q_cur){
      #print(q_selected)
      
      #output_re=Re(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
      #output_im=Im(I_q_cur[q_ori_ring_loc_unique_index[[i_q_selected]],])/(sz)
      
      q_selected=q[i_q_selected]
      
      sigma_2=A_cur[i_q_selected]/4
      
      acf0 = sigma_2*exp(-q_selected^2*MSD/4) ##assume 2d
      acf=acf0
      acf[1] = acf[1]+eta ##for grad this is probably no adding
      acf=as.numeric(acf)
      
      Tz <- Toeplitz$new(len_t,acf=acf) ##maybe create it outside?
      Hessian=matrix(NA,p+1,p+1) ##last one is
      acf_grad=matrix(NA,len_t,p+1)
      for(i_p in 1:p){
        acf_grad[,i_p]=-acf0*q_selected^2/4* MSD_grad[,i_p]*grad_trans[i_p]
      }
      #acf_grad[,p+1]=(-acf0*0.5/sigma_2)
      acf_grad[,p+1]=(-acf0*0.5/sigma_2)*sign(I_o_q_2_ori[i_q_selected] - B_cur/2)
      acf_grad[1,p+1]= acf_grad[1,p+1]+0.5
      acf_grad[,p+1]= acf_grad[,p+1]*sigma_2_0_hat
      
      for(i_p in 1:(p+1) ){
        for(j_p in 1:(p+1) ){
          Hessian[i_p,j_p]=Tz$trace_hess(as.numeric(acf_grad[,i_p]), as.numeric(acf_grad[,j_p]) )
        }
      }
      
      #Hessian_list[[i_q_selected]]=Hessian
      Hessian_sum=Hessian_sum+Hessian*length(q_ori_ring_loc_unique_index[[i_q_selected]])
      ###a litte more conservation is to say they are perfectly correlated in a ring
      #Hessian_sum=Hessian_sum+Hessian
      
      
      #Hessian_inv_list[[i_q_selected]]=solve(Hessian)
      #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]*length(q_ori_ring_loc_unique_index[[i_q_selected]])
      ##a litte more conservation is to say they are perfectly correlated in a ring
      #Hessian_inv_sum=Hessian_inv_sum+ Hessian_inv_list[[i_q_selected]]
      
      
      #cor(Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][2],]),Re(I_q_cur[q_ori_ring_loc_unique_index[[1]][1],]))
    }
    
    Hessian_sum=Hessian_sum*M/sum(lengths(q_ori_ring_loc_unique_index[1:num_q_cur]))
    sd_theta_B=sqrt(diag(solve(Hessian_sum)))
    #sd_theta_B=sqrt(diag(Hessian_inv_sum/sum(lengths(q_ori_ring_loc_unique_index))^2 ))
    
    param_range[1,]=param_range[1,]-sd_theta_B*qnorm(0.975) ##this is 10 times larger to account for not estimating A and B correct and other misspecification
    param_range[2,]=param_range[2,]+sd_theta_B*qnorm(0.975)
    half_length_param_range_est=sd_theta_B*qnorm(0.975)
    
  }
  
  #result = list()
  #result$param_range = param_range
  #result$sd_theta_B = sd_theta_B
  #return(result)
  return(param_range)
}

  
###yue he



# ## simulation 
bm_particle_intensity <- function(pos0,N,len,mu,sigma){
  pos=matrix(NA,N*len,2)
  pos[1:N,]=pos0
  for(i in 1:(len-1)){
    pos[i*N+(1:N),] = pos[(i-1)*N+(1:N),]+matrix(rnorm(2*N,mean=mu,sd=sigma),N,2)
  }
  return(pos)
}
bm_particle_anisotropic_intensity <- function(pos0,N,len,mu,sigma){
  pos=matrix(NA,N*len,2)
  pos[1:N,]=pos0
  for(i in 1:(len-1)){
    pos[i*N+(1:N),1] = pos[(i-1)*N+(1:N),1]+matrix(rnorm(N,mean=mu,sd=sigma[1]),N,1)
    pos[i*N+(1:N),2] = pos[(i-1)*N+(1:N),2]+matrix(rnorm(N,mean=mu,sd=sigma[2]),N,1)
  }
  return(pos)
}

ou_particle_intensity <- function(pos0,N,len,mu,sigma,rho,drift_dir){
  if(drift_dir==0){
    pos=matrix(NA,N*len,2)
    #pos[1:N,]=pos0
    pos[1:N,]=pos0+matrix(rnorm(2*N, sd=sigma),N,2)
    
    sd_innovation_OU=sqrt(sigma^2*(1-rho^2))
    #center=rep(floor(len/2),2) ##need to center ?
    for(i in 1:(len-1)){
      pos[i*N+(1:N),] = rho*(pos[(i-1)*N+(1:N),]-pos0-(i-1)*mu)+(i-1)*mu+pos0+sd_innovation_OU*matrix(rnorm(2*N),N,2)
    }
  }else if(drift_dir==1){
    theta = 2*pi*runif(N)
    pos=matrix(NA,N*len,2)
    #pos[1:N,]=pos0
    pos[1:N,]=pos0+matrix(rnorm(2*N, sd=sigma),N,2)
    
    sd_innovation_OU=sqrt(sigma^2*(1-rho^2))
    for(i in 1:(len-1)){
      pos[i*N+(1:N),1] = rho*(pos[(i-1)*N+(1:N),1]-pos0[,1]-(i-1)*mu*cos(theta))+(i-1)*mu*cos(theta)+pos0[,1]+sd_innovation_OU*rnorm(N)
      pos[i*N+(1:N),2] = rho*(pos[(i-1)*N+(1:N),1]-pos0[,1]-(i-1)*mu*sin(theta))+(i-1)*mu*sin(theta)+pos0[,2]+sd_innovation_OU*rnorm(N)
    }
  }
  return(pos)
}

fbm_particle_intensity <- function(pos0,N,len,sigma,H){
  pos=matrix(NA,N*len,2)
  pos[,1] = rep(pos0[,1],len)
  pos[,2] = rep(pos0[,2],len)
  fBM_cov = cov_fBM(len,H)
  L = t(chol(sigma^2*fBM_cov))
  increments1 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  increments2 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  pos[(N+1):(N*len),1] =pos[(N+1):(N*len),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos[(N+1):(N*len),2] =pos[(N+1):(N*len),2]+as.numeric(t(apply(increments2,2,cumsum)))
  return(pos)
}
fbm_particle_anisotropic_intensity <- function(pos0,N,len,sigma,H){
  pos=matrix(NA,N*len,2)
  pos[,1] = rep(pos0[,1],len)
  pos[,2] = rep(pos0[,2],len)
  fBM_cov1 = cov_fBM(len,H[1])
  L1 = t(chol(sigma[1]^2*fBM_cov1))
  fBM_cov2 = cov_fBM(len,H[2])
  L2 = t(chol(sigma[2]^2*fBM_cov2))
  increments1 = L1%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  increments2 = L2%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  pos[(N+1):(N*len),1] =pos[(N+1):(N*len),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos[(N+1):(N*len),2] =pos[(N+1):(N*len),2]+as.numeric(t(apply(increments2,2,cumsum)))
  return(pos)
}

fbm_ou_particle_intensity <- function(pos0,N,len,sigma_fbm=1,sigma_ou=1,H,rho){
  pos1=matrix(NA,N*len,2)
  pos2=matrix(NA,N*len,2)
  pos=matrix(NA,N*len,2)
  pos0 = pos0 + matrix(rnorm(2*N, sd=sigma_ou),N,2)
  
  pos[1:N,] = pos0
  pos1[1:N,] = pos0
  
  sd_innovation_OU=sqrt(sigma_ou^2*(1-rho^2))
  #center=rep(floor(len/2),2) ##need to center ?
  for(i in 1:(len-1)){
    pos1[i*N+(1:N),] = rho*(pos1[(i-1)*N+(1:N),]-pos0)+pos0+sd_innovation_OU*matrix(rnorm(2*N),N,2)
  }
  
  pos2[,1] = rep(pos0[,1],len)
  pos2[,2] = rep(pos0[,2],len)
  fBM_cov = cov_fBM(len,H)
  L = t(chol(sigma_fbm^2*fBM_cov))
  increments1 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  increments2 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  pos2[(N+1):(N*len),1] = pos2[(N+1):(N*len),1]+as.numeric(t(apply(increments1,2,cumsum)))
  pos2[(N+1):(N*len),2] = pos2[(N+1):(N*len),2]+as.numeric(t(apply(increments2,2,cumsum)))
  pos = pos1+pos2 - cbind(rep(pos0[,1],len), rep(pos0[,2],len))
  
  # pos=matrix(NA,N*len,2)
  # pos[1:N,] = pos0
  # 
  # fBM_cov = cov_fBM(len,H)
  # L = t(chol(sigma_fbm^2*fBM_cov))
  # increments1 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  # increments2 = L%*%matrix(rnorm((len-1)*N),nrow=len-1,ncol=N)
  # sd_innovation_OU=sqrt(sigma^2*(1-rho^2))
  # 
  # for(i in 1:(len-1)){
  #   pos[i*N+(1:N),1] = rho*(pos[(i-1)*N+(1:N),1]-pos0[,1])+pos0[,1]+sd_innovation_OU*rnorm(N)+increments1[i,]
  #   pos[i*N+(1:N),2] = rho*(pos[(i-1)*N+(1:N),2]-pos0[,2])+pos0[,2]+sd_innovation_OU*rnorm(N)+increments2[i,]
  # }
  # 
  return(pos)
}


#bm_ou_particle_intensity <- function(pos0,N,len,sigma){
#
#}

simulation <- function(frame_size=480,len=500,noise="gaussian",I0=20, Imax=255,M=50,SP="BM", pos0=matrix(nrow=N,ncol=2),mu=0,rho=0.95, drift_dir=0,H=0.3,
                       sigma_bm=2, sigma_p=2, sigma_ou=2, sigma_fbm=2){
  if(sum(is.na(pos0))>=1){
    #pos0 = matrix(frame_size*runif(N*2),nrow=N,ncol=2)  ## does it make it move out too easily?
    #pos0 = matrix(0.25*frame_size+0.5*frame_size*runif(N*2),nrow=N,ncol=2)
    pos0 = matrix(frame_size/8+0.75*frame_size*runif(M*2),nrow=M,ncol=2)
    #lhs_sample=as.vector(maximinLHS(n=N,k=2))
    #pos0 = matrix(0.25*frame_size+0.5*frame_size*lhs_sample,nrow=N,ncol=2)
    
  }
  if(SP == "BM"){
    #pos = bm(pos0,len=len,sigma=sigma, mu=mu)
    ##let's try R, this is not slow
    pos = bm_particle_intensity(pos0=pos0,N=M,len=len,mu=mu,sigma=sigma_bm)
  }else if(SP == "OU"){
    #pos = ou_same_dir(pos0,len=len,sigma=sigma,rho=rho)
    pos = ou_particle_intensity(pos0=pos0,N=M,len=len,mu=mu,sigma=sigma_ou,rho=rho,drift_dir=drift_dir)
  }else if(SP == "FBM"){
    pos = fbm_particle_intensity(pos0=pos0,N=M,len=len,sigma=sigma_fbm,H=H)
  }else if(SP == "OU+FBM"){
    pos = fbm_ou_particle_intensity(pos0=pos0,N=M,len=len,H=H, rho=rho, sigma_ou = sigma_ou, sigma_fbm = sigma_fbm)
  }else if(SP=='BM_anisotropic'){
    pos = bm_particle_anisotropic_intensity(pos0=pos0,N=M,len=len,mu=mu,sigma=sigma_bm)
  }else if(SP=='FBM_anisotropic'){
    pos = fbm_particle_anisotropic_intensity(pos0=pos0,N=M,len=len,sigma=sigma_fbm,H=H)
  }
  
  
  if(length(I0) == len){
    if(noise == "uniform"){
      I = matrix(runif(frame_size*frame_size*len)-0.5, nrow=len,ncol = frame_size*frame_size)
      I = I*I0
    }else if(noise == "gaussian"){
      I = matrix(rnorm(frame_size*frame_size*len), nrow=len,ncol = frame_size*frame_size)
      I = I*sqrt(I0)
    }
  }else if(length(I0) == 1){
    if(noise == "uniform"){
      I = matrix(I0*(runif(frame_size*frame_size*len)-0.5), nrow=len,ncol = frame_size*frame_size)
    } else if(noise == "gaussian"){
      I = matrix(sqrt(I0)*rnorm(frame_size*frame_size*len), nrow=len,ncol = frame_size*frame_size)
    }
  }
  
  sim_list=list()
  
  Ic = rep(Imax,M)
  if(length(Imax)==1){
    sim_list$intensity = fill_intensity(len=len,N=M,I=I,pos=pos,Ic=Ic,frame_size=frame_size, sigma_p=sigma_p)
  }
  sim_list$pos=pos
  colnames(sim_list$pos)=c('x','y')
  sim_list$M=M
  sim_list$frame_size=frame_size
  sim_list$len=len
  sim_list$SP=SP
  
  return(sim_list)
}





particle_trajectory <- function(pos, N, len, frame_size,SP,title=NA,fame_name = NA, save = 0){
  traj1 = pos[seq(1,dim(pos)[1],by=N),]
  if(is.na(title)==T){title=SP}
  if(save ==1){
    pdf(file=file_name,height=4.5,width=5.5)
    plot(traj1[,1],traj1[,2],ylim=c(0,frame_size),xlim=c(0,frame_size),
         type="l",col=1, xlab = "frame size", ylab="frame size", main=title)
    for (i in 2:N){
      v = pos[seq(i,dim(pos)[1],by=N),]
      lines(v[,1],v[,2],ylim=c(0,frame_size),xlim=c(0,frame_size),type="l", col=i)
    }
    dev.off()
  }else{
    plot(traj1[,1],traj1[,2],ylim=c(0,frame_size),xlim=c(0,frame_size),
         type="l",col=1, xlab = "frame size", ylab="frame size", main=title)
    for (i in 2:N){
      v = pos[seq(i,dim(pos)[1],by=N),]
      lines(v[,1],v[,2],ylim=c(0,frame_size),xlim=c(0,frame_size),type="l", col=i)
    }
  }
}

numerical_msd <- function(pos, N,len){
  pos_msd = array(pos, dim=c(N, len,2))
  msd_i = matrix(nrow=N,ncol=len-1)
  for(dt in 1:(len-1)){
    ndt = len-dt
    xdiff = pos_msd[,1:ndt,1]-pos_msd[,(1+dt):(ndt+dt),1]
    ydiff = pos_msd[,1:ndt,2]-pos_msd[,(1+dt):(ndt+dt),2]
    mean_square =xdiff^2+ydiff^2
    if (length(dim(mean_square))>1){
      msd_i[,dt] = apply(mean_square,1,function(x){mean(x,na.rm=T)}) 
    }else
      msd_i[,dt] = mean_square 
  }
  result_list=list()
  num_msd_mean = apply(msd_i,2,function(x){mean(x,na.rm=T)})
  result_list$num_msd_mean = num_msd_mean
  result_list$num_msd = msd_i
  return(result_list)
}


# Aug 23, 2025
get_Gp_Gpp_GSER <- function(MSD_est, d_input, input_training, simul_model_name, sigma_bm, sigma_fbm, H, sigma_ou, rho,
                            kB=1, Temp=1, a=1) {
  # d_input starts from 0
  if(simul_model_name == 'BM'){
    beta = 2*sigma_bm^2
    
    MSD_truth = beta*d_input
    MSD_ln_dev_truth = 1  # d lnMSD/d ln d_input
  } else if(simul_model_name == 'FBM'){
    beta = 2*sigma_fbm^2
    alpha = 2*H
    
    MSD_truth = beta*d_input^alpha
    MSD_ln_dev_truth = alpha
  } else if(simul_model_name == 'OU'){
    amplitude = 4*sigma_ou^2
    
    MSD_truth = (amplitude*(1-rho^d_input))
    #dMSD_d_input = amplitude*(-log(rho))*(rho^d_input) 
    MSD_ln_dev_truth = d_input*(-log(rho))*(rho^d_input)/(1-rho^d_input)
    MSD_ln_dev_truth = MSD_ln_dev_truth[-1]
  } else if(simul_model_name == 'OU+FBM'){
    amplitude = 4*sigma_ou^2
    beta = 2*sigma_fbm^2
    alpha = 2*H 
    
    MSD_truth = beta*d_input^alpha + (amplitude*(1-rho^d_input)) 
    dMSD_d_input = beta*alpha*(d_input^(alpha-1)) + amplitude*(-log(rho))*(rho^d_input) 
    MSD_ln_dev_truth = (d_input[-1]/MSD_truth[-1])*dMSD_d_input[-1]
  }

  Gcomplex_truth = 2*kB*Temp/(3*pi*a*MSD_truth[-1]*gamma(1+MSD_ln_dev_truth))
  Gprime_truth = Gcomplex_truth*cos(pi*MSD_ln_dev_truth/2)
  Gdoubleprime_truth = Gcomplex_truth*sin(pi*MSD_ln_dev_truth/2)

  # estimated Gp and Gpp
  # grad = numeric(length(d_input))
  # for (i in 1:length(d_input)) {
  #   step = numeric(length(d_input))
  #   step[i] = 1e-3  # Small step in the ith direction
  #   MSD_est_step = Get_MSD_nonparametric(theta=theta_est,d_input=d_input+step,input_training=input_training,model_name="direct_nonparametric")
  #   grad[i] = (MSD_est_step[i] - MSD_est[i]) / step[i]
  # }
  # 
  # MSD_ln_dev_est = (d_input[-1]/(MSD_est[-1])) * grad[-1]
  
  n = length(d_input)-1
  MSD_ln_dev_est = numeric(n)
  ln_MSD_est = log(MSD_est[-1])
  ln_time = log(d_input[-1])
  if (n >= 3) {
    MSD_ln_dev_est[2:(n-1)] = (ln_MSD_est[3:n] - ln_MSD_est[1:(n-2)]) / (ln_time[3:n] - ln_time[1:(n-2)])
  }
  # ends
  MSD_ln_dev_est[1] = (ln_MSD_est[2] - ln_MSD_est[1]) / (ln_time[2] - ln_time[1])
  MSD_ln_dev_est[n] = (ln_MSD_est[n] - ln_MSD_est[n-1]) / (ln_time[n] - ln_time[n-1])
  
  Gcomplex_est= 2*kB*Temp/(3*pi*a*MSD_est[-1]*gamma(1+MSD_ln_dev_est))
  Gprime_est = Gcomplex_est*cos(pi*MSD_ln_dev_est/2)
  Gdoubleprime_est = Gcomplex_est*sin(pi*MSD_ln_dev_est/2)
  omega = 2*pi/d_input[-1]
  
  if(simul_model_name == 'BM') {
    Gprime_truth = rep(0, n)
    Gprime_est = rep(0, n)
  }


  return_list=list()
  return_list$omega=omega
  return_list$Gprime_truth=Gprime_truth
  return_list$Gdoubleprime_truth=Gdoubleprime_truth
  return_list$Gprime_est=Gprime_est
  return_list$Gdoubleprime_est=Gdoubleprime_est
  return_list
}



# for experimental data 
get_Gp_Gpp_GSER_real_data <- function(MSD_est, d_input, kB=1, Temp=1, a=1, polynomial = F) {
   
   #grad = numeric(length(d_input))
   #for (i in 1:length(d_input)) {
   #  step = numeric(length(d_input))
   #  step[i] = 1e-3  # Small step in the ith direction
   #  MSD_est_step = Get_MSD_nonparametric(theta=theta_est,d_input=d_input+step,input_training=input_training,model_name="direct_nonparametric")
   #  grad[i] = (MSD_est_step[i] - MSD_est[i]) / 1e-3
   #}
   #
   #MSD_ln_dev_est = (d_input[-1]/(MSD_est[-1])) * grad[-1]
  
  idx = which(MSD_est > 0)
  #n = length(d_input)-1
  n = length(idx)
  MSD_ln_dev_est = numeric(n)
  ln_MSD_est = log(MSD_est[idx])
  ln_time = log(d_input[idx])
  
  if(polynomial == T) {
    ##4th order
    lm_ln_MSD_est=lm(ln_MSD_est ~ poly(ln_time, 4, raw=TRUE))
    MSD_ln_dev_est=4*lm_ln_MSD_est$coefficients[5]*ln_time^3+3*lm_ln_MSD_est$coefficients[4]*ln_time^2+2*lm_ln_MSD_est$coefficients[3]*ln_time+lm_ln_MSD_est$coefficients[2]
    
  } else {
  
    if (n >= 3) {
      MSD_ln_dev_est[2:(n-1)] = (ln_MSD_est[3:n] - ln_MSD_est[1:(n-2)]) / (ln_time[3:n] - ln_time[1:(n-2)])
    }
    # ends
    MSD_ln_dev_est[1] = (ln_MSD_est[2] - ln_MSD_est[1]) / (ln_time[2] - ln_time[1])
    MSD_ln_dev_est[n] = (ln_MSD_est[n] - ln_MSD_est[n-1]) / (ln_time[n] - ln_time[n-1])
  }
  
  Gcomplex_est= 2*kB*Temp/(3*pi*a*MSD_est[idx]*gamma(1+MSD_ln_dev_est))
  Gprime_est = Gcomplex_est*cos(pi*MSD_ln_dev_est/2)  
  Gdoubleprime_est = Gcomplex_est*sin(pi*MSD_ln_dev_est/2) 
  omega = 2*pi/d_input[idx]
  
  return_list=list()
  return_list$omega=omega 
  return_list$Gprime_est=Gprime_est
  return_list$Gdoubleprime_est=Gdoubleprime_est
  return_list
}


