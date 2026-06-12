
###this allows fixed variance
###if fixed_noise_variance=T, then the output_weights is the inverse variance 
rcalibration_no_discrepancy_local<-function(design, observations, 
             p_theta,#simul_type=0,
             simul_nug, ##with a nugget in the emulator 
             input_simul,output_simul,
             theta_range,output_weights=NULL,
             loc_index_emulator,
             S,S_0,#discrepancy_type='no-discrepancy',
             sd_proposal, fixed_noise_variance,
             initial_values=NULL,
             thinning=1){ ##if for fixed noise, then the output_weights are inverse variance
  
  ##step 1 build emulator
  if(dim(output_simul)[2]==1){ ##vector output
    emulator=rgasp(design=input_simul, response=output_simul,nugget.est=simul_nug)
    emulator_type='rgasp'
    #model@emulator_rgasp=emulator
    
  }else{ ##matrix output, ppgasp
    emulator=ppgasp(design=input_simul, response=output_simul,nugget.est=simul_nug)
    emulator_type='ppgasp'
    #model@emulator_ppgasp=emulator
    
    #if(!is.null(loc_index_emulator)){
      #model@loc_index_emulator=loc_index_emulator ##only useful for ppgasp emulator
    #}else{
     # model@loc_index_emulator=1:dim(output_simul)[2]
    #}
  }
  
  ##step 2, initialize 
  p_x=dim(design)[2]
  num_obs=length(observations)
  
  ##X is the additional mean basis
  X=matrix(0,dim(as.matrix(design))[1],1)
  
  if(is.null(initial_values)){
    par_cur=rowMeans(theta_range)  ###current par
  }else{
    par_cur=initial_values
  }
  
  ##if there is variance, then we have a variance par
  #if(fixed_noise_variance==F){
    par_cur= c(par_cur,1) #var par always 1 if fixed variance
  #}
  

  #post_sample(design,observations,emulator)
  if(!is.null(output_weights)){
     inv_output_weights=1/output_weights;
  }
  if(fixed_noise_variance==T){
    record_par=matrix(0,ceiling(S/thinning),p_theta+1);
    record_par[,p_theta+1]=rep(1,ceiling(S/thinning))
  }else{
    record_par=matrix(0,S,p_theta+1);
  }
  record_post=rep(0,ceiling(S/thinning))
  record_S=length(1:ceiling(S/thinning))
  
  accept_S_theta=0;
  accept_S_dec=0; # this is to count how many sample points are outside the boundary and get rejected
  
  c_prop=1/4
  
  sd_theta=sd_proposal[1:p_theta];
  
  
  post_cur=0;
  
  theta_cur=par_cur[1:p_theta];
  theta_sample=par_cur[1:p_theta];
  r_ratio=0;
  
  
  
  ##current variance of computer model
  cm_obs_cur=mathematical_model_eval(input=design,par_cur[1:p_theta],simul_type=0,emulator=emulator,emulator_type=emulator_type,
                                     loc_index_emulator=loc_index_emulator,math_model=NULL);
  
  post_cur=Log_marginal_post_no_discrepancy(par_cur,observations,p_theta,X,have_mean=F, inv_output_weights,cm_obs_cur,S_2_f=0,num_obs_all=num_obs);
  
  ##step 3, MCMC
  for (i_S in 1:S){
    if(i_S==floor(S*c_prop)){
      #cat(post_cur,'\n')
      cat(c_prop*S, ' of ', S, ' posterior samples are drawn \n')
      c_prop=c_prop+1/4
    }
    
    ##sample a noisevariance parameter then sample it, this is from RobustCalibration package 
    if(fixed_noise_variance==F){
      par_cur[(p_theta+1)]=Sample_sigma_2_theta_m_no_discrepancy(par_cur,observations,p_theta,X,have_trend=F, 
                                                                 inv_output_weights,cm_obs_cur,S_2_f,num_obs_all=num_obs)
    }
    
    #sample theta
    theta_cur=par_cur[1:p_theta];
    decision_0=F;
    
    for(i_theta in 1:p_theta){
      theta_sample[i_theta]=rnorm(1,mean=theta_cur[i_theta],sd=sd_theta[i_theta]*(theta_range[i_theta,2]- theta_range[i_theta,1]) )  ##here maybe truncated
      if((theta_sample[i_theta]>theta_range[i_theta,2])| (theta_sample[i_theta]<theta_range[i_theta,1]) ){
        decision_0=T  ##reject directly
        accept_S_dec=accept_S_dec+1
        break;
      }
    }
    
    
    if(decision_0==F){
      #theta_sample= c(0.72, 0.286, 4.5)
      
      param_propose=par_cur;
      param_propose[1:p_theta]=theta_sample;
      
      cm_obs_propose=mathematical_model_eval(design,theta_sample,simul_type=0,emulator,emulator_type,loc_index_emulator,math_model=NULL);
      
      post_propose=Log_marginal_post_no_discrepancy(param_propose,observations,p_theta,X,have_mean=F, inv_output_weights,cm_obs_propose,S_2_f=0,num_obs_all=num_obs);
      
      
      #plot(as.numeric(cm_obs_propose_save),col='blue',type='l')
      #lines(as.numeric(cm_obs_propose),col='red',type='l')
      #lines(observations,type='l')
      
      
      #  Log_marginal_post(param_propose,L,output,p_theta,p_x,X,have_trend,CL,a,b, cm_obs_propose);
      r_ratio=exp(post_propose-post_cur);
      decision=Accept_proposal(r_ratio);
      if(decision){
        par_cur=param_propose;
        post_cur=post_propose;
        cm_obs_cur=cm_obs_propose;
        accept_S_theta=accept_S_theta+1;
      }
      
    }
    
    ##record the thinning sample
    if(i_S%%thinning==0){
      
      record_par[floor(i_S/thinning),]=par_cur;
      record_post[floor(i_S/thinning)]=post_cur;
    }
    
    
  }
  
  return_list=list()
  return_list$post_sample=record_par
  return_list$record_post=record_post
  return_list$acceptance_theta=accept_S_theta
  return_list$emulator=emulator
  return_list$emulator_type=emulator_type
  return_list$p_theta=p_theta
  cat(accept_S_theta, ' of ', S, ' posterior samples are accepted \n')
  
  return(return_list)
}

predict_rcalibration_no_discrepancy_local<-function(object, testing_input, X_testing=NULL,
                               n_thinning=10, testing_output_weights=NULL, 
                               interval_est=NULL,interval_data=F,
                               math_model=NULL,test_loc_index_emulator=NULL,...){
  if(is.null(X_testing)){
    X_testing=matrix(0,dim(testing_input)[1],1)
  }
  if(is.null(testing_output_weights)){
    testing_output_weights=rep(1,dim(testing_input)[1])
  }
  
  #if(object$simul_type==0){
    if(is.null(test_loc_index_emulator)){
      if(object$emulator_type=='ppgasp'){
        test_loc_index_emulator=1:dim(object$emulator@output)[2]
      }else{
        test_loc_index_emulator=c(1)
      }
    }
    emulator=object$emulator

  #}
  
  ##Pred
  record_cm_pred=0
  record_cm_pred_no_trend=0
  
  SS=floor(dim(object$post_sample)[1]/n_thinning)
  c_prop=1/4
  
  if(!is.null(interval_est)){
    record_interval=matrix(0,dim(testing_input)[1],length(interval_est));
    if(interval_data!=T){ ##interval for mean 
      record_sample_interval_mean=matrix(NA,dim(testing_input)[1],SS)
    }
  }
  count_i_S=0;
  for(i_S in (1: SS)*n_thinning ){
    count_i_S=count_i_S+1
    #print(i_S)
    
    if(i_S==floor(SS*n_thinning*c_prop)){
      cat(c_prop*100, 'percent is completed \n')
      c_prop=c_prop+1/4
    }
    
    # print(i_S)
    
    theta=object$post_sample[i_S,1:object$p_theta]
    sigma_2_0=object$post_sample[i_S,object$p_theta+1]
    
    mean_cm_test=mathematical_model_eval(testing_input,theta,simul_type=0,emulator,object$emulator_type,
                                         loc_index_emulator=test_loc_index_emulator,math_model);
    
    if(object$emulator_type=='ppgasp'){
      
      #if(dim(mean_cm_test)[1]==1){
      mean_cm_test=as.vector(mean_cm_test)
      #}
    }
    
     mean_cm_test_no_trend=mean_cm_test
    # if(object$have_trend){
    #   theta_m=object$post_sample[i_S,(object@p_theta+2):(object@p_theta+1+object@q)]
    #   mean_cm_test=mean_cm_test+X_testing%*%theta_m
    # }
    # 
    
    record_cm_pred=record_cm_pred+mean_cm_test
    
    record_cm_pred_no_trend=record_cm_pred_no_trend+mean_cm_test_no_trend
    if(!is.null(interval_est)){
      if(interval_data==T){
        qnorm_all=qnorm(interval_est);
        for(i_int in 1:length(interval_est) ){
          record_interval[,i_int]=record_interval[,i_int]+mean_cm_test+qnorm_all[i_int]*sqrt(sigma_2_0/testing_output_weights)
          
        }
        #if(interval_data==T){
        # for(i_int in 1:length(interval_est) ){
        #   record_interval[,i_int]=record_interval[,i_int]+qnorm_all[i_int]*sqrt(sigma_2_0/testing_output_weights)
        # }
        #}
      }else{##for mean 
        record_sample_interval_mean[,count_i_S]=mean_cm_test
      }
      
    }
  }
  record_cm_pred=record_cm_pred/floor(SS)
  record_cm_pred_no_trend=record_cm_pred_no_trend/floor(SS)
  #output.list <- list()
  
  #output.list$math_model_mean=record_cm_pred
  predictobj=list()
  predictobj$math_model_mean=record_cm_pred
  predictobj$math_model_mean_no_trend=record_cm_pred_no_trend
  predictobj$mean=predictobj$math_model_mean
  
  if(!is.null(interval_est)){
    if(interval_data==T){
      record_interval=record_interval/floor(SS)
      #output.list$interval=record_interval
      predictobj$interval=record_interval
      #ans.list[[2]]=record_interval
    }else{ ###interval for mean
      for(i_test in 1:dim(testing_input)[1]){
        record_interval[i_test,]=quantile(record_sample_interval_mean[i_test,],interval_est,na.rm = T)
        predictobj$interval=record_interval
        
      }
    }
  }
  #matplot(as.matrix((record_sample_interval_mean[,seq(1000,5000,5)])),type='l')
  return(predictobj)
}
  