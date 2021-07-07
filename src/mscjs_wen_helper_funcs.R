


#function that prints fixed and random effects and their standard deviations and correlations in a nice table
print_out<-function(mscjs_fit){
  
  
  #function that prints fixed and random effects and their standard deviations and correlations in a nice table for a given model (i.e. phi, p, or psi)
  print_cov<-function(param){
    mat_out<-matrix(NA,nrow=0,ncol=4+100,dimnames = list(c(),c("Groups","Name","SD","SE(SD)","Corr",rep("",99))))
    ind<-0
    for(i in 1:length(mod_rep[[paste0("sd_",tolower(param))]])){
      ind2<-ind+length(mod_rep[[paste0("sd_",tolower(param))]][[i]])
      
      ##RE variance
      value<-ad_rep_vals[names(ad_rep_vals)==paste0("exp(theta_",tolower(param),")")][(ind+1):ind2]
      
      ###SE
      sd_val<-ad_rep_sds[names(ad_rep_vals)==paste0("exp(theta_",tolower(param),")")][(ind+1):ind2]
      ##lower tri corr
      cor_mat<-mod_rep[[paste0("corr_",tolower(param))]][[i]]
      cor_mat[upper.tri(cor_mat,TRUE)]<-NA
      
      ##combine
      mat_i<-matrix(NA,nrow=length(value),ncol=4+100,dimnames = list(rep("",length(value)),c("Groups","Name","SD","SE(SD)","Corr",rep("",99))))
      mat_i[1,1]<-mscjs_fit[[paste0(param,".design.glmmTMB")]]$grpVar[i]
      mat_i[,2]<-mscjs_fit[[paste0(param,".design.glmmTMB")]]$condList$reTrms$cnms[i][[1]]
      mat_i[,3]<-round(value,4)
      mat_i[,4]<-round(sd_val,4)
      if(ncol(cor_mat)>0){ mat_i[,5:(4+ncol(cor_mat))]<-round(cor_mat,4)
      }
      mat_out<-rbind(mat_out,mat_i)
      
      ind<-ind2
      
    }
     out<-mat_out[,apply(mat_out,2,function(x)sum(!is.na(x)))>0]
    
    print(out,na.print = "",quote = FALSE)
    
    
    out[,1:2]<-t(apply(out[,2,drop=FALSE],1,function(x){
      y<-unlist(strsplit(x,":"))
      if(length(y)==1)y<-c(y,"All")
      return(y)
    }))
      
    out <- as_tibble(out) %>% 
      rename("Occasion"=1,"LH"=2) %>%  
      mutate(Occasion=substr(Occasion,nchar(Occasion),nchar(Occasion)), 
             LH = case_when(LH=="LHsummer" ~"Sum.0", 
                            LH=="LHfall" ~ "Fal.0", 
                            LH=="LHsmolt" ~ "Spr.1",
                            LH=="age_classsmolt" ~ "Age.1", 
                            LH=="age_classsub" ~ "Age.0", 
                            TRUE ~ LH))
  
    
    return(out)
  }
  
  
  
  #mscjs_fit$fit$time_for_run
  mscjs_fit$fit$time_for_MLE
  mscjs_fit$fit$time_for_sdreport
  print(mscjs_fit$fit$time_for_run)
  print(paste("AIC=",mscjs_fit$fit$AIC))
  #sd report
  #sd_rep<-sdreport(mod)
  #sd_rep
  fixed_est<-mscjs_fit$fit$SD$par.fixed
  fixed_sd<-sqrt(diag(mscjs_fit$fit$SD$cov.fixed))
  
  #report
  mod_rep<-mscjs_fit$mod$report(mscjs_fit$last_par_best)
  #ad report
  ad_rep_vals<-mscjs_fit$fit$SD$value
  ad_rep_sds<-mscjs_fit$fit$SD$sd
  
  #par_out<-mscjs_fit$mod$env$parList()
  
  print("Phi")
  print("fixed")
  tibble(par_name=colnames(mscjs_fit$Phi.design.glmmTMB$data.tmb$X),estimate=ad_rep_vals[names(ad_rep_vals)=="beta_phi"],sd=ad_rep_sds[names(ad_rep_vals)=="beta_phi"]) %>%  mutate(p_value= round(2*pnorm(abs(estimate/sd), lower.tail = FALSE),3)) %>% print(n=100)
 print("Random covariance")
  try(phi_rand_cov<- print_cov("Phi"))
  
  
  
  print("p")
  print("fixed")
  tibble(par_name=colnames(mscjs_fit$p.design.glmmTMB$data.tmb$X),estimate=ad_rep_vals[names(ad_rep_vals)=="beta_p"],sd=ad_rep_sds[names(ad_rep_vals)=="beta_p"]) %>%  mutate(p_value= round(2*pnorm(abs(estimate/sd), lower.tail = FALSE),3)) %>% print(n=100)
  print("Random covariance")
  try(p_rand_cov<-print_cov("p"))
  
  
  print("Psi")
  print("fixed")
  tibble(par_name=colnames(mscjs_fit$Psi.design.glmmTMB$data.tmb$X),estimate=ad_rep_vals[names(ad_rep_vals)=="beta_psi"],sd=ad_rep_sds[names(ad_rep_vals)=="beta_psi"]) %>%  mutate(p_value= round(2*pnorm(abs(estimate/sd), lower.tail = FALSE),3)) %>% print(n=100)
  try(psi_rand_cov<-print_cov("Psi"))
  
  return(list(phi_rand_cov,p_rand_cov,psi_rand_cov))
}




#function to make a table of coeficients for environmental covariates for plotting
env_cov_tab_func<-function(covs=c("time1:win_air:LHfall",
                                   "time1:win_air:LHsummer" , 
                                   
                                   "time1:LHfall:win_flow"   ,                    
                                   "time1:LHsummer:win_flow" ,
                                   
                                   "time4:age_classsmolt:ersstWAcoast.sum"     ,  
                                   "time4:age_classsub:ersstWAcoast.sum"   ,      
                                   "time4:age_classsmolt:cui.spr"          ,      
                                   "time4:age_classsub:cui.spr" ),fit_obj=mscjs_fit){
  
  ad_rep_vals<-fit_obj$fit$SD$value
  ad_rep_sds<-fit_obj$fit$SD$sd

mat_vars<-strsplit(covs,":") %>% unlist() %>% matrix(ncol=3,byrow = TRUE)

tab_coef<-tibble(time=substr(mat_vars[,1],5,5),
       LH=apply(mat_vars,1,function(x)x[which(substr(x,1,2)%in%c("LH","ag"))]),
       var=apply(mat_vars[,-1],1,function(x)x[which(!(substr(x,1,2)%in%c("LH","ag")))]),
       est=(ad_rep_vals[names(ad_rep_vals)=="beta_phi"][match(covs,colnames(fit_obj$Phi.design.glmmTMB$data.tmb$X))]),
       sd=(ad_rep_sds[names(ad_rep_vals)=="beta_phi"][match(covs,colnames(fit_obj$Phi.design.glmmTMB$data.tmb$X))])
) %>% mutate( LH=case_when(LH=="LHfall"~"Fal-0",
                               LH=="LHsummer"~"Sum-0",
                               LH%in%c("LHsmolt","age_classsmolt")~"Spr-1",
                               TRUE~"Age-0"),
                  LH=fct_relevel(LH,"Sum-0","Fal-0","Age-0","Spr-1"),
              var=case_when(var=="win_air"~"win.air.temp",
                            var=="win_flow"~"win.flow",
                            var=="RIS_flow_juv"~"flow",
                            var=="RIS_temp_juv"~"temp",
                            var%in%c("MC_spill_pct_juv",
                                     "UC_spill_pct_juv")~"spill",
                            var=="McN_flow"~"flow",
                            var=="McN_temp"~"temp",
                            var=="ersstWAcoast.sum"~"SST.sum",
                            var=="ersstArc.win"~"SST.win",
                            var=="cui.spr"~"CUI.spr"
                                        ))


ggplot(data=tab_coef,aes(x=var,y=est,color=LH))+facet_grid(cols=vars(time),labeller =  labeller(time=c(`1`="Release to Lower Wenatchee",`2`="Lower Wenatchee to McNary", `3` = "McNary to Bonneville ",`4` = "Marine survival")),scales = "free")+geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .75),size=5)+geom_linerange(aes(ymin=est-1.96*sd,ymax=est+1.96*sd),position=position_dodge(width = .75),lwd=1.25)+xlab("")+ylab(expression(beta))+ggthemes::scale_colour_colorblind()+labs(color = "Life\nhistory")
}

# function to plot coefficient effects of redds
plot_redd_effect_fun<-function(covs, 
      fit_obj=mscjs_fit){
  
  ad_rep_vals<-fit_obj$fit$SD$value
  ad_rep_sds<-fit_obj$fit$SD$sd
  
  mat_vars<-strsplit(covs,":") %>% unlist() %>% matrix(ncol=4,byrow = TRUE)
  
  mat_vars<-strsplit(covs,":") %>% unlist() %>% matrix(ncol=4,byrow = TRUE) %>% as_tibble() %>% `colnames<-`(c("time","LH" ,"stream","redd")) %>% 
    mutate(time=substr(time,5,5),
           LH=apply(mat_vars,1,function(x)x[which(substr(x,1,2)%in%c("LH","ag"))]),
           LH=case_when(LH=="LHfall"~"Fal-0",
                        LH=="LHsummer"~"Sum-0",
                        LH%in%c("LHsmolt","age_classsmolt")~"Spr-1",
                        TRUE~"Age-0"),
           LH=fct_relevel(LH,"Sum-0","Fal-0","Age-0","Spr-1"),
           stream=apply(mat_vars,1,function(x)x[which(substr(x,1,2)=="st")]),
           stream= substr(stream, 7,nchar(stream)),
                   est=(ad_rep_vals[names(ad_rep_vals)=="beta_phi"][match(covs,colnames(fit_obj$Phi.design.glmmTMB$data.tmb$X))]),
                   sd=(ad_rep_sds[names(ad_rep_vals)=="beta_phi"][match(covs,colnames(fit_obj$Phi.design.glmmTMB$data.tmb$X))])
  ) 
  
  
  ggplot(data=mat_vars,aes(x=stream,y=est,color=LH))+facet_grid(cols=vars(time),scales = "free",labeller =  labeller(time=c(`1`="Release to Lower Wenatchee",`2`="Lower Wenatchee to McNary", `3` = "McNary to Bonneville ",`4` = "Marine survival")))+geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width = .75),size=5)+geom_linerange(aes(ymin=est-1.96*sd,ymax=est+1.96*sd),position=position_dodge(width = .75),lwd=1.25)+xlab("")+ylab(expression(beta))+ggthemes::scale_colour_colorblind()+labs(color = "Life\nhistory")
}


#function for simulation from a multivariate normal, taken fron J. Thorson's Fish Utils. Package
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}


# function to plot downstream survival probs
plot_func<-function(){
  sd_rep<-mscjs_fit$fit$SD
  mod_rep<-mscjs_fit$mod$report()
  
  Phi.design.dat2 <- mscjs_dat$Phi.design.dat%>% mutate(
    eta_phi= sd_rep$value[names(sd_rep$value)=="eta_phi"],
    eta_phi_sd=sd_rep$sd[names(sd_rep$value)=="eta_phi"],
    phi_fit=plogis(eta_phi),
    lcl_phi=plogis(qnorm(.025,eta_phi,eta_phi_sd)),
    ucl_phi=plogis(qnorm(.975,eta_phi,eta_phi_sd))) %>% 
    mutate(LH=fct_relevel(LH,"summer","fall","smolt")) %>%
    group_by(LH,mig_year,stream,time) %>%
    mutate(cohort_freq=sum(freq)) %>% 
    mutate(across(phi_fit:ucl_phi,~sum(.*freq/cohort_freq))) %>% 
    filter(!(time=="1"&LH=="Unk"))#%>% 
  # summarize(sum(length_bin*freq/cohort_freq))
  
  # Phi.design.dat2<-Phi.design.dat2 %>% select(time ,stratum , LH , stream ,mig_year, phi_fit)
  
  p.design.dat <- mscjs_dat$p.design.dat %>% mutate(#p_fit=mod_rep$p[-length(mod_rep$p)],
    eta_p=sd_rep$value[names(sd_rep$value)=="eta_p"],
    eta_p_sd=sd_rep$sd[names(sd_rep$value)=="eta_p"],
    p_fit=plogis(eta_p),
    lcl_p=plogis(qnorm(.025,eta_p,eta_p_sd)),
    ucl_p=plogis(qnorm(.975,eta_p,eta_p_sd)))%>% 
    filter(!(time=="1"&LH=="Unk"))
  
  
  Psi.design.dat<-mscjs_dat$Psi.design.dat %>% filter(tostratum==2) %>% select(LH,mig_year,stream,LWe_J,McN_J) %>% cbind(mod_rep$psi) %>% filter(LH!="summer",stream=="Chiwawa",McN_J==0,LWe_J==0) %>% pivot_longer(`1`:`3`,names_to="years",values_to="prop") %>% mutate(LH=ifelse(LH=="fall","Subyearling","Smolt"))
  
  
  #downstream  survival
  # dev.new()
  ggplot(data= as_tibble(Phi.design.dat2) %>% filter(Time<=2),aes(x=(mig_year) ,y=phi_fit,color=stream)) + geom_point(position=position_dodge(width = .75))+geom_linerange(aes(ymin=lcl_phi,ymax=ucl_phi),position=position_dodge(width = .75))+facet_grid(rows=vars(time),cols=vars(LH),labeller =  labeller(time=c(`1`="Release to Lower Wenatchee",`2`="Lower Wenatchee to McNary", `3` = "McNary to Bonneville ")),scales = "free")+ylim(0,1)+xlab("Smolt year")+ylab("Survival")
}                

#function to plot histogranm of detection DOY
plot_doy<-function(site,mark_file_CH,main){
  dat<-mark_file_CH %>% group_by(LH) %>% rename(x=all_of(site)) %>% mutate(med.doy=median(x,na.rm=T)) %>% mutate(LH=as.factor(LH),LH=fct_relevel(LH,"summer","fall","Unk","smolt")) %>% ungroup()
  
  ggplot(dat=dat , aes(x=x)) +geom_histogram(fill="grey")+geom_vline(dat=dat , aes(xintercept=med.doy))+
    facet_wrap(~LH,ncol=1,scales="free_y")+ggtitle(main)
}

#-----------------------------------------------------------------------------------------------
#                      Goodness of Fit
#-----------------------------------------------------------------------------------------------


# function to calculate posterior predictive p values with Freeman Tukey discrepency function
Freem_Tuk_P<-function(obs_dat_long,  # observed data
                      mscjs_fit,     # model object used for generating simulations
                      sim_posterior, # posterior samples
                      mscjs_dat,     #object containing some information on releases (stream, LH, etc. for housekeeping)
                      sim_rand, nsamps=250){     #should random year effects be samples from the hypderdistribution
  
  # tell model whether to sample random year effects from the hypderdistribution
  mscjs_fit$mod$env$data$sim_rand=sim_rand
  
  # set up empty vectors to holde Freeman-Tukey statistics
  FT_ref_vec<-FT_sim_vec<-numeric(nsamps)
  #set up empty matrix to hold posterior predictive simulations of data
  post_pred<-matrix(NA,nrow=nrow(obs_dat_long),ncol=nsamps)
  
  #turn off dplyer warnings
  options(dplyr.summarise.inform = FALSE)
  gc() #cleanup memory
  
  # loop through posterior samples
  for(i in 1:nsamps){
    # simulate data
    sim<-mscjs_fit$mod$simulate(par=sim_posterior[,i])
    
    # sumarize expected observations based on paramater set
    exp_det<-cbind(mscjs_dat$releases %>% ungroup() %>% select(LH,stream,sea_Year_p) , sim$det_1,sim$det_2,sim$det_3) %>% `colnames<-`(colnames(obs_dat %>% select(LH  :Tum_A_3))) %>% 
      #sum detection by stream, year, and life history
      group_by(LH,stream,sea_Year_p,) %>% summarise(across(LWe_J :Tum_A_3, sum)) %>% ungroup() %>% pivot_longer(LWe_J :Tum_A_3) %>% 
      filter(!(LH=="Unk" & name==("LWe_J")))
    # filter(!(LH=="Unk" & !name%in%c("Bon_J","McN_J")))
    
    # summarize simulated data based on parameter set
    sim_obs<-cbind(mscjs_dat$releases %>% ungroup() %>% select(LH,stream,sea_Year_p) , sim$sim_det_1,sim$sim_det_2,sim$sim_det_3) %>% `colnames<-`(colnames(obs_dat %>% select(LH  :Tum_A_3))) %>% 
      #sum detection by stream, year, and life history
      group_by(LH,stream,sea_Year_p) %>% summarise(across(LWe_J :Tum_A_3, sum)) %>% ungroup() %>% pivot_longer(LWe_J :Tum_A_3) %>% 
      filter(!(LH=="Unk" & name==("LWe_J")))
    # filter(!(LH=="Unk" &  !name%in%c("Bon_J","McN_J"))) 
    
    # save simulated data
    post_pred[,i]<-sim_obs$value
    
    # calculate Freeman-Tukey statistics for simulated and observed data
    FT_ref_vec[i]<- sum((sqrt(obs_dat_long$value)-sqrt(exp_det$value))^2) # observed
    FT_sim_vec[i]<- sum((sqrt(sim_obs$value)-sqrt(exp_det$value))^2)      # simulated
    if((i/50)%%1==0){gc()} #clear memory every 50 iterations
  }
  
  #calculate Bayesian P-value
  p<-sum(na.exclude(FT_sim_vec>FT_ref_vec))/length(na.exclude(FT_sim_vec>FT_ref_vec))
  
  return(list(p=p,post_pred=post_pred))
  
}


#--------------------------------------------------------------------------------------------------------

#function to calculate standardized quantile residuals
make_stan_res<-function(mscjs_fit,mscjs_dat,obs_dat_long,sim_rand,post_pred=NULL,n_samps=250){
  last_best<-mscjs_fit$last_par_best
   mscjs_fit$mod$env$data$sim_rand<-sim_rand
  
  # summarize the expected number of detections in each cell based on the fitted model parameters
  ## extract the "simulation object" from the fitted, which has values of expected detections.
  sim<-mscjs_fit$mod$simulate(par=last_best)
  ## summarize expected number of detections
  exp_det_MLE<-cbind(mscjs_dat$releases %>% ungroup() %>% select(LH,stream,sea_Year_p) , sim$det_1,sim$det_2,sim$det_3) %>% `colnames<-`(colnames(obs_dat %>% select(LH  :Tum_A_3))) %>% 
    #sum detection by stream, year, and life history
    group_by(LH,stream,sea_Year_p) %>% summarise(across(LWe_J :Tum_A_3, sum)) %>% ungroup() %>% pivot_longer(LWe_J :Tum_A_3) %>% 
    filter(!(LH=="Unk" & name==("LWe_J")))
  
  # #residuals
  # test<-(cbind(obs_dat_long %>% select(LH:name),value=((obs_dat_long$value)-(exp_det_MLE$value))))
  # #sum residuals across years
  # test %>% group_by(LH,stream,name) %>% summarise(across(value, sum)) %>% pivot_wider()
  # #sum residuals acorss years, streams, and LHs
  # test %>%  group_by(name) %>% summarise(across(value, sum))
  
  
  if(is.null(post_pred)){
  gc()
  mscjs_fit$mod$env$data$sim_rand<-0
  post_pred<-replicate(n_samps,{
    sim_data<-mscjs_fit$mod$simulate(par=last_best)
    # * simulated detection from model *
    sim_obs<-cbind(mscjs_dat$releases %>% ungroup() %>% select(LH,stream,sea_Year_p) , sim_data$sim_det_1,sim_data$sim_det_2,sim_data$sim_det_3) %>% `colnames<-`(colnames(obs_dat %>% select(LH  :Tum_A_3))) %>% 
      #sum detection by stream, year, and life history
      group_by(LH,stream,sea_Year_p) %>% summarise(across(LWe_J :Tum_A_3, sum)) %>% ungroup() %>% pivot_longer(LWe_J :Tum_A_3) %>% 
      filter(!(LH=="Unk" & name==("LWe_J")))
    
    # filter(!(LH=="Unk"  &  !name%in%c("Bon_J","McN_J")))
    sim_obs$value
  })
  }
  
  dharm_sim<-createDHARMa(simulatedResponse = post_pred,
                          observedResponse = obs_dat_long$value,
                          fittedPredictedResponse = exp_det_MLE$value,
                          integerResponse = TRUE)
  
  return(dharm_sim)
}
