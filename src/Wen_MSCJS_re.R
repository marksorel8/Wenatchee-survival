library(here)
library(tidyverse)
library(TMB)
if(file.exists(here("src","mark_fil_ch_12032020.Rdata"))){
  load(here("src","mark_fil_ch_12032020.Rdata"))
}else{
source(here("src","data_proc.r"))
}


all_bio_data<-read_csv(here("data","all_bio_data.csv"))

make_dat<-function(mark_file_CH=mark_file_CH,sites=c("McN_J","JDD_J","Bon_J","Est_J","Bon_A","McJ_A","PRa_A","RIs_A","Tum_A")){

#subset columns needed for analysis
dat_out<-select(mark_file_CH,sea_Year_p,LH,stream, #grouping variables
                all_of(sites)) %>% #sites/occasions to include in model
  #first year with all stream data through last year where data on all three return ages is available (because it is 2020)
  filter(sea_Year_p>=2007&sea_Year_p<=2016) %>%
  #make grouping variables factors
  mutate_at(vars(sea_Year_p:stream),~as.factor(as.character(.x))) %>% 
  #create multistate capture histories
  mutate(ch=select(., sites[1]: sites[length(sites)]) %>%  reduce(paste0)) %>%
  #mutate(ch=select(., McN_J,Bon_J,Est_J) %>% reduce(paste0)) %>%
  mutate(ch=paste0("1",ch)) %>% 
  #reduce data to unqiue capture history/ groups combos and counts
  group_by_all() %>% summarise(freq=n()) %>% as.data.frame()


#Occasion sites
occasion_sites<-colnames(select(dat_out,sites[1]: sites[length(sites)]))
#get number of site/occasions
nOCC<-nchar(dat_out$ch[1])-1
#get number of downstream sites/occasiosn
nDS_OCC<-sum(substr(occasion_sites,5,5)=="J")
#numbr of unique capture histories
n_unique_CH<-nrow(dat_out)
#number of states
n_states<-3
# if Est_J included. This is the towed array. Will fix survival at 1 between Bonneville and Estuary Towed array. Therefore mortality in that period will end up in the ocean period. 
# Est_J_time<-ifelse(max(occasion_sites=="Est_J"), which(occasion_sites=="Est_J"),100)
# if JDD_J included, fix survival at 1 between Bonneville and Estuary Towed array. Therefore mortality in that period will end up between JDD and Bonneville_J
# JDD_J_time<-ifelse(max(occasion_sites=="JDD_J"), which(occasion_sites=="JDD_J"),100)
# if PRa_A or RIs_A included, fix survival at 1 from McN_A to these sites
# PRa_RIs_A_times<-match(c("PRa_A","RIs_A"),(occasion_sites))
# PRa_RIs_A_times[is.na(PRa_RIs_A_times)]<-100


#process data using RMark function. Specifies grouping variables for parameters, and type of model and hence parameters. "Multistrate" used S(Phi), p, and psi
wenatchee.processed<-RMark::process.data(dat_out,model="Multistrata",groups=c("LH","stream","sea_Year_p"))

#Make design data using RMark function. Sets up matrices for each parameter where there is a row for each combination of the grouping variables, occasion, age, etc. Specifiying the pim.type can reduce some of the combinations/# of rows. See?RMark::make.design.data() for more info. 
wenatchee.ddl<-RMark::make.design.data(wenatchee.processed,parameters=list(S=list(pim.type="time"), #survival changes by occasion and potentially stratum(i.e. fish age)
                                                                           p=list(pim.type="time"),#detection changes by occasion and potentially stratum(i.e. fish age)
                                                                           Psi=list(pim.type="constant")), # state (age at return) transitions only occur at one time period.
                                       remove.unused = TRUE)

#~~~~
#extract the individual design data for each parameter
##phi/s survival
Phi.design.dat<-wenatchee.ddl$S %>% 
  filter(Time>(nDS_OCC)|stratum==1) %>%  #cant be in strata (fish age) other than 1 on downstream, or ocean for survival. Note "Time" column indexing starts at 0, so occasion/Time (nDS_OCC-1) is the last juvenile detection occasion
  # filter(!Time%in%(c(JDD_J_time,Est_J_time,PRa_RIs_A_times)-1)) %>% 
   # filter(!Time%in%(c(JDD_J_time,Est_J_time)-1)) %>% 
  # if Est_J of JDD included, I will fix survival to be 1 between McN_J and JDD_J and Bon_J and Est_J within the cpp code, but want to pull it out of design matrix here
  mutate(mig_year=as.factor(ifelse(Time<=(nDS_OCC),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum)))),
         mig_year_num=as.numeric(mig_year)) %>%   #add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults
cbind(.,model.matrix(~time-1+stream-1+LH-1+stratum-1,data=.))
  
  
  # add column for first time
#downstream time
#ocean time
#upstream times
  
  
  
##p detection, same as Phi
p.design.dat<-wenatchee.ddl$p %>% 
  filter(Time>(nDS_OCC-1)|stratum==1) %>% # can be in multiple state for detection at time (nDS_OCC) but nor for survival.
  filter(Time<(nOCC-1)) %>%  # assuming detection at last time is 1, so not including this time in the detection design data
  mutate(mig_year=as.factor(ifelse(Time<=(nDS_OCC-1),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum))))) %>%   #add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults
  cbind(.,model.matrix(~time-1+stream-1+LH-1+stratum-1,data=.))

##Psi transition. 
Psi.design.dat<-wenatchee.ddl$Psi %>% 
  filter(stratum==1) %>%  #Can only transition from state 1 (Juvenile entering ocean)
  arrange(tostratum,group) %>%    # sort by stratum and group so the first half of rows represents the alr probs of transitioning to age 2 and the second half of the rows the alr probs of transitioning to state 3. This is neccesary for the way I am coding this in TMB, to take use the two halves of the vectors when doing the backtransformation from alr to simplex. 
  cbind(.,model.matrix(~stream-1+LH-1+tostratum-1,data=.)) %>% 
  rename("mig_year"="sea_Year_p")

#~~~~
#make PIMs, which are matrices (nCH by nTimes) that give the index of the paramater wiithin the vector of phi or p parameters
##Phi pim. list, where each element is a PIM for a different state
###list of empty matrices (nCH x nOCC)
Phi_pim<-rep(list(matrix(NA,n_unique_CH,nOCC)),n_states)
p_pim<-rep(list(matrix(NA,n_unique_CH,(nOCC-1))),n_states) #dont include last tiem
####fill in pim matrix for Phi (survival)
for ( i in 1:(nDS_OCC + 1)){#loop over downstream "times" (e.g. McN_j, Bon_j, Est_j)and ocean
  #fill in matrix for state 1 (only possible state for downstream migration)
  Phi_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(Phi.design.dat$group,"time",Phi.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
  if(i<=nDS_OCC){
  p_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i+1), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
  }
}

for ( i in (nDS_OCC+2):nOCC){ #loop over upstream "times"/occasions 
  for(j in 1:n_states){ #loop over "states" (i.e. age of fish/years at sea)
    Phi_pim[[j]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                            paste0(Phi.design.dat$group,"time",Phi.design.dat$time,"stratum",Phi.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)

  p_pim[[j]][,i-1]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time,"stratum",p.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)

}
}

##Psi pim. (nCH length vector) for each CH prob of row of ALR prob matrix of transition to states 2 or 3 (columns 1 or 2) corresponding to returning to Bonneville as adults after 2 or 3 years. The design data has been ordered such that the first half is for transition to state 2 and the second hallf for transition to state 3, so I only need the index corrsponding with the transitin to state 2 in the design data matrix
Psi_pim<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"tostratum",2),
               paste0(Psi.design.dat$group,"tostratum",Psi.design.dat$tostratum))-1 #subtract 1 becauuse TMB indexing starts at 0

#number of groups in the psi deisgn matrix
n_groups<-nrow(Psi.design.dat)/2


return(list(dat_out=dat_out,
            sites=sites,
            Phi.design.dat=Phi.design.dat,
            p.design.dat=p.design.dat,
            Psi.design.dat=Psi.design.dat,
            Phi_pim=Phi_pim,
            p_pim=p_pim,
            Psi_pim=Psi_pim,
            n_groups=n_groups,
            occasion_sites=occasion_sites,
            nOCC=nOCC,
            nDS_OCC=nDS_OCC,
            n_unique_CH=n_unique_CH,
            n_states=n_states
            # Est_J_time=Est_J_time,
            # JDD_J_time=JDD_J_time,
            # PRa_RIs_A_times=PRa_RIs_A_times
            ))

}


fit_wen_mscjs<-function(x,phi_formula, p_formula, psi_formula,doFit=TRUE,silent=FALSE){
#~~~~
#glmmTMB objects to get design matrices etc. for each parameter
## phi
Phi.design.glmmTMB<-glmmTMB::glmmTMB(formula(phi_formula), data=x$Phi.design.dat,dispformula = ~0,doFit=FALSE)
#time+time:LH+time:stream+diag(0+time|stream:LH:mig_year)
## p
p.design.glmmTMB<-glmmTMB::glmmTMB(formula(p_formula), data=x$p.design.dat,dispformula = ~0,doFit=FALSE)
#par.index~time+time:LH+time:stream
## psi
Psi.design.glmmTMB<-glmmTMB::glmmTMB(formula(psi_formula), data=x$Psi.design.dat,dispformula = ~0,doFit=FALSE)
#par.index~tostratum+tostratum:LH,

#+diag(0+time|stream:LH:mig_year)
#+us(tostratum|stream:LH:sea_Year_p)


#set detection at the last site to be very high. I will fix this value, because otherwise it is unidentifiable, and other analyses  (i.e. looking at detection of adult fish that were picked up on instream arrays and whether they were detected at Tumwate) suggests detestion at Tumwater Dam IS 100%

#Make data for TMB
dat_TMB<-with(x,list(
  n_OCC=nOCC,
  nDS_OCC=nDS_OCC,
  n_states=n_states,
  n_groups=n_groups,
  n_unique_CH=n_unique_CH,
  # Est_J_time=Est_J_time-1,
  # JDD_J_time=JDD_J_time-1,
  # PRa_RIs_A_times=PRa_RIs_A_times-1,
  CH=select(dat_out,sites[1]:sites[length(sites)]) %>% as.matrix(),
  freq=dat_out$freq,
  X_phi=Phi.design.glmmTMB$data.tmb$X,
  X_p=p.design.glmmTMB$data.tmb$X,
  X_psi=Psi.design.glmmTMB$data.tmb$X,
  Z_phi=Phi.design.glmmTMB$data.tmb$Z,
  Z_p=p.design.glmmTMB$data.tmb$Z,
  Z_psi=Psi.design.glmmTMB$data.tmb$Z,
  Phi_pim=Phi_pim, #indexing starts at 0
  p_pim=p_pim,
  Psi_pim=Psi_pim,
  phi_terms= Phi.design.glmmTMB$data.tmb$terms,
  p_terms= p.design.glmmTMB$data.tmb$terms,
  psi_terms= Psi.design.glmmTMB$data.tmb$terms
))


#make param inits for TMB
par_TMB<-list(
  beta_phi=Phi.design.glmmTMB$parameters$beta,
  beta_p=p.design.glmmTMB$parameters$beta,
  beta_psi=Psi.design.glmmTMB$parameters$beta,
  b_phi=Phi.design.glmmTMB$parameters$b,
  b_p=p.design.glmmTMB$parameters$b,
  b_psi=Psi.design.glmmTMB$parameters$b,
  theta_phi=Phi.design.glmmTMB$parameters$theta,
  theta_p=p.design.glmmTMB$parameters$theta,
  theta_psi=Psi.design.glmmTMB$parameters$theta
)  

if(doFit){
  #~~~~
  setwd(here("Src"))
  #compile and load TMB model
  TMB::compile("wen_mscjs_re.cpp")
  dyn.load(dynlib("wen_mscjs_re"))
#initialize model
mod<-TMB::MakeADFun(data=dat_TMB,parameters = par_TMB,random=c("b_phi","b_p","b_psi"),DLL ="wen_mscjs_re", silent = silent)

fit<-TMBhelper::fit_tmb(mod,newtonsteps = 1,getsd = TRUE)
}else{
  fit<-NA
  mod<-NA
}



return(list(fit=fit,
            mod=mod,
            Phi.design.glmmTMB=Phi.design.glmmTMB,
            p.design.glmmTMB=p.design.glmmTMB,
            Psi.design.glmmTMB=Psi.design.glmmTMB,
            par_TMB=par_TMB,
            dat_TMB=dat_TMB))
}




