library(here)
library(tidyverse)
library(TMB)
source(here("src","data_proc.r"))

#subset columns needed for analysis
dat_out<-select(mark_file_CH,sea_Year_p,LH,stream, #grouping variables
                McN_J,Bon_J,Est_J,Bon_A,Tum_A) %>% #sites/occasions to include in model
  #first year with all stream data through last year where data on all three return ages is available (because it is 2020)
  filter(sea_Year_p>=2007&sea_Year_p<=2017) %>%
  #make grouping variables factors
  mutate_at(vars(sea_Year_p:stream),~as.factor(as.character(.x))) %>% 
  #create multistate capture histories
  mutate(ch=select(., McN_J:Tum_A) %>%  reduce(paste0)) %>%
  #mutate(ch=select(., McN_J,Bon_J,Est_J) %>% reduce(paste0)) %>%
  mutate(ch=paste0("1",ch)) %>% 
  #reduce data to unqiue capture history/ groups combos and counts
  group_by_all() %>% summarise(freq=n()) %>% as.data.frame()


#Occasion sites
occasion_sites<-colnames(select(dat_out,McN_J:Tum_A))
#get number of site/occasions
nOCC<-nchar(dat_out$ch[1])-1
#get number of downstream sites/occasiosn
nDS_OCC<-sum(substr(occasion_sites,5,5)=="J")
#numbr of unique capture histories
n_unique_CH<-nrow(dat_out)
#number of states
n_states<-3
# is Est_J included. This is the towed array. Will fix survival at 1 between Bonneville and Estuary Towed array. Therefore mortality in that period will end up in the ocean period. 
Est_J_time<-ifelse(max(occasion_sites=="Est_J"), which(occasion_sites=="Est_J"),100)



#process data using RMark function. Specifies grouping variables for parameters, and type of model and hence parameters. "Multistrate" used S(Phi), p, and psi
wenatchee.processed<-RMark::process.data(dat_out,model="Multistrata",groups=c("LH","stream","sea_Year_p"))

#Make design data using RMark function. Sets up matrices for each parameter where there is a row for each combination of the grouping variables, occasion, age, etc. Specifiying the pim.type can reduce some of the combinations/# of rows. See?RMark::make.design.data() for more info. 
wenatchee.ddl<-RMark::make.design.data(wenatchee.processed,parameters=list(S=list(pim.type="time"), #survival changes by occasion and potentially stratum(i.e. fish age)
                                                                           p=list(pim.type="time"),#detection changes by occasion and potentially stratum(i.e. fish age)
                                                                           Psi=list(pim.type="constant")), # state (age at return) transitions only occur at one time period.
                                       remove.unused = TRUE)


#extract the individual design data for each parameter
##phi/s survival
Phi.design.dat<-wenatchee.ddl$S %>% 
  filter(Time>(nDS_OCC-1)|stratum==1) %>%  #cant be in strata (fish age) other than 1 on downstream, or ocean for survival. Note "Time" column indexing starts at 0, so occasion/Time (nDS_OCC-1) is the last juvenile detection occasion
  filter(Time!=(Est_J_time-1)) %>% 
  # if Est_J included, I will fix survival to be 1 between Bon_J and Est_J within the cpp code, but want to pull it out of design matrix here
  mutate(mig_year=as.factor(ifelse(Time<=(nDS_OCC-1),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum)))))#add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults

##p detection, same as Phi
p.design.dat<-wenatchee.ddl$p %>% 
  filter(Time>(nDS_OCC-1)|stratum==1) %>% # can be in multiple state for detection at time (nDS_OCC) but nor for survival.
  filter(Time<(nOCC-1)) # assuming detection at last time is 1, so not including this time in the detection design data
 


  
  
##Psi transition. 
Psi.design.dat<-wenatchee.ddl$Psi %>% 
  filter(stratum==1) %>%  #Can only transition from state 1 (Juvenile entering ocean)
  arrange(tostratum,group) # sort by stratum and group so the first half of rows represents the alr probs of transitioning to age 2 and the second half of the rows the alr probs of transitioning to state 3. This is neccesary for the way I am coding this in TMB, to take use the two halves of the vectors when doing the backtransformation from alr to simplex. 

#make PIMs, which are matrices (nCH by nTimes) that give the index of the paramater wiithin the vector of phi or p parameters

##Phi pim. list, where each element is a PIM for a different state
###list of empty matrices (nCH x nOCC)
Phi_pim<-rep(list(matrix(NA,n_unique_CH,nOCC)),n_states)
p_pim<-rep(list(matrix(NA,n_unique_CH,(nOCC-1))),n_states) #dont include last tiem
####fill in pim matrix for Phi (survival)
for ( i in 1:nDS_OCC){#loop over downstream "times" (e.g. McN_j, Bon_j, Est_j)
  #fill in matrix for state 1 (only possible state for downstream migration)
  Phi_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(Phi.design.dat$group,"time",Phi.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
  
  p_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i+1), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
}

for ( i in (nDS_OCC+1):nOCC){ #loop over upstream "times"/occasions
  for(j in 1:n_states){ #loop over "states" (i.e. age of fish/years at sea)
    Phi_pim[[j]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                            paste0(Phi.design.dat$group,"time",Phi.design.dat$time,"stratum",Phi.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)
  if(i<nOCC){ # dont include last time because assuming perfecft detection
  p_pim[[j]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i+1,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time,"stratum",p.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)
  }
}
}


##Psi pim. (nCH length vector) for each CH prob of row of ALR prob matrix of transition to states 2 or 3 (columns 1 or 2) corresponding to returning to Bonneville as adults after 2 or 3 years. The design data has been ordered such that the first half is for transition to state 2 and the second hallf for transition to state 3, so I only need the index corrsponding with the transitin to state 2 in the design data matrix
Psi_pim<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"tostratum",2),
               paste0(Psi.design.dat$group,"tostratum",Psi.design.dat$tostratum))-1 #subtract 1 becauuse TMB indexing starts at 0

#number of groups in the psi deisgn matrix
n_groups<-nrow(Psi.design.dat)/2

#glmmTMB objects to get design matrices etc. for each parameter
## phi
Phi.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+time+LH+diag(LH:time|mig_year),data=Phi.design.dat,dispformula = ~0,doFit=FALSE)
## p
p.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+time,data=p.design.dat,dispformula = ~0,doFit=FALSE)

## psi
Psi.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+tostratum,data=Psi.design.dat,dispformula = ~0,doFit=FALSE)

setwd(here("Src"))
#compile and load TMB model
TMB::compile("wen_mscjs_re.cpp")
dyn.load(dynlib("wen_mscjs_re"))

#set detection at the last site to be very high. I will fix this value, because otherwise it is unidentifiable, and other analyses  (i.e. looking at detection of adult fish that were picked up on instream arrays and whether they were detected at Tumwate) suggests detestion at Tumwater Dam IS 100%

#Make data for TMB
dat_TMB<-list(
  n_OCC=nOCC,
  nDS_OCC=nDS_OCC,
  n_states=n_states,
  n_groups=n_groups,
  n_unique_CH=n_unique_CH,
  Est_J_time=Est_J_time-1,
  CH=select(dat_out,McN_J:Tum_A) %>% as.matrix(),
  freq=dat_out$freq,
  X_phi=Phi.design.glmmTMB$data.tmb$X,
  X_p=p.design.glmmTMB$data.tmb$X,
  X_psi=Psi.design.glmmTMB$data.tmb$X,
  Z_phi=Phi.design.glmmTMB$data.tmb$Z,
  Phi_pim=Phi_pim, #indexing starts at 0
  p_pim=p_pim,
  Psi_pim=Psi_pim,
  fix_p_last=numeric(0),
  phi_terms= Phi.design.glmmTMB$data.tmb$terms
)


#make param inits for TMB
par_TMB<-list(
  beta_phi=Phi.design.glmmTMB$parameters$beta,
  beta_p=p.design.glmmTMB$parameters$beta,
  beta_psi=Psi.design.glmmTMB$parameters$beta,
  b_phi=Phi.design.glmmTMB$parameters$b,
  theta_phi=Phi.design.glmmTMB$parameters$theta
)  


#initialize model
mod<-TMB::MakeADFun(data=dat_TMB,parameters = par_TMB,random="b_phi",DLL ="wen_mscjs_re")

TMBhelper::fit_tmb(mod,getsd = FALSE)
#fit model
for ( i in 1:3){
  fit<-nlminb(mod$env$last.par.best,mod$fn,mod$gr)
  }
fit
#sd report
sd_rep<-sdreport(mod)
sd_rep
#report
mod_rep<-mod$report()
head(mod_rep$psi)
mod_rep$phi[1:20]
mod_rep$p[1:3]

Phi.design.dat <- Phi.design.dat %>% mutate(phi_fit=mod_rep$phi)

dev.new()
ggplot(data= as_tibble(Phi.design.dat),aes(x=as.numeric(mig_year) ,y=phi_fit,color=LH)) + geom_point()+geom_line()+facet_wrap(~time,scales="free")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mark model for comparison
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_out<-select(mark_file_CH,sea_Year_p,LH,stream, #grouping variables
                McN_J,Bon_J,Bon_A,Tum_A) %>% #sites/occasions to include in model
  #first year with all stream data through last year where data on all three return ages is available (because it is 2020)
  filter(sea_Year_p>=2007&sea_Year_p<=2017) %>%
  #make grouping variables factors
  mutate_at(vars(sea_Year_p:stream),~as.factor(as.character(.x))) %>% 
  #create multistate capture histories
  mutate(ch=select(., McN_J:Tum_A) %>%  reduce(paste0)) %>%
  #mutate(ch=select(., McN_J,Bon_J,Est_J) %>% reduce(paste0)) %>%
  mutate(ch=paste0("1",ch)) 

wenatchee.processed=process.data(dat_out,model="Multistrata",groups=c("LH","stream"))
#wenatchee.processed=process.data(test_dat,model="MSCJS",strata.labels=c("A","B","C"))
wenatchee.ddl=make.design.data(wenatchee.processed)



#
# p
#
#
# States B&C dont exist for juveniles (times 2-5)

#wenatchee.ddl$p<-wenatchee.ddl$p[-wenatchee.ddl$p$Time<=2&wenatchee.ddl$p$stratum!="A",]
#wenatchee.ddl$S<-wenatchee.ddl$S[-wenatchee.ddl$S$Time<=2&wenatchee.ddl$S$stratum!="A",]



#column for migration year where onupstream state A is seaward year +1, b +2, and c +3
seaward_year_vec<-as.numeric(substr(wenatchee.ddl$p$sea_Year_p,2,5))
wenatchee.ddl$p$det_year<-factor(ifelse(wenatchee.ddl$p$Time<2,seaward_year_vec,ifelse(wenatchee.ddl$p$stratum==1,seaward_year_vec+1,ifelse(wenatchee.ddl$p$stratum==2,seaward_year_vec+2,seaward_year_vec+3))))



#wenatchee.ddl$p$fix=ifelse((wenatchee.ddl$p$stratum!="A"&wenatchee.ddl$p$Time<=4),0,NA)
# set detection at Tumwater adult at 100% based on other work
wenatchee.ddl$p$fix[wenatchee.ddl$p$time==5]=1


#
# Psi
#
# can only change state (other than death) at time 4
# rest are fixed values
wenatchee.ddl$Psi$fix=NA
# stay in current state (year) other than time 4 (ocean)
wenatchee.ddl$Psi$fix[wenatchee.ddl$Psi$stratum==
                        wenatchee.ddl$Psi$tostratum & wenatchee.ddl$Psi$time!=(3)]=1

wenatchee.ddl$Psi$fix[wenatchee.ddl$Psi$stratum!=
                        wenatchee.ddl$Psi$tostratum &wenatchee.ddl$Psi$time!=(3)]=0

wenatchee.ddl$Psi$fix[wenatchee.ddl$Psi$stratum!=1 &
                        wenatchee.ddl$Psi$time==3]=0



# 
# S
#
#Add a columna for whether first occasion, which is necessarily  different for subyearlings and yearlings
wenatchee.ddl$S$time_1<-ifelse(wenatchee.ddl$S$time==1,1,0)
#
#column for yearlings vs subs
wenatchee.ddl$S$Yrlng<-ifelse(wenatchee.ddl$S$LH=="smolt",1,0)
#
#columns for whether white river, which might be different because above lake Wenatchee
wenatchee.ddl$S$White<-ifelse(wenatchee.ddl$S$stream=="White",1,0)
#
#
#column for migration year where onupstream state A is seaward year +1, b +2, and c +3
seaward_year_vec<-as.numeric(substr(wenatchee.ddl$S$sea_Year_p,2,5))
wenatchee.ddl$S$det_year<-factor(ifelse(wenatchee.ddl$S$Time<2,seaward_year_vec,ifelse(wenatchee.ddl$S$stratum=="A",seaward_year_vec+1,ifelse(wenatchee.ddl$S$stratum=="B",seaward_year_vec+2,seaward_year_vec+3))))


wenatchee.ddl$S$upstream_time<-ifelse(wenatchee.ddl$S$Time>2,1,0)


# fixing S to 1 in strata where fish cant be
wenatchee.ddl$S$fix<-NA
wenatchee.ddl$S$fix[wenatchee.ddl$S$stratum!=1 &
                      wenatchee.ddl$S$Time<3]=1 #using numeric "Time" which starts at 0 hence 3 is 4





# fit model
p.timexstratum.tag=list(formula=~-1+time:LH)
Psi.sxtime=list(formula=~-1+tostratum)
S.stratumxtime=list(formula=~-1+time:LH)



test<-RMark::mark(wenatchee.processed,wenatchee.ddl,
                  model.parameters=list(S=S.stratumxtime,p= p.timexstratum.tag,Psi=Psi.sxtime),run=TRUE)

real_test<-get.real(test,par="p")
real_test$`Group:LHsmolt.streamChiwawa.sea_Year_py2007`

test$results$beta

sd_rep$value

#compare parameters
sd_rep$par.fixed - test$results$beta[,1]
#compare standard errors
sqrt(diag(sd_rep$cov.fixed)) - test$results$beta[,2]
