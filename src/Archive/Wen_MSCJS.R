library(here)
library(tidyverse)
library(TMB)
source(here("src","data_proc.r"))

#subset columns needed for analysis
dat_out<-select(mark_file_CH,sea_Year_p,LH,stream, #grouping variables
                McN_J,Bon_J,Bon_A,Tum_A) %>% #sites/occasions to include in model
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


#get number of site/occasions
nOCC<-nchar(dat_out$ch[1])-1
#get number of downstream sites/occasiosn
nDS_OCC<-sum(substr(colnames(dat_out %>% select( McN_J:Tum_A)),5,5)=="J")
#numbr of unique capture histories
n_unique_CH<-nrow(dat_out)
#number of states
n_states<-3

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
  filter(Time>(nDS_OCC-1)|stratum==1) #cant be in strata (fish age) other than 1 on downstream, or ocean for survival. Note "Time" column indexing starts at 0, so occasion/Time (nDS_OCC-1) is the last juvenile detection occasion

##p detection, same as Phi
p.design.dat<-wenatchee.ddl$p %>% 
  filter(Time>(nDS_OCC-1)|stratum==1)# can be in multiple strate for detection at time (nDS_OCC) but nor for survival. Will fix detection at last time to 1 down below when nitializing model.

##Psi transition. 
Psi.design.dat<-wenatchee.ddl$Psi %>% 
  filter(stratum==1) %>%  #Can only transition from state 1 (Juvenile entering ocean)
  arrange(tostratum,group) # sort by stratum and group so the first half of rows represents the alr probs of transitioning to age 2 and the second half of the rows the alr probs of transitioning to state 3. This is neccesary for the way I am coding this in TMB, to take use the two halves of the vectors when doing the backtransformation from alr to simplex. 

#make PIMs, which are matrices (nCH by nTimes) that give the index of the paramater wiithin the vector of phi or p parameters

##Phi pim. list, where each element is a PIM for a different state
###list of empty matrices (nCH x nOCC)
Phi_pim<-rep(list(matrix(NA,n_unique_CH,nOCC)),n_states)
####fill in pim matrix for Phi (survival)
for ( i in 1:nDS_OCC){#loop over downstream "times" (e.g. McN_j, Bon_j, Est_j)
  #fill in matrix for state 1 (only possible state for downstream migration)
  Phi_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(Phi.design.dat$group,"time",Phi.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
}

for ( i in (nDS_OCC+1):nOCC){ #loop over upstream "times"/occasions
  for(j in 1:n_states) #loop over "states" (i.e. age of fish/years at sea)
    Phi_pim[[j]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                            paste0(Phi.design.dat$group,"time",Phi.design.dat$time,"stratum",Phi.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)
}

##phi_pim and p_pim are the same. There is a detection and survival for each occasion, although I will fix detection at the last occasion below.
p_pim<-Phi_pim

##Psi pim. (nCH length vector) for each CH prob of row of ALR prob matrix of transition to states 2 or 3 (columns 1 or 2) corresponding to returning to Bonneville as adults after 2 or 3 years. The design data has been ordered such that the first half is for transition to state 2 and the second hallf for transition to state 3, so I only need the index corrsponding with the transitin to state 2 in the design data matrix
Psi_pim<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p) %>% reduce(paste0),"tostratum",2),
               paste0(Psi.design.dat$group,"tostratum",Psi.design.dat$tostratum))-1 #subtract 1 becauuse TMB indexing starts at 0

#number of groups in the psi deisgn matrix
n_groups<-nrow(Psi.design.dat)/2

#glmmTMB objects to get design matrices etc. for each parameter
## phi
Phi.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+time,data=Phi.design.dat,dispformula = ~0,doFit=FALSE)
## p
p.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+time,data=p.design.dat,dispformula = ~0,doFit=FALSE)

## psi
Psi.design.glmmTMB<-glmmTMB::glmmTMB(par.index~-1+tostratum,data=Psi.design.dat,dispformula = ~0,doFit=FALSE)

setwd(here("Src"))
#compile and load TMB model
TMB::compile("wen_mscjs.cpp")
dyn.load(dynlib("wen_mscjs"))

#set detection at the last site to be very high. I will fix this value, because otherwise it is unidentifiable, and other analyses  (i.e. looking at detection of adult fish that were picked up on instream arrays and whether they were detected at Tumwate) suggests detestion at Tumwater Dam IS 100%

#Make data for TMB
dat_TMB<-list(
  n_OCC=nOCC,
  nDS_OCC=nDS_OCC,
  n_states=n_states,
  n_groups=n_groups,
  n_unique_CH=n_unique_CH,
  CH=select(dat_out,McN_J:Tum_A) %>% as.matrix(),
  freq=dat_out$freq,
  X_phi=Phi.design.glmmTMB$data.tmb$X,
  X_p=p.design.glmmTMB$data.tmb$X,
  X_psi=Psi.design.glmmTMB$data.tmb$X,
  Phi_pim=Phi_pim, #indexing starts at 0
  p_pim=p_pim,
  Psi_pim=Psi_pim,
  fix_p_last=numeric(0)
)


#make param inits for TMB
par_TMB<-list(
  beta_phi=Phi.design.glmmTMB$parameters$beta,
  beta_p=p.design.glmmTMB$parameters$beta,
  beta_psi=Psi.design.glmmTMB$parameters$beta
)  


#map to fix parameter
map<-list(
  beta_p=factor(seq(length(par_TMB$beta_p )))
  #  beta_psi=factor(c(NA,NA))
)
#fix detection prob at last time at 1
map$beta_p[substr(colnames(dat_TMB$X_p),5,5)==(nOCC+1)]<-NA #last occastion detection is mapped
map$beta_p<-map$beta_p %>% droplevels()
fix_p_last<-which(p.design.dat$time==(nOCC+1))-1 # fix lest occasion detection at 1 (see TMB model .cpp file for details)
dat_TMB$fix_p_last<-fix_p_last

#initialize model
mod<-TMB::MakeADFun(data=dat_TMB,parameters = par_TMB,map=map,DLL ="wen_mscjs")
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
mod_rep$phi[1:4]
mod_rep$p[1:3]

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
p.timexstratum.tag=list(formula=~-1+time)
Psi.sxtime=list(formula=~-1+tostratum)
S.stratumxtime=list(formula=~-1+time)



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
