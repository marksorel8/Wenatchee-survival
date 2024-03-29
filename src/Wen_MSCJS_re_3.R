library(here)
library(tidyverse)
library(TMB)

if(file.exists(here("Data","mark_file_CH.csv"))){ # if processed data file exists in directory
  mark_file_CH<-read.csv(here("Data","mark_file_CH.csv")) #read it
}else{
  #otherwise try to read it from the  where I stashed a copy on Zenodo
  try(mark_file_CH<-read_csv("https://zenodo.org/record/6216373/files/mark_file_CH.csv?download=1"))
  
}

if(!exists("mark_file_CH")){ #if cant get processed from project directory or off web, try processing raw data. Will only work if raw data files are in project directory
  source(here("src","data_proc.r"))
  save(mark_file_CH,file=here("Data","mark_file_CH.csv")) #save processed file
  
}

if(!exists("mark_file_CH")){ # send helpful warning if capture history data cant be loaded
  stop("Can't load processed capture history data in Wen_MSCJS_re_3.R. Connect to internet and try again or contact authors.")
}






#load environmental data covariates
if(file.exists(here("Data","env_dat.csv"))){
  env_dat<-read.csv(here("Data","env_dat.csv"))
}else{
  source(here("src","Env data funcs.r"))
  env_dat<-get_dat()
  write.csv(env_dat,file=here("Data","env_dat.csv"),row.names = FALSE)
}

make_dat<-function(mark_file_CH=mark_file_CH,sites=c("LWe_J","McN_J","JDD_J","Bon_J","Est_J","Bon_A","McN_A","PRa_A","RIs_A","Tum_A"),start_year=2006, end_year=2017,cont_cov,length_bin=5,doy_bin=10,inc_unk=FALSE,exc_unk=FALSE){


  
  #drop lower trap releases if not ussing
   if(exc_unk){mark_file_CH <-mark_file_CH %>% filter(LH!="Unk") %>% droplevels()}
  
  if(is.null(all_of(cont_cov))){
    dat_out<- mark_file_CH %>%  
      #subset columns needed for analysis
      select(sea_Year_p,LH,stream, #grouping variables
             all_of(sites))
  }else{if(length(all_of(cont_cov))==2){
    dat_out<- mark_file_CH %>%  
      #add grouped length and release day columns
      mutate(length_bin=ceiling(Length.mm/length_bin)*length_bin-(length_bin/2),
             rel_DOY_bin=ceiling((Mark.Day.Number+ifelse(LH=="smolt",365,0))/doy_bin)*doy_bin-(doy_bin/2) ,
             rel_DOY_bin=ifelse(LH=="Unk",rel_DOY_bin[1],rel_DOY_bin))  %>%
      #subset some very small or large length
      filter(length_bin>=55 &length_bin<=200 & rel_DOY_bin>10) %>%
      #subset columns needed for analysis
      select(sea_Year_p,LH,stream, #grouping variables
             all_of(sites),all_of(cont_cov))
  }else{if(all_of(cont_cov)=="rel_DOY_bin"){
    dat_out<- mark_file_CH %>%  
      #add grouped length and release day columns
      mutate(rel_DOY_bin=ceiling((Mark.Day.Number+ifelse(LH=="smolt",365,0))/doy_bin)*doy_bin-(doy_bin/2),
             rel_DOY_bin=ifelse(LH=="Unk",rel_DOY_bin[1],rel_DOY_bin)) %>%
      #subset some very small or large length
      filter(rel_DOY_bin>10) %>%

      #subset columns needed for analysis
      select(sea_Year_p,LH,stream, #grouping variables
             all_of(sites),all_of(cont_cov))
    
  }else{
    dat_out<- mark_file_CH %>%  
      #add grouped length and release day columns
      mutate(length_bin=ceiling(Length.mm/length_bin)*length_bin-(length_bin/2),
      ) %>%
      #subset some very small or large length
      filter(length_bin>=55 &length_bin<=200 ) %>%
      mutate(across(c(length_bin),scale)) %>% 
      #subset columns needed for analysis
      select(sea_Year_p,LH,stream, #grouping variables
             all_of(sites),all_of(cont_cov))
  }}}
  
  
  dat_out<- dat_out %>%  
    #subset columns needed for analysis
    select(sea_Year_p,LH,stream, #grouping variables
           all_of(sites),cont_cov) %>% 
    #sites/occasions to include in model
    #first year with all stream data through last year where data on all three return ages is available (because it is 2020)
    filter(sea_Year_p>=start_year & sea_Year_p<=end_year) %>%
    #make grouping variables factors
    mutate_at(vars(sea_Year_p:stream),~as.factor(as.character(.x))) %>% 
    #create multistate capture histories
    mutate(ch=select(., sites[1]: sites[length(sites)]) %>%  reduce(paste0)) %>%
    #mutate(ch=select(., McN_J,Bon_J,Est_J) %>% reduce(paste0)) %>%
    mutate(ch=paste0("1",ch)) %>% 
    #reduce data to unqiue capture history/ groups combos and counts
    group_by_all() %>% summarise(freq=n()) %>% as.data.frame() %>% 
    arrange(LH) # arrange so unknown LH (marked at lower wenatchee trap) comes last
  
#Occasion sites
occasion_sites<-colnames(select(dat_out,sites[1]: sites[length(sites)]))
#get number of site/occasions
nOCC<-nchar(dat_out$ch[1])-1
#get number of downstream sites/occasiosn
nDS_OCC<-sum(substr(occasion_sites,5,5)=="J")
#numbr of unique capture histories
n_unique_CH<-nrow(dat_out)
#number of unique CHs that are unknown
n_known_CH<-sum(dat_out$LH!="Unk")
#number of states
n_states<-3
#occasions corresponding to lower wenatchee and mcnary juveniles (for trap dependence)
# trap_dep<-which(sites%in%c("LWe_J", "McN_J"))


#~~~~
#extract the individual design data for each parameter
##phi/s survival
Phi.design.dat<-
  dat_out %>% 
  #release group
  select(sea_Year_p:stream,all_of(cont_cov)) %>% distinct() %>%
  # expand(sea_Year_p,LH,stream,get(cont_cov)) %>% #include combinations with no data
  # rename("rel_DOY_bin"=4) %>% 
  #add time
  arrange(LH) %>% 
  full_join(tibble(time=1:length(sites)),by=character()) %>% mutate(Time=time-1) %>% 
  #add trap dependency
  # full_join(tibble(LWe_J=0:1),by=character()) %>% 
  # full_join(tibble(McN_J=0:1),by=character()) %>% 
  #add stratum
  full_join(tibble(stratum=1:3),by=character()) %>% 
  #make time and stratum factor variables
  mutate(time=as.factor(time),stratum=as.factor(stratum)) %>% 
  #make group column
  mutate(group=select(., LH,stream,sea_Year_p,all_of(cont_cov)) %>%  reduce(paste0)) %>% 
  mutate(par.index=1:nrow(.)) %>% 
  filter(Time>(nDS_OCC)|stratum==1) %>%  #cant be in strata (fish age) other than 1 on downstream, or ocean for survival. Note "Time" column indexing starts at 0, so occasion/Time (nDS_OCC-1) is the last juvenile detection occasion
   mutate(mig_year=as.factor(ifelse(Time<=(nDS_OCC),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum)))),
         mig_year_num=as.numeric(mig_year)) %>%   #add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults
cbind(.,model.matrix(~time+stream+LH+stratum-1,data=.)) %>% 
 mutate(streamChiwawa=as.numeric(streamNason+streamWhite==0),
        age_class=as.factor(ifelse(LH%in%c("fall","summer"),"sub",as.character(LH))), 
        LHfall=as.numeric(LH=="fall"),age_0=as.numeric(LH!="smolt"), 
        LWe_new=as.numeric(as.numeric(as.character(sea_Year_p)>2011)),
        LH_t1=ifelse(LH=="Unk","smolt",as.character(LH)),
        stream_t1=ifelse(stream=="LWE","Chiwawa",as.character(stream))) %>% 
  #make length bin and release DOY numeric
  # mutate(across(.cols=all_of(cont_cov),.fns=function(x)scale(as.numeric(as.character(x))))) %>% 
  #make par index a sequence
  mutate(par.index=(1:nrow(.))-1) %>% 
  #add environmental covariates
  left_join(env_dat %>% mutate(mig_year=as.factor(mig_year)),by="mig_year") %>% 
  mutate(across(win_flow:last_col(),scale)) %>% 
  # mutate(across(sum_flow:pdo.aut,scale)) %>% 
 replace(is.na(.), 0) %>% # note these do not affect the likelihood. The only missing values are for marine covariates in recent years when adults returned but no juveniles went to sea (i.e. years when those covariates are used). Still, they end up in design matrix so cause things to crash if na.
 
  ungroup() %>%
  #add a column of 1's to use as a goruping variable when specifying penalized cemplexity priors
  mutate(one="1")
  # add column for first time
#downstream time
#ocean time
#upstream times
  
  
##need to set lower weantchee detection to 0 in 2011-2012 because no trap. Find out which occ is LWe_J
occ_LWe_J<-which(sites=="LWe_J")-1; if(length(occ_LWe_J)==0){occ_LWe_J<- -5}

##p detection, same as Phi
p.design.dat<-
  dat_out %>% 
  #release group
  select(sea_Year_p:stream,all_of(cont_cov)) %>% distinct() %>% 
  #add time
  full_join(tibble(time=2:length(sites)),by=character()) %>% mutate(Time=time-2) %>% 
  #add trap dependency
  # full_join(tibble(LWe_J=0:1),by=character()) %>% 
  # full_join(tibble(McN_J=0:1),by=character()) %>% 
  #add stratum
  full_join(tibble(stratum=1:3),by=character()) %>% 
  #make time and stratum factor variables
  mutate(time=as.factor(time),stratum=as.factor(stratum)) %>% 
  #make group column
  mutate(group=select(., LH,stream,sea_Year_p,all_of(cont_cov)) %>%  reduce(paste0)) %>% 
  mutate(par.index=1:nrow(.)) %>% 
  filter(Time>(nDS_OCC-1)|stratum==1) %>% # can be in multiple state for detection at time (nDS_OCC) but nor for survival.
  filter(Time<(nOCC-1)) %>%  # assuming detection at last time is 1, so not including this time in the detection design data
  mutate(mig_year=as.factor(ifelse(Time<=(nDS_OCC-1),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum))))) %>%   #add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults
  cbind(.,model.matrix(~time+stream+LH+stratum-1,data=.)) %>% 
  mutate(streamChiwawa=as.numeric(stream=="Chiwawa")) %>%
  #make length bin and release DOY numeric
  # mutate(across(.cols=all_of(cont_cov),.fns=function(x)scale(as.numeric(as.character(x))))) %>% 

  mutate( age_class=as.factor(ifelse(LH%in%c("fall","summer"),"sub",as.character(LH))),
          LHfall=as.numeric(LH=="fall"),age_0=as.numeric(LH!="smolt"),
          LWe_new=as.numeric(as.numeric(as.character(sea_Year_p)>2011)),
          LH_t1=as.factor(ifelse(LH=="Unk","smolt",as.character(LH))),
          stream_t1=as.factor(ifelse(stream=="LWE","Chiwawa",as.character(stream)))) %>% 
  dplyr::filter( !(as.numeric(as.character(sea_Year_p))%in%(2011:2012) & Time==occ_LWe_J )) %>%   
  
  filter((LH=="Unk" &stream=="LWE" )|(LH!="Unk" &stream!="LWE")) %>% 
  filter(!((LH=="Unk" &stream=="LWE")&(as.numeric(as.character(sea_Year_p))%in%(2011:2012)))) %>% #remove years with no releases at lower trap, so don't estimate proportions below

  #make par index a sequence
  mutate(par.index=(1:nrow(.))-1)%>% 
  #add environmental covariates
  left_join(env_dat %>% mutate(mig_year=as.factor(mig_year)),by="mig_year") %>% 
  mutate(across(sum_flow:last_col(),scale)) %>% 
  replace(is.na(.), 0) %>% # note these do not affect the likelihood. The only missing values are for marine covariates in recent years when adults returned but no juveniles went to sea (i.e. years when those covariates are used). Still, they end up in design matrix so cause things to crash if na.

  ungroup() %>% 
  #add a column of 1's to use as a goruping variable when specifying penalized cemplexity priors
  mutate(one="1")

# try(p.design.dat<-p.design.dat %>% mutate(LWe_J=as.numeric(as.character(LWe_J))))
# try(p.design.dat<-p.design.dat %>% mutate(McN_J=as.numeric(as.character(McN_J))))
# try(p.design.dat<-p.design.dat %>% mutate(LWe_J=as.numeric(as.character(LWe_J))))

##Psi transition. 
Psi.design.dat<-
  dat_out %>% 
  #release group
  select(sea_Year_p:stream,all_of(cont_cov)) %>% distinct() %>% 
  #add trap dependency
  # full_join(tibble(LWe_J=0:1),by=character()) %>% 
  # full_join(tibble(McN_J=0:1),by=character()) %>% 
  #add stratum
  full_join(tibble(stratum=1),by=character()) %>% 
  full_join(tibble(tostratum=as.factor(2:3)),by=character()) %>% 
  #make group column
  mutate(group=select(., LH,stream,sea_Year_p,all_of(cont_cov)) %>%  reduce(paste0)) %>% 
  mutate(par.index=1:nrow(.)) %>% 
  filter(stratum==1) %>%  #Can only transition from state 1 (Juvenile entering ocean)
  arrange(tostratum,group) %>%    # sort by stratum and group so the first half of rows represents the alr probs of transitioning to age 2 and the second half of the rows the alr probs of transitioning to state 3. This is neccesary for the way I am coding this in TMB, to take use the two halves of the vectors when doing the backtransformation from alr to simplex. 
  cbind(.,model.matrix(~stream-1+LH-1+tostratum-1,data=.)) %>% 
  mutate(streamChiwawa=as.numeric(streamNason+streamWhite==0), age_class=as.factor(ifelse(LH!="smolt","sub","yrlng")), LHfall=as.numeric(LH=="fall"),age_0=as.numeric(LH!="smolt")) %>%
  rename("mig_year"="sea_Year_p") %>% 
  #make length bin and release DOY numeric
  # mutate(across(.cols=all_of(cont_cov),.fns=function(x)scale(as.numeric(as.character(x))))) %>% 
  #make par index a sequence
  mutate(par.index=(1:nrow(.))-1)%>% 
  #add a column of 1's to use as a goruping variable when specifying penalized cemplexity priors
  mutate(one="1")


#~~~~
#make PIMs, which are matrices (nCH by nTimes) that give the index of the paramater wiithin the vector of phi or p parameters
##Phi pim. list, where each element is a PIM for a different state
###list of empty matrices (nCH x nOCC)
Phi_pim<-rep(list(matrix(NA,n_unique_CH,nOCC)),n_states)
p_pim<-rep(list(matrix(NA,n_unique_CH,(nOCC-1))),n_states) #dont include last tiem
####fill in pim matrix for Phi (survival)
for ( i in 1:(nDS_OCC + 1)){#loop over downstream "times" (e.g. McN_j, Bon_j, Est_j)and ocean
  #fill in matrix for state 1 (only possible state for downstream migration)
  Phi_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p,all_of(cont_cov)) %>% reduce(paste0),"time",i), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(Phi.design.dat$group,"time",Phi.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
  if(i<=nDS_OCC){
  p_pim[[1]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p,all_of(cont_cov)) %>% reduce(paste0),"time",i+1), #match group (LH,stream, seaward migration year) and occasion (time) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time))-1 #with corresponding row in design data (which will become corresponding element of parameter vector). subtract 1 becauuse TMB indexing starts at 0
  }
}

for ( i in (nDS_OCC+2):nOCC){ #loop over upstream "times"/occasions 
  for(j in 1:n_states){ #loop over "states" (i.e. age of fish/years at sea)
    Phi_pim[[j]][,i]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p,all_of(cont_cov)) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                            paste0(Phi.design.dat$group,"time",Phi.design.dat$time,"stratum",Phi.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)

  p_pim[[j]][,i-1]<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p,all_of(cont_cov)) %>% reduce(paste0),"time",i,"stratum",j), #match group (LH,stream, seaward migration year), occasion (time), and stratum (fish age/return year) for each CH
                          paste0(p.design.dat$group,"time",p.design.dat$time,"stratum",p.design.dat$stratum))-1#with corresponding row in design data (which will become corresponding element of parameter vector)

}
}

### replace NAs in p_pim (where detection fixed at 0) with value after largest value in p_pim
### becaue going to add a value of 0 on to the vector of detection probs. 
p_pim[[1]][is.na(p_pim[[1]])]<-nrow(p.design.dat)

##Psi pim. (nCH length vector) for each CH prob of row of ALR prob matrix of transition to states 2 or 3 (columns 1 or 2) corresponding to returning to Bonneville as adults after 2 or 3 years. The design data has been ordered such that the first half is for transition to state 2 and the second hallf for transition to state 3, so I only need the index corrsponding with the transitin to state 2 in the design data matrix
Psi_pim<-match(paste0(dat_out %>% select(LH,stream,sea_Year_p,all_of(cont_cov)) %>% reduce(paste0),"tostratum",2),
               paste0(Psi.design.dat$group,"tostratum",Psi.design.dat$tostratum))-1 #subtract 1 becauuse TMB indexing starts at 0

## *** Make data for simulating capture histories and calculating expected detections

#### number released per year, LH, stream, and continuous covariate bin
releases<-dat_out %>% select(LH,stream,(cont_cov),sea_Year_p,freq) %>% group_by(across(LH:sea_Year_p))%>% summarise(freq =sum(freq))

#cohorts thaty aren't releases form lower trap
n_known_CH_sim<-sum(releases$LH!="Unk")

### phi pim for simulation
phi_pim_sim<-inner_join(releases,Phi.design.dat %>% select(par.index,time,stratum,c("LH","stream","sea_Year_p",cont_cov))) %>% arrange(time,stratum) %>% # combine releases with design data
  pivot_wider(values_from=par.index,names_from=c(time,stratum,)) %>% ungroup()%>% select(paste(1:length(sites),1,sep="_"),paste((nOCC-nDS_OCC+2) :length(sites),2,sep="_"),paste((nOCC-nDS_OCC+2):length(sites),3,sep="_")) %>% as.matrix()




#### p pim for simulation
p_pim_sim<-inner_join(releases,p.design.dat %>% select(par.index,time,stratum,c("LH","stream","sea_Year_p",cont_cov))) %>% # combine releases with design data
  pivot_wider(values_from=par.index,names_from=c(time,stratum)) %>% ungroup() %>%  select(paste(2:length(sites),1,sep="_"), paste((nOCC-nDS_OCC+2) :length(sites),2,sep="_"), paste((nOCC-nDS_OCC+2):length(sites),3,sep="_")) %>%
  #fill in NAs (years to be set to detection of 0) with the correct index
  replace(is.na(.), nrow(p.design.dat)) %>%  as.matrix()


#### psi pim for simulation
psi_pim_sim<-left_join(releases %>% rename(mig_year=sea_Year_p)   %>% ungroup(),
  Psi.design.dat %>% select(par.index,stratum,tostratum,c("LH","stream",cont_cov,"mig_year",)) %>% filter(stratum==1&tostratum==2) %>% as_tibble() %>% distinct(across(stratum:mig_year),.keep_all=TRUE) )




## *** PIMS for calculating weighted averages for unknown LH fish
#placeholders
Phi.pim_unk<-NA
p.pim_unk<-NA
Phi.pim_unk_years <- NA
p.pim_unk_years<- NA
n_unk_LH_phi<-0
n_unk_LH_p<-0
n_known_LH_phi<-nrow(Phi.design.dat)
n_known_LH_p<-nrow(p.design.dat)

if(inc_unk){
####Subset design data to just unknown LH fish (marked in tributaries)
Phi.design.dat_unk <- Phi.design.dat %>% filter(LH=="Unk" &stream=="LWE") %>% droplevels()
#### number of rows or phi design data for known LH parameters
n_unk_LH_phi<-nrow(Phi.design.dat_unk)

p.design.dat_unk <- p.design.dat %>% filter(LH=="Unk" &stream=="LWE") %>% droplevels()
#### number of rows or phi design data for known LH parameters
n_unk_LH_p<-nrow(p.design.dat_unk)

# Psi.design.dat_unk <- Psi.design.dat %>% filter(LH=="Unk" &stream=="LWE") %>% droplevels()

####Subset design data to just known LH fish (marked in tributaries)
Phi.design.dat <- Phi.design.dat %>% filter(LH!="Unk" &stream!="LWE") %>% droplevels()
#### number of rows or phi design data for known LH parameters
n_known_LH_phi<-nrow(Phi.design.dat)

p.design.dat <- p.design.dat %>% filter(LH!="Unk" &stream!="LWE") %>% droplevels()
#### number of rows or phi design data for known LH parameters
n_known_LH_p<-nrow(p.design.dat)

# Psi.design.dat <- Psi.design.dat %>% filter(LH!="Unk" &stream!="LWE") %>% droplevels()
#### pims


#make real pims

Phi.pim_unk <- inner_join(Phi.design.dat_unk %>% select(sea_Year_p,time,stratum)%>% droplevels() %>% ungroup %>% distinct(), Phi.design.dat %>% filter(stream=="Chiwawa"&LH!="summer") %>% select(par.index,sea_Year_p,LH,time,stratum) %>% droplevels() %>% ungroup() %>% distinct(across(sea_Year_p:McN_J),.keep_all = TRUE),by=c("sea_Year_p","time","stratum")) %>% pivot_wider(values_from=par.index,names_from=c(LH)) 

p.pim_unk <- inner_join(p.design.dat_unk %>% select(sea_Year_p,time,stratum,)%>% droplevels() %>% ungroup %>% distinct(), p.design.dat %>% filter(stream=="Chiwawa"&LH!="summer") %>% select(par.index,sea_Year_p,LH,time,stratum) %>% droplevels() %>% ungroup() %>% distinct(across(sea_Year_p:McN_J),.keep_all = TRUE),by=c("sea_Year_p","time","stratum")) %>% pivot_wider(values_from=par.index,names_from=c(LH)) 

#no subyearlings released in 2006 so make subyearling parameters smolts parameters

Phi.pim_unk[Phi.pim_unk$sea_Year_p==2006,"fall"]<-Phi.pim_unk[Phi.pim_unk$sea_Year_p==2006,"smolt"]
p.pim_unk[p.pim_unk$sea_Year_p==2006,"fall"]<-p.pim_unk[p.pim_unk$sea_Year_p==2006,"smolt"]

Phi.pim_unk_years <- Phi.pim_unk%>% pull(sea_Year_p ) %>% as.numeric() - 1
p.pim_unk_years<- p.pim_unk %>% pull(sea_Year_p ) %>% as.numeric() - 1

Phi.pim_unk <- Phi.pim_unk %>% select(fall:smolt) %>% as.matrix()
p.pim_unk <- p.pim_unk %>% select(fall:smolt) %>% as.matrix()
}
# Psi.pim_unk <- inner_join(Psi.design.dat_unk %>% select(mig_year) %>% distinct(), Psi.design.dat %>% filter(stratum==1&tostratum==2&stream=="Chiwawa"&LH!="summer") %>% select(par.index,stratum,tostratum,mig_year,LH )%>% distinct(across(stratum:mig_year),.keep_all=TRUE) %>% droplevels()) %>% pivot_wider(values_from=par.index,names_from=c(LH))

#number of groups in the psi deisgn matrix
n_groups<-nrow(Psi.design.dat)/2

#occasions with trap dependent detection (effect of detection/non-detection at previous occasion)

#add number of fish in each group (freq) to design data for calculating averages
Phi.design.dat<-left_join(Phi.design.dat,releases) 
p.design.dat<-left_join(p.design.dat,releases)
Psi.design.dat<-left_join(Psi.design.dat,releases %>% rename(mig_year=sea_Year_p)   %>% ungroup())

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
            n_known_CH=n_known_CH,
            n_known_CH_sim=n_known_CH_sim,
            n_states=n_states,
            releases=releases,
            n_released=releases$freq,
            f_rel=ifelse(releases$LH=="Unk",1,0),
            phi_pim_sim=phi_pim_sim,
            p_pim_sim=p_pim_sim ,
            psi_pim_sim=psi_pim_sim$par.index,
            # TD_occ=TD_occ,
            # TD_i= TD_i,
            Phi.pim_unk=Phi.pim_unk,
            p.pim_unk=p.pim_unk ,
            Nyears=length(start_year:end_year),
            Phi.pim_unk_years = Phi.pim_unk_years,
            p.pim_unk_years= p.pim_unk_years ,
            n_unk_LH_phi=n_unk_LH_phi,
            n_unk_LH_p=n_unk_LH_p,
            n_known_LH_phi=n_known_LH_phi,
            n_known_LH_p=n_known_LH_p,
            f=ifelse(dat_out$LH=="Unk",1,0),
            inc_unk=inc_unk
            ))

}



fit_wen_mscjs<-function(x,phi_formula, p_formula, psi_formula,doFit=TRUE,silent=FALSE,sd_rep=TRUE,sim_rand=1,REML=FALSE,hypersd=1,map_hypers=c(FALSE,FALSE),pen=c(1,1),start_par=NULL){

#~~~~
#glmmTMB objects to get design matrices etc. for each parameter
## phi
Phi.design.glmmTMB<-glmmTMB::glmmTMB(formula(phi_formula), data=x$Phi.design.dat,dispformula = ~0,doFit=FALSE,contrasts = list(LH="contr.sum",stream="contr.sum",age_class="contr.sum"))
#,contrasts = list(LH="contr.sum",time="contr.sum",stream="contr.sum",age_class="contr.sum")
#time+time:LH+time:stream+diag(0+time|stream:LH:mig_year)
## p
p.design.glmmTMB<-glmmTMB::glmmTMB(formula(p_formula), data=x$p.design.dat,dispformula = ~0,doFit=FALSE,contrasts = list(LH_t1="contr.sum",stream_t1="contr.sum",LH="contr.sum",stream="contr.sum",age_class="contr.sum"))
#par.index~time+time:LH+time:stream
## psi
Psi.design.glmmTMB<-glmmTMB::glmmTMB(formula(psi_formula), data=x$Psi.design.dat,dispformula = ~0,doFit=FALSE,contrasts = list(age_class="contr.sum"))
#par.index~tostratum+tostratum:LH,




#set detection at the last site to be very high. I will fix this value, because otherwise it is unidentifiable, and other analyses  (i.e. looking at detection of adult fish that were picked up on instream arrays and whether they were detected at Tumwate) suggests detestion at Tumwater Dam IS 100%

#Make data for TMB
dat_TMB<-with(x,list(
  n_OCC=nOCC,
  nDS_OCC=nDS_OCC,
  n_states=n_states,
  n_groups=n_groups,
  n_unique_CH=n_unique_CH,
  # n_known_CH_sim=n_known_CH_sim,
  # n_known_CH=n_known_CH,
  CH=select(dat_out,sites[1]:sites[length(sites)]) %>% as.matrix(),
  freq=dat_out$freq,
  X_phi= Phi.design.glmmTMB$data.tmb$X,#[,-c(1:3)],
  X_p=  p.design.glmmTMB$data.tmb$X,
  X_psi= Psi.design.glmmTMB$data.tmb$X,
  Z_phi=Phi.design.glmmTMB$data.tmb$Z,
  Z_p=p.design.glmmTMB$data.tmb$Z,
  Z_psi=Psi.design.glmmTMB$data.tmb$Z,
  Phi_pim=Phi_pim, #indexing starts at 0
  p_pim=p_pim,
  Psi_pim=Psi_pim,
  phi_terms= Phi.design.glmmTMB$data.tmb$terms,
  p_terms= p.design.glmmTMB$data.tmb$terms,
  psi_terms= Psi.design.glmmTMB$data.tmb$terms,
  pen_prior=c(0,50), #halfnormal prior on penalized complexity rate parameters
  n_released=n_released,
  f_rel=f_rel,
  phi_pim_sim=phi_pim_sim,
  p_pim_sim=p_pim_sim ,
  psi_pim_sim=psi_pim_sim,
  f=f,
  hyper_mean=0,
  hyper_SD=hypersd,
  beta_phi_pen_ind=Phi.design.glmmTMB$data.tmb$X[1,-(1:(x$nOCC))] %>% names %>% substr(5,5) %>% as.factor(),
  beta_p_pen_ind=p.design.glmmTMB$data.tmb$X[1,-(1:(x$nOCC-1))] %>% names %>% substr(5,5) %>% as.factor() ,
  sim_rand = sim_rand #draw random effects from hyperdistribution in simulation rather than sampling from posterior.
))


#make param inits for TMB
if(!is.null(start_par)){
  par_TMB<-start_par  
}else

par_TMB<-list(
  beta_phi_ints=Phi.design.glmmTMB$parameters$beta[(1:(x$nOCC))], #intercept for each site and  LH 
  beta_phi_pen=Phi.design.glmmTMB$parameters$beta[-(1:(x$nOCC))],
  beta_p_ints=p.design.glmmTMB$parameters$beta[1:(x$nOCC-1)],#intercept for each site except the last and unique LH 
  beta_p_pen=p.design.glmmTMB$parameters$beta[-(1:(x$nOCC-1))],
  beta_psi_ints=Psi.design.glmmTMB$parameters$beta[1:2], #intercept for each strata
  beta_psi_pen=Psi.design.glmmTMB$parameters$beta[-(1:2)],
  log_pen_sds_phi=rep(0,length(Phi.design.glmmTMB$parameters$beta[-(1:(x$nOCC))])), # standard deviations of penalty priors
  log_pen_sds_p=rep(0,length(p.design.glmmTMB$parameters$beta[-(1:(x$nOCC-1))])),
  log_pen_sds_psi=rep(0,length(Psi.design.glmmTMB$parameters$beta[-(1:2)])),
  b_phi=Phi.design.glmmTMB$parameters$b,
  b_p=p.design.glmmTMB$parameters$b,
  b_psi=Psi.design.glmmTMB$parameters$b,
  theta_phi=Phi.design.glmmTMB$parameters$theta,
  theta_p=p.design.glmmTMB$parameters$theta,
  theta_psi=Psi.design.glmmTMB$parameters$theta,
  logit_p_subs=rep(0,(x$Nyears-2)),
  hyper_SD=hypersd,
  pen_phi=numeric(length(unique(dat_TMB$beta_phi_pen_ind)))+2,
  pen_p=numeric(length(unique(dat_TMB$beta_p_pen_ind)))+2,
  pen_psi=numeric(1)+2,
  pen=pen,
  pen_rand_phi=numeric(length(dat_TMB$phi_terms))+2,
  pen_rand_p=numeric(length(dat_TMB$p_terms))+2,
  pen_rand_psi=numeric(length(dat_TMB$psi_terms))+2,
  hyper_mean=0
)  
  fit<-NA
  mod<-NA

  #~~~~
  setwd(here("Src"))
  
  random<-c("b_phi","b_p","b_psi","beta_phi_pen","beta_p_pen","beta_psi_pen")
  if(REML){random<-c("beta_phi_ints","beta_p_ints","beta_psi_ints",random)}
  
  # set hyper parameters of distribution of proportion of subyearlings at trap as fixed.

  

#initialize model
  if(!x$inc_unk){  #model excluding unknown LH stream fish released at LWe_J (this is what is used in the paper)
    
    #compile TMB model if necessary
      TMB::compile("wen_mscjs_re_4.cpp")
    # load TMB model
    dyn.load(dynlib("wen_mscjs_re_4"))
    random=c(random,"pen_phi","pen_p","pen_psi","pen_rand_phi","pen_rand_p","pen_rand_psi")
mod<-TMB::MakeADFun(data=dat_TMB,parameters = par_TMB,random=random,DLL ="wen_mscjs_re_4", silent = silent)

  }else{ #model including unknown LH stream fish released at LWe_J
    
    #compile and load TMB model
      map<-list()
  if(map_hypers[1]){map$hyper_mean=factor(NA)}
  if(map_hypers[2]){map$hyper_SD=factor(NA)}
      random<-c(random,"logit_p_subs")
    TMB::compile("wen_mscjs_re_2.cpp")
    dyn.load(dynlib("wen_mscjs_re_2"))
    mod<-TMB::MakeADFun(data=dat_TMB,parameters = par_TMB,random=random,map=map,DLL ="wen_mscjs_re_2", silent = silent)
  }
  
  
  if(doFit){ 
    upper<-rep(Inf,length(mod$par))
try(fit<-TMBhelper::fit_tmb(mod,newtonsteps = 1,getsd = sd_rep,getJointPrecision = sd_rep))
}



return(list(fit=fit,
            last_par_best=mod$env$last.par.best,
            mod=mod,
            Phi.design.glmmTMB=Phi.design.glmmTMB,
            p.design.glmmTMB=p.design.glmmTMB,
            Psi.design.glmmTMB=Psi.design.glmmTMB,
            par_TMB=par_TMB,
            dat_TMB=dat_TMB))
}




