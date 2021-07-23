library(here)
library(tidyverse)


#Wild or natural origin spring Chinook captured by screw trap in the Chiwawa, Nason, or White
mark_file<-read_csv(here("Data","ptagis","Tagging detail.csv")) %>% 
  filter(`Rear Type Code`=="W",
         `Run Name`=="Spring",
         `Species Name`=="Chinook",
         `Capture Method Code`=="SCREWT"
         )%>% mutate(mark_time=0)

#Wild or natural origin spring Chinook captured by screw trap in the Lower Wenatchee Trap
LWe_mark_file<-read_csv(here("Data","ptagis","Lower Wen Tagging detail.csv")) %>% 
  filter(`Rear Type Code`=="W",
         `Run Name`=="Spring",
         `Species Name`=="Chinook",
         `Capture Method Code`=="SCREWT"
  ) %>% mutate(mark_time=1)



mark_file<-rbind(mark_file,LWe_mark_file)
rm(LWe_mark_file)

#code for getting day of year from ptagis date column
#test<-lower_mid_col_dams %>% mutate(`First Date MMDDYYYY2`=as.Date(`First Date MMDDYYYY`,format="%m/%d/%Y"), doy=lubridate::yday(`First Date MMDDYYYY2`))

# obs_dat<-read_csv(here("Data","ptagis","ptagis_obs_data_wen_trib_screwt.csv")) %>% 
#   mutate(date_det=lubridate::date(lubridate::mdy_hms(obs_date))) %>% 
# group_by(tag_id,obs_site) %>% mutate(first_date=min(date_det),last_det=max(date_det)) %>% 
#   select(!obs_date:date_det) %>% droplevels() %>% 
#  distinct() %>% ungroup() %>% 
#   mutate("First Year YYYY"=lubridate::year(first_date),"Last Year YYYY"=lubridate::year(last_det)) %>% 
#   mutate("First DOY"=lubridate::yday(first_date),"Last DOY"=lubridate::yday(last_det)) %>% 
#   rename("Site Code Value"=obs_site,"Tag Code"=tag_id)

#recaptures in lower Wenatchee dams
lower_wen_traps_recaps<-read_csv(here("Data","ptagis","ptagis_recap_data_wen_trib_screwt.csv")) %>% mutate(recap_date=lubridate::mdy_hms(recap_date),"First Year YYYY"=lubridate::year(recap_date),"First Day Num"=lubridate::yday(recap_date)) %>% rename("Tag Code" = tag_id ) %>% arrange(recap_date)

#Lower and mid Columbia mainstem dam detections
lower_mid_col_dams<-read_csv(here("Data","ptagis","mid_lower_Col_Interrogation_Summary.csv"))

#upper Columbia mainstem Dam and Tumwater Dam detections
upper_col_dams_tum<-read_csv(here("Data","ptagis","upper_Col_Wen_Interrogation_summary.csv"))

#detections on instream arrays
inst_array_det<-read_csv(here("Data","ptagis","Interrogation Summary wen trib tag adult array det.csv"))



#detection data for fish released in lower Wenatchee trap from Dan W. 
obs_dat_LWe_releases<-read_csv(here("Data","ptagis","ptagis_obs_data_lower_wen_screwt.csv")) %>%
  mutate(date_det=lubridate::date(lubridate::mdy_hms(obs_date))) %>%
  group_by(tag_id,obs_site) %>% mutate(first_date=min(date_det),last_det=max(date_det)) %>%
  select(!obs_date:date_det) %>% droplevels() %>%
  distinct() %>% ungroup() %>%
  mutate("First Year YYYY"=lubridate::year(first_date),"Last Year YYYY"=lubridate::year(last_det)) %>%
  mutate("First Day Num"=lubridate::yday(first_date),"Last Day Num"=lubridate::yday(last_det)) %>%
  rename("Site Code Value"=obs_site,"Tag Code"=tag_id)

#add data for detections of fish released at lower wenatchee trap
lower_mid_col_dams<-bind_rows(lower_mid_col_dams,obs_dat_LWe_releases %>% filter(`Site Code Value`%in% c("MCJ","JDJ","B1J","B2J","BCC","TWX","BO1", "BO2", "BO3", "BO", "MC1", "MC2","JO1", "JO2","TD1", "TD2")))

upper_col_dams_tum<-bind_rows(upper_col_dams_tum,obs_dat_LWe_releases %>% filter(`Site Code Value`%in% c(  "PRA","RIA","TUF")))


#load length cutoffs for subyearling/yearling delineation for all days <=179
load(here("Data","cutoffs_and_props.Rdata"))


#add capture histories to mark file
mark_file_CH<-mark_file %>% 

#Lower Wenatchee Juvenile recapture
 left_join(lower_wen_traps_recaps %>%
              filter(`recap_site` %in% c("WENA4T","WENATT")) %>% #subset McNary juvenile detection sites
              select(c(`Tag Code`,`First Year YYYY`,`First Day Num`)) %>% #take columns first detection year and pit code
             # arrange(recap_date) %>% #arrange in order of increasing recap_date
             distinct(`Tag Code`,.keep_all=TRUE) %>% #remove duplicate detections at multiple "sites" at this trap 
              rename("LWe_J"="First Year YYYY",
                     "LWe_J_doy"="First Day Num"),by="Tag Code") %>%   #rename to represent site for merging
  
  #add release doy for LWe_J
  mutate(LWe_J_doy=ifelse(mark_time==1,`Release Day Number`,LWe_J_doy)) %>% 

#McNary Juvenile detection
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` =="MCJ") %>% #subset McNary juvenile detection sites 
            select(c(`Tag Code`,`First Year YYYY`,`First Day Num`)) %>% #take columns first detection year and pit code 
            rename("McN_J"="First Year YYYY",
                   "McN_J_doy"="First Day Num"),by="Tag Code") %>%   #rename to represent site for merging

#John Day Juvenle detection
  left_join(lower_mid_col_dams %>% 
              filter(`Site Code Value` =="JDJ") %>% 
              select(c(`Tag Code`,`First Year YYYY`,`First Day Num`)) %>%
              rename("JDD_J"="First Year YYYY",
                     "JDD_J_doy"="First Day Num"),by="Tag Code") %>%

#Bonneville Juvenile detection
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("B1J","B2J","BCC")) %>% 
            select(c(`Tag Code`,`First Year YYYY`,`First Day Num`)) %>% 
            arrange(`First Year YYYY`) %>% #arrange in order of increasing first year
            distinct(`Tag Code`,.keep_all=TRUE) %>% #remove duplicate detections at multiple "sites" at this dam 
            rename("Bon_J"="First Year YYYY",
                   "Bon_J_doy"="First Day Num"),by="Tag Code") %>%   

#estuary towed array
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` =="TWX") %>% 
            select(c(`Tag Code`,`First Year YYYY`,`First Day Num`)) %>%
            rename("Est_J"="First Year YYYY",
                   "Est_J_doy"="First Day Num"),by="Tag Code") %>% 

#Bonneville Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("BO1", "BO2",  "BO3", "BO")) %>% 
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>% 
            arrange(desc(`Last Year YYYY`)) %>% #latest year detected first, so when we take distinct below that is the one we keep
            distinct(`Tag Code`,.keep_all=TRUE) %>% 
            rename("Bon_A"="Last Year YYYY",
                   "Bon_A_doy"="Last Day Num"),by="Tag Code") %>% 

#The Dalles Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("TD1", "TD2")) %>% 
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>%
            arrange(desc(`Last Year YYYY`)) %>%
            distinct(`Tag Code`,.keep_all=TRUE) %>%
            rename("TDD_A"="Last Year YYYY",
                   "TDD_A_doy"="Last Day Num"),by="Tag Code") %>%

#John Day Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("JO1", "JO2")) %>% 
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>%
            arrange(desc(`Last Year YYYY`)) %>%
            distinct(`Tag Code`,.keep_all=TRUE) %>%
            rename("JDD_A"="Last Year YYYY",
                   "JDD_A_doy"="Last Day Num"),by="Tag Code") %>% 

#McNary Adult
left_join(lower_mid_col_dams %>% 
          filter(`Site Code Value` %in% c("MC1", "MC2")) %>% 
          select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>% 
          arrange(desc(`Last Year YYYY`)) %>%
          distinct(`Tag Code`,.keep_all=TRUE) %>%
          rename("McN_A"="Last Year YYYY",
                 "McN_A_doy"="Last Day Num"),by="Tag Code") %>% 
  
#Priest Rapids Adult
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("PRA")) %>% 
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>% 
            rename("PRa_A"="Last Year YYYY",
                   "PRa_A_doy"="Last Day Num"),by="Tag Code") %>% 
  
#Rock Island Adult
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("RIA")) %>% 
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>% 
            rename("RIs_A"="Last Year YYYY",
                   "RIs_A_doy"="Last Day Num"),by="Tag Code") %>% 
  
#Tumwater adult #TODO check if juveniles/ ghost tags detected
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("TUF")) %>%
            select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>%
            rename("Tum_A"="Last Year YYYY",
                   "Tum_A_doy"="Last Day Num"),by="Tag Code") %>% 

  #instream arrays above Tumwater day
  left_join(inst_array_det %>% 
              select(c(`Tag Code`,`Last Year YYYY`,`Last Day Num`)) %>%
              arrange(desc(`Last Year YYYY`),desc(`Last Day Num`)) %>% #arrange in order of increasing first year
              distinct(`Tag Code`,.keep_all=TRUE) %>% #remove duplicate detections at multiple "sites" at this dam 
              rename("instr_array"="Last Year YYYY",
                     "instr_array_A_doy"="Last Day Num"),by="Tag Code") %>%   
  
#move all the detection day of year columns to after all the detections year columns
  relocate( ends_with("doy"),.after=instr_array) %>% 
  #Add juvenile life history
  
  #add age
  mutate( age = case_when(mark_time==1~ "Unk",
                          (mark_file$`Mark Day Number`>179)~   "sub", #if DOY > 179 then subyearling
                          is.na(mark_file$`Length mm`)~       NA_character_,
                          (mark_file$`Length mm`>=cutoffs_and_props[[1]]$y[ifelse((mark_file$`Mark Day Number`-49)>0,
                                                                                  (mark_file$`Mark Day Number`-49),1)]) ~     "YCW",#assign age based on cutoff rule
                          TRUE~                                 "sub"
  )
  ) %>% 
  filter(!is.na(age)) %>% # drop 329 fish captured before doy 179 but without length data
  
  
  #get rid of fish whose tag ages are inconsistent with observed juvenile migration year
  mutate(seaward_year_obs= select(., LWe_J:Est_J) %>% reduce(pmin,na.rm=TRUE)) %>%  # add year of seaward migration
  mutate(mark_to_seward=`seaward_year_obs`-`Mark Year YYYY`) %>% # add years between tagging and seaward migration (for those fish detected migrating downstream as juveniles)
  filter(
    is.na(`mark_to_seward`)|(age=="sub"&`mark_to_seward`==1)|
      (age=="YCW"&`mark_to_seward`==0)|(age=="Unk"&`mark_to_seward`==0) 
  ) %>%#get rid of fish known to have migrated at an age not consistant with their assigned life history (14/7469 subs, 5/13630 yrlngs) 
  
  
  #assign life history (LH)
  mutate(LH= case_when(
    age=="Unk" ~ "Unk",
    age=="YCW" ~ "smolt",
    `Mark Day Number`<= 139 ~"fry",
    `Mark Day Number`<=262 ~"summer",
    TRUE ~"fall",
    
  )) %>% 
  
  
  # create stream column that is more readable
  mutate(stream = case_when(
    `Mark Site Info Code`=="CHIWAT" ~ "Chiwawa",
    `Mark Site Info Code`=="NASONC" ~ "Nason",
    `Mark Site Info Code`=="WHITER" ~ "White",
    TRUE ~ "LWE"
  )) %>% 
  #add predicted seaward migration year
  mutate(sea_Year_p = ifelse(age%in%c("YCW","Unk"),`Mark Year YYYY`,`Mark Year YYYY`+1)) %>% 
  
  #mutate detections to be years between seaward migration year and detection year 
  mutate_at(vars(LWe_J:instr_array), ~.x-sea_Year_p) %>%  
  #make juvenile detection "1" and non-detections "0"
  mutate_at(vars(LWe_J:Est_J), ~case_when(is.na(.x)~ 0,
                                          .x==0~1)) %>%
  # #make adult detections after 1 year "1", after 2 year "2" after three years "3" otherwise "0"
  mutate_at(vars(Bon_A:instr_array), ~case_when(is.na(.x)~ 0,
                                          .x==0~0,
                                          .x==1~1,
                                          .x==2~2,
                                          .x==3~3,
                                          TRUE ~0)) %>% 
  #if only adult detection is at tumwater dam presume "ghost tag" (fish died and some years later tag washed downstream over dam) and change detection to 0
  mutate(Tum_A=ifelse( select(., Bon_A:RIs_A)%>% rowSums() ==0& Tum_A !=0 ,0,Tum_A)) %>% 
  
  #zero out instream array adult detections unless they occurred after a valid ddownstream adult detection
  mutate(instr_array=ifelse(
    select(., Bon_A:RIs_A) %>%  reduce(pmax.int,na.rm=T)!=instr_array | #adult mainstem year is not equal to instream array year
      select(., Bon_A_doy:RIs_A_doy) %>%  reduce(pmax.int,na.rm=T)>instr_array_A_doy ,#adult mainstem day is greater than instream array day
    0, instr_array )) %>% 

  #make detection of LWe_J 1 if that is release location (for trap dependence)
  mutate(LWe_J=if_else(stream=="LWE",1,LWe_J)) %>%
  #subset columns of interest
  select(sea_Year_p,LH,stream,LWe_J:instr_array,`Length mm`,`Mark Day Number`,LWe_J_doy:instr_array_A_doy)


write.csv(mark_file_CH,file=here("Data","mark_file_CH.csv"))

