library(here)
library(tidyverse)


#Wild or natural origin spring Chinook captured by screw trap in the Chiwawa, Nason, or White
mark_file<-read_csv(here("Data","ptagis","Tagging detail.csv")) %>% 
  filter(`Rear Type Code`=="W",
         `Run Name`=="Spring",
         `Species Name`=="Chinook",
         `Capture Method Code`=="SCREWT"
         )

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
lower_wen_traps_recaps<-read_csv(here("Data","ptagis","ptagis_recap_data_wen_trib_screwt.csv")) %>% mutate(recap_date=lubridate::mdy_hms(recap_date),"First Year YYYY"=lubridate::year(recap_date)) %>% rename("Tag Code" = tag_id ) %>% arrange(recap_date)

#Lower and mid Columbia mainstem dam detections
lower_mid_col_dams<-read_csv(here("Data","ptagis","mid_lower_Col_Interrogation_Summary.csv"))

#upper Columbia mainstem Dam and Tumwater Dam detections
upper_col_dams_tum<-read_csv(here("Data","ptagis","upper_Col_Wen_Interrogation_summary.csv"))

#load length cutoffs for subyearling/yearling delineation for all days <=179
load(here("Data","cutoffs_and_props.Rdata"))

#add capture histories to mark file
mark_file_CH<-mark_file %>% 

#Lower Wenatchee Juvenile recapture
 left_join(lower_wen_traps_recaps %>%
              filter(`recap_site` %in% c("WENA4T","WENATT")) %>% #subset McNary juvenile detection sites
              select(c(`Tag Code`,`First Year YYYY`)) %>% #take columns first detection year and pit code
             # arrange(recap_date) %>% #arrange in order of increasing recap_date
             distinct(`Tag Code`,.keep_all=TRUE) %>% #remove duplicate detections at multiple "sites" at this trap 
              rename("LWe_J"="First Year YYYY"),by="Tag Code") %>%   #rename to represent site for merging

#McNary Juvenile detection
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` =="MCJ") %>% #subset McNary juvenile detection sites 
            select(c(`Tag Code`,`First Year YYYY`)) %>% #take columns first detection year and pit code 
            rename("McN_J"="First Year YYYY"),by="Tag Code") %>%   #rename to represent site for merging

#John Day Juvenle detection
  left_join(lower_mid_col_dams %>% 
              filter(`Site Code Value` =="JDJ") %>% 
              select(c(`Tag Code`,`First Year YYYY`)) %>%
              rename("JDD_J"="First Year YYYY"),by="Tag Code") %>%

#Bonneville Juvenile detection
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("B1J","B2J","BCC")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>% 
            arrange(`First Year YYYY`) %>% #arrange in order of increasing first year
            distinct(`Tag Code`,.keep_all=TRUE) %>% #remove duplicate detections at multiple "sites" at this dam 
            rename("Bon_J"="First Year YYYY"),by="Tag Code") %>%   

#estuary towed array
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` =="TWX") %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>%
            rename("Est_J"="First Year YYYY"),by="Tag Code") %>% 

#Bonneville Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("BO1", "BO2",  "BO3", "BO")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>% 
            arrange(desc(`First Year YYYY`)) %>% #latest year detected first, so when we take distinct below that is the one we keep
            distinct(`Tag Code`,.keep_all=TRUE) %>% 
            rename("Bon_A"="First Year YYYY"),by="Tag Code") %>% 

#The Dalles Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("TD1", "TD2")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>%
            arrange(desc(`First Year YYYY`)) %>%
            distinct(`Tag Code`,.keep_all=TRUE) %>%
            rename("TDD_A"="First Year YYYY"),by="Tag Code") %>%

#John Day Adult
  left_join(lower_mid_col_dams %>% 
            filter(`Site Code Value` %in% c("JO1", "JO2")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>%
            arrange(desc(`First Year YYYY`)) %>%
            distinct(`Tag Code`,.keep_all=TRUE) %>%
            rename("JDD_A"="First Year YYYY"),by="Tag Code") %>% 

#McNary Adult
left_join(lower_mid_col_dams %>% 
          filter(`Site Code Value` %in% c("MC1", "MC2")) %>% 
          select(c(`Tag Code`,`First Year YYYY`)) %>% 
          arrange(desc(`First Year YYYY`)) %>%
          distinct(`Tag Code`,.keep_all=TRUE) %>%
          rename("McN_A"="First Year YYYY"),by="Tag Code") %>% 
  
#Priest Rapids Adult
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("PRA")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>% 
            rename("PRa_A"="First Year YYYY"),by="Tag Code") %>% 
  
#Rock Island Adult
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("RIA")) %>% 
            select(c(`Tag Code`,`First Year YYYY`)) %>% 
            rename("RIs_A"="First Year YYYY"),by="Tag Code") %>% 
  
#Tumwater adult #TODO check if juveniles/ ghost tags detected
  left_join(upper_col_dams_tum %>% 
            filter(`Site Code Value` == c("TUF")) %>%
            select(c(`Tag Code`,`First Year YYYY`)) %>%
            rename("Tum_A"="First Year YYYY"),by="Tag Code") %>% 

#Add juvenile life history
  
#add age
  mutate( age = case_when(
    (mark_file$`Mark Day Number`>179)~   "sub", #if DOY > 179 then subyearling
   is.na(mark_file$`Length mm`)~       NA_character_,
   (mark_file$`Length mm`>=cutoffs_and_props[[1]]$y[(mark_file$`Mark Day Number`-49)])~     "YCW",#assign age based on cutoff rule
  TRUE~                                 "sub"
  )
  ) %>% 
  filter(!is.na(age)) %>% # drop 329 fish captured before doy 179 but without length data
  

#get rid of fish whose tag ages are consistent with observed juvenile migration year
 mutate(seaward_year_obs= select(., LWe_J:Est_J) %>% reduce(pmin,na.rm=TRUE)) %>%  # add year of seaward migration
mutate(mark_to_seward=`seaward_year_obs`-`Mark Year YYYY`) %>% # add years between tagging and seaward migration (for those fish detected migrating downstream as juveniles)
filter(
  is.na(`mark_to_seward`)|(age=="sub"&`mark_to_seward`==1)|
    (age=="YCW"&`mark_to_seward`==0) 
         ) %>%#get rid of fish known to have migrated at an age not consistant with their assigned life history (14/7469 subs, 5/13630 yrlngs) 
  
  
#assign life history (LH)
  mutate(LH= case_when(
    age=="YCW" ~ "smolt",
    `Mark Day Number`<= 139 ~"fry",
    `Mark Day Number`<=262 ~"summer",
    TRUE ~"fall",
    
  )) %>% 


# create stream column that is more readable
mutate(stream = case_when(
  `Mark Site Info Code`=="CHIWAT" ~ "Chiwawa",
  `Mark Site Info Code`=="NASONC" ~ "Nason",
  TRUE ~ "White"
)) %>% 
#add predicted seaward migration year
mutate(sea_Year_p = ifelse(age=="YCW",`Mark Year YYYY`,`Mark Year YYYY`+1)) %>% 

#mutate detections to be years between seaward migration year and detection year 
mutate_at(vars(LWe_J:Tum_A), ~.x-sea_Year_p) %>%  
   #make juvenile detection "A" and non-detections "0"
  mutate_at(vars(LWe_J:Est_J), ~case_when(is.na(.x)~ 0,
                                          .x==0~1)) %>%
  # #make adult detections after 1 year "1", after 2 year "2" after three years "3" otherwise "0"
  mutate_at(vars(Bon_A:Tum_A), ~case_when(is.na(.x)~ 0,
                                          .x==0~0,
                                          .x==1~1,
                                          .x==2~2,
                                          .x==3~3,
                                          TRUE ~0)) %>% 
  #drop rows with capture history where only adult detection is at tumwater dam (because of concern about "ghost tags)
  filter(!( select(., Bon_A:RIs_A)%>% rowSums() ==0& Tum_A !=0))

