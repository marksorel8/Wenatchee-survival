
get_dat<-function(){
 
#function to load and download packages in necessary
pkgTest <- function(x){
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function 
pkgTest("dataRetrieval")
require("dataRetrieval")
pkgTest("here")
require("here")
#function to load Chiwawa discharge data
discharge_func<-function(){
  
  
  
  Wen.Peshastin.flowDV<- readNWISdv(siteNumber="12459000",
                           "00060", "1990-01-01", "2018-12-31",statCd="00003")
  
  #munging
  colnames(Wen.Peshastin.flowDV)[4]<-"y"
  Wen.Peshastin.flowDV$date<-as.Date(Wen.Peshastin.flowDV$Date,format="%Y-%m-%d")
  Wen.Peshastin.flowDV$Year<-as.integer(format(Wen.Peshastin.flowDV$date,"%Y"))
  Wen.Peshastin.flowDV$week<-as.integer(format(Wen.Peshastin.flowDV$date,"%W"))
  Wen.Peshastin.flowDV$month<-as.integer(format(Wen.Peshastin.flowDV$date,"%m"))
  
  Wen.Peshastin.flowDV$day<-as.integer(format(Wen.Peshastin.flowDV$date,"%d"))
  Wen.Peshastin.flowDV$doy<-as.integer(format(Wen.Peshastin.flowDV$date,"%j"))
  return(Wen.Peshastin.flowDV)
}


    
## air temperature data. 
###I hate that I can;t figure out a better way to download these data, but for now I'm going to the NWS website (https://w2.weather.gov/climate/xmacis.php?wfo=otx) and using the "NOWData" tool to retrieve a table of monthly summarized average temperatures at Wenatachee Pangborn airpoirt from 2000-2021. I'm then copying the table into excel and saving it as a .csv (I know...uncool). 

Wen_air<-read.csv(here("Data","env cov","wen prang monthly air temp.csv")) %>% pivot_longer(cols=2:13,names_to="month",names_prefix = "X",values_to="temp") %>% select(-2) %>% rename(Year=year,y=temp)
 
# Chiw_temp<-readxl::read_xlsx(here("Data","env cov","Temp Request 3-9-21.xlsx"))
# dev.new()
# ggplot(Chiw_temp,aes(x=Date,y=`AvgOfTemp, °C`))+geom_line()


#dowload discharage data
Dis<-discharge_func()

# *download mainstem discharge and temperature data  *
library(vroom)
#temp(WGM), discharge, spill at McNary, and Rock Island
path<-"http://www.cbr.washington.edu/dart/cs/php/rpt/mg.php?sc=1&mgconfig=river&outputFormat=csvSingle&year%5B%5D=2020&year%5B%5D=2019&year%5B%5D=2018&year%5B%5D=2017&year%5B%5D=2016&year%5B%5D=2015&year%5B%5D=2014&year%5B%5D=2013&year%5B%5D=2012&year%5B%5D=2011&year%5B%5D=2010&year%5B%5D=2009&year%5B%5D=2008&year%5B%5D=2007&year%5B%5D=2006&loc%5B%5D=BON&loc%5B%5D=MCN&loc%5B%5D=RIS&data%5B%5D=Outflow&data%5B%5D=Temp+%28WQM%29&startdate=1%2F1&enddate=12%2F31&avgyear=0&consolidate=1&grid=1&y1min=0&y1max=&y2min=&y2max=&size=medium"

mainstem_flow_temp<-vroom(path) %>% mutate(date2=lubridate::mdy(paste0(`mm-dd`,"-",year)),
                                           month=lubridate::month(date2),doy=lubridate::yday(date2)) %>% rename(y=value)  


path<-"http://www.cbr.washington.edu/dart/cs/php/rpt/mg.php?sc=1&mgconfig=river&outputFormat=csvSingle&year%5B%5D=2021&year%5B%5D=2020&year%5B%5D=2019&year%5B%5D=2018&year%5B%5D=2017&year%5B%5D=2016&year%5B%5D=2015&year%5B%5D=2014&year%5B%5D=2013&year%5B%5D=2012&year%5B%5D=2011&year%5B%5D=2010&year%5B%5D=2009&year%5B%5D=2008&year%5B%5D=2007&year%5B%5D=2006&year%5B%5D=2005&loc%5B%5D=PRD&loc%5B%5D=RIS&loc%5B%5D=WAN&data%5B%5D=Spill&data%5B%5D=Spill+Percent&startdate=1%2F1&enddate=12%2F31&avgyear=0&consolidate=1&grid=1&y1min=0&y1max=&y2min=&y2max=&size=medium"

UC_spill <-vroom(path) %>% mutate(date2=lubridate::mdy(paste0(`mm-dd`,"-",year)),
                                          month=lubridate::month(date2),doy=lubridate::yday(date2)) %>% rename(y=value) 


path<-"http://www.cbr.washington.edu/dart/cs/php/rpt/mg.php?sc=1&mgconfig=river&outputFormat=csvSingle&year%5B%5D=2021&year%5B%5D=2020&year%5B%5D=2019&year%5B%5D=2018&year%5B%5D=2017&year%5B%5D=2016&year%5B%5D=2015&year%5B%5D=2014&year%5B%5D=2013&year%5B%5D=2012&year%5B%5D=2011&year%5B%5D=2010&year%5B%5D=2009&year%5B%5D=2008&year%5B%5D=2007&year%5B%5D=2006&year%5B%5D=2005&loc%5B%5D=BON&loc%5B%5D=JDA&loc%5B%5D=MCN&loc%5B%5D=TDA&data%5B%5D=Spill&data%5B%5D=Spill+Percent&startdate=1%2F1&enddate=12%2F31&avgyear=0&consolidate=1&grid=1&y1min=0&y1max=&y2min=&y2max=&size=medium"

MC_spill <-vroom(path) %>% mutate(date2=lubridate::mdy(paste0(`mm-dd`,"-",year)),
                                  month=lubridate::month(date2),doy=lubridate::yday(date2)) %>% rename(y=value) 

#functions to calculat winter maximum flow and summer minimum flow, code copied from M. Sheuerell. https://github.com/mdscheuerell/skagit_sthd/blob/master/analysis/App_1_Retrieve_covariates.pdf

#winter maximum (in year when fish would be migrating to the ocean (mig_year))
winter_mean_func<-function(dat,resp_name){
  ## autumn flows in year t
  aut <- subset(dat, (month>=10 & month<=12))
  ## spring flows in year t+1
  spr <- subset(dat,
                     (month>=1 & month<=2))
  ## change spr year index to match aut
  aut[,"Year"] <- as.integer(aut[,"Year"]) + 1
  ## combine flows indexed to winter start year & calculate max flow
  dat_wtr <- aggregate(y ~ Year, data = rbind(aut,spr), mean) %>% rename(mig_year=Year)
  dat_wtr[,"y"] <- round(dat_wtr[,"y"], 1) 
  
  return(dat_wtr) 
}

#winter variability 
winter_CV_func<-function(dat,resp_name){
  ## autumn flows in year t
  aut <- subset(dat, (month>=10 & month<=12))
  ## spring flows in year t+1
  spr <- subset(dat,
                (month>=1 & month<=2))
  ## change spr year index to match aut
  aut[,"Year"] <- as.integer(aut[,"Year"]) + 1
  ## combine flows indexed to winter start year & calculate max flow
  dat_wtr <- aggregate(y ~ Year, data = rbind(aut,spr), function(x){sd(x)/mean(x)}) %>% rename(mig_year=Year)
  dat_wtr[,"y"] <- round(dat_wtr[,"y"], 1) 
  
  return(dat_wtr) 
}



#summer mean 
summer_mean_func<-function(dat){
  ## autumn flows in year t
  summer <- subset(dat, (month>=6 & month<=9))
  ## calculate min flow
  dat_sum <- aggregate(y ~ Year, data = summer, mean) %>% 
    mutate(mig_year=Year+1)

  return(dat_sum) 
}


#spring mean 
spring_mean_func<-function(dat,first_month,last_month){
  ## autumn flows in year t
  summer <- subset(dat, (month>=first_month & month<=last_month))
  ## calculate min flow
  dat_sum <- aggregate(y ~ year, data = summer, mean) %>% 
    rename(mig_year=year)
  
  return(dat_sum) 
}


env_dat_out<-
#mainstem wenatchee flow and air temp
summer_mean_func(as.data.frame(Dis)) %>% rename(sum_flow=y) %>% 
full_join(winter_mean_func(as.data.frame(Dis)) %>% rename(win_flow=y)) %>% 
full_join(winter_mean_func(as.data.frame(Dis)) %>% rename(win_flow_CV=y)) %>% 
full_join(summer_mean_func(as.data.frame(Wen_air)) %>% rename(sum_air=y)) %>% 
full_join(winter_mean_func(as.data.frame(Wen_air)) %>% rename(win_air=y)) %>% 
full_join(spring_mean_func(as.data.frame(Wen_air)%>%rename(year=Year),3,4) %>% rename(spr_air=y)) %>% 
full_join(spring_mean_func(as.data.frame(Dis)%>%rename(year=Year),3,4) %>% rename(spr_dis=y)) %>%


#mainstem Columbia flow and water temp
##Rock Island
####juvenile
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="RIS",parameter=="outflow"),first_month = 4,last_month = 5) %>% rename(RIS_flow_juv=y)) %>% 
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="RIS",parameter=="tempc"),first_month = 4,last_month = 5) %>% rename(RIS_temp_juv=y)) %>% 
## upper columbia spill %
  full_join(spring_mean_func(UC_spill %>% filter(location%in%c("PRD","RIS","WAN"),parameter=="spillpct") %>% select(-location) ,first_month = 4,last_month = 5) %>% rename(UC_spill_pct_juv=y)) %>% 
#### adult
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="RIS",parameter=="outflow"),first_month = 5,last_month = 6) %>% rename(RIS_flow_ad=y)) %>% 
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="RIS",parameter=="tempc"),first_month = 5,last_month = 6) %>% rename(RIS_temp_ad=y)) %>% 
##McNary
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="MCN",parameter=="outflow"),first_month = 5,last_month = 6) %>% rename(McN_flow=y)) %>% 
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="MCN",parameter=="tempc"),first_month = 5,last_month = 6) %>% rename(McN_temp=y)) %>% 
  full_join(spring_mean_func(MC_spill %>% filter(location=="MCN",parameter=="spillpct"),first_month = 5,last_month = 6) %>% rename(McN_spill=y)) %>% 
  ##Bonneville
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="BON",parameter=="outflow"),first_month = 5,last_month = 6) %>% rename(Bon_flow=y)) %>% 
  full_join(spring_mean_func(mainstem_flow_temp %>% filter(location=="BON",parameter=="tempc"),first_month = 5,last_month = 6) %>% rename(Bon_temp=y)) %>% 
  full_join(spring_mean_func(MC_spill %>% filter(location=="BON",parameter=="spillpct"),first_month = 5,last_month = 6) %>% rename(Bon_spill=y)) %>% 
  ## mid columbia spill %
  full_join(spring_mean_func(MC_spill %>% filter(location%in%c("JDA","TDA","MCN" ), parameter=="spillpct") %>% select(-location) ,first_month = 4,last_month = 5) %>% rename(MC_spill_pct_juv=y)) %>% 
  filter(mig_year>=2006&mig_year<=2020) 
  
# Covariates that I recieved from Brian Burke for SAR model
envDir<-here("Data","LCM environmental data")

envFiles<-list.files(path=envDir,pattern="*.csv")

envdata<-tibble(year=2006:2020)

for ( i in envFiles){y<-read.csv(paste(envDir,i,sep="/"))
envdata<-envdata %>% left_join(y)}

env_dat_out<-full_join(env_dat_out,envdata %>% rename(mig_year=year))

return(env_dat_out)
}





