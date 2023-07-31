###############################################################################
###                                                                         ###
###                 Resp R script - Hertiability                            ###
###                        RespRv2.3.1                                      ### 
###                                                                         ### 
############################################################################### 

#--- download necessary libraries ---# 
library("respR")
library("tidyverse") 
library("data.table") 
library("hms")

#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter3_Heritability")

#--- fish data ---# 
FISH_ID = "133.1"
mass = 0.0007676 #in kilograms 
chamber = "ch4" 
temperature = 28.5 
salinity = 35 
chamber_vol = 0.04860
system1 = "DELL" 
Date1 = "01/07/2023"

#--- import data ---# 
firesting <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 56AM/Oxygen data raw/firesting.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Time (HH:MM:SS)` = col_time(format = "%H:%M:%S")), 
                             trim_ws = TRUE, skip = 19)

#--- data manipulation ---#
firesting2 <- firesting %>% 
  select(c(1:3,5:8)) %>% 
  rename(TIME = `Time (HH:MM:SS)`, 
         dTIME = `Time (s)`, 
         DATE = Date, 
         ch1 = Ch1...5, 
         ch2 = Ch2...6,
         ch3 = Ch3...7, 
         ch4 = Ch4...8) %>% 
  select(c("dTIME", all_of(chamber),"TIME")) %>% 
  rename(Oxygen = chamber) %>% 
  mutate(dTIME = as.numeric(dTIME), 
         Oxygen = as.numeric(Oxygen), 
         TIME = as.character(TIME))

#--- inspect file ---#
ap <- inspect(firesting2)

#--- import cycle data ---# 
cycle1 <- read.csv("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 56AM/All slopes/Cycle_1.txt", 
                          sep = ";")

#- find out when the cycle started -#
cycle1.start = cycle1[1,1] 
#- find out when the cycle ended - #
cycle1.end = tail(cycle1, n=1)[1,1]

#- find out which row in the firesting file has the same time as the cycle start time -# 
cycle1.start.row <- which(firesting2$TIME == cycle1.start); cycle1.start.row
#- find out which row in the firesting file has the same time as the cycle end time -#
cycle1.end.row <- which(firesting2$TIME == cycle1.end); cycle1.end.row

#- subset the first cycle from the firesting file -#
cycle1_data <- firesting2 %>% 
  subset_data(from = cycle1.start.row, 
              to = cycle1.end.row, 
              by = "row") %>% 
  inspect() 

#--- end experiment cycle period (all measurements in seconds) ---# 
buffer = 0 
measure = 300
flush = 195 
wait = 45 

cycle.length = measure + flush + wait
reps <- seq(cycle1.start.row, cycle1.end.row, cycle.length) 

starts = reps + buffer 
ends = reps + buffer + measure


mmr <- list()

mmr <- subset_data(cycle1_data, from = starts, to = ends, by = "TIME") %>% 
  auto_rate(method = "rolling", plot = T, width = 30, by = "time")  

#--- background rates ---# 
pre1 <- read.csv("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 15AM/All slopes/Cycle_2.txt", 
                   sep = ";")
pre2 <- read.csv("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 15AM/All slopes/Cycle_3.txt", 
                 sep = ";")
pre3 <- read.csv("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 15AM/All slopes/Cycle_4.txt", 
                 sep = ";")

firesting_pre <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 15AM/Oxygen data raw/firesting.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, skip = 19) 

firesting_pre2 <-  firesting %>% 
  select(c(1:3,5:8)) %>% 
  rename(TIME = `Time (HH:MM:SS)`, 
         dTIME = `Time (s)`, 
         DATE = Date, 
         ch1 = Ch1...5, 
         ch2 = Ch2...6,
         ch3 = Ch3...7, 
         ch4 = Ch4...8) %>% 
  select(c("dTIME", all_of(chamber),"TIME")) %>% 
  rename(Oxygen = chamber) %>% 
  mutate(dTIME = as.numeric(dTIME), 
         Oxygen = as.numeric(Oxygen), 
         TIME = as.character(TIME)) %>% 
  inspect()

cycle_list <- list(pre1,pre2,pre3)
bg <- list()

for (i in 1:3) {
  bg[[i]] <- cycle_list[[i]] %>% 
    select(c(1:2,4:7))  %>% 
    rename(TIME = Time, 
           dTIME = Seconds.from.start.for.linreg,  
           ch1 = ch1.po2, 
           ch2 = ch2.po2,
           ch3 = ch3.po2, 
           ch4 = ch4.po2) %>% 
    select(c("dTIME", all_of(chamber),"TIME")) %>% 
    rename(Oxygen = chamber) %>% 
    mutate(dTIME = as.numeric(dTIME), 
           Oxygen = as.numeric(Oxygen), 
           TIME = as.character(TIME)) %>% 
    inspect() %>% 
    calc_rate.bg()
}


pre_final[[1]]$rate.bg;pre_final[[2]]$rate.bg;pre_final[[3]]$rate.bg

pre1_final <- pre1 %>% 
  select(c(1:2,4:7)) %>% 
  rename(TIME = Time, 
         dTIME = Seconds.from.start.for.linreg,  
         ch1 = ch1.po2, 
         ch2 = ch2.po2,
         ch3 = ch3.po2, 
         ch4 = ch4.po2) %>% 
  select(c("dTIME", all_of(chamber),"TIME")) %>% 
  rename(Oxygen = chamber) %>% 
  mutate(dTIME = as.numeric(dTIME), 
         Oxygen = as.numeric(Oxygen), 
         TIME = as.character(TIME)) %>% 
  inspect() %>% 
  calc_rate.bg()


bg_pre <- inspect(pre1, dtime = 1, oxygen=2) |> 
  calc_rate.bg() %>% 
  summary()


mmr_adj <- adjust_rate(mmr, by=c(bg[[1]]$rate.bg,bg[[2]]$rate.bg,bg[[3]]$rate.bg), method = "mean")
mmr_adj

mmr_adj2 <- mmr_adj %>% 
  convert_rate(oxy.unit = "%Air", 
               time.unit = "secs", 
               output.unit = "mg/h/kg", 
               volume = chamber_vol,
               mass = mass,
               S = salinity, 
               t = temperature, 
               P = 1.013253)

summary(mmr_adj2)

mmr_final <- mmr_adj2 %>% 
  select_rate(method = "rsq", n = c(0.95,1)) %>%
  select_rate(method = "highest", n = 1) %>% 
  plot(type = "full") %>%
  summary(export = TRUE)

FISH_ID;mmr_final$rate.output;mmr_final$rate.output*mass

`mmr(mg/h/kg)` <- mmr_final$rate.output *-1
mmr <- mmr_final$rate.output*mass * -1

results <- data.frame(FISH_ID = FISH_ID, 
                      Mass = mass, 
                      Chamber = chamber, 
                      System = system1,
                      Volume = chamber_vol, 
                      Date = Date1,
                      Salinity = salinity, 
                      Temperature = temperature, 
                      `mmr(mg/h/kg)` = `mmr(mg/h/kg)`, 
                      mmr = mmr) %>% 
  view()
