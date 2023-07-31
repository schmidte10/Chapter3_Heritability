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

#--- set working directory ---# 
#setwd(C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter3_Heritability)

#--- fish data ---# 
FISH_ID = "example1"
mass = 0.0007676 #in kilograms 
chamber = "ch1" 
tempearture = 28.5 

#--- import data ---# 
firesting <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/2023/Resp_backup/2023_Resp/Dell7440/Experiment_ 01 July 2023 11 15AM/Oxygen data raw/firesting.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Time (HH:MM:SS)` = col_time(format = "%H:%M:%S")), 
                             trim_ws = TRUE, skip = 19)
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
  rename(Oxygen = chamber) 



