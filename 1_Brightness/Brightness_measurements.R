# Clean up
rm(list=ls())

# The easiest way to get tidyr is to install the whole tidyverse:
#install.packages("tidyverse")
#install.packages("gghighlight")
#install.packages("drc")
#install.packages("plyr") 
#install.packages("quantreg")
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)	
library(drc)
library(plyr)


## General input for the data
list_of_files <- list.files(path = "~your folder path",
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = TRUE)
#Remove the first output file (remove this line if you run the code the first time)
list_of_files <-list_of_files[-31]
df <- list_of_files %>%
  set_names() %>% 
  map_df(read_csv, .id = "file_name")

#change column name
df <-dplyr::rename(df, index = ...1)

#only look at the mean values
df2 <- df%>% dplyr::select(file_name, index, starts_with("Mean"))
ncol <-ncol(df2)
#Change layout of the data and change names around. 
df3 <- pivot_longer(df2,cols=3:ncol, names_to ="number", values_to = "Mean", values_drop_na = TRUE)
df4 <- pivot_wider (df3, names_from = index, values_from = Mean)
df4 <-dplyr::rename(df4, farred = 3)
df4 <-dplyr::rename(df4, green = 4)

#make two individual data frames for Green and farred
df_green <-df4%>% dplyr::select(file_name, number, green)
df_farred <-df4%>% dplyr::select(file_name, number, farred)

#distribute by the experiment so I can perform background correction with first value
df_green2 <- pivot_wider (df_green, names_from = file_name, values_from = green)
df_farred2 <- pivot_wider (df_farred, names_from = file_name, values_from = farred)
df_mean_g <- df_green2%>% dplyr::select(number)
df_mean_fr <- df_farred2%>% dplyr::select(number)


#background correction
bg_green <-filter(df_green2, number == "Mean1")
bg_farred <-filter(df_farred2, number == "Mean1")

df_green3 <- {df_green2[,-1] - bg_green[rep(1, nrow(df_green2)), -1]}
df_farred3 <- {df_farred2[,-1] - bg_farred[rep(1, nrow(df_farred2)), -1]}

#bind the tags together again
df_green4 <- bind_cols(df_mean_g,df_green3)
df_farred4 <- bind_cols(df_mean_fr,df_farred3)
ncol2 <-ncol(df_green4)

#make longer to join green and farred
df_green5 <- pivot_longer(df_green4,cols=2:ncol2, names_to ="file_name", values_to = "green", values_drop_na = TRUE)
df_farred5 <- pivot_longer(df_farred4,cols=2:ncol2, names_to ="file_name", values_to = "farred", values_drop_na = TRUE)
df_both <-full_join(df_green5,df_farred5)
df_both2 <-arrange(df_both, file_name)

#remove the background one
df_both3 <-filter(df_both2, number != "Mean1")


df5 <-mutate(df_both3, ratio = farred/green)
df6 <-separate(df5, file_name, sep= "HeLA_", into = c("base", "end"))
df7 <-separate(df6, end, sep= "/", into = c("sample", "end"))

#select the data I want 
df8 <-df7%>% dplyr::select(sample, number, end, ratio)

#combine all from 1 replicate
df9 <- pivot_wider (df8, names_from = sample, values_from = ratio)

#choose the name here according to the type of analysis
out <- "Analysis_BasalBrightness"
fullpath = getwd()
directoryname = basename(fullpath)
ext <- ".csv"
filename7 <- paste(directoryname,out,ext)
write.table(df9, file =filename7, sep = ",")

