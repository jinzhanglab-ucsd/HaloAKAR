# Clean up
rm(list=ls())

# The easiest way to get tidyr is to install the whole tidyverse:
#install.packages("tidyverse")
#install.packages("gghighlight")
#install.packages("drc")
#nleinstall.packages("plyr") 

#install.packages("quantreg")
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)	
library(drc)
library(plyr)


## General input for the data
file.list <- list.files(pattern='*.xlsx')
nleng <- length(file.list)
m<-1
ratio_all<-data.frame()
ratio_all2<-data.frame()
ratio_time<-data.frame()
#define how many channels you measured
channels <- 2
#define the channel you are interested in
ch_int <-2
#define which time point you are interested in, if its the beginning define here if it is the endpoint define in the for loop 
#time_p <-1

# read all files after each other, calculate the ratio and append data to a file
for (i in 1:nleng) {
	#reads in the file
	filename <-file.list[m]
	n <- read_excel(file.list[m])
	#finds how many cell there are in total 
	l_tot<- nrow(n)
	#if endpoint then time_p = l_tot
	time_p <-11
	time_p2 <-23
	w_tot<- ncol(n)
	#define how many cells were measured
	n_cell<-(w_tot-1)/channels
	#define colum range in file that is relevant
	start<-n_cell*(ch_int-1)+2
	end <-n_cell*ch_int+1
	#Gets the name of the file
	index <- as.data.frame(filename, stringsAsFactors=FALSE)
	#gets the specific time point, 
	n_time <-n[time_p,]
	n_time2 <-n[time_p2,]
	n_time_cells <-n_time[,start:end]
	n_time_cells2 <-n_time2[,start:end]
		
	#binds the two values together and appends to the previous values
		values<-tibble(index,n_time_cells)
	values2 <-tibble(index, n_time_cells2)
	#ratios <-values2/values
	ratio_all<- rbind.fill(ratio_all, values)
	ratio_all2<-rbind.fill(ratio_all2, values2)
	#ratio_time <- rbind.fill(ratio_time, ratios)
	m<-m+1
}

ratio_t <- t(ratio_all)
ratio_t2 <-t(ratio_all2)
#ratio_time_t <-t(ratio_time)

# Save data
#choose the name here according to the type of analysis
	out <- "Analysis_3min"
	out2 <- "Analysis_end"
	fullpath = getwd()
	directoryname = basename(fullpath)
	ext <- ".csv"
	filename7 <- paste(directoryname,out,ext)
	filename8 <- paste(directoryname,out2,ext)
	write.table(ratio_t, file =filename7, sep = ",")
	write.table(ratio_t2, file =filename8, sep = ",")
	


