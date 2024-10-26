# Clean up the work space 
rm(list=ls())

#install libraries 
library(tidyverse)

#Read in the 5 files prepared by the Bash script (Input library and the four sorted libraries)
inputL <-as_tibble(read.table("PKA-left_G_S99_L002_R1.CPM", header=FALSE, sep= "\t")) 
inputL <- rename(inputL,input=V2)

notreat_H <-as_tibble(read.table("PKA-left-notreat-High_S102_L002_R1.CPM", header=FALSE, sep= "\t"))
notreat_H <- rename(notreat_H,ntH=V2)

notreat_M <-as_tibble(read.table("PKA-left-notreat-Med_S103_L002_R1.CPM", header=FALSE, sep= "\t"))
notreat_M <- rename(notreat_M,ntM=V2)

treat_H <-as_tibble(read.table("PKA-left-treat-High_S100_L002_R1.CPM", header=FALSE, sep= "\t"))
treat_H  <- rename(treat_H,tH=V2)

treat_M <-as_tibble(read.table("PKA-left-treat-Med_S101_L002_R1.CPM", header=FALSE, sep= "\t"))
treat_M <- rename(treat_M,tM=V2)

#combine all of them together
data_1 <- full_join(inputL,notreat_H, by = "V1")
data_1 <- full_join(data_1,notreat_M, by = "V1")
data_1 <- full_join(data_1,treat_H, by = "V1")
data_1 <- full_join(data_1,treat_M, by = "V1")

#calculate the Enrichment
data_2 <- mutate(data_1, ntH / input)
data_2 <- mutate(data_2, ntM / input)
data_2 <- mutate(data_2, tH / input)
data_2 <- mutate(data_2, tM / input)

# remove the STOP codons
noSTOP <- filter(data_2, !grepl('Stop', V1)) 

#calculate the conditions and the enrichment score with the frequencies with the STOPs in there 
data_3 <-mutate(noSTOP, cond1=ifelse(ntH/input<1, 1, 0))
data_3 <-mutate(data_3, cond2=ifelse(ntM/input>1, 1, 0))
data_3 <-mutate(data_3, cond3=ifelse(tH/input>1, 1, 0))
data_3 <-mutate(data_3, cond4=ifelse(tM/input<1, 1, 0))
data_3 <-mutate(data_3, condtot=cond1+cond2+cond3+cond4)
data_3 <-mutate(data_3, score = tH/input*ntM/input)

#Parental FAR = Phe|Ala|Arg
value_row <- filter(data_3, V1 == "Phe|Ala|Arg")$score
value_rowf <- filter(data_3, V1 == "Phe|Ala|Arg")

#Data that fulfills all 4 criteria
data_four <- filter(data_3, condtot==4) 
data_four <-add_row(data_four,value_rowf)
write.csv(data_four, "R1_L_Data_four.csv")

#Data that fulfills all 4 criteria and has a high score (higher than the parental sequence)
data_four_high <- filter(data_3, condtot==4 & score > value_row)
data_four_high <-add_row(data_four_high,value_rowf)
write.csv(data_four_high, "R1_L_Data_four_and_high.csv")

#Data that fulfills all 4 criteria or has a high score
data_four_or_high <- filter(data_3, condtot==4 | score > value_row)
data_four_or_high <-add_row(data_four_or_high,value_rowf)
write.csv(data_four_or_high, "R1_L_Data_four_or_high.csv")

#Data that has a high score
data_high <- filter(data_3, score > value_row )
data_high <-add_row(data_high,value_rowf)
write.csv(data_high, "R1_L_Data_high.csv")


#### calculate everything without the STOP for comparison
# calculate the new normalized fractions without the STOPS 
data_1b <- mutate(noSTOP, input_NnoSTOP = input/sum(input)*1000000)
data_1b <- mutate(data_1b, ntH_NnoSTOP = ntH/sum(ntH)*1000000)
data_1b <- mutate(data_1b, ntM_NnoSTOP = ntM/sum(ntM)*1000000)
data_1b <- mutate(data_1b, tH_NnoSTOP = tH/sum(tH)*1000000)
data_1b <- mutate(data_1b, tM_NnoSTOP = tM/sum(tM)*1000000)

#calculate the enrichment without the STOPS 
data_2b <- mutate(data_1b, ntH_NnoSTOP / input_NnoSTOP)
data_2b <- mutate(data_2b, ntM_NnoSTOP / input_NnoSTOP)
data_2b <- mutate(data_2b, tH_NnoSTOP / input_NnoSTOP)
data_2b <- mutate(data_2b, tM_NnoSTOP / input_NnoSTOP)

# calculate the new normalized fractions
data_3b <-mutate(data_2b, cond1=ifelse(ntH_NnoSTOP/input_NnoSTOP<1, 1, 0))
data_3b <-mutate(data_3b, cond2=ifelse(ntM_NnoSTOP/input_NnoSTOP>1, 1, 0))
data_3b <-mutate(data_3b, cond3=ifelse(tH_NnoSTOP/input_NnoSTOP>1, 1, 0))
data_3b <-mutate(data_3b, cond4=ifelse(tM_NnoSTOP/input_NnoSTOP<1, 1, 0))
data_3b <-mutate(data_3b, condtot2=cond1+cond2+cond3+cond4)
data_3b <-mutate(data_3b, score2 = tH_NnoSTOP/input_NnoSTOP*ntM_NnoSTOP/input_NnoSTOP)

#Parental FAR = Phe|Ala|Arg
value_row_b <- filter(data_3b, V1 == "Phe|Ala|Arg")$score2
value_rowf_b <- filter(data_3b, V1 == "Phe|Ala|Arg")

#Data that fulfills all 4 criteria
data_four_b <- filter(data_3b, condtot2==4) 
data_four_b <-add_row(data_four_b,value_rowf_b)
write.csv(data_four_b, "R1_L_Data_four_noSTOP.csv")

#Data that fulfills all 4 criteria and has a high score
#based on this data set cutt ofs were defined
data_four_high_b <- filter(data_3b, condtot2==4 & score2 > value_row_b)
data_four_high_b <-add_row(data_four_high_b,value_rowf_b)
write.csv(data_four_high_b, "R1_L_Data_four_and_high_noSTOP.csv")

#Data that fulfills all 4 criteria or has a high score
data_four_or_high_b <- filter(data_3b, condtot2==4 | score2 > value_row_b)
data_four_or_high_b <-add_row(data_four_or_high_b,value_rowf_b)
write.csv(data_four_or_high_b, "R1_L_Data_four_or_high_noSTOP.csv")

#Data that has a high score
data_high_b <- filter(data_3b, score2 > value_row_b)
data_high_b <-add_row(data_high_b,value_rowf_b)
write.csv(data_high_b, "R1_L_Data_high_noSTOP.csv")

#Data combined
data_four_all <-full_join(data_four, data_four_b)
write.csv(data_four_all, "R1_L_Data_high_both.csv")