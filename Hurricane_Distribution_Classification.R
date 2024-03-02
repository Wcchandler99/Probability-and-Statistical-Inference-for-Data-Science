#Prepared by:
# Will Chandler, wcc44
# Vignesh Krishna, vk505
# Reuben S Varghese, rsv39

install.packages("ggplot2")
install.packages("fitdistrplus")
install.packages("tidyverse")
library("ggplot2")
library("fitdistrplus")
library("tidyverse")
library(dplyr)

# Download Data
column_names <- c("Year", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

df <- read.csv2("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/climate_data.txt", sep = " ", col.names = column_names)

cat_4_1 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_4_data_1.csv", header = TRUE)
cat_4_2 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_4_data_2.csv", header = TRUE)
cat_4_3 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_4_data_3.csv", header = TRUE)
cat_4_4 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_4_data_4.csv", header = TRUE)
cat_4_5 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_4_data_5.csv", header = TRUE)
cat_5 <- read.csv("C:/Users/wccha/Documents/Rutgers/Prob_&_Stat/R/cat_5_data.csv", header = TRUE)

# Data Munging Wikipedia Data
for (i in (1:length(cat_4_1$Duration))){
  if ("," %in% strsplit(cat_4_1$Duration[i], split = "")[[1]]){
    cat_4_1$Year[i] <- strsplit(cat_4_1$Duration[i], ", ")[[1]][2]
  }
  else{
    cat_4_1$Year[i] <- strsplit(cat_4_1$Duration[[i]], " ")[[1]][2]
  }
}

for (i in (1:length(cat_4_1$Year))){
  if (length(strsplit(cat_4_1$Year[[i]], split = " ")[[1]]) > 1){
    cat_4_1$Year[i] <- strsplit(cat_4_1$Year[[i]], split = " ")[[1]][2]
  }
}

cat_4_2$Year <- cat_4_2$Season
cat_4_3$Year <- cat_4_3$Season
cat_4_4$Year <- cat_4_4$Season

for (i in (1:length(cat_5$Dates.at.Category.5.intensity))){
  if ("," %in% strsplit(cat_5$Dates.at.Category.5.intensity[i], split = "")[[1]]){
    cat_5$Year[i] <- strsplit(cat_5$Dates.at.Category.5.intensity[i], ", ")[[1]][2]
  }
  else{
    cat_5$Year[i] <- strsplit(cat_5$Dates.at.Category.5.intensity[[i]], " ")[[1]][2]
  }
}

decade_names <- c("1850", "1860", "1870", "1880", "1890", "1900", "1910", "1920", "1930", "1940", "1950", "1960", "1970", "1980", "1990", "2000", "2010", "2020")

decade_counts <- c(2, 1, 1, 4, 4, 2, 6, (6+2), (10+6), (8+1), (15+2), (10+5), (6+3), (7+3), (12+2), (15+8), (12+6), (10+2))

decades_count_df <- data.frame(decade_names, decade_counts)

# Data Munging AMO Data
df <- data.frame(df)
df <- data.frame(lapply(df, as.numeric))
for (i in (1:length(df$Year))){
  for (j in (0:17)){
    if (df$Year[i] >= (1850 + (j)*10)){
      if (df$Year[i] < (1860 + (j)*10)){
        df$Decade[i] = 1850 + (j)*10
      }
    }
  }
}

df$totals <- rowSums(df[2:13])
decade_totals <- df %>% group_by(Decade) %>% summarise(total = sum(totals))
decade_totals$temp_change_per_decade <- decade_totals$total/120

decade_totals$decade_counts <- decade_counts[1:16]

# Plots
qqplot(qpois(c(1:18)/19, mean(decade_counts)), decade_counts, 
       xlab = "Expected Poisson Values", 
       ylab = "Actual Hurricane Counts", 
       main = "QQplot Poisson") 

qqplot(qnorm(c(1:18)/19, mean(decade_counts)), decade_counts, 
       xlab = "Expected Normal Values", 
       ylab = "Actual Hurricane Counts", 
       main = "QQplot Normal")

ggplot(decade_totals[1:14,]) +
  geom_line(aes(x = Decade, y = total, color = "Temperature Change")) + 
  geom_line(aes(x = Decade, y = decade_counts, color = "Number of Hurricanes")) +
  ggtitle("Temperature Change per Decade and Number of Hurricanes") +
  ylab("Number of Hurricanes") +
  theme_gray(base_size = 24)

pairs(decade_totals[(1:14),c(1:4)])

ggplot(decade_totals) + 
  geom_col(aes(x = Decade, y = decade_counts)) +
  ylab("Hurricane Counts") + 
  ggtitle("Hurricane Count Bar Chart") +
  theme_gray(base_size = 24)

View(decade_totals)
help(geom_col)
