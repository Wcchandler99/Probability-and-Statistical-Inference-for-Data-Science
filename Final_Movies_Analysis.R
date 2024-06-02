install.packages("mclust")
library("mclust")

install.packages("tidyverse") 
library(tidyverse)

data <- read.csv("C:\\Users\\wccha\\Documents\\Rutgers\\Prob_&_Stat\\Final\\archive\\movies_genres.csv")

# Combining Genres
data$T_H_C_M <- rowSums(data[, c("Thriller", "Crime", "Mystery", "Horror")])
data$A_A_F_S <- rowSums(data[, c("Action", "Adventure", "Fantasy", "Science.Fiction")])

data <- data[, !(names(data) %in% c("Thriller", "Crime", "Mystery", "Action", "Horror", "Adventure", "Fantasy", "Science.Fiction"))]

clusters = 5

# EM Algorithm
result <- Mclust(data, G = clusters)  # G is the number of clusters

cat("Cluster assignments:", result$classification, "\n")

print(result$parameters)

####################################################################
data$cluster <- result$classification

cluster_to_plot <- 3
filtered_data <- data[data$cluster == cluster_to_plot, ]

sum_by_cluster <- filtered_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), sum))

melted_data <- sum_by_cluster %>%
  pivot_longer(cols = -cluster, names_to = "genre", values_to = "sum_instances")

colors <- c("Drama" = "#A46422", "Comedy" = "#800080", "Romance" = "#E06F8B", "Action" = "#BE2633",
            "Thriller" = "#493C2B", "Crime" = "#ffcccb", "Horror" = "#1B2632", "Adventure" = "#EB8931",
            "Mystery" = "#005784", "Family" = "#31A2F2", "Fantasy" = "#44891A", "Music" = "#D2B48C",
            "Sci-Fi" = "#A3CE27", "Documentary" = "#2F484E", "History" = "#B2DCEF", "War" = "#000000", 
            "Animation" = "#CBC3E3", "Western" = "#918151", "Sport" = "#FFFF00", "Film-Noir" = "#A9A9A9",
            "NA" = "#FFFFFF")

genre_colors <- colors[melted_data$genre]

melted_data <- melted_data[order(-melted_data$sum_instances),]

pie(melted_data$sum_instances, labels = melted_data$genre, col = genre_colors, main = paste("Genre Distribution for Action, Horror (Cluster 4)")) 
####################################################################
data$cluster <- result$classification

sum_by_cluster <- data %>%
  group_by(cluster) %>%
  summarise(across(everything(), sum))

melted_data <- sum_by_cluster %>%
  pivot_longer(cols = -cluster, names_to = "genre", values_to = "sum_instances")

colors <- c("Drama" = "#A46422", "Comedy" = "#800080", "Romance" = "#E06F8B", "Action" = "#BE2633",
            "Thriller" = "#493C2B", "Crime" = "#ffcccb", "Horror" = "#1B2632", "Adventure" = "#EB8931",
            "Mystery" = "#005784", "Family" = "#31A2F2", "Fantasy" = "#44891A", "Music" = "#D2B48C",
            "Sci-Fi" = "#A3CE27", "Biography" = "#2F484E", "History" = "#B2DCEF", "War" = "#000000", 
            "Animation" = "#CBC3E3", "Western" = "#918151", "Sport" = "#FFFF00", "Film-Noir" = "#A9A9A9",
            "NA" = "#FFFFFF")

cluster_colors <- c("1" = "#A46422", "2" = "#BE2633", "3" = "#493C2B", "4" = "#800080", "5" = "#2F484E", "6" = "#E06F8B", "7" = "#800080", "8" = "#800080", "9" = "lightblue", "10" = "lightgreen" )

ggplot(melted_data, aes(x = genre, y = sum_instances, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Sum of Each Genre by Cluster",
       x = "Genre",
       y = "Sum of Instances") +
  scale_fill_manual(values = cluster_colors, labels = c("Drama", "Action, Adventure, Sci Fi, Fantasy", "Thriller, Horror, Crime, Mystery", "Action, Adventure, Sci Fi, Fantasy", "Documentary", "Romance", "Romance, Comedy")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#########################################################
#Pie chart of world genres by count
count.genre1 <- movies2 %>% filter(genre1 != "Adult") %>% filter(genre1 != "Reality-TV") %>% filter(genre1 != "News") %>% count(genre1, name = "count.genre1") %>% rename(genre = genre1)
count.genre2 <- movies2 %>% filter(genre1 != "Adult") %>% filter(genre1 != "Reality-TV") %>% filter(genre1 != "News") %>% count(genre2, name = "count.genre2") %>% rename(genre = genre2)
count.genre3 <- movies2 %>% filter(genre1 != "Adult") %>% filter(genre1 != "Reality-TV") %>% filter(genre1 != "News") %>% count(genre3, name = "count.genre3") %>% rename(genre = genre3)
count.genre1 <- count.genre1[!(is.na(count.genre1$genre) | count.genre1$genre==""), ]#https://stackoverflow.com/questions/9126840/delete-rows-with-blank-values-in-one-particular-column
count.genre2 <- count.genre2[!(is.na(count.genre2$genre) | count.genre2$genre==""), ]
count.genre3 <- count.genre3[!(is.na(count.genre3$genre) | count.genre3$genre==""), ]
count.genre <- merge(count.genre1, count.genre2, by = "genre", all = TRUE)
count.genre <- merge(count.genre, count.genre3, by = "genre", all = TRUE)
count.genre[is.na(count.genre)] <- 0
count.genre <- count.genre %>% mutate(count.genre = count.genre1 + count.genre2 + count.genre3)
count.genre <- count.genre[order(-count.genre$count.genre),]
count.genre$genre <- factor(count.genre$genre, count.genre$genre)
pie(count.genre$count.genre, labels = count.genre$genre, col = colors)