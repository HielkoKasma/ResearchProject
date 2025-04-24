#Meeting 1 - analyzing metadata for the ATLANTIS study
RP_meta <- metatable_ATL_repeated_subjects
RP_meta
RP_meta$meth_file_id == "208107860041_R01C01"
table(RP_meta$gender)
table(RP_meta$gender,RP_meta$asthma.status)
table(RP_meta$gender, RP_meta$smoking.status)

table(RP_meta$gender, RP_meta$age)

table(RP_meta$gender, RP_meta$asthma.status, RP_meta$smoking.status)
fac <- factor(RP_meta$RNA_seq_id)
RP_meta$gender[,fac]
c(RP_meta$gender)
table(c(RP_meta$gender, RP_meta$age))
RP_meta


female <- c(13, 11, 2, 2, 11)
male <- c(11, 9, 2, 2, 9)
df <- data.frame(female, male)
rownames(df) <- c("number of patients", "asthma", "healthy", "ex-smoker", "smoker")
df
table(RP_meta$gender, RP_meta$astma.status)


#Meeting 2 - PCA plot - X and Y coordinates are stored in variable pca_data
#Importing the randomly selected 10% of the dataset
RP_data <- RP_tenper_random_namesfix 
#names of SubjectX_VisitX is assigned in csv file
head(RP_data)
RP_data[1,]


#packages required
install.packages("tidyverse")
install.packages("ggplot2")
library(tidyverse)
library(ggplot2)
library(ggrepel)

#principle component analysis 
#This code is used to find the PCs for data set RP_data

RP_numeric <- RP_data %>%
  select(where(is.numeric)) 
head(RP_numeric)

metatable <- `metatable_ATL_repeated_subjects.(1)` %>%
  mutate(pca_ids = colnames(RP_numeric))

prcomp(t(RP_numeric), scale. = T) #repeated in next line without saving in variable
pca <- prcomp(t(RP_numeric), scale. = T) #Scale = T --> data is measured in amount of standard deviations from mean
#this is used when comparing two variables that are different, and differently scaled (eg hight vs age)
#our data is already scaled, therefore we should run the code again with scale = F
head(pca)
plot(pca$x[,1], pca$x[,2]) #PCA plot of RP_data
head(RP_numeric)

#Scree plot to visualize impact of each PC
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var/sum(pca_var)*100, 1)
pca_var_percent 
barplot(pca_var_percent, main="Scree Plot", xlab="Principal Component (1 --> 48)", 
        ylab="Percent Variation", col = "dodgerblue2")
#PC 1 explains 29.1% of variation
#PC 2 explains 7.4% of variation
#PC1 and 2 explain the majority of the variation in the dataset

#pca_data contains the the values from RP_data, projected on PC1 and PC2 
#if the values withing [] are changed, a different PC is selected
#pca$x[,1] corresponds to PC1, pca$x[,48] corresponds to PC48
pca_data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1], 
                       Y=pca$x[,2]) %>%
  mutate(group_id = str_remove(Sample, "_.*"))  %>%
  left_join(metatable, by = c("Sample" = "pca_ids"))
print(pca_data)

#PCA plot with labeled data points
library(ggrepel)
ggplot(data = pca_data, aes(x = X, y = Y, label = Sample)) +
  geom_point(aes(color = group_id )) +
  geom_text_repel(size = 2.5, max.overlaps = Inf, 
                  segment.color = "grey0") +
  xlab(paste("PC1 - ", pca_var_percent[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_percent[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("PCA Plot")

#to view a labeled PCA plot for different PCs: change values in [] for ggplot and pca_data like described before

aes(colour = gender) #put this in geom_point() to assign colors based on groups
head(pca_data)

  
#to view a labeled PCA plot for different PCs: change values in [] for ggplot and pca_data like described before


#meeting 3
#code for making PC matrix
head(pca_data)
PC1_values <- pca_data$X
head(PC1_values)
distance_PC1 <- dist(PC1_values)
matrix_PC1 <- as.matrix(distance_PC1)
view(matrix_PC1)

# Get indices of the paired rows
self_indices <- seq(1, nrow(matrix_PC1), by = 2)

# Get distances between each subject’s two visits
self_distances <- mapply(function(i) matrix_PC1[i, i+1], self_indices)

# Mean self-stability
mean_self_distances <- mean(self_distances)
print(mean_self_distances)


non_self_distances <- c()

# Loop over all rows (samples)
for (i in 1:nrow(matrix_PC1)) {
  
  # Paired index: if i is odd, its pair is i+1; if even, pair is i-1
  pair_idx <- ifelse(i %% 2 == 1, i + 1, i - 1)
  
  # All other sample indices *excluding* the paired one
  compare_to <- setdiff(1:nrow(matrix_PC1), pair_idx)
  
  # Get distances from sample i to all other (non-pair) samples
  dists <- matrix_PC1[i, compare_to]
  
  # Add those distances to the vector
  non_self_distances <- c(non_self_distances, dists)
}

# Now calculate the mean non-self stability
non_self_distances_no_zeros <- non_self_distances[non_self_distances != 0]
mean_non_self_distances <- mean(non_self_distances_no_zeros)
print(mean_non_self_distances)

view(non_self_distances)
view(non_self_distances_no_zeros)

view(self_stability)

stability_PC1 <- mean_self_distances - mean_non_self_distances
stability_PC1
sd_self_distances <- sd(self_distances)
sd_non_self_distances <- sd(non_self_distances_no_zeros)


#this code is used to create results dataframe, only run the FIRST TIME
results_df <- data.frame(Value = numeric(5))
rownames(results_df) <- c(
  "Mean_Self_Distance",
  "SD_Self_Distance",
  "Mean_NonSelf_Distance",
  "SD_NonSelf_Distance",
  "Stability")

colnames(results_df) <- "PC1"
results_df["Mean_Self_Distance", "Value"] <- mean_self_distances
results_df["SD_Self_Distance", "Value"] <- sd_self_distances
results_df["Mean_NonSelf_Distance", "Value"] <- mean_non_self_distances
results_df["SD_NonSelf_Distance", "Value"] <- sd(non_self_distances_no_zeros)
results_df["Stability", "Value"] <- mean_self_distances - mean_non_self_distances
results_df

#this code adds new columns to the results dataframe, ONLY RUN FROM PC2 AND ONWARDS!!!

new_column <- c(
  mean_self_distances,
  sd_self_distances,
  mean_non_self_distances,
  sd(non_self_distances_no_zeros),
  mean_self_distances - mean_non_self_distances)

results_df[["PC6"]] <- new_column
results_df
view(results_df)

rownames(results_df) <- c(
  "Mean Self Distance",
  "SD Self Distance",
  "Mean Non-Self Distance",
  "SD Non-Self Distance",
  "Stability")
results_df
view(results_df)








