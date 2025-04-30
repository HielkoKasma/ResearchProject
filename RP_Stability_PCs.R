
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

view(restls_df)

#calculating ICC: 

library(dplyr)
library(tidyr)
install.packages("psych")
library(psych)

#needed to have the same sample names in the metatable as int he methylation data (s1_v1), etc... this took me ages. some of the code is redundant
sample_names <- colnames(RP_data)
sample_names
sample_names <- sample_names[!sample_names %in% "X"]
metatable$pca_ids <- sample_names
all(metatable$pca_ids%in% colnames(RP_data))
library(psych)
install.packages("tidyverse")
library(tidyverse)
install.packages("tidyr")
library(tidyr)
RP_data[] <- lapply(RP_data, as.numeric)
RP_data_ICC <- RP_data %>%
  pivot_longer(cols = everything(), names_to = "SampleID_Visit", values_to = "Methylation_Value")
View(RP_data_ICC)
head(RP_data_ICC)
metatable <- metatable %>%
  rename(SampleID_Visit = pca_ids)

#i managed to make the merged data table: 

merged_data <- left_join(RP_data_ICC, metatable, by = "SampleID_Visit")
View(merged_data)
library(tidyverse)
library(tibble)
RP_data <- RP_data %>%
  rownames_to_column(var = "CpG_Site")
View

#made a copy of the methylation data with CPG sites as a column title instead of: 

NEW_RP<-RP_tenper_random_namesfix
View(RP_tenper_random_namesfix)
View(NEW_RP)
head(rownames(NEW_RP))
head(NEW_RP[, 1:5])
colnames(NEW_RP)[1] <- "CpG_Site" 

#i then chnaged the data into longform for the methylation 

if (!"CpG_Site" %in% colnames(NEW_RP)) {
  RP_data_long <- NEW_RP %>%
    rownames_to_column(var = "CpG_Site") %>%
    pivot_longer(-CpG_Site, names_to = "SampleID_Visit", values_to = "Methylation_Value")
} else {
  RP_data_long <- NEW_RP %>%
    pivot_longer(-CpG_Site, names_to = "SampleID_Visit", values_to = "Methylation_Value")
}

# merged data: 
  
View(RP_data_long)
merged_data <- left_join(RP_data_long, metatable, by = "SampleID_Visit")
View(merged_data)
head(merged_data)
nrow(merged_data)
sum(is.na(merged_data))

#made it wide again:

data_wide <- merged_data %>%
  select(CpG_Site, SampleID_Visit, Methylation_Value) %>%
  pivot_wider(names_from = SampleID_Visit, values_from = Methylation_Value)
View(data_wide)
install.packages("lme4")
library(lme4)
library(psych)

#PSYCH NOT WORKING

# Transpose the numeric data: rows become subjects, columns become CpG sites
icc_results <- psych::ICC(t(data_wide[,-1]))  # Remove CpG_Site column before transpose

# View results
View(icc_results)
print(icc_results$results)
colnames(data_wide)
data_wide <- data_wide %>%
  select(-CpG_Site)
icc_results <- psych::ICC(t(data_wide), model = "twoway", type = "consistency", unit = "average")

library(irr)

# Assuming your data is in long format where rows are methylation sites and columns are sample_visit
icc_data <- data_wide  # Replace with your data

# Calculate ICC using the 'icc()' function from 'irr' package
icc_results <- icc(icc_data, model = "twoway", type = "consistency", unit = "average")

# Print the results
print(icc_results)


#TRYING IRR

install.packages("irr")
library(irr)

# Assuming your data is in long format where rows are methylation sites and columns are sample_visit
icc_data <- data_wide  # Replace with your data

# Calculate ICC using the 'icc()' function from 'irr' package
icc_results <- icc(icc_data, model = "twoway", type = "consistency", unit = "average")

# Print the results
print(icc_results)

CpG_rand<-NEW_RP[9110,]
head(CpG_rand)
View(CpG_rand)
CpG_rand_clean <- CpG_rand[,-1]

library(tidyverse)
long_df <- CpG_rand_clean %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Value") %>%
  separate(Sample, into = c("SampleID", "Visit"), sep = "_") %>%
  mutate(
    SampleID = factor(SampleID, levels = paste0("S", 1:24)),
    Visit = factor(Visit, levels = c("V1", "V2"))
  ) %>%
  arrange(SampleID) %>%
  pivot_wider(names_from = Visit, values_from = Value) %>%
  select(V1, V2)

View(long_df)

library(psych)
icc_result_rand<-ICC(long_df)
print(icc_result_rand)

library(irr)
icc_result_rand_irr<-icc(long_df, model = "twoway", type = "consistency", unit = "single")
print(icc_result_rand_irr)
