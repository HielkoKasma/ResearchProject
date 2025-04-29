
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
