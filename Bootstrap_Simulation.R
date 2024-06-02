library(ggplot2)
library(reshape2)

# This is the function for the bootstrap estimator. 
# Argument 1: Sampled data from lognormal distribution
# Argument 2: The number of bootstrap samples we generate
bootstrap_estimator <- function(x1, nboot) {
  n <- length(x1) # Length of the samples vector
  bootstrap_x1 = c()
  
  # Creating 'nboot' bootstrap samples and calculating mean and sd for each sample
  bootstrap_x1 <- replicate(nboot, {
    sample_data <- sample(x1, length(x1), replace = TRUE)
    c(mean = mean(sample_data), sd = sd(sample_data))
  })
  
  if (n==3){
  # Convert to matrix and filter out rows where sd is zero
  bootstrap_x1_matrix <- t(bootstrap_x1)  # Transpose to get the correct structure
  bootstrap_x1 <- t(bootstrap_x1_matrix[bootstrap_x1_matrix[, "sd"] != 0, ])
  }
  
  
  
  bootstrap_means <- bootstrap_x1[1, ]  # Extracting bootstrap means
  bootstrap_sds <- bootstrap_x1[2, ]    # Extracting bootstrap standard deviations
  
  # Calculating the standard deviation for the means of boostrap samples
  standard_deviation_bootstrap_x1 <- sd(bootstrap_means)
  
  # Calculating the bias
  #bias_bootstrap_x1 <- mean(bootstrap_x1) - mean(x1)
  bias_bootstrap_x1 <- mean(bootstrap_means) - mean(x1)
  
  # Correct the sample mean using the bootstrap bias
  bias_corrected_mean = mean(x1) - bias_bootstrap_x1
  
  # Calculating the 95% confidence interval
  confidence_interval_bootstrap_x1 <- quantile(bootstrap_means, probs = c(.025, .975))
  
  # Returns a list of the original mean, bias corrected mean, bias, standard deviation and the bootstrap 95% confidence interval
  list(
    bootstrap_x1 = bootstrap_means,
    sd_for_pivotal = bootstrap_sds,
    original_mean = mean(x1),
    bias_corrected_mean = bias_corrected_mean,
    bias = bias_bootstrap_x1,
    sd = standard_deviation_bootstrap_x1,
    ci = confidence_interval_bootstrap_x1
  )
}

# This is the function for the jackknife estimator
# Argument : Sampled data from the lognormal distribution
jackknife_estimator <- function(x2){
  n <- length(x2) # Length of the samples vector
  jackknife_x2 <- c()
 
  # Creating n jackknife samples and calculating the mean of each
  for (i in 1:n) {
    jackknife_x2[i] <- mean(x2[-i])  # Exclude the ith observation
  }
  
  # Compute the bias-corrected estimate of the mean
  bias_corrected_mean_jackknife <- n * mean(x2) - (n-1) * mean(jackknife_x2)
  # Compute the bias
  bias_jackknife_x2 <- mean(jackknife_x2) - mean(x2)
  # Compute the standard deviation
  standard_deviation_jackknife_x2 <- sd(jackknife_x2)
  # Compute the 95% confidence interval based on bias corrected mean
  confidence_interval_jackknife_x2 <- c(bias_corrected_mean_jackknife - 1.96 * sd(jackknife_x2) / sqrt(n), bias_corrected_mean_jackknife + 1.96 * sd(jackknife_x2) / sqrt(n))
  
  # Returns a list of bias, standard deviation, and 95% confidence interval
  list(
    bias = bias_jackknife_x2,
    sd = standard_deviation_jackknife_x2,
    ci = confidence_interval_jackknife_x2
  )
}

# This is the function for the simulator
# Argument 1 : Number of samples to be drawn from the lognormal distribution
# Argument 2 : Value of alpha to be used in the confidence interval
simulation <- function(n, alpha){
  
  # Sampling n values from a lognormal distribution
  x3 <- rlnorm(n)
  
  mean_x3 <- mean(x3) # Mean of the samples
  sd_x3 <- sd(x3) # Standard deviation of the samples
  
  # Calculating confidence interval based on central limit theorem
  confidence_interval_clt <- c(mean_x3 - qnorm(1-alpha/2) * sd_x3/sqrt(n), mean_x3 + qnorm(1-alpha/2) * sd_x3/sqrt(n))
  
  # Results returned from bootstrap estimator
  bootstrap_results <- bootstrap_estimator(x3, 10000)
  
  # Bootstrap percentile confidence interval
  bootstrap_results <- bootstrap_estimator(x3, 10000)
  confidence_interval_bootstrap_percentile <- bootstrap_results$ci
  
  # Bootstrap pivotal confidence interval
  pivot_adjust <- c()
  nboot <- 10000
  count_sd <- 1
  for (i in bootstrap_results$bootstrap_x1){
    pivot_adjust <- append(pivot_adjust, (i - mean_x3)/(bootstrap_results$sd_for_pivotal[count_sd]/sqrt(n)))
    count_sd <- count_sd+1
  }
  pivot_raw <- quantile(pivot_adjust, probs = c(.025, .975))
  m1 <- mean(bootstrap_results$bootstrap_x1)
  s1 <- bootstrap_results$sd
  c1 <- pivot_raw[1]
  c2 <- pivot_raw[2]
  
  
  confidence_interval_bootstrap_pivotal <- c(mean_x3-c2*(sd_x3/sqrt(n)) , mean_x3-c1*(sd_x3/sqrt(n)))
  
  # Results returned from jackknife estimator
  jackknife_results <- jackknife_estimator(x3)
  
  # Jackknife confidence interval
  confidence_interval_jackknife <- jackknife_results$ci
  
  # List of CI from central limit theorem, Bootstrap percentile and pivotal CIs
  list(
    ci_clt = confidence_interval_clt,
    ci_bootstrap_percentile = confidence_interval_bootstrap_percentile,
    ci_bootstrap_pivotal = confidence_interval_bootstrap_pivotal
    #ci_jackknife = confidence_interval_jackknife
  )
}

# This is the function to compare the coverage rates of the different CIs
# Argument 1 : Vector of the samples sizes for which we need to compare (3,10,30,100)
# Argument 2 : Number of simulations that need to be done
# Argument 3 : The true mean of the lognormal distribution
#compare_coverage_rates <- function(sample_sizes, n_simulations=1000, true_mean=exp(0.5 + 0.5^2/2)) {
compare_coverage_rates <- function(sample_sizes, n_simulations=1000, true_mean=exp(0.5)) {
  
  coverage_rates <- matrix(0, nrow=length(sample_sizes), ncol=3)
  
  for(i in 1:length(sample_sizes)) {
    n <- sample_sizes[i] # Taking each sample size at a time
    
    coverage_counts <- rep(0, 3)
    
    for(j in 1:n_simulations) {
      results <- simulation(n, .05)
      
      CIs <- list(results$ci_clt, results$ci_bootstrap_percentile, results$ci_bootstrap_pivotal)
      
      # Checking if the true mean falls within the range of each of the CIs. If so, the count for that CI gets increased
      for(k in 1:3) {
        if(CIs[[k]][1] <= true_mean){
          if(CIs[[k]][2] >= true_mean) {
            coverage_counts[k] <- coverage_counts[k] + 1
          }
        }
      }
    }
    
    coverage_rates[i, ] <- coverage_counts / n_simulations # To get the ratio of the simulations in which true mean was present in the CIs
  }
  
  rownames(coverage_rates) <- sample_sizes
  colnames(coverage_rates) <- c("CLT", "Bootstrap Percentile", "Bootstrap Pivotal")
  
  return(coverage_rates) # For each of the sample sizes, the ratios for the three CIs is returned
}

# To print the coverage rates for sample sizes 3,10,30,100 for the CLT CI, Bootstrap percentile and pivotal CIs
coverage_rates <- compare_coverage_rates(c(3, 10, 30, 100))

print(coverage_rates)

# Transforming the matrix into a data frame for ggplot
coverage_rates_df <- as.data.frame(coverage_rates)
coverage_rates_df$SampleSize <- rownames(coverage_rates_df)
rownames(coverage_rates_df) <- NULL

# Melting the data frame into a long format suitable for ggplot
long_coverage_rates_df <- melt(coverage_rates_df, id.vars = 'SampleSize', variable.name = 'CI_Method', value.name = 'Coverage_Rate')
long_coverage_rates_df$SampleSize <- as.numeric(as.character(long_coverage_rates_df$SampleSize))

# Plotting the data
ggplot(long_coverage_rates_df, aes(x = SampleSize, y = Coverage_Rate, group = CI_Method, color = CI_Method)) +
  geom_line() +
  geom_point() +
  labs(title = "Coverage Rates for Different CI Methods Across Sample Sizes", x = "Sample Size", y = "Coverage Rate") +
  theme_minimal()

# This function is to estimate the bias in standard deviation when divided by n and compare for bootstrap and jackknife estimators
# Argument 1 : Number of elements to be sampled from lognormal distribution
# Argument 2 : Number of simulations to be performed
estimate_sd_bias <- function(n, n_simulations = 1000) {
  # Stores data for 1000 simulations
  biases <- replicate(n_simulations, {
    sample_data <- rnorm(n) # Sampling n elements from a lognormal distribution
    
    # Standard deviation when dividing by n
    sample_sd_n <- sqrt(sum((sample_data - mean(sample_data))^2)/n)
    
    # Standard deviation when dividing by n-1
    sample_sd_n_minus_1 <- sd(sample_data)
    
    # Bias returned from bootstrap estimator function
    bias_bootstrap <- bootstrap_estimator(sample_data, 10000)$bias
    
    # Bias returned from jackknife estimator function
    bias_jackknife <- jackknife_estimator(sample_data)$bias
    
    # A vector of standard deviations when taking n and n-1, bootstrap bias and jackknife bias
    c(n = sample_sd_n, n_minus_1 = sample_sd_n_minus_1, bootstrap = bias_bootstrap, jackknife = bias_jackknife)
  })
  
  # Returns the means of each variable over 1000 simulations
  rowMeans(biases)
}

# Prints the comparison for SDs and bias
#print(estimate_sd_bias(100,1000))

# Getting biases for n = 100 and n_simulations = 1000
biases <- estimate_sd_bias(100, 1000)
print(biases)
