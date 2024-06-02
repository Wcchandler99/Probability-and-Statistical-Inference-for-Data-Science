#BINOMIAL

# Sample data
sample_data <- rbinom(100, size = 100, prob = 0.5)  

sum_data <- sum(sample_data)
sum_data_squared <- sum(sample_data^2)

# Calculating sample mean and variance
sample_mean <- mean(sample_data)

# Method of moments estimates
p_hat <- 1+sample_mean-(sum_data_squared/sum_data)
n_hat <- sample_mean/(1+sample_mean-(sum_data_squared/sum_data))

# Output the estimates
cat("Estimated n through method of moments is:", n_hat, "\n")
cat("Estimated p through method of moments is:", p_hat, "\n")

n <- 100
successes <- sample_data

# MLE for p
p_hat <- mean(successes) / n

# Output results
cat("MLE for p:", p_hat, "\n")

########################################

#EXPONENTIAL

# Sample data
sample_data <- rexp(100, rate = 0.5)

# Calculating the sample mean
sample_mean <- mean(sample_data)

# Method of moments estimate for lambda
lambda_hat <- 1 / sample_mean

# Output the estimate
cat("Estimated lambda from method of moments is :", lambda_hat, "\n")

durations <- sample_data

# Calculate MLE for lambda (λ)
mle_lambda = 1 / mean(durations)

# Output the results
cat("MLE of lambda:", mle_lambda, "\n")

######
# BETA

n_samples <- 100
alpha_true <- 2
beta_true <- 5
samples <- rbeta(n_samples, alpha_true, beta_true)

# Calculate sample mean and variance
sample_mean <- mean(samples)
sample_variance <- var(samples)

# Method of Moments estimates for alpha and beta
alpha_hat <- ((1 - sample_mean) / sample_variance - 1 / sample_mean) * sample_mean^2
beta_hat <- alpha_hat * (1 / sample_mean - 1)

# Output the estimates
cat("Estimated Alpha from method of moments is :", alpha_hat, "\n")
cat("Estimated Beta from method of moments is :", beta_hat, "\n")

library(stats4)

data <- samples

# Log-likelihood function for the Beta distribution using digamma function
logLikBeta <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  sum_log_x <- sum(log(data))
  sum_log_one_minus_x <- sum(log(1 - data))
  n <- length(data)
  logLik <- (alpha - 1) * sum_log_x +
    (beta - 1) * sum_log_one_minus_x -
    n * lbeta(alpha, beta)
  return(-logLik) # negative because 'optim' minimizes
}

# Score function for the Beta distribution that includes digamma
scoreFunction <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  n <- length(data)
  sum_log_x <- sum(log(data))
  sum_log_one_minus_x <- sum(log(1 - data))
  score_alpha <- n * (digamma(alpha + beta) - digamma(alpha)) + sum_log_x
  score_beta <- n * (digamma(alpha + beta) - digamma(beta)) + sum_log_one_minus_x
  return(c(score_alpha, score_beta))
}

# Initial parameter guesses for alpha and beta
start_params <- c(alpha = 2, beta = 5)

# Optimize the parameters
mle <- optim(start_params, logLikBeta, gr = scoreFunction, method = "BFGS", control = list(fnscale = -1))

# Extract the MLE for alpha and beta
alpha_mle <- mle$par[1]
beta_mle <- mle$par[2]

# Output the results
cat("MLE for alpha:", alpha_mle, "\n")
cat("MLE for beta:", beta_mle, "\n")

####
# GAMMA

n_samples <- 100
shape_true <- 2 # alpha
scale_true <- 3 # beta
samples <- rgamma(n_samples, shape = shape_true, scale = scale_true)

# Calculate sample mean and variance
sample_mean <- mean(samples)
sample_variance <- var(samples)

# Method of Moments estimates for shape (alpha) and scale (beta)
alpha_hat <- sample_mean^2 / sample_variance
beta_hat <- sample_variance / sample_mean

# Output the estimates
cat("Estimated Shape (Alpha):", alpha_hat, "\n")
cat("Estimated Scale (Beta):", beta_hat, "\n")



# Load required libraries
library(stats)

data <- samples

# Your data and initial estimates
initial_alpha <- alpha_hat  # Your initial estimate for alpha
initial_beta <- beta_hat   # Your initial estimate for beta

# Ensure data only contains positive values
data <- data[data > 0]

# Log-likelihood function for the gamma distribution
gamma_log_likelihood <- function(params) {
  # Ensure alpha and beta are positive
  alpha <- abs(params[1])
  beta <- abs(params[2])
  n <- length(data)
  
  # Log-likelihood calculation with safeguards
  ll <- ifelse(alpha > 0 & beta > 0,
               n * alpha * log(beta) + (alpha - 1) * sum(log(data)) - beta * sum(data) - n * lgamma(alpha),
               -Inf)  # Return negative infinity if alpha or beta are non-positive to avoid invalid values
  return(-ll)  # Return negative log-likelihood for minimization
}

# Optimization using optim with method L-BFGS-B to add bounds
optim_result <- optim(
  par = c(initial_alpha, initial_beta),  # Initial estimates for alpha and beta
  fn = gamma_log_likelihood,             # Objective function to minimize
  method = "L-BFGS-B",                   # Optimization method that supports bounds
  lower = c(1e-6, 1e-6),                 # Lower bounds close to zero for alpha and beta
  upper = c(Inf, Inf)                    # Upper bounds set to infinity for alpha and beta
)

# Retrieve the MLE estimates for alpha and beta
mle_alpha <- optim_result$par[1]
mle_beta <- optim_result$par[2]

# Print the MLE estimates
cat("MLE for alpha:", mle_alpha, "\n")
cat("MLE  for beta:", mle_beta, "\n")

####
# GEOMETRIC

n_samples <- 100
p_true <- 0.3
samples <- rgeom(n_samples, prob = p_true) + 1 # rgeom returns number of failures before the first success, hence add 1

# Calculate sample mean
sample_mean <- mean(samples)

# Method of Moments estimate for p
p_hat <- 1 / sample_mean

# Output the estimate
cat("Estimated Probability (p):", p_hat, "\n")

trials <- samples

# Calculate MLE for p
mle_p = 1 / mean(trials)

# Output the results
cat("MLE of p:", mle_p, "\n")

#####
# POISSON

n_samples <- 100
lambda_true <- 4.5  
samples <- rpois(n_samples, lambda = lambda_true)

# Calculate sample mean
sample_mean <- mean(samples)

# Method of Moments estimate for lambda
lambda_hat <- sample_mean

# Output the estimate
cat("Estimated Lambda from method of moments is :", lambda_hat, "\n")

mle_lambda <- sample_mean

cat("Estimated Lambda from MLE:", mle_lambda, "\n")

####
# UNIFORM

n_samples <- 100
a_true <- 2
b_true <- 5
samples <- runif(n_samples, min = a_true, max = b_true)

# Calculate sample mean and variance
sample_mean <- mean(samples)
sample_variance <- var(samples)

# Method of Moments estimates for a and b
a_hat <- sample_mean - sqrt(3 * sample_variance)
b_hat <- sample_mean + sqrt(3 * sample_variance)

# Output the estimates
cat("Estimated a from method of moments is :", a_hat, "\n")
cat("Estimated b from method of moments is:", b_hat, "\n")

cat("Estimated a from MLE:", min(samples), "\n")
cat("Estimated b from MLE:", max(samples), "\n") 

# Creating the covariance matrix

a<-a_true
b<-b_true

denom_cov_uniform <- (n/(b-a)-n/((b-a)^2)) - ((2*n/((b-a)^3))^2)

h_matrix_uniform <- matrix(
  data = c(n/((b-a)^2) ,  (-2)*n/((b-a)^3) ,  (-2)*n/((b-a)^3) ,  n/((b-a)^2)), 
  nrow = 2,                             
  ncol = 2,                             
  byrow = TRUE                          
)

cov_matrix = 1/denom_cov_uniform*h_matrix_uniform

print("PRINTING COVARIANCE MATRIX FOR UNIFORM DISTRIBUTION")
print(cov_matrix)


#####
# NORMAL

n_samples <- 100
mu_true <- 0    
sigma_true <- 1 
samples <- rnorm(n_samples, mean = mu_true, sd = sigma_true)

# Calculate sample mean and standard deviation
sample_mean <- mean(samples)

# Method of Moments estimates for mu and sigma
mu_hat <- sample_mean

sigma_squared_hat <- (-1)*sample_mean^2 + 1/n_samples*sum(samples^2)

# Output the estimates
cat("Estimated Mean (mu) from method of moments is :", mu_hat, "\n")
cat("Estimated Sigma Squared from method of moments is :", sigma_squared_hat, "\n")

mle_mean <- sample_mean
mle_variance <- sum((samples-sample_mean)^2)/n_samples

# Output the estimates
cat("Estimated Mean from MLE :", mle_mean, "\n")
cat("Estimated Variance from MLE :", mle_variance, "\n")

# Getting the covariance matrix

mu <- mu_true
ss <- sigma_true^2 # true variance
n <- n_samples
ss2 <- ss^2 # variance squared
ss3 <- ss^3 # variance cubed

denom_cov_part1 <- (-1)*n/ss*(n/(2*ss2) - 1/ss3*(sum(samples)-n*mu))
denom_cov_part2 <- ((-1)/ss2*(sum(samples)-n*mu)) ^ 2

denom_cov_normal <- denom_cov_part1-denom_cov_part2

matrix_1_1 <- n/(2*ss2)-1/ss3*((sum(samples)-mu) ^ 2)
matrix_1_2 <- 1/ss2*(sum(samples)-n*mu)
matrix_2_1 <- 1/ss2*(sum(samples)-n*mu)
matrix_2_2 <- (-1)*n/ss

h_matrix_normal <- matrix(
  data = c(matrix_1_1, matrix_1_2, matrix_2_1, matrix_2_2), 
  nrow = 2,                             
  ncol = 2,                             
  byrow = TRUE                          
)

cov_matrix = 1/denom_cov_normal*h_matrix_normal

print("PRINTING COVARIANCE MATRIX FOR NORMAL DISTRIBUTION")
print(cov_matrix)
