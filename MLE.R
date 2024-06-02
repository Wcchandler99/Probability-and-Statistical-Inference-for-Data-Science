library(stats)
library(stats4)

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

gamma_mle <- function(samples) {
  
  # Calculate sample mean and variance
  sample_mean <- mean(samples)
  sample_variance <- var(samples)
  
  # Method of Moments estimates for shape (alpha) and scale (beta)
  alpha_hat <- sample_mean^2 / sample_variance
  beta_hat <- sample_variance / sample_mean
  
  data <- samples
  initial_alpha <- alpha_hat  # Your initial estimate for alpha
  initial_beta <- beta_hat   # Your initial estimate for beta
  
  # Ensure data only contains positive values
  data <- data[data > 0]
  

  
  # Optimization using optim with method L-BFGS-B to add bounds
  optim_result <- optim(
    par = c(initial_alpha, initial_beta),  # Initial estimates for alpha and beta
    fn = gamma_log_likelihood,             # Objective function to minimize
    method = "L-BFGS-B",                   # Optimization method that supports bounds
    lower = c(1e-6, 1e-6),                 # Lower bounds close to zero for alpha and beta
    upper = c(Inf, Inf),                    # Upper bounds set to infinity for alpha and beta
    data = data
  )
  
  # Retrieve the MLE estimates for alpha and beta
  mle_alpha <- optim_result$par[1]
  mle_beta <- optim_result$par[2]
  
  # Print the MLE estimates
  cat("MLE for alpha:", mle_alpha, "\n")
  cat("MLE  for beta:", mle_beta, "\n")
  
  return(list(mle_alpha=mle_alpha, mle_beta=mle_beta))
}

n_samples <- 100
shape_true <- 2 # alpha
scale_true <- 3 # beta
samples <- rgamma(n_samples, shape = shape_true, scale = scale_true)
samples
gamma_mle(samples)


# Log-likelihood function for the Beta distribution using digamma function
logLikBeta <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  epsilon <- 1e-8
  adjusted_data <- pmin(pmax(data, epsilon), 1 - epsilon) # Adjust data to avoid log(0)
  sum_log_x <- sum(log(adjusted_data))
  sum_log_one_minus_x <- sum(log(1 - adjusted_data))
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

beta_mle <- function(samples) {
  # Calculate sample mean and variance
  sample_mean <- mean(samples)
  sample_variance <- var(samples)
  
  # Method of Moments estimates for alpha and beta
  alpha_hat <- ((1 - sample_mean) / sample_variance - 1 / sample_mean) * sample_mean^2
  beta_hat <- alpha_hat * (1 / sample_mean - 1)
  
  # Output the estimates
  #cat("Estimated Alpha from method of moments is :", alpha_hat, "\n")
  #cat("Estimated Beta from method of moments is :", beta_hat, "\n")
  
  data <- samples
  
  # Initial parameter guesses for alpha and beta
  start_params <- c(alpha = 2, beta = 5)
  
  # Optimize the parameters
  mle <- optim(start_params, logLikBeta, gr = scoreFunction, method = "BFGS", control = list(fnscale = -1), data = data)
  
  # Extract the MLE for alpha and beta
  alpha_mle <- mle$par[1]
  beta_mle <- mle$par[2]
  
  # Output the results
  cat("MLE for alpha:", alpha_mle, "\n")
  cat("MLE for beta:", beta_mle, "\n")
  
  return(list(mle_alpha=alpha_mle, mle_beta=beta_mle))
}

n_samples <- 100
alpha_true <- 2
beta_true <- 5
samples <- rbeta(n_samples, alpha_true, beta_true)

beta_mle(samples)
