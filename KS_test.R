library(stats)
library(stats4)

# Log-likelihood function for the gamma distribution
gamma_log_likelihood <- function(params, data) {
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
  
  return(list(mle_alpha=mle_alpha, mle_beta=mle_beta))
}

# Log-likelihood function for the Beta distribution using digamma function
logLikBeta <- function(params, data) {
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
scoreFunction <- function(params, data) {
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
  
  return(list(mle_alpha=alpha_mle, mle_beta=beta_mle))
}

goodness_of_fit <- function(x, nboot, test_function){
  
  if (test_function == 'binomial'){
    x_test <- rbinom(n = length(x), prob = mean(x)/length(x), size = 1)
    
    ks <- ks.test(x, 'pbinom', prob = mean(x)/length(x), size = length(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'pbinom', prob = mean(column_data)/length(x_test), size = length(x))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'exponential'){
    x_test <- rexp(n = length(x), rate = 1/mean(x))
    
    ks <- ks.test(x, 'pexp', rate = 1/mean(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'pexp', rate = 1/mean(column_data))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'beta'){
    alpha <- beta_mle(x)$mle_alpha
    beta <- beta_mle(x)$mle_beta
    
    x_test <- rbeta(n = length(x), shape1 = alpha, shape2 = beta)
    
    ks <- ks.test(x, 'pbeta', shape1 = alpha, shape2 = beta)
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    
    ks_test_column <- function(column_data){
      alpha_test <- beta_mle(column_data)$mle_alpha
      beta_test <- beta_mle(column_data)$mle_beta
      ks_results <- ks.test(column_data, 'pbeta', shape1 = alpha_test, shape2 = beta_test)
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'gamma'){
    alpha <- gamma_mle(x)$mle_alpha
    beta <- gamma_mle(x)$mle_beta
    
    x_test <- rgamma(n = length(x), shape = alpha, scale  = beta)
    
    ks <- ks.test(x, 'pgamma', shape = alpha, scale  = beta)
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      alpha_test <- beta_mle(column_data)$mle_alpha
      beta_test <- beta_mle(column_data)$mle_beta
      ks_results <- ks.test(column_data, 'pgamma', shape = alpha_test, scale = beta_test)
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'geometric'){
    x_test <- rgeom(n = length(x), 1/mean(x))
    
    ks <- ks.test(x, 'pgeom', p = 1/mean(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    alpha_test <- beta_mle(x_test)[1]
    beta_test <- beta_mle(x_test)[2]
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'pgeom', p = 1/mean(x))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'poisson'){
    x_test <- rpois(n = length(x), mean(x))
    
    ks <- ks.test(x, 'ppois', mean(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'ppois', mean(x))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'uniform'){
    x_test <- runif(n = length(x), min(x), max(x))
    
    ks <- ks.test(x, 'punif', min(x), max(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x_test), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'punif', min(x), max(x))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  if (test_function == 'normal'){
    x_test <- rnorm(length(x))
    
    ks <- ks.test(x, 'pnorm', mean(x), sd(x))
    
    bootstrap_x <- replicate(nboot, {
      sample_data <- sample(x_test, length(x), replace = TRUE)
    })
    
    ks_test_column <- function(column_data){
      ks_results <- ks.test(column_data, 'pnorm', mean(column_data), sd(column_data))
      return(ks_results$p.value)
    }
    
    ks_bootstrap <- apply(bootstrap_x, 2, ks_test_column)
  }
  
  p_value_comparisons <- table(ks['p.value'] < ks_bootstrap)
  
  p_value_proportions <- prop.table(p_value_comparisons)
  
  return(p_value_proportions)
}


x <- rgamma(100, 3, 2)

goodness_of_fit(x, 10000, 'uniform')

