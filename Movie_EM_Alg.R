library("mnorm")

check_significant_change <- function(t_values, prev_t_values) {
  max_change <- max(abs(t_values - prev_t_values))
  return(max_change > 0.00000005)
}

identify_cluster <- function(t_values) {
  clus_num <- apply(t_values, 1, which.max)
  return(clus_num)
}

create_pair_plot <- function(data, out, clus_col) {
  pairs(data, pch = 16, col = clus_col, main = "Pairs Plot")
  for (i in seq_along(out$mus)) {
    points(out$mus[[i]], pch = 3, col = i)
  }
  cat("\n\nPress Enter to Continue\n\n")
  readline()
  return(clus_col)
}

create_diff_plot <- function(data, out) {
  d_values <- lapply(1:(length(out$mus) - 1), function(i) out$mus[[i + 1]] - out$mus[[i]])
  
  x_s <- c()
  y_s <- c()
  
  for (i in seq_along(data[, 1])) {
    x_s <- append(x_s, sapply(d_values, function(d) (unlist(data[i, ]) %*% t(unlist(d))) / (sum(unlist(d)^2))^0.5 ))
  }
  
  two_d_dat <- data.frame(x = x_s, clus = identify_cluster(out$t_values))
  
  colors <- rainbow(length(out$mus))[two_d_dat$clus]
  plot(two_d_dat$x, col = colors, pch = 16, xlab = 'X', ylab = 'Y', main = "Plot of Orthogonalized Values")
  legend("topright", legend = seq_along(out$mus), col = rainbow(length(out$mus)), pch = 16, title = "Cluster")
  cat("\n\nPress Enter to Continue\n\n")
  readline()
  return(two_d_dat)
}

print_final_values <- function(out, two_d_dat) {
  for (i in seq_along(out$mus)) {
    cat(paste("Cluster", i, "mean :\n"))
    print(out$mus[[i]])
    cat("\nCluster", i, "sigma :\n")
    print(out$sigmas[[i]])
  }
  cat("\n\n", two_d_dat, "\n")
}

# Function for the EM Algorithm
library(MASS)

library(mvtnorm)

EM_alg <- function(data, num_clusters = 3, max_iterations = 1000, convergence_threshold = 0.00000005) {
  mus <- lapply(1:num_clusters, function(i) as.matrix(data[sample(nrow(data), 1), , drop = FALSE]))
  sigmas <- lapply(1:num_clusters, function(i) diag(ncol(data)))
  
  t_values <- matrix(1, nrow = nrow(data), ncol = num_clusters)
  
  prev_t_values <- numeric(length = nrow(data))
  
  kk <- 0
  
  while (kk < max_iterations && check_significant_change(t_values, prev_t_values)) {
    print(kk)
    kk <- kk + 1
    
    prev_t_values <- t_values
    
    p_values <- rep(1/num_clusters, num_clusters)
    
    t_values <- matrix(0, nrow = nrow(data), ncol = num_clusters)
    
    for (i in seq_along(mus)) {
      for (x in seq_len(nrow(data))) {
        t_values[x, i] <- p_values[i] * dmvnorm(data[x, , drop = FALSE], mean = as.vector(mus[[i]]), sigma = sigmas[[i]])
      }
    }
    
    #t_values <- t(t_values / rowSums(t_values))
    
    mus <- lapply(seq_along(mus), function(i) {
      numerator_values <- lapply(1:nrow(data), function(x) t_values[x, i] * as.matrix(data[x, , drop = FALSE]))
      numerator <- Reduce('+', numerator_values)
      numerator / sum(t_values[, i])
    })
    
    sigmas <- lapply(seq_along(sigmas), function(i) {
      numerator_values <- lapply(1:nrow(data), function(x) t_values[x, i] * tcrossprod(as.matrix(data[x, , drop = FALSE] - mus[[i]])))
      numerator <- Reduce('+', numerator_values)
      numerator / sum(t_values[, i])
    })
  }
  
  return(list(mus = mus, sigmas = sigmas, t_values = t_values))
}
# Perform all the functions and produce output
data <- read.csv("C:\\Users\\wccha\\Documents\\Rutgers\\Prob_&_Stat\\Final\\archive\\user_genre_ratings_no_0s.csv")
data <- data[, -1]
out <- EM_alg(data, num_clusters = 10)
clus_col <- identify_cluster(out$t_values)
create_pair_plot(data, out, clus_col)
two_d_dat <- create_diff_plot(data, out)
print_final_values(out, two_d_dat)