#Importing the mnorm library which is used for obtaining density function of Multivariate Normal
library("mnorm")

check_significant_change <- function(t1, t2, t3, pt1, pt2, pt3) {
  max_change <- max(c(max(abs(t1 - pt1)), max(abs(t2 - pt2)), max(abs(t3 - pt3))))
  return(max_change > 0.00000005)
}

identify_cluster <- function(t1, t2, t3) {
  clus_num <- apply(cbind(t1, t2, t3), 1, function(row) which.max(row))
  return(clus_num)
}

create_pair_plot <- function(data, out, clus_col) {
  pairs(data, pch = 16, col = clus_col, main = "Pairs Plot")
  points(rbind(out$mu1, out$mu2, out$mu3), pch = 3, col = 1:3, cex = 2)
  cat("\n\nPress Enter to Continue\n\n")
  readline()
  return(clus_col)
}

create_diff_plot <- function(data, out) {
  d1 <- out$mu2 - out$mu3
  d2 <- out$mu3 - out$mu1
  
  x_s <- c()
  y_s <- c()
  
  for(i in 1:length(data[,1])) {
    #The below are the formulas used to calculate projected values using Gram-Schmidt Method
    x_s <- append(x_s, (unlist(data[i,])%*%t(unlist(d1))) / ((unlist(d1)%*%t(unlist(d1)))^0.5) )
    y_s <- append(y_s, (unlist(data[i,])%*%t(unlist(d2))) / ((unlist(d2)%*%t(unlist(d2)))^0.5) )
  }
  
  two_d_dat <- data.frame(x=x_s, y=y_s, clus=clus_col) 
  
  colors <- c("#e41a1c","#377eb8","#4daf4a")[two_d_dat$clus]
  plot(two_d_dat$x, two_d_dat$y, col = colors, pch = 16, xlab = 'X', ylab = 'Y', main="Plot of Gram-Schmidt Orthogonalized Values")
  legend("topright", legend = levels(factor(two_d_dat$clus)), col = c("#e41a1c","#377eb8","#4daf4a"), pch = 16, title = "Cluster")
  cat("\n\nPress Enter to Continue\n\n")
  readline()
  return(two_d_dat)
}

print_final_values <- function(em_out, two_d_dat) {
  cat("Cluster 1 mean :\n")
  print(em_out$mu1)
  cat("\nCluster 2 mean :\n")
  print(em_out$mu2)
  cat("\nCluster 3 mean :\n")
  print(em_out$mu3)
  cat("\n\nCluster 1 sigma :\n")
  print(em_out$sigma1)
  cat("\nCluster 2 sigma :\n")
  print(em_out$sigma2)
  cat("\nCluster 3 sigma :\n")
  print(em_out$sigma3)
  cat("\n\n", two_d_dat, "\n")
}


# Function for the EM Algorithm
EM_alg <- function(data) {
  
  #Initializing Cluster Means to be a random data point
  mu1 <- as.matrix(data[1,])
  mu2 <- as.matrix(data[2,])
  mu3 <- as.matrix(data[3,])
  
  #Initializing the Cluster Covariance Matrices to be Diagonal Matrices
  sigma1 <- diag(5)
  sigma2 <- diag(5)
  sigma3 <- diag(5)
  
  
  #Initializing previous Tau values to 0s
  prev_tau1 <- numeric(length(data[,1]))
  prev_tau2 <- numeric(length(data[,1]))
  prev_tau3 <- numeric(length(data[,1]))
  
  #Initializing current Tau values to 1s
  tau1 <- seq(1,1,length.out = length(data[,1]))
  tau2 <- seq(1,1,length.out = length(data[,1]))
  tau3 <- seq(1,1,length.out = length(data[,1]))
  
  #Counter to keep track of number of iterations
  kk<-0
  
  # A while loop which runs until there is no significant change in parameters
  while(check_significant_change(tau1,tau2,tau3,prev_tau1,prev_tau2,prev_tau3)) { 
    
    print(kk)
    kk <- kk + 1
    
    #Updating Previous Tau Values
    prev_tau1 <- tau1
    prev_tau2 <- tau2
    prev_tau3 <- tau3
    
    #Assigning Current Tau values to be empty matrices
    tau1 <- c()
    tau2 <- c()
    tau3 <- c()
    
    #Taking cluster probability values to be random but their sum must be 1 
    p1 <- 0.2
    p2 <- 0.1
    p3 <- 0.7
    
    # Calculating Tau value for each data point using Multivariate Normal Density Function
    for (x in (1:length(data[,1]))) {
      t1 <- p1*(dmnorm(as.matrix(data[x, ]), mu1, sigma1)$den)
      t2 <-  p2*(dmnorm(as.matrix(data[x, ]), mu2, sigma2)$den)
      t3 <-  p3*(dmnorm(as.matrix(data[x, ]), mu3, sigma3)$den)
      
      #This is a condition to prevent 0 division error as there is a chance of very small Tau for which system may assume value to be 0
      if(t1+t2+t3 != 0){ 
        tau1 <-  append(tau1, t1/(t1+t2+t3))
        tau2 <-  append(tau2, t2/(t1+t2+t3))
        tau3 <-  append(tau3, t3/(t1+t2+t3))
      }
      
      else{
        tau1 <-  append(tau1, 0)
        tau2 <-  append(tau2, 0)
        tau3 <-  append(tau3, 0)
      }
      
    }
    
    
    #Finding the sum of Tau values
    sum_tau1 <- sum(tau1)
    sum_tau2 <- sum(tau2)
    sum_tau3 <- sum(tau3)
    
    #Updating P values
    p1 <- sum_tau1/(sum_tau1+sum_tau2+sum_tau3)
    p2 <- sum_tau2/(sum_tau1+sum_tau2+sum_tau3)
    p3 <- sum_tau3/(sum_tau1+sum_tau2+sum_tau3)
    
    #Take a list for storing mean values
    mu_numerator_values1 <- list()
    mu_numerator_values2 <- list()
    mu_numerator_values3 <- list()
  
    
    #Use matrix operations to avoid for loop and at once calculate new mean values
    for (x in (1:length(data[,1]))){
      mu_numerator_values1[[x]] <- as.matrix(tau1[x]*as.matrix(data[x, ]))
      mu_numerator_values2[[x]] <- as.matrix(tau2[x]*as.matrix(data[x, ]))
      mu_numerator_values3[[x]] <- as.matrix(tau3[x]*as.matrix(data[x, ]))
    }
    
    #Use Reduce operation to find sum of mean values for each data point to get sum of means
    mu_numerator1 <- Reduce('+', mu_numerator_values1)
    mu_numerator2 <- Reduce('+', mu_numerator_values2)
    mu_numerator3 <- Reduce('+', mu_numerator_values3)
    
    #Find the final new mean values by dividing by sum of that particular cluster Tau value
    mu1 <- mu_numerator1/sum_tau1
    mu2 <- mu_numerator2/sum_tau2
    mu3 <- mu_numerator3/sum_tau3
    
    #Take a list for storing covariance values
    sigma_numerator_values1 <- list()
    sigma_numerator_values2 <- list()
    sigma_numerator_values3 <- list()
    
    
    #Use matrix operations to avoid for loop and at once calculate new covariance values
    for (x in (1:length(data[,1]))){
      sigma_numerator_values1[[x]] <-  as.matrix(tau1[x]*( t(as.matrix(data[x,]) - mu1)%*%(as.matrix(data[x,]) - mu1) ))
      sigma_numerator_values2[[x]] <-  as.matrix(tau2[x]*( t(as.matrix(data[x,]) - mu2)%*%(as.matrix(data[x,]) - mu2) ))
      sigma_numerator_values3[[x]] <-  as.matrix(tau3[x]*( t(as.matrix(data[x,]) - mu3)%*%(as.matrix(data[x,]) - mu3) ))
    }
    
    #Use Reduce operation to find sum of covariance values for each data point to get sum of covariances
    sigma_numerator1 = Reduce('+', sigma_numerator_values1)
    sigma_numerator2 = Reduce('+', sigma_numerator_values2)
    sigma_numerator3 = Reduce('+', sigma_numerator_values3)
    
    
    #Find the final new covariance values by dividing by sum of that particular cluster Tau value
    sigma1 <- sigma_numerator1/sum_tau1
    sigma2 <- sigma_numerator2/sum_tau2
    sigma3 <- sigma_numerator3/sum_tau3
  
  }
  
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,tau1=tau1,tau2=tau2,tau3=tau3))
}


#Perform all the functions and produce output

data <- read.csv("C:\\Users\\wccha\\Documents\\Rutgers\\Prob_&_Stat\\Final\\archive\\user_genre_ratings.csv")
data <- data[ , -1]
em_out <- EM_alg(data)
clus_col <- identify_cluster(em_out$tau1,em_out$tau2,em_out$tau3)
create_pair_plot(data, em_out, clus_col)
two_d_dat <- create_diff_plot(data, em_out)
print_final_values(em_out,two_d_dat)
