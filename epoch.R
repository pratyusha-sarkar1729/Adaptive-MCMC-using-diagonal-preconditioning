library(mvtnorm)
library(coda)
library(ggplot2)
library(reshape2) 


CCIPCA <- function(X, l, k) {
  # X is an N x d data matrix.
  # l is the so called `amnesic` parameter; see Weng et al 2003 equation (8)
  # k is the number of e.vecs/vals output
  
  # Outputs an estimate of the first k leading eigenvalue + vector pairs of the
  # covariance of X (where it exists).
  
  N <- nrow(X)
  d <- ncol(X)
  
  # storage
  e_vecs <- matrix(nrow = k, ncol = d)
  e_vals <- vector(length = k)
  # surrogate data along orthogonal directions to save orthonormalising at each step:
  us <- matrix(nrow = k, ncol = d)
  for (n in 1:N) {
    us[1, ] <- X[n, ]
    for (i in 1:min(k, n)) {
      if (i == n) {
        e_vecs[i, ] <- us[i, ]
      } else {
        inner_prod <- sum(us[i, ] * e_vecs[i, ])
        euclidean_norm <- sqrt(sum(e_vecs[i, ] ^ 2))
        e_vecs[i, ] <- ((n - 1 - l) / n) * e_vecs[i, ] + inner_prod * ((1 + l) / (n * euclidean_norm)) * us[i, ]
        
        inner_prod <- sum(us[i, ] * e_vecs[i, ])
        euclidean_norm <- sqrt(sum(e_vecs[i, ] ^ 2))
        if (i < k) {
          us[i + 1, ] <- us[i, ] - (inner_prod / (euclidean_norm ^ 2)) * e_vecs[i, ]
          # # Sometimes us[i + 1, ] is the zero vector and sometimes it is close to
          # # us[i, ], so we have to correct for this:
          # size <- sqrt(sum(us[i + 1, ] * us[i + 1, ]))
          # absolute_angle <- abs(sum(us[i + 1, ] * us[i, ])) / (sqrt(sum(us[i + 1, ] * us[i + 1, ])) * sqrt(sum(us[i, ] * us[i, ])))
          # if (size < 10 ^ -5 | absolute_angle > 0.5){
          #   random_index <- sample(1:N, size = 1)
          #   random_sample <- X[random_index, ]
          #   us[i + 1, ] <- random_sample
          #   for (j in 1:i) {
          #     us[i + 1, ] <- us[i + 1, ] - ((sum(random_sample * us[j, ])) / sum(us[j, ] * us[j, ])) * us[j, ]
          #   }
          # }
        }
      }
    }
  }
  for (i in 1:k) {
    e_vals[i] <- sqrt(sum(e_vecs[i, ] ^ 2))
    e_vecs[i, ] <- (1 / e_vals[i]) * e_vecs[i, ]
  }
  return(list(e_vecs = e_vecs, e_vals = e_vals))
}

rwm_normal <- function(x, sigma, D, mu_pi, Sigma_pi, nits) {
  # Diagonally preconditioned Random Walk Metropolis on a normal target with
  # mean mu_pi and covariance Sigma_pi.
  d <- length(x)
  pi_curr <- dmvnorm(x, mean = mu_pi, sigma = Sigma_pi)
  
  # storage
  chain <- matrix(nrow = nits + 1, ncol = d)
  chain[1, ] <- x
  
  for (i in 1:nits) {
    # proposal
    x_prop <- x + sigma * D * rnorm(d) ##  # Propose a candidate state by adding a random step from a Gaussian proposal
    pi_prop <- dmvnorm(x_prop, mean = mu_pi, sigma = Sigma_pi)
    if(pi_prop == 0){
      alpha = 1
    }else{
      alpha <- min(1, pi_prop / pi_curr) 
    }
    if (runif(1) < alpha) {
      # accept
      x <- x_prop
      pi_curr <- pi_prop
    }
    # storage
    chain[i + 1, ] <- x
  }
  return(list(chain = chain))
}

rwm_normal_diagonal_adaptation <- function(x, sigma, D, mu_pi, Sigma_pi, nits) {
  # Diagonal adaptive RWM on a normal target with mean mu_pi and
  # covariance Sigma_pi
  d <- length(x)
  pi_curr <- dmvnorm(x, mean = mu_pi, sigma = Sigma_pi)
  mu <- x
  
  # storage
  chain <- matrix(nrow = nits + 1, ncol = d)
  chain[1, ] <- x
  diagonals <- matrix(nrow = nits + 1, ncol = d)
  diagonals[1, ] <- D
  
  for (i in 1:nits) {
    # proposal
    x_prop <- x + sigma * D * rnorm(d)
    pi_prop <- dmvnorm(x_prop, mean = mu_pi, sigma = Sigma_pi)
    alpha <- min(1, pi_prop / pi_curr)
    if (runif(1) < alpha) {
      # accept
      x <- x_prop
      pi_curr <- pi_prop
    }
    # adaptive step
    mu <- mu + (1 / (i + 1)) * (x - mu)
    D <- D + (1 / (i + 1)) * ((x - mu) ^ 2 - D)
    
    # storage
    chain[i + 1, ] <- x
    diagonals[i + 1, ] <- D
  }
  return(list(chain = chain, diagonals = diagonals))
}

construct_L <- function(lambdas, vs) {
  # vs is a k x d matrix whose rows are orthonormal eigenvector
  # Constructs a preconditioner based on eigeninformation
  k <- length(lambdas)
  d <- ncol(vs)
  # check whether there are the same amount of vs as lambdas
  if (nrow(vs) != k) {
    return(NaN)
  }
  # pad D out with 1s
  D <- diag(c(lambdas ^ (-0.5), rep(1, d - k)))
  vs <- t(vs)
  # Construct Q using SVD
  Q <- svd(vs, nu = d)$u
  return(Q %*% D %*% t(Q))
}

condition_number <- function(Sigma) {
  # outputs the condition number of a positive definite matrix
  # Sigma
  e_vals <- eigen(Sigma)$values
  kappa <- max(e_vals) / min(e_vals)
  return(kappa)
}

leading_eigen_info_normal <- function(x, sigma, epoch_lengths, mu_pi, Sigma_pi) {
  # x: the starting value of the chain
  # sigma: the global proposal scale
  # epoch_lengths: a vector of the lengths of the epochs
  # mu_pi: target mean
  # Sigma_pi: target variance matrix
  
  # Outputs an adaptive MCMC chain targeting a normal distribution. Adaptation
  # occurs at the end of each epoch, and infers eigeninformation and scale
  # information to use in the next epoc
  d <- length(x)
  no_epochs <- length(epoch_lengths)
  epoch_ends <- cumsum(epoch_lengths) + 1
  # storage
  chain <- matrix(nrow = sum(epoch_lengths) + 1, ncol = d)
  ESSs <- matrix(nrow = no_epochs, ncol = d)
  Sigmas <- array(dim = c(d, d, no_epochs + 1))
  Sigmas[,, 1] <- Sigma_pi
  kappas <- vector(length = no_epochs + 1)
  kappas[1] <- condition_number(Sigma_pi)
  epoch = 1
  for (epoch in 1:no_epochs) {
    if (epoch == 1) {
      # nonadaptive MCMC on pi ^ (1) (original target)
      nits <- epoch_lengths[epoch]; D <- rep(1, d)
      X <- rwm_normal(x, sigma, D, mu_pi, Sigma_pi, nits)$chain
      # storage
      chain[1:epoch_ends[epoch], ] <- X
      ESSs[epoch, ] <- effectiveSize(as.mcmc(X))
      
      # Adaptive stuff:
      # IMPORTANT: we must subtract the mean first before performing the CCIPCA
      X_shifted <- scale(X, scale = F)
      
      # get the leading eigeninformation
      l <- 2; k <- 1
      leading_eigen_info <- CCIPCA(X_shifted, l, k)
      lambdas <- leading_eigen_info$e_vals; vs <- leading_eigen_info$e_vecs
      
      # construct an appropriate transformation
      L <- construct_L(lambdas, vs)
      
      # transform the chain: Y \sim \pi ^ (2)
      Y <- t(L %*% t(X_shifted))
      
      # extract diagonal info
      D <- sqrt(diag(var(Y)))
      
      # storage
      new_Sigma <- diag(D ^ -1) %*% L %*% Sigma_pi %*% t(L) %*% diag(D ^ -1)
      Sigmas[, , epoch + 1] <- new_Sigma
      kappas[epoch + 1] <- condition_number(new_Sigma)
    } else {
      # nonadaptive MCMC
      nits <- epoch_lengths[epoch]
      # nonadaptive RWM on \pi ^ (epoch)
      X <- rwm_normal(Y[nrow(Y), ], sigma, D, mu_pi, L %*% Sigma_pi %*% t(L), nits)$chain
      # storage
      chain[(epoch_ends[epoch - 1] + 1):epoch_ends[epoch], ] <- X[2:nrow(X), ]
      ESSs[epoch, ] <- effectiveSize(as.mcmc(X))
      
      # Adaptive stuff:
      # IMPORTANT: we must subtract the mean first before performing the CCIPCA
      X_shifted <- scale(X, scale = F)
      
      # get the leading eigeninformation
      l <- 2; k <- 1
      leading_eigen_info <- CCIPCA(X_shifted, l, k)
      lambdas <- leading_eigen_info$e_vals; vs <- leading_eigen_info$e_vecs
      
      # construct an appropriate transformation from pi ^ (epoch) to
      # pi ^ (epoch + 1)
      new_L <- construct_L(lambdas, vs)
      
      # transform the chain: Y \sim \pi ^ (epoch + 1)
      Y <- t(new_L %*% t(rbind(Y, X_shifted)))
      # Y <- t(new_L %*% t(X_shifted))
      
      # extract diagonal info
      D <- sqrt(diag(var(Y)))
      
      # construct an appropriate transformation from pi ^ (1) to
      # pi ^ (epoch + 1)
      L <- new_L %*% L
      
      # storage
      new_Sigma <- diag(D ^ -1) %*% L %*% Sigma_pi %*% t(L) %*% diag(D ^ -1)
      Sigmas[, , epoch + 1] <- new_Sigma
      kappas[epoch + 1] <- condition_number(new_Sigma)
    }
  }
  return(list(chain = chain, ESSs = ESSs, Sigmas = Sigmas, kappas = kappas))
}

####################################################################################
####################################################################################

# Define a custom color palette
my_colors <- c("red", "blue", "green", "purple")  # You can customize the colors

# Define the number of iterations and the dimensions of the result array
num_iterations <- 5
num_rows <- length(D)
num_columns <- length(epoch_lengths) + 1

# Initialize a vector to store the results
results_array_1 <- numeric(3 * num_iterations * length(D))

# Loop through each iteration
for (i in 1:num_iterations) {
  rho_1 <- 0.1
  rho_2 <- 0.5
  rho_3 <- 0.9
  
  out_1 <- matrix(0, length(D), length(epoch_lengths) + 1)
  out_2 <- matrix(0, length(D), length(epoch_lengths) + 1)
  out_3 <- matrix(0, length(D), length(epoch_lengths) + 1)
  
  for (k in 1:length(D)) {
    Sigma_pi <- matrix(rep(rho_1, 2 * D[k]), nrow = D[k], ncol = D[k])
    for (j in 1:D[k]) {
      Sigma_pi[j, j] <- 1
    }
    mu_pi <- rep(0, D[k])
    x <- rnorm(D[k])
    sigma <- 2.38 / sqrt(length(x))
    epoch_lengths <- c(1000, 1000, 1000, 1000, 1000, 1000)
    
    out_1[k, ] = leading_eigen_info_normal(x, sigma, epoch_lengths, mu_pi, Sigma_pi)$kappas
    
    # Store the value of out_1[, 7] in the results_array_1
    results_array_1[(i - 1) * length(D) + k] <- out_1[k, 7]
  }
  
  for (k in 1:length(D)) {
    Sigma_pi <- matrix(rep(rho_2, 2 * D[k]), nrow = D[k], ncol = D[k])
    for (j in 1:D[k]) {
      Sigma_pi[j, j] <- 1
    }
    mu_pi <- rep(0, D[k])
    x <- rnorm(D[k])
    sigma <- 2.38 / sqrt(length(x))
    epoch_lengths <- c(1000, 1000, 1000, 1000, 1000, 1000)
    
    out_2[k, ] = leading_eigen_info_normal(x, sigma, epoch_lengths, mu_pi, Sigma_pi)$kappas
    
    # Store the value of out_2[, 7] in the results_array_1
    results_array_1[(i - 1) * length(D) + k + num_iterations * length(D)] <- out_2[k, 7]
  }
  
  for (k in 1:length(D)) {
    Sigma_pi <- matrix(rep(rho_3, 2 * D[k]), nrow = D[k], ncol = D[k])
    for (j in 1:D[k]) {
      Sigma_pi[j, j] <- 1
    }
    mu_pi <- rep(0, D[k])
    x <- rnorm(D[k])
    sigma <- 2.38 / sqrt(length(x))
    epoch_lengths <- c(1000, 1000, 1000, 1000, 1000, 1000)
    
    out_3[k, ] = leading_eigen_info_normal(x, sigma, epoch_lengths, mu_pi, Sigma_pi)$kappas
    
    # Store the value of out_3[, 7] in the results_array_1
    results_array_1[(i - 1) * length(D) + k + 2 * num_iterations * length(D)] <- out_3[k, 7]
  }
}

# Create a data frame with the combined results
df_combined <- data.frame(D = rep(D, 3 * num_iterations),
                          rho = rep(c(rho_1, rho_2, rho_3), each = length(D) * num_iterations),
                          X7 = results_array_1)


plot <- ggplot(df_combined, aes(x = factor(D), y = X7, fill = factor(D))) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none") + 
  facet_wrap(~rho)+xlab("dimensions")+ylab("values")

plot

# Define the number of iterations and the dimensions of the result array
num_iterations <- 5
num_rows <- length(D)
num_columns <- length(epoch_lengths) + 1

# Initialize a vector to store the results
results_array_2 <- numeric(2 * num_iterations * length(D))

# Initialize an index to keep track of the position in results_array_1
index <- 1

# Loop through each iteration
for (i in 1:num_iterations) {
  rho_1 <- 0.1
  rho_2 <- 0.5
  rho_3 <- 0.9
  
  out_1 <- matrix(0, length(D), length(epoch_lengths) + 1)
  out_2 <- matrix(0, length(D), length(epoch_lengths) + 1)
  
  for (k in 1:length(D)) {
    Sigma_pi <- matrix(rep(rho_1, 2 * D[k]), nrow = D[k], ncol = D[k])
    for (j in 1:D[k]) {
      Sigma_pi[j, j] <- 1
    }
    mu_pi <- rep(0, D[k])
    x <- rnorm(D[k])
    sigma <- 2.38 / sqrt(length(x))
    epoch_lengths <- c(1000, 1000, 1000, 1000, 1000, 1000)
    
    out_1[k, ] = leading_eigen_info_normal(x, sigma, epoch_lengths, mu_pi, Sigma_pi)$kappas
    
    # Store the value of out_1[, 7] in the results_array_1
    results_array_2[index] <- out_1[k, 7]
    index <- index + 1
  }
  
  for (k in 1:length(D)) {
    Sigma_pi <- matrix(rep(rho_2, 2 * D[k]), nrow = D[k], ncol = D[k])
    for (j in 1:D[k]) {
      Sigma_pi[j, j] <- 1
    }
    mu_pi <- rep(0, D[k])
    x <- rnorm(D[k])
    sigma <- 2.38 / sqrt(length(x))
    epoch_lengths <- c(1000, 1000, 1000, 1000, 1000, 1000)
    
    out_2[k, ] = leading_eigen_info_normal(x, sigma, epoch_lengths, mu_pi, Sigma_pi)$kappas
    
    # Store the value of out_2[, 7] in the results_array_1
    results_array_2[index] <- out_2[k, 7]
    index <- index + 1
  }
  
  # Create a data frame with the combined results
  df_combined_1 <- data.frame(D = rep(D, 2*num_iterations), rho = rep(c(rho_1,rho_2), each = length(D)*num_iterations), X7 = results_array_2)
}

## plot for 5 iterations
plot_1_5 <- ggplot(df_combined_1, aes(x = factor(D), y = X7, fill = factor(D))) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none") + 
  facet_wrap(~rho, ncol = 3)+xlab("dimension")+ylab("values")+ggtitle("Plot for 5 iterations")

# plot for 10 iterations
plot_1_10 <- ggplot(df_combined_1, aes(x = factor(D), y = X7, fill = factor(D))) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none") + 
  facet_wrap(~rho, ncol = 3)+xlab("iterations")+ylab("values")+ggtitle("Plot for 10 iterations")
