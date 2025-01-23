#' Summarize Outbreak Reconstruction
#'
#' This function summarizes the MCMC output.
#'
#' @param results List returned by run_mcmc().
#' @param burnin Proportion of MCMC iterations to discard as burnin. Defaults to 0.2.
#' @return A list consisting of a matrix of direct transmissions and their posterior probabilities, a matrix of indirect transmissions and their posterior probabilities, posterior samples of various parameters, and the time of the MRCA (if the tree is unrooted).
#' @export
summarize <- function(results, burnin = 0.2){
  n_reps <- length(results[[1]])
  names <- results[[3]] # Sequence names
  rooted <- results[[4]] # Is the tree rooted?
  s_max <- results[[5]]
  n_obs <- length(names) # Number of observed hosts
  indirect <- matrix(0, ncol = n_obs, nrow = n_obs)
  direct <- matrix(0, ncol = n_obs, nrow = n_obs)
  burnin <- n_reps * burnin + 1

  mus <- c()
  N_effs <- c()
  pis <- c()
  Rs <- c()

  # Time of infection for each observed host
  ts <- matrix(ncol = n_obs - 1, nrow = 0)

  if(!rooted){
    tmrca <- c()
  }
  for (i in burnin:n_reps) {
    mus <- c(mus, results[[2]][[i]]$mu)
    N_effs <- c(N_effs, results[[2]][[i]]$N_eff)
    pis <- c(pis, results[[2]][[i]]$pi)
    Rs <- c(Rs, results[[2]][[i]]$R)

    h <- results[[2]][[i]]$h
    n <- results[[2]][[i]]$n
    w <- sapply(results[[2]][[i]]$seq, length) - 1

    ts <- rbind(ts, sapply(results[[2]][[i]]$seq[2:n_obs], function(v){v[1]}))

    if(!rooted){
      tmrca <- c(tmrca, (results[[2]][[i]]$seq[[1]])[1])
    }

    # Most recent observed ancestor
    h_obs <- c()
    for (j in 2:n_obs) {
      h_obs[j] <- h[j]
      while (h_obs[j] > n_obs) {
        h_obs[j] <- h[h_obs[j]]
      }
    }

    trans <- cbind(h[2:n], 2:n)
    direct_trans <- trans[trans[,1] <= n_obs & trans[,2] <= n_obs & w[2:n] == 0, ]

    indirect_trans <- cbind(h_obs[2:n_obs], 2:n_obs)

    direct[direct_trans] <- direct[direct_trans] + 1
    indirect[indirect_trans] <- indirect[indirect_trans] + 1
  }
  indirect <- indirect / (n_reps - burnin + 1)
  direct <- direct / (n_reps - burnin + 1)

  # Name the matrices
  names[1] <- paste(names[1], "(root)")
  rownames(direct) <- names
  colnames(direct) <- names
  rownames(indirect) <- names
  colnames(indirect) <- names

  if(!rooted){
    direct <- direct[2:n_obs, 2:n_obs]
    indirect <- indirect[2:n_obs, 2:n_obs]
  }

  # Process times of infection
  mean_ts <- colMeans(ts) + s_max
  lower_ts <- apply(ts, 2, function(v){quantile(v, 0.025)}) + s_max
  upper_ts <- apply(ts, 2, function(v){quantile(v, 0.975)}) + s_max

  names(mean_ts) <- names[2:n_obs]
  names(lower_ts) <- names[2:n_obs]
  names(upper_ts) <- names[2:n_obs]



  out <- list(
    log_likelihood = results[[1]],
    direct_transmissions = direct,
    indirect_transmissions = indirect,
    mu = mus,
    N_effs = N_effs,
    pi = pis,
    R = Rs,
    mean_t_inf = mean_ts,
    lower_95_hpd_t_inf = lower_ts,
    upper_95_hpd_t_inf = upper_ts
  )

  # Visualization
  plot(out$log_likelihood, type = "l", main = "Log-Likelihood", xlab = "Iteration", ylab = "Log-Likelihood")
  hist(out$mu, xlab = "Value", main = "Evolution Rate (subs/site/day)")
  hist(out$N_eff, xlab = "Value", main = "Within-host Effective Population Size")
  hist(out$pi, xlab = "Value", main = "Sampling Rate")
  hist(out$R, xlab = "Value", main = "Reproductive Number")

  # If unrooted, return samples of time of MRCA
  if(!rooted){
    tmrca <- s_max + tmrca
    hist(tmrca, xlab = "Date", main = "Epidemic Start Date", breaks = "days")
    out <- c(out, list(time_of_MRCA = tmrca))
  }

  return(out)
}
