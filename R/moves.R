# MIT License
#
# Copyright (c) 2023 Ivan Specht
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### MCMC moves

moves <- list()

## Update one of the w_i's by adding or subtracting either:
# rounded N(0, sqrt(delta_t * lambda_g / a_g))  (strategic) OR
# N(0, 3) (random)
moves$w <- function(mcmc, data){
  # Choose random host with ancestor
  if(data$rooted){
    i <- sample(2:mcmc$n, 1)
  }else{
    i <- sample(which(mcmc$h != 1), 1)
  }

  h <- mcmc$h[i]

  delta_t <- mcmc$t[i] - mcmc$t[h]

  # Proposal
  prop <- mcmc

  if(delta_t > 50){
    change <- round(rnorm(1, 0, sqrt(delta_t * mcmc$lambda_g / mcmc$a_g)))
  }else{
    change <- round(rnorm(1, 0, 3))
  }

  prop$w[i] <- mcmc$w[i] + change

  if(prop$w[i] < 0){
    return(mcmc)
  }

  ## Experimental version...
  if(data$experimental){
    prop$seq[[i]] <- c(
      prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
    )
  }


  if(data$experimental){
    hastings <-
      # P(new to old): draw the seq values in mcmc
      lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
      # P(old to new): draw the seq values in prop
      lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))
  }else{
    # With new coalescent:
    hastings <- 0
  }

  return(accept_or_reject(prop, mcmc, data, i, hastings))

}

## Update time of infection for a host on an edge leading to i
moves$seq <- function(mcmc, data){
  # Choose random host with ancestor
  if(data$rooted){
    i <- sample(2:mcmc$n, 1)
  }else{
    i <- sample(which(mcmc$h != 1), 1)
  }

  # If no intermediate hosts, nothing to do
  if(mcmc$w[i] == 0){
    return(mcmc)
  }else{

    # Proposal: resample all intermediate hosts' times of infection
    prop <- mcmc
    prop$seq[[i]] <- c(
      prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
    )

    return(accept_or_reject(prop, mcmc, data, i))

  }



}

## Update one of the t_i's using a N(0,1) proposal density if observed; N(0, 10) if not
moves$t <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(setdiff(2:mcmc$n, mcmc$external_roots), 1)
  # Proposal
  prop <- mcmc
  # Wider variance when it's unobserved, or it's the root of an unrooted tree
  prop$t[i] <- rnorm(1, mcmc$t[i], ifelse((i > data$n_obs) | (!data$rooted & mcmc$h[i] == 1), 10, 1))

  if(prop$t[i] < prop$t[prop$h[i]]){
    return(mcmc)
  }

  # Resample seq
  if(data$experimental){
    prop$seq[[i]] <- c(
      prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
    )

    hastings <-
      # P(new to old): draw the seq values in mcmc
      lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
      # P(old to new): draw the seq values in prop
      lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))

  }else{
    hastings <- 0
  }

  update <- c(i, which(mcmc$h ==i))


  return(accept_or_reject(prop, mcmc, data, update, hastings))
}

## Update t_i and w_i and w_j's simultaneously, where j is the child of i
# If recursive = T, also update time of infection for all ancestors of i
moves$w_t <- function(mcmc, data, recursive = F){
  # Choose random host with ancestor
  choices <- setdiff(2:mcmc$n, mcmc$external_roots)

  if(length(choices) == 0){
    return(mcmc)
  }else{
    i <- ifelse(length(choices) == 1, choices, sample(choices, 1))
    js <- which(mcmc$h == i)
    h <- mcmc$h[i]

    # Maximum time at which i can be infected
    if(i <= data$n_obs){
      max_t <- min(c(mcmc$t[js], data$s[i]))
    }else{
      max_t <- min(mcmc$t[js])
    }

    # Minimum time at which i can be infected
    min_t <- mcmc$t[h]

    # If recursive, SD is proportional to max_t - t[i]
    if(recursive){
      sd <- (max_t - mcmc$t[i]) / 5
    }else{
      sd <- (max_t - mcmc$t[i]) / 10
    }

    # Make the proposal of the change in time and generation
    delta_t <- rnorm(1, 0, sd)
    delta_w <- round(delta_t / (mcmc$a_g / mcmc$lambda_g))

    # If move results in negative evolutionary time...reject immediately
    if(mcmc$t[i] + delta_t > max_t){
      return(mcmc)
    }else{
      prop <- mcmc

      if(recursive){
        # All ancestors, not including root (case 1)
        is <- ancestry(mcmc$h, i)[-1]

        prop$t[is] <- mcmc$t[is] + delta_t
        prop$w[is[1]] <- mcmc$w[is[1]] + delta_w

        # All kids of all i's, excluding i's
        all_js <- setdiff(
          which(mcmc$h %in% is),
          is
        )

        prop$w[all_js] <- mcmc$w[all_js] - delta_w

        new_sd <- (max_t - prop$t[i]) / 5
        hastings <- dnorm(delta_t, 0, new_sd, log = T) - dnorm(delta_t, 0, sd, log = T)

        update <- c(is, all_js)

      }else{
        prop$t[i] <- mcmc$t[i] + delta_t
        prop$w[i] <- mcmc$w[i] + delta_w
        prop$w[js] <- mcmc$w[js] - delta_w
        hastings <- 0
        update <- c(i, js)
      }

      if(any(prop$w[update] < 0) | any(prop$t[update] <= prop$t[prop$h[update]])){
        return(mcmc)
      }

      # For each case in update: resample seq
      for (i in update) {
        prop$seq[[i]] <- c(
          prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
        )

        hastings <- hastings +
          # P(new to old): draw the seq values in mcmc
          lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
          # P(old to new): draw the seq values in prop
          lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))
    }
  }
}

## Update b using a N(0,0.01) proposal density
moves$b <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$b <- rnorm(1, mcmc$b, 0.1)
  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update b using a N(0,0.01) proposal density
moves$pi <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$pi <- rnorm(1, mcmc$pi, 0.05)
  if(prop$pi <= 0 | prop$pi >= 1){
    return(mcmc)
  }else{
    update <- 2:mcmc$n
    return(accept_or_reject(prop, mcmc, data, update))
  }

}

## Update lambda using a N(0,0.5^2) proposal density
moves$lambda <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$lambda <- rnorm(1, mcmc$lambda, 0.5)
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update a_g using a N(0,1) proposal density
moves$a_g <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$a_g <- rnorm(1, mcmc$a_g, 0.5)

  # Also update reproductive number and hence psi to maintain constant growth rate
  prop$psi <- prop$rho / (exp((prop$a_g / prop$lambda_g) * data$growth) + prop$rho) # second parameter, NBin offspring distribution (computed in terms of R0)

  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update a_s using a N(0,1) proposal density
moves$a_s <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$a_s <- rnorm(1, mcmc$a_s, 1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update mu using a N(0,1e-7) proposal density
moves$mu <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$mu <- rnorm(1, mcmc$mu, data$init_mu / 5)

  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update p using a N(0,1e-7) proposal density
moves$p <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$p <- rnorm(1, mcmc$p, data$init_mu / 10)
  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update v using a N(0,100) proposal density (rounded)
moves$v <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$v <- round(rnorm(1, mcmc$v, 1000))
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update rho using a N(0,0.1) proposal density
moves$rho <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$rho <- rnorm(1, mcmc$rho, 0.1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update R
moves$R <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$R <- rnorm(1, mcmc$R, 0.1)
  if(prop$R <= 0){
    return(mcmc)
  }else{
    prop$psi <- prop$rho / (prop$R + prop$rho)
    prop$e_lik <- e_lik(prop, data)
    prop$prior <- prior(prop)

    if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
      return(prop)
    }else{
      return(mcmc)
    }
  }
}

## Update psi using a N(0,0.1) proposal density
moves$psi <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$psi <- rnorm(1, mcmc$psi, 0.1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update genotype at (a) missing sites in observed host, or (b) all sites in unobserved host
#### May need to check hastings ratio here...
moves$genotype <- function(mcmc, data){

  # Choose random host with ancestor
  i <- sample(setdiff(2:mcmc$n, mcmc$external_roots), 1)
  js <- which(mcmc$h == i) # Children
  # Let h denote the ancestor of i; never used in computations

  # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
  interest <- unique(c(
    mcmc$m01[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m10[js]),
    unlist(mcmc$m1y[js]),
    mcmc$m10[[i]],
    mcmc$mx0[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m0y[js])
  ))

  # If i is observed, we can only change sites with missing data
  if(i <= data$n_obs){
    interest <- intersect(interest, data$snvs[[i]]$missing$call)
  }

  if(length(interest) == 0){
    return(mcmc)
  }else{

    # Proposal
    prop <- mcmc

    # Pick one SNV to update. We switch whether it exists or not in i.
    snv <- ifelse(length(interest) == 1, interest, sample(interest, 1))

    prop <- flip_genotype(prop, mcmc, i, js, snv)

    update <- c(i, js) # For which hosts must we update the genomic likelihood?
    return(accept_or_reject(prop, mcmc, data, update, hastings))

  }
}

### Topological moves

## Move the ancestor of a node one step upstream (towards tips) or one step downstream (towards root) onto next/previous tracked host
moves$h_step <- function(mcmc, data, upstream = TRUE, resample_t = FALSE, resample_w = FALSE){
  # Choose random host with ancestor
  if(resample_t){
    i <- sample(setdiff(2:mcmc$n, mcmc$external_roots), 1)
  }else{
    i <- sample(2:mcmc$n, 1)
  }

  h_old <- mcmc$h[i]

  # Proposal
  prop <- mcmc

  # Are we going upstream or downstream?
  #upstream <- runif(1) < 1/2

  if(upstream){
    # Who are the other children of h_old?
    children <- setdiff(which(mcmc$h == h_old), i)

    # What's the maximum time at which i can be infected?
    max_t <- get_max_t(mcmc, data, i)

    # Which ones have a compatible time of infection?
    children <- children[mcmc$t[children] < max_t]

    # Children not allowed to be frozen
    children <- setdiff(children, mcmc$external_roots)

    # If no valid children, reject
    # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
    if(length(children) == 0 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      return(mcmc)
    }else{

      # Pick one
      h_new <- ifelse(length(children) == 1, children, sample(children, 1))

      prop <- shift_upstream(prop, data, i, h_old, h_new, resample_t, resample_w)

      # What's the change in edge weight for i?
      change <- prop$w[i] - mcmc$w[i]



      # Update seq for i
      if(prop$w[i] < 0 | prop$t[i] <= prop$t[prop$h[i]]){
        return(mcmc)
      }

      #print(prop$w[i])

      if(data$experimental){
        prop$seq[[i]] <- c(
          prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
        )
      }

      update <- i
      # If updating t, need to also change genomic likelihood of children of i
      if(resample_t){
        update <- c(update, which(mcmc$h == i))
      }

      hastings <- log(length(children)) # P(new -> old): 1; P(old -> new): choose from among #[children] people to be h_new

      if(resample_t){
        hastings <- hastings - log(max_t - mcmc$t[h_old]) + # P(new -> old): uniform draw of time of infection
          log(max_t - prop$t[h_new]) # P(old -> new): uniform draw of time of infection
      }
      if(resample_w){ # Poisson draw for edge weights
        hastings <- hastings + dpois(mcmc$w[i], (mcmc$t[i] - mcmc$t[h_old]) * mcmc$lambda_g / mcmc$a_g, log = T) -
          dpois(prop$w[i], (prop$t[i] - prop$t[h_new]) * mcmc$lambda_g / mcmc$a_g, log = T)
      }

      if(data$experimental){
        hastings <- hastings +
          # P(new to old): draw the seq values in mcmc
          lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
          # P(old to new): draw the seq values in prop
          lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))

    }
  }else{
    if(
      h_old == 1 |
      (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2) |
      (!data$rooted & mcmc$h[h_old] == 1)
    ){
      # If no downstream move, reject
      # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
      # Also reject if tree is unrooted and we're moving a new node onto the root
      return(mcmc)
    }else{

      # New ancestor of i is ancestor's ancestor
      h_new <- mcmc$h[h_old]

      prop <- shift_downstream(prop, data, i, h_old, h_new, resample_t, resample_w)

      # What's the change in edge weight for i? (Positive)
      change <- mcmc$w[h_old] + 1

      ## Compute the number of possible children who could be chosen by i in the new config
      # Who are the other children of h_old?
      children <- setdiff(which(prop$h == h_new), i)

      # What's the maximum time at which i can be infected?
      max_t <- get_max_t(mcmc, data, i)

      # Which ones have a lesser time of infection than max_t?
      children <- children[prop$t[children] < max_t]

      # Children can't be frozen
      children <- setdiff(children, mcmc$external_roots)

      # What's the change in edge weight for i?
      change <- prop$w[i] - mcmc$w[i]

      # Update seq for i
      if(prop$w[i] < 0 | prop$t[i] <= prop$t[prop$h[i]]){
        return(mcmc)
      }
      if(data$experimental){
        prop$seq[[i]] <- c(
          prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
        )
      }

      update <- i
      # If updating t, need to also change genomic likelihood of children of i
      if(resample_t){
        update <- c(update, which(mcmc$h == i))
      }

      hastings <- -log(length(children)) # P(new -> old): choose from among #[children] people to be h_new; P(old -> new): 1

      if(resample_t){
        hastings <- hastings - log(max_t - mcmc$t[h_old]) + # P(new -> old): uniform draw of time of infection
          log(max_t - prop$t[h_new]) # P(old -> new): uniform draw of time of infection
      }
      if(resample_w){ # Poisson draw for edge weights
        hastings <- hastings + dpois(mcmc$w[i], (mcmc$t[i] - mcmc$t[h_old]) * mcmc$lambda_g / mcmc$a_g, log = T) -
          dpois(prop$w[i], (prop$t[i] - prop$t[h_new]) * mcmc$lambda_g / mcmc$a_g, log = T)
      }

      if(data$experimental){
        hastings <- hastings +
          # P(new to old): draw the seq values in mcmc
          lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
          # P(old to new): draw the seq values in prop
          lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))
    }
  }
}

## Global change in ancestor
# Importance sampling based on other nodes with similar additions / deletions
moves$h_global <- function(mcmc, data){
  # Choose random host with ancestor
  if(data$rooted){
    i <- sample(2:mcmc$n, 1)
  }else{
    i <- sample(which(mcmc$h != 1), 1)
  }
  h_old <- mcmc$h[i]

  # Nodes which are infected earlier than i
  choices <- which(mcmc$t < mcmc$t[i])

  choices <- setdiff(choices, mcmc$external_roots)

  if(!data$rooted){
    choices <- setdiff(choices, 1)
  }

  if(length(choices) == 0 | (h_old > data$n_obs & mcmc$d[h_old] <= 2)){
    return(mcmc)
  }else{

    # "Score" the choices: shared iSNV = +1
    scores <- softmax(sapply(choices, score, mcmc=mcmc, i=i), data$tau)

    h_new <- ifelse(length(choices) == 1, choices, sample(choices, 1, prob = scores))

    if(mcmc$t[i] < mcmc$t[h_new]){
      return(mcmc)
    }

    # Find the path from h_old to h_new
    route <- paths(mcmc$h, h_old, h_new)
    down <- route[[1]]
    up <- route[[2]]

    prop <- mcmc

    # If length of down < 2, don't need to do anything
    if(length(down) >= 2){
      for (j in 2:length(down)) {
        prop <- shift_downstream(prop, data, i, down[j-1], down[j])
      }
    }
    if(length(up) >= 2){
      for (j in 2:length(up)) {
        prop <- shift_upstream(prop, data, i, up[j-1], up[j])
      }
    }

    if(prop$w[i] < 0){
      return(mcmc)
    }



    update <- i

    rev_scores <- softmax(sapply(choices, score, mcmc=prop, i=i), data$tau)
    hastings <- log(rev_scores[which(choices == h_old)]) - log(scores[which(choices == h_new)])

    # Update seq
    prop$seq[[i]] <- c(
      prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
    )



    hastings <- hastings +
      # P(new to old): draw the seq values in mcmc
      lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
      # P(old to new): draw the seq values in prop
      lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))

    return(accept_or_reject(prop, mcmc, data, update, hastings))
  }
}

## The swap
## Switch h -> i -> j to
## h -> j -> i
moves$swap <- function(mcmc, data, exchange_children = FALSE){
  # Choose host with a parent and a grandparent
  choices <- which(mcmc$h != 1)
  choices <- setdiff(choices, mcmc$external_roots)

  if(length(choices) == 0){
    return(mcmc)
  }else{
    # Pick j
    j <- ifelse(length(choices) == 1, choices, sample(choices, 1))

    # Pick i
    i <- mcmc$h[j]

    # For exchange_children = FALSE
    # If i unobserved, i must have at least 3 children, because losing one
    # For exchange_children = TRUE
    # If i unobserved, j must have at least 2 children
    # If j unobserved, i must have at least 2 children (including j)

    if(
      (exchange_children == F & i > data$n_obs & mcmc$d[i] < 3) |
      (exchange_children == T & i > data$n_obs & mcmc$d[j] < 2) |
      (exchange_children == T & j > data$n_obs & mcmc$d[i] < 2)
    ){
      return(mcmc)
    }else{
      # Pick h
      h <- mcmc$h[i]

      # Children of each
      children_i <- setdiff(which(mcmc$h == i), j)
      children_j <- which(mcmc$h == j)

      # Update the state
      prop <- mcmc
      prop <- shift_downstream(prop, data, j, i, h) # Shift j from i onto h
      prop <- shift_upstream(prop, data, i, h, j) # Shift i from h onto j
      prop$w[j] <- mcmc$w[i] # Swapping edge weights
      prop$w[i] <- mcmc$w[j]
      prop$t[j] <- mcmc$t[i] # Swapping time of infection
      prop$t[i] <- mcmc$t[j]

      if(exchange_children){
        for (k in children_i) {
          prop <- shift_downstream(prop, data, k, i, j)
          prop$w[k] <- mcmc$w[k] # Keep edge weight the same
        }
        for (k in children_j) {
          prop <- shift_upstream(prop, data, k, j, i)
          prop$w[k] <- mcmc$w[k] # Keep edge weight the same
        }
      }


      update <- c(i, j, children_i, children_j) # For which hosts must we update the genomic likelihood?
      hastings <- 0

      if(any(prop$t[update] < prop$t[prop$h[update]]) | any(prop$w[update] < 0)){
        return(mcmc)
      }

      # For each case in update: resample seq
      for (i in update) {
        prop$seq[[i]] <- c(
          prop$t[i], sort(runif(prop$w[i], prop$t[prop$h[i]], prop$t[i]), decreasing = T)
        )

        hastings <- hastings +
          # P(new to old): draw the seq values in mcmc
          lfactorial(mcmc$w[i]) + mcmc$w[i] * log(1 / (mcmc$t[i] - mcmc$t[mcmc$h[i]])) -
          # P(old to new): draw the seq values in prop
          lfactorial(prop$w[i]) - prop$w[i] * log(1 / (prop$t[i] - prop$t[prop$h[i]]))
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))
    }
  }
}


## Create / remove a node
moves$create <- function(mcmc, data, create = T, upstream = T){
  # Are we creating or deleting an unobserved node?
  # if(runif(1) < 1/2){
  #   create <- T
  # }else{
  #   create <- F
  # }
  #
  # # Are we moving nodes onto the new node upstream or downstream?
  # if(runif(1) < 1/2){
  #   upstream <- T
  # }else{
  #   upstream <- F
  # }

  if(create){
    # Pick any node with an ancestor. (Note, some choices impossible, but this is okay!)
    j1 <- sample(2:mcmc$n, 1)
    h <- mcmc$h[j1]

    if(mcmc$w[j1] == 0 | (upstream & mcmc$d[h] == 1) | (!upstream & mcmc$d[j1] == 0)){
      return(mcmc)
    }else{

      # Who else are we attaching to i?
      if(upstream){
        kids <- setdiff(which(mcmc$h == h), j1)
      }else{
        kids <- which(mcmc$h == j1)
      }

      j2s <- kids[runif(length(kids)) < data$p_move]

      # If moving nobody, or moving everyone upstream off an unobserved node, or moving all but 0 or 1 downstream off an unobserved node, reject
      if(length(j2s) == 0 | (upstream & h > data$n_obs & length(j2s) == length(kids)) | (!upstream & j1 > data$n_obs & length(j2s) >= length(kids) - 1)){
        return(mcmc)
      }else{
        js <- c(j1, j2s)

        # How far upstream from h is the new node?
        # Maximum is min(w[js]) - 1 to preserve sum of all edge weights
        if(upstream){
          max_dist <- min(mcmc$w[js]) - 1
        }else{
          max_dist <- mcmc$w[j1] - 1
        }


        # If min is 0, can't make the move
        if(max_dist < 0){
          return(mcmc)
        }else{

          i <- mcmc$n + 1
          # Proposal
          prop <- mcmc

          # New edge weight coming into i, the new host
          if(!data$rooted & h == 1){
            dist <- rpois(1, mcmc$a_g / mcmc$lambda_g) # Distance in THIS CASE ONLY is generations back from j1
            prop$w[i] <- Inf # Initially, its ancestor is 1
          }else{
            dist <- sample(0:max_dist, 1)
            prop$w[i] <- dist # Distance here is generations after h
          }

          # Maximum time i could be infected
          max_t <- min(mcmc$t[js])



          ## Stick i onto h
          prop$n <- mcmc$n + 1
          prop$h[i] <- h

          if(!data$rooted & h == 1){
            prop$t[i] <- max_t - rgamma(1, (dist + 1) * mcmc$a_g / mcmc$lambda_g)
          }else{
            prop$t[i] <- mcmc$t[h] + (max_t - mcmc$t[h]) * rbeta(1, dist + 1, max_dist - dist + 1) # Weighted average
          }


          prop$d[i] <- 0
          prop$d[h] <- mcmc$d[h] + 1

          ## Initialize genotype for i. This is all changing, so we initialize as i loses all iSNVs to 0. Everything else stays the same
          prop$mx0[[i]] <- unique(c(
            mcmc$mx0[[j1]], mcmc$mxy[[j1]], mcmc$mx1[[j1]]
          ))
          prop$m01[[i]] <- character(0)
          prop$m10[[i]] <- character(0)
          prop$m0y[[i]] <- character(0)
          prop$m1y[[i]] <- character(0)
          prop$mx1[[i]] <- character(0)
          prop$mxy[[i]] <- character(0)

          ## Move all js onto i
          if(upstream){
            for (j in js) {
              prop <- shift_upstream(prop, data, j, h, i)
            }
          }else{
            prop <- shift_upstream(prop, data, j1, h, i)
            prop$w[j1] <- dist
            for (j2 in j2s) {
              prop <- shift_downstream(prop, data, j2, j1, i)
            }
          }

          ## Create new genotype for i
          geno <- genotype(prop, i, js, data$eps)
          prop <- geno[[1]]
          log_p <- geno[[2]]

          # Create seq for i and each j in js
          if(data$experimental){
            for (j in c(i, js)) {
              prop$seq[[j]] <- c(prop$t[j], sort(runif(prop$w[j], prop$t[prop$h[j]], prop$t[j]), decreasing = T))
            }
          }


          ## Hastings ratio
          hastings <- -length(j2s) * log(data$p_move) - (length(kids) - length(j2s)) * log(1 - data$p_move) - # P(old -> new): Probability of choosing kids to move onto i
            ifelse(
              (!data$rooted & h == 1),
              dpois(dist, mcmc$a_g / mcmc$lambda_g, log = T),
              -log(max_dist + 1) # P(old -> new): choose how far upstream
            ) -
            ifelse(
              (!data$rooted & h == 1),
              dgamma(max_t - prop$t[i], (dist + 1) * mcmc$a_g / mcmc$lambda_g, log = T),
              dbeta((prop$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), dist + 1, max_dist - dist + 1, log = T) # P(old -> new): beta density for t_i
            ) -
            log_p - # P(old -> new): probability of newly-created genotype for i
            log(sum(prop$h > data$n_obs, na.rm = T)) + # P(new -> old): pick a host with an unobserved ancestor
            log(mcmc$n - 1) # P(old -> new): pick a host with an ancestor

          if(data$experimental){
            # P(new to old): sample seq times for each j
            for (j in js) {
              hastings <- hastings + lfactorial(mcmc$w[j]) + mcmc$w[j]*log(1 / (mcmc$t[j] - mcmc$t[mcmc$h[j]]))
            }
            # P(old to new): sample seq times for each j and i
            for (j in c(i, js)) {
              hastings <- hastings - lfactorial(prop$w[j]) - prop$w[j]*log(1 / (prop$t[j] - prop$t[prop$h[j]]))
            }
          }

          update <- c(i, js)
          return(accept_or_reject(prop, mcmc, data, update, hastings))

        }
      }
    }
  }else{
    ## Delete a node by tucking it back inside its parent / child
    choices <- which(mcmc$h > data$n_obs)
    if(length(choices) == 0){
      return(mcmc)
    }else{
      j1 <- ifelse(length(choices) == 1, choices, sample(choices, 1))
      i <- mcmc$h[j1]
      h <- mcmc$h[i]
      js <- which(mcmc$h == i)
      j2s <- setdiff(js, j1)

      # upstream = T REVERSES upstream create move
      # upstream = F REVERSES !upstream create move

      if(
        (!upstream & any(mcmc$w[j2s] <= mcmc$w[j1])) |
        (!upstream & j1 %in% mcmc$external_roots) |
        (!data$rooted & h == 1 & upstream & length(j2s) > 0) # Because host 1 can't have multiple children
      ){
        return(mcmc)
      }else{

        prop <- mcmc

        # Put all children of i onto h or j1
        if(upstream){
          for (j in js) {
            prop <- shift_downstream(prop, data, j, i, h)
          }
        }else{
          if(any(mcmc$t[j2s] < mcmc$t[j1])){
            #print("bloop")
            return(mcmc)
          }
          for (j2 in j2s) {
            prop <- shift_upstream(prop, data, j2, i, j1)
          }
          # Put j1 onto h
          prop <- shift_downstream(prop, data, j1, i, h)
        }
        # Degree of h still needs to decrease by 1, because of deletion of i
        prop$d[h] <- prop$d[h] - 1

        ## For calculation of Hastings ratio:

        # How far upstream from h is the new node?
        # Maximum is min(w[js]) - 1 to preserve sum of all edge weights
        if(upstream){
          max_dist <- min(prop$w[js]) - 1
        }else{
          max_dist <- prop$w[j1] - 1
        }

        # Maximum time i could be infected
        max_t <- min(prop$t[js])

        # Who else are we attaching to i?
        if(upstream){
          kids <- setdiff(which(prop$h == h), j1)
        }else{
          kids <- which(prop$h == j1)
        }

        # Probability that genotype() returns the genotype of i in mcmc
        log_p <- genotype(mcmc, i, js, data$eps, comparison = T)


        # Create seq for each j in js
        if(data$experimental){
          for (j in js) {
            prop$seq[[j]] <- c(prop$t[j], sort(runif(prop$w[j], prop$t[prop$h[j]], prop$t[j]), decreasing = T))
          }
        }

        ## Hastings ratio

        # In the case of unrooted trees and h==1, hastings ratio computed based on generations and time between i and j1
        dist <- mcmc$w[i]

        hastings <-
          length(j2s) * log(data$p_move) + (length(kids) - length(j2s)) * log(1 - data$p_move) + # P(new -> old): Probability of choosing kids to move onto i
          ifelse(
            (!data$rooted & h == 1),
            dpois(dist, mcmc$a_g / mcmc$lambda_g, log = T),
            -log(max_dist + 1) # P(new -> old): choose how far upstream
          ) +
          ifelse(
            (!data$rooted & h == 1),
            dgamma(max_t - mcmc$t[i], (dist + 1) * mcmc$a_g / mcmc$lambda_g, log = T),
            dbeta((mcmc$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), dist + 1, max_dist - dist + 1, log = T) # P(new -> old): beta density for t_i
          ) +
          log_p - # P(new -> old): probability newly-created genotype for i equals genotype for i in "mcmc"
          log(prop$n - 2) + # P(new -> old): pick a host with an ancestor (-2 because haven't yet updated n)
          log(sum(mcmc$h > data$n_obs, na.rm = T)) # P(old -> new): pick a host with an unobserved ancestor

        if(data$experimental){
          # P(new to old): sample seq times for each j and i
          for (j in c(i, js)) {
            hastings <- hastings + lfactorial(mcmc$w[j]) + mcmc$w[j]*log(1 / (mcmc$t[j] - mcmc$t[mcmc$h[j]]))
          }
          # P(old to new): sample seq times for each j
          for (j in js) {
            hastings <- hastings - lfactorial(prop$w[j]) - prop$w[j]*log(1 / (prop$t[j] - prop$t[prop$h[j]]))
          }
        }

        ## Re-indexing: everyone above i steps down 1, i gets deleted
        prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1
        prop$n <- prop$n - 1
        prop$h <- prop$h[-i]
        prop$w <- prop$w[-i]
        if(data$experimental){
          prop$seq <- prop$seq[-i]
        }
        prop$t <- prop$t[-i]
        prop$m01 <- prop$m01[-i]
        prop$m10 <- prop$m10[-i]
        prop$m0y <- prop$m0y[-i]
        prop$m1y <- prop$m1y[-i]
        prop$mx0 <- prop$mx0[-i]
        prop$mx1 <- prop$mx1[-i]
        prop$mxy <- prop$mxy[-i]
        prop$d <- prop$d[-i]
        prop$g_lik <- prop$g_lik[-i]

        prop$external_roots[which(prop$external_roots > i)] <- prop$external_roots[which(prop$external_roots > i)] - 1

        # Update the js
        js[js > i] <- js[js > i] - 1

        update <- js
        return(accept_or_reject(prop, mcmc, data, update, hastings))
      }
    }
  }
}
