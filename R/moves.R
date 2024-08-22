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

#moves <- list()

## Update time of infection for a host on an edge leading to i
move_seq <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)

  if(i %in% mcmc$external_roots){
    fix_latest_host <- TRUE # If i is the root of a downstream cluster, can't resample time of infection of i
    update <- i
  }else{
    fix_latest_host <- FALSE
    update <- c(i, which(mcmc$h == i))
  }

  prop <- resample_seq(mcmc, data, i, fix_latest_host)

  hastings <- prop[[2]] - prop[[3]]
  prop <- prop[[1]]



  return(accept_or_reject(prop, mcmc, data, update, hastings))

}



## Update t_i and w_i and w_j's simultaneously, where j is the child of i
# If recursive = T, also update time of infection for all ancestors of i
move_w_t <- function(mcmc, data, recursive = F){
  # Choose random host with ancestor
  choices <- setdiff(2:mcmc$n, mcmc$external_roots)

  if(length(choices) == 0){
    return(mcmc)
  }

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

  # In recursive mode:
  if(recursive){

    # All ancestors, not including root (case 1)
    is <- ancestry(mcmc$h, i)[-1]

    max_delta <- Inf
    for (i in is) {
      # Maximum time at which i can be infected
      max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)

      # Maximum positive change for i
      max_delta <- min(max_delta, max_t - mcmc$seq[[i]][1])
    }

    # Children of is
    js <- setdiff(which(mcmc$h %in% is), is)

    # Earliest ancestor of i
    h <- is[1]

    # Maximum negative change for i
    min_delta <- mcmc$seq[[1]][1] - mcmc$seq[[h]][length(mcmc$seq[[h]])]

  }else{
    # Which i do we update
    is <- i

    # Maximum time at which i can be infected
    max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)

    # Maximum positive change for i
    max_delta <- max_t - mcmc$seq[[i]][1]

    # Children of is
    js <- which(mcmc$h == i)

    # Ancestor of i
    h <- mcmc$h[i]

    # Maximum negative change for i
    min_delta <- mcmc$seq[[h]][1] - mcmc$seq[[i]][length(mcmc$seq[[i]])]

  }

  # Change in time of infection for is
  delta <- runif(1, min_delta, max_delta)

  # Make proposal
  prop <- mcmc

  # Update seq
  for (i in is) {
    prop$seq[[i]] <- prop$seq[[i]] + delta
  }

  hastings <- 0
  update <- c(is, js)

  return(accept_or_reject(prop, mcmc, data, update, hastings))
}

## Update b using a N(0,0.01) proposal density
move_b <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$b <- rnorm(1, mcmc$b, 0.1)
  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update b using a N(0,0.01) proposal density
move_pi <- function(mcmc, data){
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
move_lambda <- function(mcmc, data){
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
move_a_g <- function(mcmc, data){
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
move_a_s <- function(mcmc, data){
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
move_mu <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$mu <- rnorm(1, mcmc$mu, data$init_mu / 5)

  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update p using a N(0,1e-7) proposal density
move_p <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$p <- rnorm(1, mcmc$p, data$init_mu / 10)
  update <- 2:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update))
}

## Update v using a N(0,100) proposal density (rounded)
move_v <- function(mcmc, data){
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
move_rho <- function(mcmc, data){
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
move_R <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$R <- rnorm(1, mcmc$R, 0.1)
  if(prop$R <= 0){
    return(mcmc)
  }else{
    #prop$psi <- prop$rho / (prop$R + prop$rho)
    prop$rho <- prop$R * prop$psi / (1 - prop$psi)
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
move_psi <- function(mcmc, data){
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

## Update genotype at (a) missing sites in observed host, or (b) all sites in unobserved host, based on parsimony
move_genotype <- function(mcmc, data){

  # Choose random host with ancestor
  i <- sample(setdiff(2:mcmc$n, mcmc$external_roots), 1)
  js <- which(mcmc$h == i) # Children

  # Proposal
  prop <- mcmc

  # Create new genotype for i
  prop <- genotype(mcmc, data, i)
  hastings <- prop[[2]] - prop[[3]]
  prop <- prop[[1]]

  update <- c(i, js) # For which hosts must we update the genomic likelihood?

  return(accept_or_reject(prop, mcmc, data, update, hastings, check_parsimony = c(mcmc$h[i], js)))
}

### Topological moves

## Move the ancestor of a node one step upstream (towards tips) or one step downstream (towards root) onto next/previous tracked host
move_h_step <- function(mcmc, data, upstream = TRUE){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)

  # Children of i
  js <- which(mcmc$h == i)

  if(i %in% mcmc$external_roots){
    fix_latest_host <- TRUE # If i is the root of a downstream cluster, can't resample time of infection of i
  }else{
    fix_latest_host <- FALSE
  }
  update <- c(i, js)

  # Previous ancestor of i
  h_old <- mcmc$h[i]

  # Proposal
  prop <- mcmc

  # Hastings ratio of this proposal
  hastings <- 0

  # Probability of proposing current genome and current seq
  hastings <- hastings +
    genotype(mcmc, data, i, output = "log_p_new_old") +
    resample_seq(mcmc, data, i, fix_latest_host, output = "log_p_new_old")

  if(upstream){
    # Who are the other children of h_old?
    children <- setdiff(which(mcmc$h == h_old), i)
    children <- setdiff(children, mcmc$external_roots)

    if(length(children) == 0 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      return(mcmc)
    }

    # Pick one
    h_new <- ifelse(length(children) == 1, children, sample(children, 1))

    # If times incompatible, return mcmc
    if(fix_latest_host){
      if(mcmc$seq[[i]][1] <= mcmc$seq[[h_new]][1]){
        return(mcmc)
      }
    }else{
      if(get_max_t(mcmc, data, i) <= mcmc$seq[[h_new]][1]){
        return(mcmc)
      }
    }

    # Update hastings ratio based on number of children to choose from. This is log P(old to new)
    hastings <- hastings - log(1 / length(children))

    # Update ancestor and initialize genetics
    prop$h[i] <- h_new
    prop <- update_genetics_upstream(prop, i, h_new)

    # Resample genotype
    prop <- genotype(prop, data, i)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    # Resample seq
    prop <- resample_seq(prop, data, i, fix_latest_host)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

  }else{

    ## Downstream move
    if(
      # If no downstream move, reject
      # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
      # Also reject if tree is unrooted and we're moving a new node onto the root
      h_old == 1 |
      (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2) |
      (!data$rooted & mcmc$h[h_old] == 1)
    ){
      return(mcmc)
    }

    # New ancestor of i is ancestor's ancestor
    h_new <- mcmc$h[h_old]

    # Update ancestor and initialize genetics
    prop$h[i] <- h_new
    prop <- update_genetics_downstream(prop, i, h_old)

    # Resample genotype
    prop <- genotype(prop, data, i)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    # Resample seq
    prop <- resample_seq(prop, data, i, fix_latest_host)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    # Who are the other children h_new? (For hastings ratio)
    children <- setdiff(which(prop$h == h_new), i)
    children <- setdiff(children, prop$external_roots)

    # P(new -> old): choose correct child of h_new
    hastings <- hastings + log(1 / length(children))
  }

  return(accept_or_reject(prop, mcmc, data, update, hastings = hastings, check_parsimony = c(h_old, h_new, js)))
}

## Global change in ancestor
# Importance sampling based on other nodes with similar additions / deletions
# If biassed = T, biassed towards places on tree with similar additions
# Else random location on tree
move_h_global <- function(mcmc, data, biassed = T){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)

  # Old ancestor
  h_old <- mcmc$h[i]

  # Children of i
  js <- which(mcmc$h == i)

  if(i %in% mcmc$external_roots){
    fix_latest_host <- TRUE # If i is the root of a downstream cluster, can't resample time of infection of i
  }else{
    fix_latest_host <- FALSE

  }

  update <- c(i, js)
  hastings <- 0

  # Probability of proposing current genome and current seq
  hastings <- hastings +
    genotype(mcmc, data, i, output = "log_p_new_old") +
    resample_seq(mcmc, data, i, fix_latest_host, output = "log_p_new_old")

  # Max time at which i can be infected
  max_t <- get_max_t(mcmc, data, i)

  # Nodes which are infected earlier than i
  ts <- sapply(mcmc$seq, function(v){v[1]})
  choices <- which(ts < max_t)
  choices <- setdiff(choices, mcmc$external_roots)

  if(!data$rooted){
    choices <- setdiff(choices, 1)
  }

  # Can't pick self
  choices <- setdiff(choices, i)

  # If no choices, or removing i from a node that's unobserved and has 2 kids, reject
  if(length(choices) == 0 | (h_old > data$n_obs & sum(mcmc$h[2:mcmc$n] == h_old) <= 2)){
    return(mcmc)
  }

  # "Score" the choices: shared iSNV = +1
  if(biassed){
    scores <- softmax(sapply(choices, score, mcmc=mcmc, i=i), data$tau)
  }else{
    scores <- rep(1/length(choices), length(choices))
  }


  h_new <- ifelse(length(choices) == 1, choices, sample(choices, 1, prob = scores))

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

  # Resample genotype
  prop <- genotype(prop, data, i)
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  # Resample seq
  prop <- resample_seq(prop, data, i, fix_latest_host)
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  new_max_t <- get_max_t(prop, data, i)
  if(new_max_t != max_t){
    stop("max_t error")
  }

  # Scores for proposing the move in the opposite direction
  if(biassed){
    rev_scores <- softmax(sapply(choices, score, mcmc=prop, i=i), data$tau)
  }else{
    rev_scores <- rep(1/length(choices), length(choices))
  }

  hastings <- hastings + log(rev_scores[which(choices == h_old)]) - log(scores[which(choices == h_new)])

  return(accept_or_reject(prop, mcmc, data, update, hastings, check_parsimony = c(h_old, h_new, js)))

}

## The swap
## Switch h -> i -> j to
## h -> j -> i
move_swap <- function(mcmc, data, exchange_children = FALSE){
  # Choose host with a parent and a grandparent
  choices <- which(mcmc$h != 1)
  choices <- setdiff(choices, mcmc$external_roots)

  if(length(choices) == 0){
    return(mcmc)
  }

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
    (exchange_children == F & i > data$n_obs & sum(mcmc$h[2:mcmc$n] == i) < 3) |
    (exchange_children == T & i > data$n_obs & sum(mcmc$h[2:mcmc$n] == j) < 2) |
    (exchange_children == T & j > data$n_obs & sum(mcmc$h[2:mcmc$n] == i) < 2)
  ){
    return(mcmc)
  }

  # Get h
  h <- mcmc$h[i]

  # Children of each
  children_i <- setdiff(which(mcmc$h == i), j)
  children_j <- which(mcmc$h == j)

  # Update the state
  prop <- mcmc
  prop <- shift_downstream(prop, data, j, i, h) # Shift j from i onto h
  prop <- shift_upstream(prop, data, i, h, j) # Shift i from h onto j
  prop$seq[[j]] <- mcmc$seq[[i]] # Swapping seq
  prop$seq[[i]] <- mcmc$seq[[j]]


  if(exchange_children){
    for (k in children_i) {
      prop <- shift_downstream(prop, data, k, i, j)
    }
    for (k in children_j) {
      prop <- shift_upstream(prop, data, k, j, i)
    }
  }


  update <- c(i, j, children_i, children_j) # For which hosts must we update the genomic likelihood?
  hastings <- 0

  prop_ts <- sapply(prop$seq, function(v){v[1]})

  if(any(prop_ts[update] <= prop_ts[prop$h[update]])){
    return(mcmc)
  }


  # Check that the move doesn't violate parsimony
  to_check <- c(h,i,j)
  if(exchange_children){
    to_check <- c(to_check, children_i, children_j)
  }

  return(accept_or_reject(prop, mcmc, data, update, hastings, check_parsimony = to_check))


}


## Create / remove a node
move_create <- function(mcmc, data, upstream = T){

  hastings <- 0

  if(upstream){

    # Probability of picking a node as ancestor to new node
    prob_h <- p_pick_h(mcmc, data)

    if(sum(prob_h) == 0){
      return(mcmc)
    }

    prob_h <- prob_h / sum(prob_h)

    # Sample h proportional to ds
    h <- sample(1:mcmc$n, 1, prob = prob_h)

    ## Term 1a
    # log P(old -> new): P(choose h)
    hastings <- hastings - log(prob_h[h])

    # Pick who moves
    # Number of people to move is uniform over valid choices
    min_n_move <- 2
    if(h > data$n_obs){
      max_n_move <- sum(mcmc$h[2:mcmc$n] == h) - 1
    }else{
      max_n_move <- sum(mcmc$h[2:mcmc$n] == h)
    }


  }else{

    # Choose any host with ancestor
    j1 <- sample(2:mcmc$n, 1)

    ## Term 1b
    hastings <- hastings - log(1 / (mcmc$n - 1)) # P(old to new): pick host j1

    # Old ancestor
    h <- mcmc$h[j1]

    # Children of j1
    ks <- which(mcmc$h == j1)

    # Minimum number of the ks to move onto the new node is always 1, since j1 itself is a child of i
    min_n_move <- 1

    # Max we can move is all of the ks if j1 is observed, or all but 2 if not
    if(j1 <= data$n_obs){
      max_n_move <- length(ks)
    }else{
      max_n_move <- length(ks) - 2
    }

    if(min_n_move > max_n_move){
      return(mcmc)
    }
  }

  n_move <- ifelse(min_n_move == max_n_move, min_n_move, sample(min_n_move:max_n_move, 1))

  ## Term 2
  hastings <- hastings - log(1 / (max_n_move - min_n_move + 1))


  if(upstream){
    # Who moves?
    js <- which(mcmc$h == h)

    ## Term 3a
    hastings <- hastings + lchoose(length(js), n_move) # P(old to new): pick who moves. (1 / length(js) choose n_move)
    js <- sample(js, n_move, replace = FALSE)
  }else{

    ## Term 3b
    hastings <- hastings + lchoose(length(ks), n_move) # P(old to new): pick who moves. (1 / length(ks) choose n_move)
    ks <- ifelse(length(ks) == 1, ks, sample(ks, n_move, replace = FALSE))
    js <- c(j1, ks) # Now, these are the children of i
  }

  # New node index
  i <- mcmc$n + 1
  # Proposal
  prop <- mcmc

  # Stick i onto h
  prop$n <- mcmc$n + 1
  prop$h[i] <- h

  ## Initialize genotype for i. This is all changing, so we initialize as i loses all iSNVs to 0. Everything else stays the same
  prop$mx0[[i]] <- unique(c(
    mcmc$mx0[[js[1]]], mcmc$mxy[[js[1]]], mcmc$mx1[[js[1]]]
  ))
  prop$m01[[i]] <- character(0)
  prop$m10[[i]] <- character(0)
  prop$m0y[[i]] <- character(0)
  prop$m1y[[i]] <- character(0)
  prop$mx1[[i]] <- character(0)
  prop$mxy[[i]] <- character(0)

  # Also initialize mcmc$isnv
  prop$isnv$call[[i]] <- character(0)
  prop$isnv$af[[i]] <- numeric(0)

  if(upstream){
    ## Move all js onto i
    for (j in js) {
      prop <- shift_upstream(prop, data, j, h, i)
    }
  }else{
    prop <- shift_upstream(prop, data, j1, h, i)
    for (k in ks) {
      prop <- shift_downstream(prop, data, k, j1, i)
    }
  }

  # Get genotype for i
  prop <- genotype(prop, data, i)
  # Term 4
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  # For the first element of prop$seq[[i]], choose uniformly random time up through prop$seq[[j]][1]
  max_t <- get_max_t(prop, data, i, fix_child_seq = FALSE)
  min_t <- prop$seq[[h]][1]

  prop$seq[[i]] <- runif(1, min_t, max_t)
  # Term 5
  hastings <- hastings - log(1 / (max_t - min_t))

  # Then, resample seq again, fixing its first element
  prop <- resample_seq(prop, data, i, fix_latest_host = TRUE)
  # Term 6
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  # Finally, resample seq for j. Also fix latest host, so we don't have to update g_lik for children of j
  for (j in js) {
    prop <- resample_seq(prop, data, j, fix_latest_host = TRUE)
    # Terms 7
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }

  ## For the move in the reverse direction: we need to pick i as the node to remove, and sample seq for js
  # Reverses Term 1 below
  hastings <- hastings + log(1 / (prop$n - data$n_obs)) # P(new to old): pick i as the host to remove

  # Reverses Terms 2 below
  for (j in js) {
    hastings <- hastings + resample_seq(mcmc, data, j, fix_latest_host = TRUE, output = "log_p_new_old")
  }

  update <- c(i, js)
  return(accept_or_reject(prop, mcmc, data, update, hastings, check_parsimony = c(h, js)))
}

# upstream = T here reverses move_create(... upstream = T) and vice versa

move_delete <- function(mcmc, data, upstream = T){

  # If no unobserved nodes, nothing to do here
  if(mcmc$n == data$n_obs){
    return(mcmc)
  }

  hastings <- 0

  # Who to delete
  choices <- (data$n_obs + 1):(mcmc$n)
  # Term 1
  hastings <- hastings - log(1 / length(choices)) # P(old to new): pick host to delete

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))
  h <- mcmc$h[i]
  js <- which(mcmc$h == i)

  if(!upstream){
    ts <- sapply(mcmc$seq[js], function(v){v[1]})
    j1 <- js[which.min(ts)]
    ks <- setdiff(js, j1)
  }

  ## Some initial updates to hastings ratio: time of mcmc$seq[[i]][1], rest of mcmc$seq[[i]], mcmc$seq[[j]], and genotype at i
  max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)
  min_t <- mcmc$seq[[h]][1]
  # Reverses Term 5 above
  hastings <- hastings + log(1 / (max_t - min_t))

  # Reverses Term 6 above
  hastings <- hastings + resample_seq(mcmc, data, i, fix_latest_host = TRUE, output = "log_p_new_old")

  # Reverse Terms 7 above
  for (j in js) {
    hastings <- hastings + resample_seq(mcmc, data, j, fix_latest_host = TRUE, output = "log_p_new_old")
  }

  # Genotype at i
  # Reverses Term 4 above
  hastings <- hastings + genotype(mcmc, data, i, output = "log_p_new_old")

  prop <- mcmc


  if(upstream){
    for (j in js) {
      prop <- shift_downstream(prop, data, j, i, h)
    }
  }else{
    for (k in ks) {
      prop <- shift_upstream(prop, data, k, i, j1)
    }
    prop <- shift_downstream(prop, data, j1, i, h)

  }


  # Resample seq for js (or j1 and ks)
  # Terms 2
  for (j in js) {
    prop <- resample_seq(prop, data, j, fix_latest_host = TRUE)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }

  ## Re-indexing: everyone above i steps down 1, i gets deleted
  prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1
  prop$n <- prop$n - 1
  prop$h <- prop$h[-i]
  prop$seq <- prop$seq[-i]

  prop$m01 <- prop$m01[-i]
  prop$m10 <- prop$m10[-i]
  prop$m0y <- prop$m0y[-i]
  prop$m1y <- prop$m1y[-i]
  prop$mx0 <- prop$mx0[-i]
  prop$mx1 <- prop$mx1[-i]
  prop$mxy <- prop$mxy[-i]

  prop$isnv$call <- prop$isnv$call[-i]
  prop$isnv$af <- prop$isnv$af[-i]

  prop$g_lik <- prop$g_lik[-i]
  prop$external_roots[which(prop$external_roots > i)] <- prop$external_roots[which(prop$external_roots > i)] - 1

  # Update the js
  js[js > i] <- js[js > i] - 1

  # Update h
  if(h > i){
    h <- h - 1
  }

  if(!upstream){
    ks[ks > i] <- ks[ks > i] - 1
    if(j1 > i){
      j1 <- j1 - 1
    }
  }

  if(upstream){
    ## For hastings ratio: probability of picking h and js
    prob_h <- p_pick_h(prop, data)

    if(sum(prob_h) == 0){
      print(mcmc$h)
      print(prop$h)
      print(prop$external_roots)
      print("Whaaaa?")
    }

    prob_h <- prob_h / sum(prob_h)

    # log P(new -> old): P(choose h)
    # Reverses Term 1a
    hastings <- hastings + log(prob_h[h])

    # Pick who moves
    # Number of people to move is uniform over valid choices
    min_n_move <- 2
    if(h > data$n_obs){
      max_n_move <- sum(prop$h[2:prop$n] == h) - 1
    }else{
      max_n_move <- sum(prop$h[2:prop$n] == h)
    }

    # Reverse Terms 2, 3
    hastings <- hastings + log(1 / (max_n_move - min_n_move + 1)) # P(new -> old): pick number of people to move
    hastings <- hastings - lchoose(sum(prop$h[2:prop$n] == h), length(js)) # P(new to old): pick who moves
  }else{


    # Probability of choosing j1
    # Reverses term 1b
    hastings <- hastings + log(1 / (prop$n - 1))

    # Minimum number of the ks to move onto the new node is always 1, since j1 itself is a child of i
    min_n_move <- 1

    # Max we can move is all of the ks if j1 is observed, or all but 2 if not
    if(j1 <= data$n_obs){
      max_n_move <- sum(prop$h[2:prop$n] == j1)
    }else{
      max_n_move <- sum(prop$h[2:prop$n] == j1) - 2
    }

    # Reverse Terms 2, 3
    hastings <- hastings + log(1 / (max_n_move - min_n_move + 1)) # P(new -> old): pick number of people to move
    hastings <- hastings - lchoose(sum(prop$h[2:prop$n] == j1), length(ks)) # P(new to old): pick who moves

  }

  update <- js

  return(accept_or_reject(prop, mcmc, data, update, hastings, check_parsimony = c(h, js)))

}
