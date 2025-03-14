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
move_seq <- function(mcmc, data, also_resample_tmu){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)
  update_e <- mcmc$h[i]
  update_g <- mcmc$h[i]
  update_m <- integer(0)

  if(data$split_bottlenecks){
    update_g <- c(update_g, i)
  }

  prop <- resample_seq(mcmc, data, i, fix_latest_host = TRUE, also_resample_tmu = also_resample_tmu)

  hastings <- prop[[2]] - prop[[3]]
  prop <- prop[[1]]

  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m, hastings))

}

move_tmu <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)

  prop <- resample_tmu(mcmc, data, i)

  hastings <- prop[[2]] - prop[[3]]
  prop <- prop[[1]]

  update_g <- mcmc$h[i]
  update_e <- integer(0)
  update_m <- integer(0)

  if(data$split_bottlenecks){
    update_g <- c(update_g, i)
  }

  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m, hastings))
}



## Update t_i and w_i and w_j's simultaneously, where j is the child of i
# If recursive = T, also update time of infection for all ancestors of i
move_w_t <- function(mcmc, data, recursive = F){

  if(data$rooted){
    # Choose random host with ancestor
    choices <- setdiff(2:mcmc$n, mcmc$external_roots)
  }else{
    # Choose random host
    choices <- setdiff(1:mcmc$n, mcmc$external_roots)
  }


  if(length(choices) == 0){
    return(mcmc)
  }

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))


  # In recursive mode:
  if(recursive){

    if(data$rooted){
      # All ancestors, not including root (case 1)
      is <- ancestry(mcmc$h, i)[-1]
    }else{
      # All ancestors, going down to the root
      is <- ancestry(mcmc$h, i)
    }

    max_delta <- Inf
    for (i in is) {
      # Maximum time at which i can be infected
      max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)

      # Maximum positive change for i
      max_delta <- min(max_delta, max_t - mcmc$seq[[i]][1])
    }

    # Children of is
    js <- setdiff(which(mcmc$h %in% is), is)

    if(data$rooted){
      # Maximum negative change for i
      min_delta <- mcmc$seq[[1]][1] - mcmc$seq[[is[1]]][1]
    }

    # Else the mechanism for picking delta is different; see below

  }else{
    # Which i do we update
    is <- i

    # Maximum time at which i can be infected
    max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)

    # Maximum positive change for i
    max_delta <- max_t - mcmc$seq[[i]][1]

    # Children of is
    js <- which(mcmc$h == i)

    # Ancestor of i; can be NA if i == 1
    h <- mcmc$h[i]

    if(data$rooted | is != 1){
      # Maximum negative change for i
      min_delta <- mcmc$seq[[h]][1] - mcmc$seq[[i]][1]
    }
  }

  # Hastings ratio for infection times of js
  hastings <- 0

  # P(new to old): pick seq for each j and oldest ancestor in is

  # If updating tmu and seq for 1, only care about seq and tmu for js
  # otherwise also must update seq and tmu for the earliest of the is
  if(!(1 %in% is)){
    who_update_seq_tmu <- c(js, is[1])
  }else{
    who_update_seq_tmu <- js
  }

  for (j in who_update_seq_tmu) {
    hastings <- hastings + resample_seq(mcmc, data, j, fix_latest_host = TRUE, output = "log_p_new_old", also_resample_tmu = TRUE)
  }

  if(!(1 %in% is)){
    # Change in time of infection for is
    delta <- runif(1, min_delta, max_delta)
  }else{
    # If delta negative, take its mean to be equal to max_delta/2
    delta <- rdelta(max_delta/2, max_delta)
    # P(old to new): sample delta (non-symmetric here)
    hastings <- hastings - ddelta(delta, max_delta/2, max_delta, log = TRUE)
    # After move: max_delta decreases by delta
    # P(new to old):
    hastings <- hastings + ddelta(-delta, (max_delta - delta)/2, max_delta - delta, log = TRUE)
  }

  # Make proposal
  prop <- mcmc

  # Update seq and tmu for is (just adding a fixed value here)
  # Doesnt affect i == 1
  for (i in is) {
    prop$seq[[i]] <- prop$seq[[i]] + delta
    prop$tmu[[i]] <- prop$tmu[[i]] + delta
  }

  # Update seq for js and (sometimes) is[1] randomly
  for (j in who_update_seq_tmu) {
    prop <- resample_seq(prop, data, j, fix_latest_host = TRUE, also_resample_tmu = TRUE)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }

  update <- is
  if(!recursive){
    if(!is.na(h)){
      update <- c(update, h)
    }
  }else{
    update <- unique(c(1, update))
  }

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update_e = update, update_g = update_g, update_m = update, hastings))
}

## Update b using a N(0,0.01) proposal density
move_pi <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$pi <- mcmc$pi * exp(rnorm(1, 0, 0.1))
  if(prop$pi <= 0 | prop$pi >= 1){
    return(mcmc)
  }else{

    if(data$ongoing){
      prop$wbar <- wbar(data$t_min, 0, prop$R * prop$psi / (1 - prop$psi), 1 - prop$psi, prop$pi, prop$a_g, 1 / prop$lambda_g, prop$a_s, 1 / prop$lambda_s, 0.1)
    }



    update_e <- 1:mcmc$n
    update_g <- integer(0)
    update_m <- integer(0)
    return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m)) # Can parallelize this because not running within a subtree
  }
}

## Update mu
move_mu <- function(mcmc, data){
  prop <- mcmc
  prop$mu <- mcmc$mu * exp(rnorm(1, 0, 0.1))

  if(prop$mu <= 0){
    return(mcmc)
  }

  update_e <- integer(0)
  update_g <- 1:mcmc$n
  update_m <- 1:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m))
}

## Update N_eff
move_N_eff <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$N_eff <- mcmc$N_eff * exp(rnorm(1, 0, 0.1))

  if(prop$N_eff <= 0){
    return(mcmc)
  }

  update_e <- integer(0)
  update_g <- 1:mcmc$n
  update_m <- integer(0)
  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m))
}

## Update both simultaneously
move_mu_N_eff <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  scale <- exp(rnorm(1, 0, 0.05))
  prop$mu <- mcmc$mu * scale
  prop$N_eff <- mcmc$N_eff * scale


  update_e <- integer(0)
  update_g <- 1:mcmc$n
  update_m <- 1:mcmc$n
  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m))
}

## Update R
move_R <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$R <- mcmc$R * exp(rnorm(1, 0, 0.1))
  if(prop$R <= 0){
    return(mcmc)
  }

  if(data$ongoing){
    prop$wbar <- wbar(data$t_min, 0, prop$R * prop$psi / (1 - prop$psi), 1 - prop$psi, prop$pi, prop$a_g, 1 / prop$lambda_g, prop$a_s, 1 / prop$lambda_s, 0.1)
  }


  update_e <- 1:mcmc$n
  update_g <- integer(0)
  update_m <- integer(0)

  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m))

}

## Update psi using a N(0,0.1) proposal density
move_psi <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$psi <- mcmc$psi * exp(rnorm(1, 0, 0.1))

  if(prop$psi < 0 | prop$psi > 1){
    return(mcmc)
  }

  update_e <- 1:mcmc$n
  update_g <- integer(0)
  update_m <- integer(0)
  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m))
}

## Update genotype at (a) missing sites in observed host, or (b) all sites in unobserved host, based on parsimony
move_genotype <- function(mcmc, data){

  ## Choose host i for which to update genotype
  if(data$rooted){
    # Choose random host with ancestor
    choices <- setdiff(2:mcmc$n, mcmc$external_roots)
  }else{
    # Choose random host
    choices <- setdiff(1:mcmc$n, mcmc$external_roots)
  }

  if(length(choices) == 0){
    return(mcmc)
  }

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

  js <- which(mcmc$h == i) # Children of i

  # Proposal
  prop <- mcmc

  # Create new genotype for i
  prop <- genotype(mcmc, data, i)
  hastings <- prop[[2]] - prop[[3]]
  prop <- prop[[1]]

  # Resample times of mutations
  if(i == 1){
    update <- i
  }else{
    update <- c(i, mcmc$h[i]) # For which hosts must we update the genomic likelihood?
  }

  if(i == 1){
    to_check <- c(i, js)
  }else{
    to_check <- c(mcmc$h[i], i, js)
  }

  update_e <- integer(0)
  update_g <- update
  update_m <- update

  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update_e, update_g, update_m, hastings, check_parsimony = to_check))
}

### Topological moves

## Move the ancestor of a node one step upstream (towards tips) or one step downstream (towards root) onto next/previous tracked host
move_h_step <- function(mcmc, data, upstream = TRUE){
  # Choose random host with ancestor
  choices <- 2:mcmc$n

  if(length(choices) == 0){
    return(mcmc)
  }

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

  # Children of i
  js <- which(mcmc$h == i)

  if(i %in% mcmc$external_roots){
    fix_latest_host <- TRUE # If i is the root of a downstream cluster, can't resample time of infection of i
  }else{
    fix_latest_host <- FALSE
  }

  # Previous ancestor of i
  h_old <- mcmc$h[i]

  # Proposal
  prop <- mcmc

  # Hastings ratio of this proposal
  hastings <- 0

  # Probability of proposing current genome and current seq
  hastings <- hastings +
    resample_seq(mcmc, data, i, fix_latest_host, output = "log_p_new_old", also_resample_tmu = (i %in% mcmc$external_roots))

  # Only resample genotype if i isn't root of another cluster
  if(!(i %in% mcmc$external_roots)){
    hastings <- hastings + genotype(mcmc, data, i, output = "log_p_new_old")
  }

  if(upstream){
    # Who are the other children of h_old?
    children <- setdiff(which(mcmc$h == h_old), i)
    children <- setdiff(children, mcmc$external_roots)

    if(
      length(children) == 0 |
      (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2) |
      (h_old == 1 & length(which(mcmc$h == h_old)) <= 2 & !data$observed_root)
    ){
      return(mcmc)
    }

    # Pick one
    h_new <- ifelse(length(children) == 1, children, sample(children, 1))

    if(h_new %in% mcmc$external_roots){
      return(mcmc)
    }

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
    prop <- shift(prop, data, i, h_old, h_new, upstream = T)

    # Resample seq
    prop <- resample_seq(prop, data, i, fix_latest_host, also_resample_tmu = (i %in% mcmc$external_roots))
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    # Resample genotype
    if(!(i %in% mcmc$external_roots)){
      prop <- genotype(prop, data, i)
      hastings <- hastings - prop[[3]]
      prop <- prop[[1]]
    }


  }else{

    ## Downstream move
    if(
      # If no downstream move, reject
      # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
      h_old == 1 |
      (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)
    ){
      return(mcmc)
    }

    # New ancestor of i is ancestor's ancestor
    h_new <- mcmc$h[h_old]

    # Update ancestor and initialize genetics
    prop <- shift(prop, data, i, h_old, h_new, upstream = F)

    # Resample seq
    prop <- resample_seq(prop, data, i, fix_latest_host, also_resample_tmu = (i %in% mcmc$external_roots))
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    # Resample genotype (which also resamples tmu)
    if(!(i %in% mcmc$external_roots)){
      prop <- genotype(prop, data, i)
      hastings <- hastings - prop[[3]]
      prop <- prop[[1]]
    }

    # Who are the other children h_new? (For hastings ratio)
    children <- setdiff(which(prop$h == h_new), i)
    children <- setdiff(children, prop$external_roots)

    # P(new -> old): choose correct child of h_new
    hastings <- hastings + log(1 / length(children))
  }

  update <- c(h_old, h_new, i)

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update, update_g, update, hastings = hastings, check_parsimony = c(h_old, h_new, i, js)))
}

## Global change in ancestor
# Importance sampling based on other nodes with similar additions / deletions
# If biassed = T, biassed towards places on tree with similar additions
# Else random location on tree
move_h_global <- function(mcmc, data, biassed = T){
  # Choose random host with ancestor
  choices <- 2:mcmc$n
  if(length(choices) == 0){
    return(mcmc)
  }
  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

  # Old ancestor
  h_old <- mcmc$h[i]

  # Children of i
  js <- which(mcmc$h == i)

  if(i %in% mcmc$external_roots){
    fix_latest_host <- TRUE # If i is the root of a downstream cluster, can't resample time of infection of i
  }else{
    fix_latest_host <- FALSE
  }

  hastings <- 0

  # Probability of proposing current genome and current seq
  hastings <- hastings +
    resample_seq(mcmc, data, i, fix_latest_host, output = "log_p_new_old", also_resample_tmu = (i %in% mcmc$external_roots))

  if(!(i %in% mcmc$external_roots)){
    hastings <- hastings + genotype(mcmc, data, i, output = "log_p_new_old")
  }

  # Max time at which i can be infected
  if(fix_latest_host){
    max_t <- mcmc$seq[[i]][1]
  }else{
    max_t <- get_max_t(mcmc, data, i)
  }


  # Nodes which are infected earlier than i
  ts <- sapply(mcmc$seq, function(v){v[1]})
  choices <- which(ts < max_t)
  choices <- setdiff(choices, mcmc$external_roots)

  # Can't pick self
  choices <- setdiff(choices, i)

  # If no choices, or removing i from a node that's unobserved and has 2 kids, reject
  if(
    length(choices) == 0 |
    (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2) |
    (h_old == 1 & length(which(mcmc$h == h_old)) <= 2 & !data$observed_root)
  ){
    return(mcmc)
  }

  # "Score" the choices: shared iSNV = +1
  if(biassed){
    lscores <- lsoftmax(sapply(choices, score, mcmc=mcmc, i=i), data$tau)

  }else{
    lscores <- rep(-log(length(choices)), length(choices))
  }

  h_new <- ifelse(length(choices) == 1, choices, sample(choices, 1, prob = exp(lscores)))

  # Find the path from h_old to h_new
  route <- paths(mcmc$h, h_old, h_new)
  down <- route[[1]]
  up <- route[[2]]

  prop <- mcmc

  # If length of down < 2, don't need to do anything
  if(length(down) >= 2){
    for (j in 2:length(down)) {
      prop <- shift(prop, data, i, down[j-1], down[j], upstream = F)
    }
  }
  if(length(up) >= 2){
    for (j in 2:length(up)) {
      prop <- shift(prop, data, i, up[j-1], up[j], upstream = T)
    }
  }

  # Resample seq
  prop <- resample_seq(prop, data, i, fix_latest_host, also_resample_tmu = (i %in% mcmc$external_roots))
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  # Resample genotype
  if(!(i %in% mcmc$external_roots)){
    prop <- genotype(prop, data, i)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }

  # Scores for proposing the move in the opposite direction
  if(biassed){
    rev_lscores <- lsoftmax(sapply(choices, score, mcmc=prop, i=i), data$tau)
  }else{
    rev_lscores <- rep(-log(length(choices)), length(choices))
  }

  hastings <- hastings + rev_lscores[which(choices == h_old)] - lscores[which(choices == h_new)]

  update <- c(h_old, h_new, i)

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update, update_g, unique(c(update, down, up, js)), hastings, check_parsimony = unique(c(down, up, i, js))))

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

  hastings <- 0
  # Resample tmu
  for (k in c(i, j, children_i, children_j)) {
    hastings <- hastings + resample_tmu(mcmc, data, k, output = "log_p_new_old")
    if(is.infinite(hastings) | is.nan(hastings)){
      stop("wut")
    }
  }

  # Update the state
  prop <- mcmc

  prop <- shift(prop, data, j, i, h, upstream = F) # Shift j from i onto h

  prop <- shift(prop, data, i, h, j, upstream = T) # Shift i from h onto j

  prop$seq[[j]] <- mcmc$seq[[i]] # Swapping seq
  prop$seq[[i]] <- mcmc$seq[[j]]



  if(exchange_children){
    for (k in children_i) {
      prop <- shift(prop, data, k, i, j, upstream = F)
    }
    for (k in children_j) {
      prop <- shift(prop, data, k, j, i, upstream = T)
    }
  }


  check_timing <- c(i, j, children_i, children_j) # For which hosts must we update the genomic likelihood?
  check_timing_obs <- check_timing[which(check_timing <= data$n_obs)]

  prop_ts <- sapply(prop$seq[check_timing], function(v){v[length(v)]})
  prop_ts_anc <- sapply(prop$seq[prop$h[check_timing]], function(v){v[1]})

  if(any(prop_ts <= prop_ts_anc)){
    return(mcmc)
  }

  if(any(sapply(prop$seq[check_timing_obs], function(v){v[1]}) >= data$s[check_timing_obs])){
    return(mcmc)
  }

  # Resample tmu
  for (k in c(i, j, children_i, children_j)) {
    prop <- resample_tmu(prop, data, k)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }


  # Check that the move doesn't violate parsimony
  to_check <- c(h,i,j)
  if(exchange_children){
    to_check <- c(to_check, children_i, children_j)
  }

  update <- c(h, i, j)

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, children_i, children_j))
  }

  return(accept_or_reject(prop, mcmc, data, update, update_g, update, hastings, check_parsimony = to_check))

}


## Create / remove a node
## If upstream = TRUE, we have a biassed option, in which we create a node that tries to destroy homoplasies

move_create <- function(mcmc, data, upstream = T, biassed = F){

  if(!upstream & biassed){
    stop("A biassed create move is only allowed when upstream = TRUE")
  }

  hastings <- 0

  if(upstream){

    if(biassed){
      # All positions that change across the whole dataset and how many times, minus 1
      n_changes <- table(unlist(mcmc$subs$pos)) - 1

      # If all 0, we can't make a biassed move, so just make the normal move where we pick two random hosts
      # Also make this move with probability 5%, to ensure reversibility
      # HASTINGS CODE: Terms B1
      if(all(n_changes == 0) | runif(1) < 0.05){
        if(!all(n_changes == 0)){
          hastings <- hastings - log(0.05)
        }

        # Sample any two js
        hastings <- hastings + lchoose(mcmc$n - 1, 2)
        js <- sample(2:mcmc$n, 2, replace = F)

      }else{

        hastings <- hastings - log(0.95)

        # Find the positions with nonzero homoplasies
        n_changes <- n_changes[n_changes > 0]

        # Position numbers
        pos_nos <- as.numeric(names(n_changes))

        # For each position number, get a list of the js that have a mutation at that position
        js_by_pos <- list()
        for (p in 1:length(pos_nos)) {
          js_by_pos[[p]] <- which(sapply(mcmc$subs$pos, function(v){pos_nos[p] %in% v}))
        }

        # Sample a position for which to create new host with that mutation
        if(length(js_by_pos) == 1){
          js <- js_by_pos[[1]]
        }else{
          js <- js_by_pos[[sample(1:length(js_by_pos), 1)]]
        }

        if(length(js) > 2){
          # Then resample to get any two of them. We will only move two people tops
          hastings <- hastings + lchoose(length(js), 2)
          js <- sample(js, 2, replace = FALSE)
        }

        # Probability of picking js
        p_pick_js <- 0
        for (p in 1:length(js_by_pos)) {
          if(all(js %in% js_by_pos[[p]])){
            p_pick_js <- p_pick_js +
              (1 / length(js_by_pos)) * # Probability of picking index p
              (1 / choose(length(js_by_pos[[p]]), 2)) # Probability of picking the specific js
          }
        }
        hastings <- hastings - log(p_pick_js)
      }



      if(length(js) == 1){
        stop("Not moving enough js...")
      }

      # ancestors of js, if unobserved, must have degree at least 3
      # or if the ancestors are the same, need degree at least 4
      if(mcmc$h[js[1]] == mcmc$h[js[2]]){
        if(mcmc$h[js[1]] > data$n_obs){
          if(length(which(mcmc$h == mcmc$h[js[1]])) < 4){
            return(mcmc)
          }
        }
      }else{
        # Ancestors are different
        if(mcmc$h[js[1]] > data$n_obs){
          if(length(which(mcmc$h == mcmc$h[js[1]])) < 3){
            return(mcmc)
          }
        }
        if(mcmc$h[js[2]] > data$n_obs){
          if(length(which(mcmc$h == mcmc$h[js[2]])) < 3){
            return(mcmc)
          }
        }
      }

      # Get MRCA of all the js
      h <- mrca(mcmc$h, js)

      # Also reject if mrca equals any of the js
      if(h %in% js){
        return(mcmc)
      }

    }else{
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

  if(!biassed){
    n_move <- ifelse(min_n_move == max_n_move, min_n_move, sample(min_n_move:max_n_move, 1))

    ## Term 2
    if(is.nan(log(1/(max_n_move - min_n_move + 1)))){
      print(min_n_move)
      print(max_n_move)
      stop("bad1")
    }

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
  }



  # New node index
  i <- mcmc$n + 1
  # Proposal
  prop <- mcmc

  # Stick i onto h
  prop$n <- mcmc$n + 1
  prop$h[i] <- h

  ## Initialize genotype for i. Initialized to same as h
  prop$subs$from[[i]] <- character(0)
  prop$subs$pos[[i]] <- integer(0)
  prop$subs$to[[i]] <- character(0)

  prop$tmu[[i]] <- numeric(0)

  prop$dropout[[i]] <- 1:data$n_bases

  prop$bot[[i]] <- list(NULL)



  if(upstream){

    if(biassed){
      # If biassed, we first need to get each j down onto h
      # We also need to update the hastings ratio based on the number of choices for where to move j at each step, for the move in the reverse direction

      # Intermediate hosts on the way down to h
      intermediates <- integer(0)
      for (j in js) {

        # HASTINGS CODE: Reverses terms B2 below

        # Probability that we stop moving j when we reach its current ancestor
        # See below comment on this probability
        n_kids <- length(setdiff(which(prop$h == prop$h[j]), c(i, js)))
        hastings <- hastings + log(1 / (n_kids + 1)) # Probability we stop moving here on the way back up

        while(prop$h[j] != h){

          intermediates <- c(intermediates, prop$h[j])
          prop <- shift(prop, data, j, prop$h[j], prop$h[prop$h[j]], F)

          # P(new -> old), so the factor is positive
          # Idea is that if you have c children, the probability of moving to any one of them is 1 / (c+1), and the probability of staying put is also 1 / (c+1)
          # Can't move one j onto another j, or onto i

          n_kids <- length(setdiff(which(prop$h == prop$h[j]), c(i, js)))
          hastings <- hastings + log(1 / (n_kids + 1))
        }
      }
    }

    ## Move all js onto i
    for (j in js) {
      prop <- shift(prop, data, j, h, i, upstream = T)
    }
  }else{
    prop <- shift(prop, data, j1, h, i, upstream = T)
    for (k in ks) {
      prop <- shift(prop, data, k, j1, i, upstream = F)
    }
  }

  # For each j, the above may introduce mutations at dropout sites.
  # Correct this here
  for (j in js) {
    keep <- which(!(prop$subs$pos[[j]] %in% prop$dropout[[j]]))
    prop$subs$from[[j]] <- prop$subs$from[[j]][keep]
    prop$subs$pos[[j]] <- prop$subs$pos[[j]][keep]
    prop$subs$to[[j]] <- prop$subs$to[[j]][keep]
  }


  # If unobserved root and has degree 1, reject
  # Can be made more efficient by rejecting earlier
  if(length(which(prop$h == 1)) < 2 & !data$observed_root){
    return(mcmc)
  }

  # For the first element of prop$seq[[i]], choose uniformly random time up through prop$seq[[j]][1]
  max_t <- get_max_t(prop, data, i, fix_child_seq = FALSE)
  min_t <- prop$seq[[h]][1]

  #print(max_t - min_t)

  prop$seq[[i]] <- runif(1, min_t, max_t)
  # Term 5
  hastings <- hastings - log(1 / (max_t - min_t))



  # Then, resample seq again, fixing its first element
  prop <- resample_seq(prop, data, i, fix_latest_host = TRUE, also_resample_tmu = FALSE)
  # Term 6
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]

  # Finally, resample seq for j. Also fix latest host, so we don't have to update g_lik for children of j
  for (j in js) {
    prop <- resample_seq(prop, data, j, fix_latest_host = TRUE, also_resample_tmu = FALSE) # The tmu for each j will be updated below
    # Terms 7
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]

    if(prop$seq[[j]][length(prop$seq[[j]])] <= prop$seq[[i]][1]){
      stop("Bad resampling of seq")
    }
  }

  # Get genotype for i
  prop <- genotype(prop, data, i)
  # Term 4
  hastings <- hastings - prop[[3]]
  prop <- prop[[1]]



  ## For the move in the reverse direction: we need to pick i as the node to remove, and sample seq for js
  # Reverses Term 1 below
  choices <- setdiff((data$n_obs + 1):(prop$n), prop$external_roots)
  if(biassed){
    # P(new to old): pick i as the host to remove, among unobserved hosts with degree 2
    ds <- sapply(choices, function(x){length(which(prop$h == x))})
    hastings <- hastings + log(1 / sum(ds == 2))
  }else{
    hastings <- hastings + log(1 / length(choices)) # P(new to old): pick i as the host to remove
  }

  # Reverses Terms 2 below
  for (j in js) {
    hastings <- hastings + resample_seq(mcmc, data, j, fix_latest_host = TRUE, output = "log_p_new_old", also_resample_tmu = TRUE)
  }

  if(upstream){
    if(biassed){
      update <- unique(c(h, i, mcmc$h[js]))
    }else{
      update <- c(h, i)
    }
  }else{
    update <- c(h, i, j1)
  }

  if(biassed){
    to_check <- unique(c(h, i, js, intermediates))
  }else{
    to_check <- c(h, i, js)
  }

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update, update_g, unique(c(update, to_check)), hastings, check_parsimony = to_check))
}

# upstream = T here reverses move_create(... upstream = T) and vice versa

move_delete <- function(mcmc, data, upstream = T, biassed = F){

  if(!upstream & biassed){
    stop("A biassed delete move is only allowed when upstream = TRUE")
  }

  # If no unobserved nodes, nothing to do here
  if(mcmc$n == data$n_obs){
    return(mcmc)
  }

  hastings <- 0

  # Who to delete
  choices <- setdiff((data$n_obs + 1):(mcmc$n), mcmc$external_roots)

  # Auto-reject if no valid choices
  if(length(choices) == 0){
    return(mcmc)
  }

  if(biassed){
    ds <- sapply(choices, function(x){length(which(mcmc$h == x))})
    choices <- choices[which(ds == 2)]

    if(length(choices) == 0){
      return(mcmc)
    }
  }

  # Term 1
  hastings <- hastings - log(1 / length(choices)) # P(old to new): pick host to delete

  i <- ifelse(length(choices) == 1, choices, sample(choices, 1))
  h <- mcmc$h[i]
  js <- which(mcmc$h == i)





  if(!upstream){
    ts <- sapply(mcmc$seq[js], function(v){v[1]})

    # If multiple times are the minimum, move is impossible, so reject
    if(sum(ts == min(ts)) > 1){
      return(mcmc)
    }
    j1 <- js[which.min(ts)]
    # print("top")
    # print(i)
    # print(js)
    # print(j1)
    # print(which(mcmc$h == j1))
    # print(mcmc$external_roots)

    # j1 will gain some kids, so cannot be an external root
    if(j1 %in% mcmc$external_roots){
      return(mcmc)
    }

    ks <- setdiff(js, j1)

  }

  ## Some initial updates to hastings ratio: time of mcmc$seq[[i]][1], rest of mcmc$seq[[i]], mcmc$seq[[j]], and genotype at i
  max_t <- get_max_t(mcmc, data, i, fix_child_seq = FALSE)
  min_t <- mcmc$seq[[h]][1]
  # Reverses Term 5 above
  hastings <- hastings + log(1 / (max_t - min_t))

  # Reverses Term 6 above
  hastings <- hastings + resample_seq(mcmc, data, i, fix_latest_host = TRUE, output = "log_p_new_old", also_resample_tmu = FALSE)

  # Reverse Terms 7 above
  for (j in js) {
    hastings <- hastings + resample_seq(mcmc, data, j, fix_latest_host = TRUE, output = "log_p_new_old", also_resample_tmu = FALSE)
  }

  # Genotype at i
  # Reverses Term 4 above
  hastings <- hastings + genotype(mcmc, data, i, output = "log_p_new_old")

  prop <- mcmc


  if(upstream){
    for (j in js) {
      prop <- shift(prop, data, j, i, h, upstream = F)
    }
  }else{
    for (k in ks) {
      prop <- shift(prop, data, k, i, j1, upstream = T)
    }
    prop <- shift(prop, data, j1, i, h, upstream = F)

  }

  # For each j, the above may introduce mutations at dropout sites.
  # Correct this here
  # Idea being we keep the genotype in each j the same, except in dropout sites, where it (of course) matches the ancestor's genotype
  for (j in js) {
    keep <- which(!(prop$subs$pos[[j]] %in% prop$dropout[[j]]))
    prop$subs$from[[j]] <- prop$subs$from[[j]][keep]
    prop$subs$pos[[j]] <- prop$subs$pos[[j]][keep]
    prop$subs$to[[j]] <- prop$subs$to[[j]][keep]
  }

  ## Re-indexing: everyone above i steps down 1, i gets deleted
  prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1
  prop$n <- prop$n - 1
  prop$h <- prop$h[-i]
  prop$seq <- prop$seq[-i]

  prop$subs$from <- prop$subs$from[-i]
  prop$subs$pos <- prop$subs$pos[-i]
  prop$subs$to <- prop$subs$to[-i]
  prop$tmu <- prop$tmu[-i]
  prop$dropout <- prop$dropout[-i]
  prop$bot <- prop$bot[-i]

  # "$bot" does not exist for unobserved hosts

  prop$e_lik <- prop$e_lik[-i]
  prop$g_lik <- prop$g_lik[-i]
  prop$m_lik <- prop$m_lik[-i]
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

  # If biassed, we need to step up each of the js
  if(biassed){

    # Intermediate hosts crossed for which we check parsimony
    intermediates <- integer(0)

    for (j in js) {
      done <- FALSE
      while (!done) {
        # Number of kids
        # Remember, we can't move one j onto another, or onto i!
        kids <- setdiff(which(prop$h == prop$h[j]), js)

        # HASTINGS CODE: Terms B2

        # With probability 1 / (length(kids) + 1), the move ends here
        if(runif(1) < 1 / (length(kids) + 1)){
          hastings <- hastings - log(1 / (length(kids) + 1))
          done <- TRUE
        }else{
          # Otherwise, pick a kid to move onto
          hastings <- hastings - log(1 / (length(kids) + 1))
          if(length(kids) == 1){
            kid <- kids
          }else{
            kid <- sample(kids, 1)
          }
          if(prop$h[kid] != prop$h[j]){
            stop("Kid error")
          }

          intermediates <- c(intermediates, kid)
          prop <- shift(prop, data, j, prop$h[j], kid, upstream = T)
        }
      }

      # If time too early, reject
      if(prop$seq[[j]][length(prop$seq[[j]])] <= prop$seq[[prop$h[j]]][1]){
        return(mcmc)
      }
    }
  }



  # Degree of h can't be too small
  if(h > data$n_obs){
    if(length(which(prop$h == h)) < 2){
      return(mcmc)
    }
  }

  # If unobserved root and has degree 1, reject
  # Can be made more efficient by rejecting earlier
  if(length(which(prop$h == 1)) < 2 & !data$observed_root){
    return(mcmc)
  }

  # If external root has degree >0, reject
  if(any(prop$external_roots %in% prop$h)){
    return(mcmc)
  }


  # Resample seq for js (or j1 and ks)
  # Terms 2
  for (j in js) {
    prop <- resample_seq(prop, data, j, fix_latest_host = TRUE, also_resample_tmu = TRUE)
    hastings <- hastings - prop[[3]]
    prop <- prop[[1]]
  }





  ## Now, to finish calculating the probability of the proposal back to the original state...

  if(upstream){
    if(biassed){
      # HASTINGS CODE: Reverses Term B1

      # All positions that change across the whole dataset and how many times, minus 1
      n_changes <- table(unlist(prop$subs$pos)) - 1

      # If all 0, move is random choice of 2 js
      if(all(n_changes == 0)){
        # P(new to old): sample js
        hastings <- hastings - lchoose(prop$n - 1, 2)
      }else{

        # Either we picked the js randomly or strategically
        # First let's get the probability of the strategic case

        # Find the positions with nonzero homoplasies
        n_changes <- n_changes[n_changes > 0]

        # Position numbers
        pos_nos <- as.numeric(names(n_changes))

        # For each position number, get a list of the js that have a mutation at that position
        js_by_pos <- list()
        for (p in 1:length(pos_nos)) {
          js_by_pos[[p]] <- which(sapply(prop$subs$pos, function(v){pos_nos[p] %in% v}))
        }

        # Probability of picking js
        p_pick_js <- 0
        for (p in 1:length(js_by_pos)) {
          if(all(js %in% js_by_pos[[p]])){
            p_pick_js <- p_pick_js +
              (1 / length(js_by_pos)) * # Probability of picking index p
              (1 / choose(length(js_by_pos[[p]]), 2)) # Probability of picking the specific js
          }
        }

        # P(new to old): with 95% probability, move according to how to delete homoplasies; otherwise random
        hastings <- hastings + log(0.95 * p_pick_js + 0.05 / lchoose(prop$n - 1, 2))
      }

    }else{
      ## Unbiassed upstream move

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


      if(is.nan(log(1/(max_n_move - min_n_move + 1)))){
        print(min_n_move)
        print(max_n_move)
        stop("bad2")
      }

      # Reverse Terms 2, 3
      hastings <- hastings + log(1 / (max_n_move - min_n_move + 1)) # P(new -> old): pick number of people to move
      hastings <- hastings - lchoose(sum(prop$h[2:prop$n] == h), length(js)) # P(new to old): pick who moves
    }
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

    if(is.nan(log(1/(max_n_move - min_n_move + 1)))){
      print(min_n_move)
      print(max_n_move)
      stop("bad3")
    }

    # Reverse Terms 2, 3
    #print(hastings)
    hastings <- hastings + log(1 / (max_n_move - min_n_move + 1)) # P(new -> old): pick number of people to move
    # print("bottom")
    # print(j1)
    # print(which(prop$h == j1))
    # print(prop$external_roots)
    # print(data$n_obs)
    # print(hastings)
    hastings <- hastings - lchoose(sum(prop$h[2:prop$n] == j1), length(ks)) # P(new to old): pick who moves
    #print(hastings)

  }

  if(upstream){
    if(biassed){
      update <- unique(c(h, prop$h[js]))
    }else{
      update <- h
    }
  }else{
    update <- c(h, j1)
  }

  if(biassed){
    to_check <- unique(c(h, js, intermediates))
  }else{
    to_check <- c(h, js)
  }

  update_g <- update
  if(data$split_bottlenecks){
    update_g <- unique(c(update_g, js))
  }

  return(accept_or_reject(prop, mcmc, data, update, update_g, unique(c(update, to_check)), hastings, check_parsimony = to_check))
}
