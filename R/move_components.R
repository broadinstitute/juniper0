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

### Building blocks of MCMC moves

## Each move returns a list of length 3: the new state of MCMC, log P(new -> old), log P(old -> new)

# Resample seq, the time of infection of a node i and the ancestors of i along the edge from h to i
# If fix_latest_host, we don't resample the time of infection for mcmc$seq[[i]][1]
# output can be "all" or "log_p_new_old"
resample_seq <- function(mcmc, data, i, fix_latest_host, output = "all", also_resample_tmu){
  # Ancestor of i
  h <- mcmc$h[i]

  # Max time of infection of i
  if(fix_latest_host){
    max_t <- mcmc$seq[[i]][1]
  }else{
    max_t <- get_max_t(mcmc, data, i)
  }

  # Min time of infection of i
  min_t <- mcmc$seq[[h]][1]

  # w, w_old, w_new always represent the number of hosts BETWEEN min_t and max_t, i.e. the ones whose times we resample, not including hosts at min_t or max_t
  # hence when fix_latest_host = TRUE, it's one less than the length of seq; otherwise length of seq

  # log P(new -> old)
  if(fix_latest_host){
    seq_old <- (mcmc$seq[[i]])[-1]
  }else{
    seq_old <- mcmc$seq[[i]]
  }

  w_old <- length(seq_old)
  log_p_new_old <- dw(w_old, mcmc, min_t, max_t, fix_latest_host) + dseq(seq_old, w_old, min_t, max_t, mcmc)

  if(also_resample_tmu){
    # Probability of sampling tmu
    log_p_new_old <- log_p_new_old + resample_tmu(mcmc, data, i, output = "log_p_new_old")
  }

  if(output == "log_p_new_old"){
    return(log_p_new_old)
  }

  # Resample number of hosts along the edge
  w_new <- rw(mcmc, min_t, max_t, fix_latest_host)
  seq_new <- rseq(w_new, min_t, max_t, mcmc)

  if(is.unsorted(rev(seq_new))){
    stop("ww")
  }

  if(length(seq_old) != w_old | length(seq_new) != w_new){
    stop("www")
  }

  log_p_old_new <- dw(w_new, mcmc, min_t, max_t, fix_latest_host) + dseq(seq_new, w_new, min_t, max_t, mcmc)

  # Resampled times of infection
  if(fix_latest_host){
    mcmc$seq[[i]] <- c(mcmc$seq[[i]][1], seq_new)
  }else{
    mcmc$seq[[i]] <- seq_new
    if(length(seq_new) == 0){
      stop("???")
    }
  }

  if(also_resample_tmu){
    # Resample tmu
    mcmc <- resample_tmu(mcmc, data, i)
    log_p_old_new <- log_p_old_new + mcmc[[3]]
    mcmc <- mcmc[[1]]
  }

  return(list(
    mcmc,
    log_p_new_old,
    log_p_old_new
  ))
}

# Resample times of mutations leading into i
resample_tmu <- function(mcmc, data, i, output = "all"){
  max_t <- mcmc$seq[[i]][1]

  h <- mcmc$h[[i]]

  # Min t is time of infection of h
  min_t <- mcmc$seq[[h]][1]

  # Number of mutations
  n_mut <- length(mcmc$subs$pos[[i]])

  ### TRYING SIMPLER VERISION: uniform draws
  if(n_mut != length(mcmc$tmu[[i]]) | max_t <= min_t){
    log_p_new_old <- -Inf
  }else{
    log_p_new_old <- n_mut * log(1 / (max_t - min_t))
  }

  mcmc$tmu[[i]] <- runif(n_mut, min_t, max_t)
  log_p_old_new <- n_mut * log(1 / (max_t - min_t))

  if(output == "log_p_new_old"){
    return(log_p_new_old)
  }else{
    return(list(
      mcmc,
      log_p_new_old,
      log_p_old_new
    ))
  }

}

## Resample the genotype for an unobserved host, or for an observed host with missing sites, based on (approximate) parsimony
genotype <- function(mcmc, data, i, output = "all", check_parsimony = F){

  js <- which(mcmc$h == i)

  # Positions on the genome where we may need to make a change to get parsimony
  pos <- c(
    mcmc$subs$pos[[i]],
    unlist(mcmc$subs$pos[js])
  )
  # The allele in i at these positions
  current <- c(
    mcmc$subs$to[[i]],
    unlist(mcmc$subs$from[js])
  )
  # The allele in neighbors at these positions
  neighbor <- c(
    mcmc$subs$from[[i]],
    unlist(mcmc$subs$to[js])
  )

  # If i is observed, the only positions that can change are those with missing data or iSNVs
  if(i <= data$n_obs){
    keep <- pos %in% data$pos[[i]]$missing | pos %in% data$pos[[i]]$isnv$pos
    pos <- pos[keep]
    current <- current[keep]
    neighbor <- neighbor[keep]
  }

  # Number of neighbors (h[i] and js)
  n_neighbors <- 1 + length(js)

  ## Loop over snvs
  log_p_new_old <- 0
  log_p_old_new <- 0

  # First get probabilities of tmu for old state
  # Unnecessary if this function is just to check parsimony
  if(!check_parsimony){
    for (j in c(i, js)) {
      log_p_new_old <- log_p_new_old + resample_tmu(mcmc, data, j, output = "log_p_new_old")
    }
  }

  for (p in unique(pos)) {

    # Counts of A, C, G, T in neighbors of i
    counts <- rep(0, 4)
    counts[1] <- sum(pos == p & neighbor == "A")
    counts[2] <- sum(pos == p & neighbor == "C")
    counts[3] <- sum(pos == p & neighbor == "G")
    counts[4] <- sum(pos == p & neighbor == "T")

    # Letter and index of current nucleotide
    letter_current <- (current[pos == p])[1]
    ind_current <- match(letter_current, c("A", "C", "G", "T"))

    # Whichever neighbors didn't exhibit a substitution have the current nucleotide
    counts[ind_current] <- counts[ind_current] + n_neighbors - sum(counts)

    if(sum(counts) != n_neighbors){
      stop("Parsimony table error")
    }

    # Parsimonious states are whichever counts equal the maximum of count
    parsimonious <- c("A", "C", "G", "T")[which(counts == max(counts))]

    # If i observed and has iSNV, parsimonious state must be one of the observed iSNV states
    if(i <= data$n_obs){
      if(p %in% data$snvs[[i]]$isnv$pos){
        # Index of iSNV in i
        ind_isnv <- match(p, data$pos[[i]]$isnv$pos)

        parsimonious <- intersect(
          parsimonious,
          c(data$snvs[[i]]$isnv$a1, data$snvs[[i]]$isnv$a2)
        )
      }
    }

    ## We now have all the info we need to make the move
    if(check_parsimony){
      if(!(letter_current %in% parsimonious)){
        return(F)
      }
    }

    if(!check_parsimony & output == "log_p_new_old"){
      if(!(letter_current %in% parsimonious)){
        stop("The current state is not parsimonious")
      }
    }

    # P(new to old): pick one of the (possibly several) parsimonious states
    log_p_new_old <- log_p_new_old + log(1 / length(parsimonious))

    if(output == "all" & !check_parsimony){

      letter_new <- sample(parsimonious, 1)
      # P(new to old): pick one of the (possibly several) parsimonious states
      log_p_old_new <- log_p_old_new + log(1 / length(parsimonious))

      ## Update genotype
      if(letter_current != letter_new){
        # Update subs[[i]]
        if(p %in% mcmc$subs$pos[[i]]){
          # What's the index of the substitution leading into i that has position "p"?
          ind <- match(p, mcmc$subs$pos[[i]])
          if(mcmc$subs$from[[i]][ind] == letter_new){
            # If the new letter at i matches the from, we're substituting from a nucleotide to itself
            mcmc$subs$from[[i]] <- mcmc$subs$from[[i]][-ind]
            mcmc$subs$pos[[i]] <- mcmc$subs$pos[[i]][-ind]
            mcmc$subs$to[[i]] <- mcmc$subs$to[[i]][-ind]
          }else{
            # Otherwise, we've created a new "to" with the same from
            mcmc$subs$to[[i]][ind] <- letter_new
          }
        }else{
          # If no mutation present, add it
          mcmc$subs$from[[i]] <- c(mcmc$subs$from[[i]], letter_current)
          mcmc$subs$pos[[i]] <- c(mcmc$subs$pos[[i]], p)
          mcmc$subs$to[[i]] <- c(mcmc$subs$to[[i]], letter_new)
        }

        # Update subs[[j]]
        for (j in js) {
          if(p %in% mcmc$subs$pos[[j]]){
            # What's the index of the substitution leading into j that has position "p"?
            ind <- match(p, mcmc$subs$pos[[j]])
            # If the new letter at i matches the to, we're substituting from a nucleotide to itself
            if(mcmc$subs$to[[j]][ind] == letter_new){
              mcmc$subs$from[[j]] <- mcmc$subs$from[[j]][-ind]
              mcmc$subs$pos[[j]] <- mcmc$subs$pos[[j]][-ind]
              mcmc$subs$to[[j]] <- mcmc$subs$to[[j]][-ind]
            }else{
              # Otherwise, we've created a new "from" with the same to
              mcmc$subs$from[[j]][ind] <- letter_new
            }
          }else{
            # If no mutation present, add it
            mcmc$subs$from[[j]] <- c(mcmc$subs$from[[j]], letter_new)
            mcmc$subs$pos[[j]] <- c(mcmc$subs$pos[[j]], p)
            mcmc$subs$to[[j]] <- c(mcmc$subs$to[[j]], letter_current)
          }
        }

        ## Update bot[[i]]
        if(i <= data$n_obs){
          if(data$vcf_present[i]){
            # What's the index of the relevant iSNV in i?
            ind <- match(p, data$snvs[[i]]$isnv$pos)
            if(letter_new == data$snvs[[i]]$isnv$a1[ind]){
              mcmc$bot[[i]][ind] <- TRUE
            }else if(letter_new == data$snvs[[i]]$isnv$a2[ind]){
              mcmc$bot[[i]][ind] <- FALSE
            }else{
              stop("The proposed allele is neither a1 nor a2")
            }
          }
        }
      }
    }
  }

  # If we got to the end of the loop and didn't violate parsimony, parsimony is TRUE
  if(check_parsimony){
    return(T)
  }

  if(output == "log_p_new_old"){
    return(log_p_new_old)
  }else{

    # Update tmu
    for (j in c(i, js)) {
      # Resample tmu
      mcmc <- resample_tmu(mcmc, data, j)
      log_p_old_new <- log_p_old_new + mcmc[[3]]
      mcmc <- mcmc[[1]]
    }

    return(list(
      mcmc,
      log_p_new_old,
      log_p_old_new
    ))
  }
}

# Get probability for different possible positions of a new node, up to constant of proportionality
p_pick_h <- function(mcmc, data){
  ## Select a node with probability proportional to its degree
  # Node degrees
  ds <- sapply(1:mcmc$n, function(x){sum(mcmc$h[2:mcmc$n] == x)})

  # Degree must be at least 2
  ds[ds < 2] <- 0

  # If unobserved, degree must be at least 3 (so that 2 can go onto i)
  if(mcmc$n > data$n_obs){
    ds[ds < 3 & 1:mcmc$n > data$n_obs] <- 0
  }

  # Probability of choosing an external root is 0
  ds[1:mcmc$n %in% mcmc$external_roots] <- 0

  return(ds)

}





