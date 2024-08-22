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
resample_seq <- function(mcmc, data, i, fix_latest_host, output = "all"){
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

  return(list(
    mcmc,
    log_p_new_old,
    log_p_old_new
  ))

}

## Resample the genotype for an unobserved host, or for an observed host with missing sites, based on (approximate) parsimony
genotype <- function(mcmc, data, i, output = "all", check_parsimony = F){

  js <- which(mcmc$h == i)

  # SNVs where we may need to make a change to get parsimony
  snvs <- unique(c(
    mcmc$m01[[i]],
    mcmc$m0y[[i]],
    mcmc$m1y[[i]],
    mcmc$m10[[i]],
    mcmc$mx0[[i]],
    mcmc$mxy[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m0y[js]),
    unlist(mcmc$m1y[js]),
    unlist(mcmc$m10[js]),
    unlist(mcmc$mx0[js]),
    unlist(mcmc$mxy[js]),
    unlist(mcmc$mx1[js])
  ))


  # If i is observed, the only positions that can change are those with missing data
  if(i <= data$n_obs){
    snvs <- intersect(snvs, data$snvs[[i]]$missing$call)
  }

  # Number of neighbors (h[i] and js)
  n_neighbors <- 1 + length(js)

  ## Loop over snvs
  log_p_new_old <- 0
  log_p_old_new <- 0

  for (snv in snvs) {

    ## Step 1: figure out if the snv is absent, isnv, or present in i
    from <- snv_status(mcmc, i, js, snv)

    ## Step 2: Figure out the "ideal" (most parsimonious) state of the snv
    if(from == "absent"){
      # How many times does snv appear as present in neighbors?
      n_present <- sum(
        c(mcmc$m10[[i]], unlist(mcmc$m01[js])) == snv
      )

      # How many times does snv appear as iSNV in neighbors?
      n_isnv <- sum(
        c(mcmc$mx0[[i]], unlist(mcmc$m0y[js])) == snv
      )

      n_absent <- n_neighbors - n_present - n_isnv
    }

    if(from == "isnv"){
      # How many times does snv appear as present in neighbors?
      n_present <- sum(
        c(mcmc$m1y[[i]], unlist(mcmc$mx1[js])) == snv
      )

      # How many times does snv appear as iSNV in neighbors?
      n_isnv <- sum(
        c(mcmc$mxy[[i]], unlist(mcmc$mxy[js])) == snv
      )

      n_absent <- n_neighbors - n_present - n_isnv
    }

    if(from == "present"){
      # How many times does snv appear as ABSENT in neighbors?
      n_absent <- sum(
        c(mcmc$m01[[i]], unlist(mcmc$m10[js])) == snv
      )

      # How many times does snv appear as iSNV in neighbors?
      n_isnv <- sum(
        c(mcmc$mx1[[i]], unlist(mcmc$m1y[js])) == snv
      )

      n_present <- n_neighbors - n_absent - n_isnv
    }

    ## In strict mode: parsimony given by whichever is more, present or absent, among neighbors.
    ## 50/50 in case of tie
    stat <- n_present + n_isnv / 2
    if(stat > n_neighbors / 2){
      ideal <- "present"
    }else if(stat < n_neighbors / 2){
      ideal <- "absent"
    }else{
      ideal <- "cointoss"
    }


    ## We now have all the info we need to make the move
    if(check_parsimony){
      if((from == "absent" & ideal == "present") | (from == "present" & ideal == "absent")){
        return(F)
      }
    }

    if(ideal == "cointoss"){
      log_p_new_old <- log_p_new_old + log(1/2) # Designator that the ideal is ambiguous
      ideal <- sample(c("present", "absent"), 1)
      log_p_old_new <- log_p_old_new + log(1/2)
    }

    if(output == "all"){
      new <- change_genotype(mcmc, data, snv, from, ideal, i, js, snv %in% observed)
      mcmc <- new[[1]]
      log_p_old_new <- log_p_old_new + new[[2]]
      if(new[[2]] != 0){
        stop("weird")
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
    return(list(
      mcmc,
      log_p_new_old,
      log_p_old_new
    ))
  }
}

# Check whether it's possible to create a new node
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





