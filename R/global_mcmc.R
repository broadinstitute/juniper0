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

# Update global parameters

global_mcmc <- function(mcmc, data){

  # Move 1
  if(!data$fixed_mu){
    mcmc <- move_mu(mcmc, data)
  }

  # Move 2
  mcmc <- move_pi(mcmc, data)

  # Move 3
  mcmc <- move_R(mcmc, data)

  #mcmc <- move_N_eff(mcmc, data)


  return(mcmc)
}

## After making global moves, chop up the tree, and make one copy of "mcmc" for each subtree
breakdown <- function(mcmc, data, old_roots){

  subtrees <- chop(mcmc, data, old_roots)
  n_subtrees <- length(subtrees[[1]]) # number of subtrees, indexed 1, 2, ..., n_subtrees
  mcmcs <- list()
  datas <- list()

  for (i in 1:n_subtrees){
    # Initialize mcmc for each subtree
    mcmcs[[i]] <- mcmc
    # Initialize a vector of nodes in subtree i that are roots of OTHER subtrees
    mcmcs[[i]]$external_roots <- integer()
    # Initialize a vector of WHICH OTHER subtree each one is a root of
    mcmcs[[i]]$external_subtrees <- integer()
  }
  for (i in 1:n_subtrees) {

    # Root of subtree i
    root <- subtrees[[1]][i]

    # People in subtree i, excluding root
    cluster <- subtrees[[2]][[i]]

    if(root == 1){
      # Record which subtree is ancestral to subtree i
      mcmcs[[i]]$anc_cluster <- NA
    }else{
      for (j in 1:n_subtrees) {
        if(root %in% c(subtrees[[1]][j], subtrees[[2]][[j]]) & j != i){

          # Record which subtree is ancestral to subtree i
          mcmcs[[i]]$anc_cluster <- j

          # Record which node in subtree j is a root, and of which subtree it's a root
          mcmcs[[j]]$external_roots <- c(mcmcs[[j]]$external_roots, root)
          mcmcs[[j]]$external_subtrees <- c(mcmcs[[j]]$external_subtrees, i)
        }
      }
    }
  }

  for (i in 1:n_subtrees) {

    # Root of subtree i
    root <- subtrees[[1]][i]

    # People in subtree i, excluding root
    cluster <- subtrees[[2]][[i]]

    joined <- c(root, cluster)

    # To save memory: extract only necessary components of MCMC and data
    mcmcs[[i]]$n <- length(joined)
    mcmcs[[i]]$h <- mcmcs[[i]]$h[joined]
    mcmcs[[i]]$h <- match(mcmcs[[i]]$h, joined)
    mcmcs[[i]]$seq <- mcmcs[[i]]$seq[joined]

    # Care only about the first element of seq[[1]], analog of epidemic start time
    mcmcs[[i]]$seq[[1]] <- mcmcs[[i]]$seq[[1]][1]

    mcmcs[[i]]$subs$from <- mcmcs[[i]]$subs$from[joined]
    mcmcs[[i]]$subs$pos <- mcmcs[[i]]$subs$pos[joined]
    mcmcs[[i]]$subs$to <- mcmcs[[i]]$subs$to[joined]
    mcmcs[[i]]$tmu <- mcmcs[[i]]$tmu[joined]

    # Clear out mutations for root
    mcmcs[[i]]$subs$from[[1]] <- character(0)
    mcmcs[[i]]$subs$pos[[1]] <- integer(0)
    mcmcs[[i]]$subs$to[[1]] <- character(0)
    mcmcs[[i]]$tmu[[1]] <- numeric(0)

    mcmcs[[i]]$dropout <- mcmcs[[i]]$dropout[joined]
    mcmcs[[i]]$bot <- mcmcs[[i]]$bot[joined] # May contain some NAs, check if this causes issues

    mcmcs[[i]]$e_lik <- mcmcs[[i]]$e_lik[joined]
    mcmcs[[i]]$g_lik <- mcmcs[[i]]$g_lik[joined]
    mcmcs[[i]]$m_lik <- mcmcs[[i]]$m_lik[joined]

    # The roots of external clusters have been relabeled
    mcmcs[[i]]$external_roots <- match(mcmcs[[i]]$external_roots, joined)

    # Each likelihood at external roots resets to 0
    mcmcs[[i]]$e_lik[mcmcs[[i]]$external_roots] <- 0
    mcmcs[[i]]$g_lik[mcmcs[[i]]$external_roots] <- 0
    mcmcs[[i]]$m_lik[mcmcs[[i]]$external_roots] <- 0

    datas[[i]] <- data
    datas[[i]]$s <- datas[[i]]$s[joined]
    datas[[i]]$n_obs <- sum(joined <= data$n_obs)
    # If the root is unobserved, we treat it as observed anyway, because it cannot be added / deleted
    if(root > data$n_obs){
      datas[[i]]$n_obs <- datas[[i]]$n_obs + 1
    }

    datas[[i]]$snvs <- datas[[i]]$snvs[joined]

    datas[[i]]$vcf_present <- (data$vcf_present[joined])[1:datas[[i]]$n_obs]
    if(is.na(datas[[i]]$vcf_present[1])){
      datas[[i]]$vcf_present[1] <- F
    }

    datas[[i]]$observed_root <- (root <= data$n_obs)

    # Root always treated as fixed within subtree, UNLESS it's the root subtree
    if(root != 1){
      datas[[i]]$rooted <- T
    }

    # No need to store names
    datas[[i]]$names <- NULL

    # e_lik needs to be re-computed based on the smaller set
    # mcmcs[[i]]$e_lik <- e_lik(mcmcs[[i]], datas[[i]])
    # if(is.infinite(mcmcs[[i]]$e_lik)){
    #   stop("e_lik error in breakdown")
    # }

    ## SAFETY MODE

    # if(any(mcmcs[[i]]$g_lik != sapply(1:mcmcs[[i]]$n, g_lik, mcmc=mcmcs[[i]], data=datas[[i]]))){
    #
    #   bad <- which(mcmcs[[i]]$g_lik != sapply(1:mcmcs[[i]]$n, g_lik, mcmc=mcmcs[[i]], data=datas[[i]]))
    #   print(joined)
    #   print(datas[[i]]$vcf_present)
    #   print(bad)
    #   print(mcmcs[[i]]$g_lik[bad])
    #   print(sapply(1:mcmcs[[i]]$n, g_lik, mcmc=mcmcs[[i]], data=datas[[i]])[bad])
    #
    #
    #   stop("g_lik error in breakdown")
    # }


    # Record only observed hosts in each cluster
    mcmcs[[i]]$cluster <- joined[joined <= data$n_obs]

  }

  return(list(mcmcs, datas, subtrees[[1]])) # Last entry is the vector of roots, used later

}
