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

# Schedule of moves within a subtree

local_mcmc <- function(mcmc, data){

  #set.seed(210)
  res <- list()
#

     # mcmc <- mcmcs[[j]]
     # data <- datas[[j]]
     # data$n_local <- 1

  for (r in 1:data$n_local) {
    # Move 11
    mcmc <- moves$w(mcmc, data)

    # Move 12
    mcmc <- moves$t(mcmc, data)

    # Move 13
    mcmc <- moves$w_t(mcmc, data)
    mcmc <- moves$w_t(mcmc, data, recursive = T)

    # Move 14
    mcmc <- moves$h_step(mcmc, data)

    # Move 15
    mcmc <- moves$h_step(mcmc, data, upstream = F)

    # Move 16
    mcmc <- moves$h_step(mcmc, data, resample_t = T)

    # Move 17
    mcmc <- moves$h_step(mcmc, data, upstream = F, resample_t = T)

    # Move 18
    mcmc <- moves$h_step(mcmc, data, resample_t = T, resample_w = T)

    # Move 19
    mcmc <- moves$h_step(mcmc, data, upstream = F, resample_t = T, resample_w = T)

    # Move 20
    mcmc <- moves$h_global(mcmc, data)

    # Move 21
    mcmc <- moves$swap(mcmc, data)

    # Move 22
    mcmc <- moves$swap(mcmc, data, exchange_children = T)

    # Move 23
    mcmc <- moves$genotype(mcmc, data)

    # Move 24
    mcmc <- moves$create(mcmc, data)

    # Move 25
    mcmc <- moves$create(mcmc, data, create = F)

    # Move 26
    mcmc <- moves$create(mcmc, data, upstream = F)

    # Move 27
    mcmc <- moves$create(mcmc, data, create = F, upstream = F)

    # Append new results
    if(r %% data$sample_every == 0){
      res <- c(res, list(mcmc))
    }
  }

  return(res)

}


## Join together results calculated in parallel across subtrees
amalgamate <- function(all_res, mcmcs, datas, mcmc, data){

  # Number of samples for each subtree
  n_samples <- length(all_res[[1]])

  # Number of subtrees
  n_subtrees <- length(mcmcs)


  # If we didn't break up the tree, nothing to do here!
  if(n_subtrees == 1){
    return(all_res[[1]])
  }else{

    ## Loop through each sample and produce an amalgamated MCMC state:

    # Create a list to store the amalgamated results
    res <- list()
    for (i in 1:n_samples) {

      # Get the root cluster of each cluster
      anc_clusters <- c()
      roots <- c()
      for (j in 1:n_subtrees) {
        roots[j] <- all_res[[j]][[i]]$root
        anc_clusters[j] <- all_res[[j]][[i]]$anc_cluster
      }

      # First determine who the unobserved hosts are in each cluster, so that they may be re-indexed
      displacement <- 0
      mappings <- list()
      for (j in 1:n_subtrees) {
        mappings[[j]] <- c(all_res[[j]][[i]]$root, all_res[[j]][[i]]$cluster)
        # Delete unobserved hosts in the initial cluster--they will be renamed
        mappings[[j]] <- mappings[[j]][mappings[[j]] <= data$n_obs]
        # How many unobserved hosts are there at this iteration?
        n_unobs <- all_res[[j]][[i]]$n - datas[[j]]$n_obs
        if(n_unobs > 0){
          mappings[[j]] <- c(mappings[[j]], (data$n_obs + displacement + 1):(data$n_obs + displacement + n_unobs))
        }
        displacement <- displacement + n_unobs
      }



      ## Now, fill an amalgamated mcmc with info from the correct subtrees

      # Initialization doesn't really matter; we will initialize to the previous amalgamated mcmc
      # Correct lengths of entries, to avoid carrying over extraneous information
      mcmc$n <- data$n_obs + displacement
      mcmc$h <- mcmc$h[1:mcmc$n]
      mcmc$w <- mcmc$w[1:mcmc$n]
      mcmc$t <- mcmc$t[1:mcmc$n]
      mcmc$m01 <- mcmc$m01[1:mcmc$n]
      mcmc$m10 <- mcmc$m10[1:mcmc$n]
      mcmc$m0y <- mcmc$m0y[1:mcmc$n]
      mcmc$m1y <- mcmc$m1y[1:mcmc$n]
      mcmc$mx0 <- mcmc$mx0[1:mcmc$n]
      mcmc$mx1 <- mcmc$mx1[1:mcmc$n]
      mcmc$mxy <- mcmc$mxy[1:mcmc$n]
      mcmc$d <- mcmc$d[1:mcmc$n]
      mcmc$g_lik <- mcmc$g_lik[1:mcmc$n]
      mcmc$root <- NULL
      mcmc$cluster <- NULL


      # Update all entries of mcmc for each cluster (plus roots of upstream clusters)
      for (j in 1:n_subtrees) {

        mcmc$h[mappings[[j]]] <- mappings[[j]][all_res[[j]][[i]]$h]

        # Update other components of mcmc (easier!)

        # mcmc$d does not update at frozen nodes
        unfrozen <- setdiff(1:all_res[[j]][[i]]$n, datas[[j]]$frozen)
        mcmc$d[mappings[[j]][unfrozen]] <- all_res[[j]][[i]]$d[unfrozen]

        mcmc$w[mappings[[j]]] <- all_res[[j]][[i]]$w
        mcmc$t[mappings[[j]]] <- all_res[[j]][[i]]$t
        mcmc$m01[mappings[[j]]] <- all_res[[j]][[i]]$m01
        mcmc$m10[mappings[[j]]] <- all_res[[j]][[i]]$m10
        mcmc$m0y[mappings[[j]]] <- all_res[[j]][[i]]$m0y
        mcmc$m1y[mappings[[j]]] <- all_res[[j]][[i]]$m1y
        mcmc$mx0[mappings[[j]]] <- all_res[[j]][[i]]$mx0
        mcmc$mx1[mappings[[j]]] <- all_res[[j]][[i]]$mx1
        mcmc$mxy[mappings[[j]]] <- all_res[[j]][[i]]$mxy
        mcmc$g_lik[mappings[[j]]] <- all_res[[j]][[i]]$g_lik

        # For e_lik, compute as differences
        mcmc$e_lik <- mcmc$e_lik + all_res[[j]][[i]]$e_lik - ifelse(i==1, mcmcs[[j]]$e_lik, all_res[[j]][[i-1]]$e_lik)
      }


      ## Correct node degrees
      # Test this; check not much changes
      #mcmc$d <- sapply(1:mcmc$n, function(x){sum(mcmc$h[2:mcmc$n] == x)}) # Node degrees
      # Can make this smarter...

      res[[i]] <- mcmc
    }

    return(res)

  }
}






