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


  #set.seed(1)
  res <- list()

  #mcmc <- mcmcs[[j]]
  #data <- datas[[j]]
  #data$n_local <- 100

  for (r in 1:data$n_local) {





    #Move 11
    mcmc <- move_seq(mcmc, data, also_resample_tmu = F)

    if(length(unlist(mcmc$tmu)) != length(unlist(mcmc$subs$from))){
      print(r)
      stop("Updated mutations wrong")
    }

    if(!all(mcmc$g_lik[2:mcmc$n] == sapply(2:mcmc$n, g_lik, mcmc=mcmc, data=data))){
      print(which(mcmc$g_lik != sapply(1:mcmc$n, g_lik, mcmc=mcmc, data=data)))
      print(r)
      stop("Genomic likelihood error")
    }

    mcmc <- move_seq(mcmc, data, also_resample_tmu = T)




    # Move 12

    # Move 13
    mcmc <- move_w_t(mcmc, data)

    mcmc <- move_w_t(mcmc, data, recursive = T)



    if(runif(1) < 1/2){
      # Move 14
      mcmc <- move_h_step(mcmc, data)
    }else{
      # Move 15
      mcmc <- move_h_step(mcmc, data, upstream = F)
    }



    # Move 20
    mcmc <- move_h_global(mcmc, data)
    mcmc <- move_h_global(mcmc, data, biassed = F)

    # Move 21
    mcmc <- move_swap(mcmc, data)

    # Move 22
    mcmc <- move_swap(mcmc, data, exchange_children = T)

    # Move 23
    mcmc <- move_genotype(mcmc, data)

    if(runif(1) < 1/2){
      # Move 24
      mcmc <- move_create(mcmc, data)
    }else{
      # Move 25
      mcmc <- move_delete(mcmc, data)
    }



    if(runif(1) < 1/2){
      # Move 26
      mcmc <- move_create(mcmc, data, upstream = F)

    }else{
      # Move 27
      mcmc <- move_delete(mcmc, data, upstream = F)
    }

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
      # anc_clusters <- c()
      # roots <- c()
      #
      # for (j in 1:n_subtrees) {
      #   roots[j] <- all_res[[j]][[i]]$root
      #   anc_clusters[j] <- all_res[[j]][[i]]$anc_cluster
      # }

      # First determine who the unobserved hosts are in each cluster, so that they may be re-indexed
      displacement <- 0
      mappings <- list()
      for (j in 1:n_subtrees) {

        # For now, mappings[[j]] is where hosts 2:n in subtree j get re-indexed to in the output mcmc

        mappings[[j]] <- all_res[[j]][[i]]$cluster

        # If root is observed, delete it; it's re-appended later
        if(datas[[j]]$observed_root){
          mappings[[j]] <- mappings[[j]][-1]
        }

        # How many unobserved hosts are there at this iteration? (Not including root, since that's updated later)
        n_unobs <- all_res[[j]][[i]]$n - datas[[j]]$n_obs

        if(n_unobs > 0){
          mappings[[j]] <- c(mappings[[j]], (data$n_obs + displacement + 1):(data$n_obs + displacement + n_unobs))
        }
        displacement <- displacement + n_unobs
      }

      # Next, append where host 1 in subtree j gets re-indexed to in output mcmc
      for (j in 1:n_subtrees) {
        k <- all_res[[j]][[i]]$anc_cluster # Which cluster is ancestral to cluster j?
        if(is.na(k)){
          # If it's the root cluster, 1 maps to 1
          mappings[[j]] <- c(1, mappings[[j]])
        }else{
          # Who is the root of cluster j in cluster k?
          who_root <- all_res[[k]][[i]]$external_roots[which(all_res[[k]][[i]]$external_subtrees == j)]
          # Append where it maps to under mappings[[k]], remembering to subtract 1 because mappings[[k]] tells us the image of 2:n, not 1:n
          mappings[[j]] <- c(
            mappings[[k]][who_root - 1],
            mappings[[j]]
          )
        }
      }



      ## Now, fill an amalgamated mcmc with info from the correct subtrees

      # Initialization doesn't really matter; we will initialize to the previous amalgamated mcmc
      # Correct lengths of entries, to avoid carrying over extraneous information
      mcmc$n <- data$n_obs + displacement
      mcmc$h <- mcmc$h[1:mcmc$n]
      mcmc$seq <- mcmc$seq[1:mcmc$n]
      mcmc$subs$from <- mcmc$subs$from[1:mcmc$n]
      mcmc$subs$pos <- mcmc$subs$pos[1:mcmc$n]
      mcmc$subs$to <- mcmc$subs$to[1:mcmc$n]
      mcmc$tmu <- mcmc$tmu[1:mcmc$n]
      mcmc$g_lik <- mcmc$g_lik[1:mcmc$n]
      mcmc$root <- NULL
      mcmc$cluster <- NULL

      # Update all entries of mcmc for each cluster (plus roots of upstream clusters)
      for (j in 1:n_subtrees) {

        # UPDATE H
        mcmc$h[mappings[[j]]] <- mappings[[j]][all_res[[j]][[i]]$h]

        # Update other components of mcmc (easier!)

        # mcmc$d does not update at frozen nodes
        unfrozen <- setdiff(1:all_res[[j]][[i]]$n, all_res[[j]][[i]]$external_roots)

        mcmc$seq[mappings[[j]]] <- all_res[[j]][[i]]$seq
        mcmc$subs[mappings[[j]]] <- all_res[[j]][[i]]$subs

        mcmc$tmu[mappings[[j]]] <- all_res[[j]][[i]]$tmu
        mcmc$bot[mappings[[j]]] <- all_res[[j]][[i]]$bot

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






