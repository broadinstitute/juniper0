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

#' Internal function to run local moves on a subtree
#'
#' @param mcmc Current state of the Markov chain.
#' @param data Fixed parameter values supplied to the Markov chain.
#' @return List of samples from the local tree.
#' @export
local_mcmc <- function(mcmc, data){


  #set.seed(1)
  res <- list()

  #mcmc <- mcmcs[[j]]
  #data <- datas[[j]]
  #data$n_local <- 100

  for (r in 1:data$n_local) {

    #print(r)

    ## SAFETY STUFF---TURN ON WHEN MAKING CHANGES

    # if(length(unlist(mcmc$tmu)) != length(unlist(mcmc$subs$from))){
    #   print(r)
    #   stop("Updated mutations wrong")
    # }
    #
    # if(!all(mcmc$g_lik == sapply(1:mcmc$n, g_lik, mcmc=mcmc, data=data))){
    #   bad <- (which(mcmc$g_lik != sapply(1:mcmc$n, g_lik, mcmc=mcmc, data=data)))
    #   print(bad)
    #   print(mcmc$external_roots)
    #   print(mcmc$g_lik[bad])
    #   print(sapply(1:mcmc$n, g_lik, mcmc=mcmc, data=data)[bad])
    #   print(r)
    #   print(data$rooted)
    #   stop("Genomic likelihood error")
    # }
    #
    #
    # if(!all(mcmc$e_lik == sapply(1:mcmc$n, e_lik_personal, mcmc=mcmc, data=data))){
    #   bad <- (which(mcmc$e_lik != sapply(1:mcmc$n, e_lik_personal, mcmc=mcmc, data=data)))
    #   print(bad)
    #   print(mcmc$external_roots)
    #   print(mcmc$e_lik[bad])
    #   print(sapply(1:mcmc$n, e_lik_personal, mcmc=mcmc, data=data)[bad])
    #   print(r)
    #   print(data$rooted)
    #   stop("e likelihood error")
    # }
    #
    # if(!all(mcmc$m_lik == sapply(1:mcmc$n, m_lik, mcmc=mcmc, data=data))){
    #   bad <- (which(mcmc$m_lik != sapply(1:mcmc$n, m_lik, mcmc=mcmc, data=data)))
    #   print(bad)
    #   print(mcmc$external_roots)
    #   print(mcmc$m_lik[bad])
    #   print(sapply(1:mcmc$n, m_lik, mcmc=mcmc, data=data)[bad])
    #   print(r)
    #   print(data$rooted)
    #   stop("m likelihood error")
    # }

    # if(abs(sum(mcmc$m_lik + mcmc$e_lik) - e_lik(mcmc, data)) > 1e-6){
    #   print(abs(sum(mcmc$m_lik + mcmc$e_lik) - e_lik(mcmc, data)))
    #   print(sum(mcmc$m_lik + mcmc$e_lik), digits = 20)
    #   print(e_lik(mcmc, data), digits = 20)
    #   print(mcmc$e_lik[mcmc$external_roots])
    #   stop("disagreement error")
    # }

    #
    # ## Check that no SNVs are listed in "dropout"
    # for (i in 1:mcmc$n) {
    #   if(any(
    #     mcmc$subs$pos[[i]] %in% mcmc$dropout[[i]]
    #   )){
    #     stop("No mutations should be listed at positions that drop out")
    #   }
    # }
    #
    # ## Check that "dropout" is always correct
    # for (i in 1:mcmc$n) {
    #   if(
    #     any(mcmc$dropout[[i]] != get_dropout(mcmc, data, i)) | (!all(mcmc$dropout[[mcmc$h[i]]] %in% mcmc$dropout[[i]]))
    #   ){
    #     stop("Dropout updated incorrectly")
    #   }
    # }
    #
    # Check that external roots have degree 0
    if(any(mcmc$external_roots %in% mcmc$h)){
      stop("External roots must have degree 0")
    }
    #
    # if(mcmc$n > data$n_obs){
    #   degs <- sapply((data$n_obs + 1):(mcmc$n), function(n){length(which(mcmc$h == n))})
    #   if(any(degs < 2 & !((data$n_obs + 1):(mcmc$n) %in% mcmc$external_roots))){
    #     print(which(degs < 2) + data$n_obs)
    #     print(mcmc$external_roots)
    #     print(r)
    #     stop("degree error")
    #   }
    # }
    #
    # if(!data$observed_root & length(which(mcmc$h == 1)) < 2){
    #   print(which(mcmc$h == 1))
    #   print(r)
    #   stop("root degree error")
    # }

    # Move 4
    mcmc <- move_seq(mcmc, data, also_resample_tmu = F)

    # Move 5
    mcmc <- move_seq(mcmc, data, also_resample_tmu = T)

    # Move 6
    mcmc <- move_w_t(mcmc, data)

    # Move 7
    mcmc <- move_w_t(mcmc, data, recursive = T)

    # Move 8
    mcmc <- move_genotype(mcmc, data)

    # Move 9
    if(runif(1) < 1/2){
      mcmc <- move_h_step(mcmc, data)
    }else{
      mcmc <- move_h_step(mcmc, data, upstream = F)
    }

    # Move 10
    mcmc <- move_h_global(mcmc, data, biassed = F)

    # Move 11
    mcmc <- move_h_global(mcmc, data)

    # Move 12
    mcmc <- move_swap(mcmc, data)

    # Move 13
    mcmc <- move_swap(mcmc, data, exchange_children = T)

    # Move 14
    if(runif(1) < 1/2){
      mcmc <- move_create(mcmc, data)
    }else{
      mcmc <- move_delete(mcmc, data)
    }

    # Move 15
    if(runif(1) < 1/2){
      mcmc <- move_create(mcmc, data, upstream = F)
    }else{
      mcmc <- move_delete(mcmc, data, upstream = F)
    }

    # Move 16
    if(runif(1) < 1/2){
      mcmc <- move_create(mcmc, data, upstream = T, biassed = T)
    }else{
      mcmc <- move_delete(mcmc, data, upstream = T, biassed = T)
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
    return(
      list(all_res[[1]], 1)
    )
  }else{

    ## Loop through each sample and produce an amalgamated MCMC state:

    # Create a list to store the amalgamated results
    res <- list()
    for (i in 1:n_samples) {

      if(i == n_samples){
        new_roots <- c()
      }

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
      mcmc$bot <- mcmc$bot[1:mcmc$n]
      mcmc$dropout <- mcmc$dropout[1:mcmc$n]

      mcmc$root <- NULL
      mcmc$cluster <- NULL

      mcmc$e_lik <- rep(0, mcmc$n)
      mcmc$g_lik <- rep(0, mcmc$n)
      mcmc$m_lik <- rep(0, mcmc$n)

      # Update all entries of mcmc for each cluster (plus roots of upstream clusters)
      for (j in 1:n_subtrees) {

        #print("hi")

        # UPDATE H
        mcmc$h[mappings[[j]]] <- mappings[[j]][all_res[[j]][[i]]$h]

        # Update other components of mcmc (easier!)
        mcmc$seq[mappings[[j]]] <- all_res[[j]][[i]]$seq

        mcmc$subs$from[mappings[[j]]] <- all_res[[j]][[i]]$subs$from
        mcmc$subs$pos[mappings[[j]]] <- all_res[[j]][[i]]$subs$pos
        mcmc$subs$to[mappings[[j]]] <- all_res[[j]][[i]]$subs$to
        mcmc$tmu[mappings[[j]]] <- all_res[[j]][[i]]$tmu



        mcmc$bot[mappings[[j]]] <- all_res[[j]][[i]]$bot
        mcmc$dropout[mappings[[j]]] <- all_res[[j]][[i]]$dropout

        # For nodes in two clusters, one likelihood is always 0, so this is OK
        mcmc$e_lik[mappings[[j]]] <- mcmc$e_lik[mappings[[j]]] + all_res[[j]][[i]]$e_lik
        mcmc$g_lik[mappings[[j]]] <- mcmc$g_lik[mappings[[j]]] + all_res[[j]][[i]]$g_lik
        mcmc$m_lik[mappings[[j]]] <- mcmc$m_lik[mappings[[j]]] + all_res[[j]][[i]]$m_lik



        if(i == n_samples){
          new_roots <- c(new_roots, mappings[[j]][1])
        }
      }

      res[[i]] <- mcmc


      ## SAFETY MODE
      if(!is.na(data$safety)){

        if(any(abs(mcmc$g_lik - sapply(1:mcmc$n, g_lik, mcmc=mcmc, data=data)) > data$safety)){
          stop("g_lik error in amalgamate")
        }

        if(any(abs(mcmc$e_lik - sapply(1:mcmc$n, e_lik_personal, mcmc=mcmc, data=data)) > data$safety)){
          stop("e_lik error in amalgamate")
        }

        if(any(abs(mcmc$m_lik - sapply(1:mcmc$n, m_lik, mcmc=mcmc, data=data)) > data$safety)){
          stop("m_lik error in amalgamate")
        }
      }
    }

    return(list(res, new_roots))

  }
}
