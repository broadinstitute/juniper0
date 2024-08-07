## Approximate initialization of the Markov Chain to parsimony

# Idea: homoplasies are rare, and in the absence of homoplasies, we get parsimony by looking at the inclusion tree of mutations
# So let's try building such an inclusion tree, skipping homoplasic mutations when they occur
# Finally, we will use some heuristic for re-incorporating homoplasies

parsimony <- function(data){

  # SNVs rel. root
  m01 <- list()
  for (i in 1:length(data$snvs)) {
    m01[[i]] <- data$snvs[[i]]$snv$call
  }

  # Assume MCMC initialized to star-shaped tree, i.e. everyone's ancestor is the root



  # Table of all SNVs, decreasing in prevalence
  snv_tab <- sort(table(unlist(m01)), decreasing = T)

  # Names of mutations
  muts <- names(snv_tab)

  # Who has each mutation?
  who <- lapply(muts, function(s){list(which(sapply(m01, function(v){s %in% v})))})

  # Construct tree of mutations
  h_mut <- list(0) # Ancestry vector of mutations. 0 means root
  skipped <- F
  for (i in 2:length(muts)) {
    j <- i-1
    done <- F

    while (j > 0 & done == F) {
      if(all(who[[i]][[1]] %in% who[[j]][[1]])){
        h_mut[[i]] <- j
        skipped[i] <- F
        done <- T
      }else if(any(who[[i]][[1]] %in% who[[j]][[1]])){
        # If some, but not all, of the people who have mutation i have mutation j, then since more or equal people have mutation j, their intersection is a proper subset of i and j
        # So, there must be a homoplasy
        skipped[i] <- T
        h_mut[[i]] <- integer(0)
        print("homoplasy!")
        done <- T
      }
      j <- j-1
    }
    if(!done){
      h_mut[[i]] <- 0
      skipped[i] <- F
    }
  }

  ## For now: just attach each input sequence to the one on our new tree that's most similar
  # Later we will properly deal with homoplasic mutations

  # Create sequences corresponding to our tree
  seqs <- list()
  for (i in 1:length(muts)) {
    if(!skipped[i]){
      seqs[[i]] <- muts[i]
      h <- h_mut[[i]]
      while (h != 0) {
        seqs[[i]] <- c(seqs[[i]], muts[h])
        h <- h_mut[[h]]
      }
    }
  }

  # Assign to each user-supplied sequence an index in "seqs", its nearest neighbor (usually perfect match)
  # Note that the mutations in m01[[i]] should always be a superset of the best fit in seqs (discrepancy due to homoplasies)
  assignments <- c()
  for (i in 1:length(data$snvs)) {
    if(length(m01[[i]]) == 0){
      assignments[i] <- 0
    }else{
      n_match <- sapply(seqs, function(s){length(setdiff(s, m01[[i]])) + length(setdiff(m01[[i]], s))})
      n_match[skipped] <- Inf
      assignments[i] <- which.min(n_match)
    }

  }

  # Next, identify the "root" of each sequence in "seqs", i.e. the earliest case to exhibit it
  # If none exists, create it
  roots <- c()
  counter <- length(data$snvs)
  h_init <- rep(NA, length(data$snvs))
  for (i in 1:length(seqs)) {
    if(!skipped[i]){
      # Who got assigned to seqs[[i]]?
      who <- which(assignments == i)

      if(length(who) == 0){
        create <- T
      }else{
        earliest <- who[which.min(data$s[who])]
        if(identical(sort(m01[[earliest]]), sort(seqs[[i]]))){
          # If some people indeed got assigned seqs[[i]], and their sequence is identical, use that as the root
          roots[i] <- earliest
          h_init[setdiff(who, roots[i])] <- roots[i]
          create <- F
        }else{
          create <- T
        }
      }

      # Otherwise create a new node with the desired sequence
      if(create){
        counter <- counter + 1
        m01[[counter]] <- seqs[[i]]
        roots[i] <- counter
      }
    }
  }

  # Update h_init for people assigned to root
  h_init[assignments == 0] <- 1


  # Update h_init for each root
  for (i in 1:length(roots)) {
    if(!skipped[i]){
      if(h_mut[[i]] == 0){
        h_init[roots[i]] <- 1
      }else{
        h_init[roots[i]] <- roots[h_mut[[i]]]
      }
    }
  }

  ## REALLY needs to be checked if this is right...
  return(h_init)






  ## Now, revisit homoplasic mutations.
  # To deal with these, we loop through all mutations again, and break off where we find the maximum number of cases we can put on one branch
  # This translates into the length of the overlap in number of mutations for i homoplasic and j
  # This is our heuristic step
  # May need revision

  # for (i in which(skipped)) {
  #
  #   people <- who[[i]][[1]] # People with the mutation to map onto tree
  #   who[[i]] <- list() # Clear the list so we can separate into groups based on evolutionary history
  #
  #   while (length(people) > 0) {
  #     n_shared <- list() # Number of others who share the mutation
  #     for (j in 1:length(muts)) {
  #       n_shared[[j]] <- list()
  #       for (k in 1:length) {
  #
  #       }
  #       n_shared[j] <- sum(people %in% who[[j]])
  #     }
  #     best <- which.max(n_shared)
  #     found <- which()
  #   }
  #
  # }


}


# For each SNV in i, how many neighbors go present (i) to absent (neighbor)?
# present_to_absent <- table(factor(c(
#   mcmc$m01[[i]],
#   unlist(mcmc$m10[js])
# ), levels = snvs))
#
# # For each SNV in i, how many neighbors go present (i) to isnv (neighbor)?
# present_to_isnv <- table(factor(c(
#   mcmc$mx1[[i]],
#   unlist(mcmc$m1y[js])
# ), levels = snvs))
#
# # For each SNV in i, how many neighbors go isnv (i) to present (neighbor)?
# isnv_to_present <- table(factor(c(
#   mcmc$m1y[[i]],
#   unlist(mcmc$mx1[js])
# ), levels = snvs))
#
# # For each SNV in i, how many neighbors go isnv (i) to absent (neighbor)?
# isnv_to_absent <- table(factor(c(
#   mcmc$m0y[[i]],
#   unlist(mcmc$mx0[js])
# ), levels = snvs))
#
# # For each SNV in i, how many neighbors go absent (i) to present (neighbor)?
# absent_to_present <- table(factor(c(
#   mcmc$m10[[i]],
#   unlist(mcmc$m01[js])
# ), levels = snvs))
#
# # For each SNV in i, how many neighbors go absent (i) to isnv (neighbor)?
# absent_to_isnv <- table(factor(c(
#   mcmc$mx0[[i]],
#   unlist(mcmc$m0y[js])
# ), levels = snvs))
#
# ## Based on the information, how do we ideally want to change the genotype at i?
# change_present_absent <- names(which(present_to_absent >= n_neighbors - 1))
#
# change_present_isnv <- setdiff(
#   names(which(present_to_absent + present_to_isnv >= 2)), # At least two neighbors have the deletion, at consensus level or in iSNV form
#   change_present_absent
# )
#
# change_isnv_present <- names(which(isnv_to_present >= n_neighbors - 1))
#
# change_isnv_absent <- names(which(isnv_to_absent >= n_neighbors - 1))
#
# change_absent_present <- names(which(absent_to_present >= n_neighbors - 1))
#
# change_absent_isnv <- setdiff(
#   names(which(absent_to_present + absent_to_isnv >= 2)), # At least two neighbors have the addition, at consensus level or in iSNV form
#   change_absent_present
# )
#
# ## Finally, we need to make some of these moves illegal based on observed data
# # Whether vcf_present or not, for an observed site, cannot ever change present to absent or vice versa
# # Ideal state is to switch to iSNV
# change_present_isnv <- union(
#   change_present_isnv,
#   intersect(change_present_absent, observed)
# )
# change_present_absent <- setdiff(change_present_absent, observed)
#
#
# change_absent_isnv <- union(
#   change_absent_isnv,
#   intersect(change_absent_present, observed)
# )
# change_absent_present <- setdiff(change_absent_present, observed)
#
# # isnv to present only valid if isnv frequency over LOD
# # isnv to absent only valid if isnv frequency under LOD
#
# observed_present <- intersect(mcmc$isnv$call[mcmc$isnv$af > 0.5], observed)
# observed_absent <- intersect(mcmc$isnv$call[mcmc$isnv$af < 0.5], observed)
#
# change_isnv_absent <- setdiff(
#   change_isnv_absent,
#   observed_present # If an iSNV has frequency above 0.5 and the site is observed, cannot change to absent
# )
# # Ideal state is to keep at iSNV; no action required
#
# change_isnv_present <- setdiff(
#   change_isnv_present,
#   observed_absent # If an iSNV has frequency below 0.5 and the site is observed, cannot change to present
# )
# # Ideal state is to keep at iSNV; no action required

# # Which ones go 0y or 1y in j, or x0 or x1 in i? (with multiplicity)
# isnvs <- c(
#   mcmc$mx0[[i]],
#   mcmc$mx1[[i]],
#   unlist(mcmc$m0y[js]),
#   unlist(mcmc$m1y[js])
# )
# # Ensure these are valid mutations to flip
# if(i <= data$n_obs){
#   isnvs <- isnvs[isnvs %in% snvs]
# }
#
# tab_isnv <- table(isnvs)
# ind_isnv <- match(names(tab_isnv), snvs) # Indices of these iSNVs in "snvs"
#
# # And the others?
# non_isnvs <- c(
#   mcmc$m10[[i]],
#   mcmc$m01[[i]],
#   unlist(mcmc$m01[js]),
#   unlist(mcmc$m10[js])
# )
# # Ensure these are valid mutations to flip
# if(i <= data$n_obs){
#   non_isnvs <- non_isnvs[non_isnvs %in% snvs]
# }
# tab_non_isnv <- table(non_isnvs)
# ind_non_isnv <- match(names(tab_non_isnv), snvs) # Indices of these SNVs in "snvs"
#
# # Probability of swapping
# probs <- rep(0, length(snvs))
# probs[ind_non_isnv] <- probs[ind_non_isnv] + unname(tab_non_isnv)
# probs[ind_isnv] <- probs[ind_isnv] + unname(tab_isnv)/2
# probs <- probs / (length(js) + 1)
#
# # By parsimony, round probability
# probs[probs < 0.5] <- 0
# probs[probs > 0.5] <- 1
#
# # Add random noise
# probs <- (data$eps/2) + (1 - data$eps)*probs
#
# # If our goal is to compute the probability that a new host i has this genotype...
# if(comparison){
#   # The probability of creating a new genotype at i equal to that in MCMC is the probability we don't swap anything in i
#   return(sum(log(1-probs)))
# }else{
#   # Which ones get swapped?
#   which_swap <- runif(length(snvs)) < probs
#   swaps <- snvs[which_swap]
#
#   # Log probability of picking this genotype
#   log_p <- sum(log(probs[which_swap])) + sum(log(1 - probs[!which_swap]))
#
#   for (snv in swaps) {
#     mcmc <- flip_genotype(mcmc, i, js, snv)
#   }
#
#   return(list(mcmc, log_p))
# }

# for (snv in snvs) {
#   is_observed <- snv %in% observed
#   if(runif(1) > data$eps){
#     # Propose ideal state
#     if(snv %in% change_present_absent){
#       new <- change_genotype(mcmc, data, snv, from = "present", to = "absent", i, js, is_observed)
#     }
#     if(snv %in% change_present_isnv){
#       new <- change_genotype(mcmc, data, snv, from = "present", to = "isnv", i, js, is_observed)
#     }
#     if(snv %in% change_isnv_present){
#       new <- change_genotype(mcmc, data, snv, from = "isnv", to = "present", i, js, is_observed)
#     }
#     if(snv %in% change_isnv_absent){
#       new <- change_genotype(mcmc, data, snv, from = "isnv", to = "absent", i, js, is_observed)
#     }
#     if(snv %in% change_absent_present){
#       new <- change_genotype(mcmc, data, snv, from = "absent", to = "present", i, js, is_observed)
#     }
#     if(snv %in% change_absent_isnv){
#       new <- change_genotype(mcmc, data, snv, from = "absent", to = "isnv", i, js, is_observed)
#     }
#     # Otherwise do nothing; we're already in the ideal state
#
#     mcmc <- new[[1]]
#     log_p <- log_p + new[[2]] + log(1 - data$eps)
#
#   }else{
#     # Make a random change
#     if(snv %in% change_present_absent){
#       # This can only happen when no data on the SNV; hence, we sample the options isnv and present with equal probabilities
#       if(runif(1) < 1/2){
#         new <- change_genotype(mcmc, data, snv, from = "present", to = "isnv", i, js, is_observed)
#       }
#       log_p <- log_p + log(1/2)
#     }
#     if(snv %in% change_present_isnv){
#       # If snv %in% observed_present, the only other option is to keep in present, so do nothing
#       # Otherwise ...
#       if(!(snv %in% observed_present)){
#         if(runif(1) < 1/2){
#           new <- change_genotype(mcmc, data, snv, from = "present", to = "absent", i, js, is_observed)
#         }
#         log_p <- log_p + log(1/2)
#       }
#     }
#     if(snv %in% change_isnv_present){
#       # If snv %in% observed_present, the only other option is to keep in isnv, so do nothing
#       # Otherwise ...
#       if(!(snv %in% observed_present)){
#         if(runif(1) < 1/2){
#           new <- change_genotype(mcmc, data, snv, from = "isnv", to = "absent", i, js, is_observed)
#         }
#         log_p <- log_p + log(1/2)
#       }
#     }
#     if(snv %in% change_isnv_absent){
#       # If snv %in% observed_absent, the only other option is to keep in isnv, so do nothing
#       # Otherwise ...
#       if(!(snv %in% observed_absent)){
#         if(runif(1) < 1/2){
#           new <- change_genotype(mcmc, data, snv, from = "isnv", to = "present", i, js, is_observed)
#         }
#         log_p <- log_p + log(1/2)
#       }
#     }
#     if(snv %in% change_absent_present){
#       # This can only happen when no data on the SNV; hence, we sample the options isnv and present with equal probabilities
#       if(runif(1) < 1/2){
#         new <- change_genotype(mcmc, data, snv, from = "absent", to = "isnv", i, js, is_observed)
#       }
#       log_p <- log_p + log(1/2)
#     }
#     if(snv %in% change_absent_isnv){
#       # If snv %in% observed_absent, the only other option is to keep in absent, so do nothing
#       # Otherwise ...
#       if(!(snv %in% observed_absent)){
#         if(runif(1) < 1/2){
#           new <- change_genotype(mcmc, data, snv, from = "absent", to = "present", i, js, is_observed)
#         }
#         log_p <- log_p + log(1/2)
#       }
#     }
#
#     mcmc <- new[[1]]
#     log_p <- log_p + new[[2]] + log(data$eps)
#   }
# }
