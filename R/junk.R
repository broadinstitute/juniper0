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


#### OLD VERSION OF G_LIK - IMPORTANT

# Compute genomic log likelihood for each person

# g_lik <- function(mcmc, data, i){
#
#   if(mcmc$v < 0 | mcmc$mu < 0 | mcmc$p < 0 | mcmc$b < 0 | mcmc$lambda < 0 | mcmc$b > 1 | mcmc$w[i] < 0){
#     return(-Inf)
#   }else{
#
#     # Ancestor of host i
#     h <- mcmc$h[i]
#
#     # Time of end of exponential growth phase in h
#     t_g <- mcmc$t[h] + log(1/sqrt(mcmc$p)) / (mcmc$mu / mcmc$p) / log(mcmc$v)
#
#     # Evolutionary time
#     delta_t <- mcmc$t[i] - mcmc$t[h]
#
#     # If i unobserved, or has no iSNV info provided, simply evolve from end of growth phase in h[i] to end of growth phase in i
#     isnv_info <- TRUE
#     if(i > data$n_obs){
#       isnv_info <- FALSE
#     }
#     if(i <= data$n_obs){
#       if(!data$vcf_present[i]){
#         isnv_info <- FALSE
#       }
#     }
#
#     # if(!isnv_info){
#     #   delta_t <- mcmc$t[i] - mcmc$t[h]
#     # }
#
#     # Evolutionary time from first downstream host of h to infection time of i, approx
#     if(is.infinite(delta_t)){
#       delta_t_prime <- Inf
#     }else{
#       if(data$experimental){
#         # Time of first transmission on the chain from h to i
#         t_1st_trans <- mcmc$seq[[i]][length(mcmc$seq[[i]])]
#         delta_t_prime <- mcmc$t[i] - t_1st_trans
#
#         if(t_1st_trans < mcmc$t[h]){
#           return(-Inf)
#         }
#       }else{
#         delta_t_prime <- delta_t * mcmc$w[i] / (mcmc$w[i] + 1)
#       }
#     }
#
#
#
#
#     if(delta_t <= 0){
#       -Inf
#     }else{
#       # log probability of SPECIFIC iSNV in expo growth phase
#       log_p_isnv <- log(
#         (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * (1 - denovo_cdf(data$filters$af, mcmc$p)) / 3
#       )
#
#       # log probability of no iSNV in expo growth phase
#       log_p_no_isnv <- log(
#         (1 - mcmc$p)^(1/sqrt(mcmc$p)) +
#           (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * denovo_cdf(data$filters$af, mcmc$p)
#       )
#
#       # Number of sites without a mutation
#       no_mut <- data$n_bases -
#         length(mcmc$m01[[i]]) -
#         length(mcmc$m10[[i]]) -
#         length(mcmc$m0y[[i]]) -
#         length(mcmc$m1y[[i]]) -
#         length(mcmc$mx0[[i]]) -
#         length(mcmc$mx1[[i]]) -
#         length(mcmc$mxy[[i]]) #-
#       #length(data$filters$common)
#
#       # print(1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t))
#       # print(evolveJC(1, mcmc$mu, delta_t))
#
#       # Likelihood from x = 0, y = 0 or x = 1, y = 1
#       out <- no_mut * (log(evolveJC(1, mcmc$mu, delta_t)) + ifelse(isnv_info, log_p_no_isnv, 0)) +
#
#         # Likelihood from x = 0, y = 1 and x = 1, y = 0
#         (length(mcmc$m01[[i]]) + length(mcmc$m10[[i]])) * (log(evolveJC(0, mcmc$mu, delta_t)) + ifelse(isnv_info, log_p_no_isnv, 0))
#
#       # If there actually are any reported iSNVs...
#       if(i <= data$n_obs){
#         if(length(data$snvs[[i]]$isnv$call) > 0){
#
#           # Frequencies of added iSNVs
#           freq_0y <- data$snvs[[i]]$isnv$af[match(mcmc$m0y[[i]], data$snvs[[i]]$isnv$call)]
#
#           # Frequencies of deleted iSNVs
#           freq_1y <- data$snvs[[i]]$isnv$af[match(mcmc$m1y[[i]], data$snvs[[i]]$isnv$call)]
#
#           out <- out +
#             # Likelihood from x = 0, 0 < y < 1
#             length(mcmc$m0y[[i]]) * log_p_isnv +
#             sum(log(
#               (evolveJC(1, mcmc$mu, delta_t)) * denovo_normed(freq_0y, mcmc$p, data$filters) +
#                 (evolveJC(0, mcmc$mu, delta_t)) * denovo_normed(1 - freq_0y, mcmc$p, data$filters)
#             )) +
#
#             # Likelihood from x = 1, 0 < y < 1
#             length(mcmc$m1y[[i]]) * log_p_isnv +
#             sum(log(
#               (evolveJC(1, mcmc$mu, delta_t))*denovo_normed(1 - freq_1y, mcmc$p, data$filters) +
#                 (evolveJC(0, mcmc$mu, delta_t))*denovo_normed(freq_1y, mcmc$p, data$filters)
#             ))
#
#           if(h <= data$n_obs){
#             if(length(data$snvs[[h]]$isnv$call) > 0){
#
#               # Frequency of shared iSNV in ancestor of case i
#               freq_xy_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[h]]$isnv$call)]
#
#               # If transmission occurred before end of exponential growth phase, need to back-mutate these frequencies.
#               if(t_1st_trans < t_g){
#                 # To do this, we need to figure out what frequency they started at at t[h], on average
#                 start_freqs_xy <- rep(0, length(mcmc$mxy[[i]]))
#
#                 # Which of mcmc$mx0 are 1y in h?
#                 start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$m1y[[h]])] <- 1
#
#                 # If xy in h: start freq is 1/2 (approximation; can try improving this)
#                 start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$mxy[[h]])] <- 1/2
#
#                 # Proportion of exponential growth phase completed
#                 prop_exp_complete <- (t_1st_trans - mcmc$t[h]) / (t_g - mcmc$t[h])
#
#                 # Linearly interpolate frequencies
#                 freq_xy_anc <- freq_xy_anc * prop_exp_complete + start_freqs_xy * (1 - prop_exp_complete)
#               }
#
#               # Frequency of shared iSNV in case i
#               freq_xy <- data$snvs[[i]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[i]]$isnv$call)]
#
#
#               out <- out +
#                 # Likelihood from 0 < x < 1, 0 < y < 1
#                 sum(log(
#                   (evolveJC(1, mcmc$mu, delta_t_prime)*(1 - freq_xy_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*freq_xy_anc) * exp(log_p_isnv) * denovo_normed(freq_xy, mcmc$p, data$filters) * (1 - p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)) +
#                     (evolveJC(1, mcmc$mu, delta_t_prime)*freq_xy_anc + evolveJC(0, mcmc$mu, delta_t_prime)*(1 - freq_xy_anc)) * exp(log_p_isnv) * denovo_normed(1 - freq_xy, mcmc$p, data$filters) * (1 - p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)) +
#                     p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)
#                 ))
#             }
#           }
#         }
#       }
#
#       if(h <= data$n_obs){
#         if(length(data$snvs[[h]]$isnv$call) > 0){
#
#           # Frequency of iSNV in ancestor of case i
#           freq_x0_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mx0[[i]], data$snvs[[h]]$isnv$call)]
#           freq_x1_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mx1[[i]], data$snvs[[h]]$isnv$call)]
#
#           # If transmission occurred before end of exponential growth phase, need to back-mutate these frequencies.
#           if(t_1st_trans < t_g){
#             # To do this, we need to figure out what frequency they started at at t[h], on average
#             start_freqs_x0 <- rep(0, length(mcmc$mx0[[i]]))
#
#             # Which of mcmc$mx0 are 1y in h?
#             start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$m1y[[h]])] <- 1
#
#             # If xy in h: start freq is 1/2 (approximation; can try improving this)
#             start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$mxy[[h]])] <- 1/2
#
#             # Repeat for mcmc$mx1
#             start_freqs_x1 <- rep(0, length(mcmc$mx1[[i]]))
#
#             # Which of mcmc$mx0 are 1y in h?
#             start_freqs_x1[which(mcmc$mx1[[i]] %in% mcmc$m1y[[h]])] <- 1
#
#             # If xy in h: start freq is 1/2 (approximation; can try improving this)
#             start_freqs_x1[which(mcmc$mx1[[i]] %in% mcmc$mxy[[h]])] <- 1/2
#
#             # Proportion of exponential growth phase completed
#             prop_exp_complete <- (t_1st_trans - mcmc$t[h]) / (t_g - mcmc$t[h])
#
#             # Linearly interpolate frequencies
#             freq_x0_anc <- freq_x0_anc * prop_exp_complete + start_freqs_x0 * (1 - prop_exp_complete)
#             freq_x1_anc <- freq_x1_anc * prop_exp_complete + start_freqs_x1 * (1 - prop_exp_complete)
#
#             if(any(freq_x1_anc < 0) | any(freq_x0_anc < 0)){
#               print("warning2")
#               print(mcmc$seq[[i]])
#               print(t_1st_trans)
#               print(mcmc$t[h])
#               print(freq_x0_anc)
#               print(freq_x1_anc)
#             }
#           }
#
#
#
#           out <- out +
#             # Likelihood from 0 < x < 1, y = 0
#             length(mcmc$mx0[[i]]) * ifelse(isnv_info, log_p_no_isnv, 0) +
#             sum(log(
#               evolveJC(1, mcmc$mu, delta_t_prime)*(1 - freq_x0_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*freq_x0_anc
#             )) + sum(log(1 - p_all_split(mcmc$b, mcmc$w[i], freq_x0_anc))) + # probability we don't transmit successive split bottlenecks
#
#             # Likelihood from 0 < x < 1, y = 1
#             length(mcmc$mx1[[i]]) * ifelse(isnv_info, log_p_no_isnv, 0) +
#             sum(log(
#               evolveJC(1, mcmc$mu, delta_t_prime)*(freq_x1_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*(1 - freq_x1_anc)
#             )) + sum(log(1 - p_all_split(mcmc$b, mcmc$w[i], freq_x1_anc))) # probability we don't transmit successive split bottlenecks
#
#         }
#       }
#
#       if(delta_t_prime < 0){
#         print("warning")
#         print(delta_t)
#         print(delta_t_prime)
#         print(mcmc$seq[[i]])
#         print(out)
#       }
#
#
#
#       return(out)
#     }
#   }
# }

