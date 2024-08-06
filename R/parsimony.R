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
