### Experimentation with parsimony initialize step

## Idea: under infinitely-many-sites model, easy to initialize to parsimony, since the set of cases containing a given mutation is always a superset/subset of the set of cases containing another mutation, or they don't intersect at all
## The idea here is we will build a graph of mutations that respect this superset/subset relation, to the extent possible.
## We will then find some method to accommodate homoplasic mutations

## In fact, let's try as the first step to identify instances of homoplasy. Then build a tree without homoplasy. Then layer back on these mutations.
## This is easy enough to do in O(n^2), n being the total number of SNVs across the whole dataset

# Helper: which sequences contain a given snv?
which_contain_snv <- function(snv, genotypes){
  which(sapply(genotypes, function(v){snv %in% v}))
}

# Build tree of SNVs
snv_tree <- function(mcmc){
  # Get list of unique sequences observed in the dataset
  unique_genotypes <- unique(mcmc$m01)

  # Get list of all SNVs
  snvs <- unique(unlist(unique_genotypes))

  homoplasic <- character(0)

  # Loop over all (unordered) pairs of snvs to check homoplasy
  for (i in 1:(length(snvs) - 1)) {
    for(j in (i+1):length(snvs)){
      which_i <- which_contain_snv(snvs[i], unique_genotypes)
      which_j <- which_contain_snv(snvs[j], unique_genotypes)
      # Check whether their relationship is homoplasic
      n_shared <- length(intersect(which_i, which_j))

      if(n_shared > 0 & n_shared < length(which_i) & n_shared < length(which_j)){
        homoplasic <- c(homoplasic, snvs[i], snvs[j])
      }
    }
  }

  # Get all non-homoplasic SNVs
  non_homoplasic_snvs <- setdiff(snvs, homoplasic)

  # How often does each of these appear in unique_genotypes?
  freqs <- sapply(non_homoplasic_snvs, function(snv){sum(unlist(unique_genotypes) == snv)})

  # Sort by frequency, decreasing
  freqs <- sort(freqs, decreasing = T)

  # Re-decleare non_homoplasic_snvs based on the new order
  non_homoplasic_snvs <- names(freqs)

  # Output is a three-column table: (1) SNV name, (2) parent SNV, (3) whether the sets of sequences containing SNV and its parent are equal (hence their order is exchangeable)
  tree <- data.frame(snv = non_homoplasic_snvs, parent = NA, equal = NA)
  tree$parent[1] <- "root"
  tree$equal[1] <- F

  for (i in 2:nrow(tree)) {
    # Which sequences contain snv i?
    which_i <- which_contain_snv(non_homoplasic_snvs[i], unique_genotypes)
    # Have we found the parent of i yet?
    found <- FALSE
    j <- i-1
    while(!found & j > 0){
      # Which sequences contain snv j?
      which_j <- which_contain_snv(non_homoplasic_snvs[j], unique_genotypes)

      # If all of the sequences containing snv i also contain snv j, i is descended directly from j (due to sort step)
      if(all(which_i %in% which_j)){
        tree$parent[i] <- non_homoplasic_snvs[j]

        # If which_i and which_j identical, order doesn't matter
        if(all(which_j %in% which_i)){
          tree$equal[i] <- T
        }else{
          tree$equal[i] <- F
        }

        # We have found the parent of i
        found <- T
      }else{
        # We have not found the parent of i
        j <- j - 1
      }
    }

    # If we never found the parent of i, it's descended from the root
    if(!found){
      tree$parent[i] <- "root"
      tree$equal[i] <- F
    }
  }



}
