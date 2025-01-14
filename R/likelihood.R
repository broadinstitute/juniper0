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

# Likelihood from mutations and prior on rho
m_lik <- function(mcmc, data, i, js = NULL){

  if(i %in% mcmc$external_roots){
    return(0)
  }

  if(is.null(js)){
    js <- which(mcmc$h == i)
  }

  if(length(js) == 0){
    return(0)
  }

  # Jukes-Cantor evolutionary probability
  prob_JC <- -mcmc$mu * sum(sapply(js, evo_time, mcmc=mcmc, data=data)) + length(unlist(mcmc$subs$from[js])) * log(mcmc$mu / 3)

  return(
    prob_JC
  )


}

# Compute e_lik for each individual
e_lik_personal <- function(mcmc, data, i, js = NULL){

  if(i %in% mcmc$external_roots){
    return(0)
  }

  if(is.null(js)){
    # Children of i
    js <- which(mcmc$h == i)
  }



  if(length(js) > 0){
    # Reverse order of mcmc$seq[js]
    seq_js <- lapply(mcmc$seq[js], rev)

    # Times of infection of js and hosts leading into js
    t_i <- unlist(seq_js)

    # Which indices have ancestor 1?
    ls <- sapply(seq_js, length)
    anc_1 <- c(0, cumsum(ls[1:length(ls) - 1])) + 1

    t_s <- rep(NA, length(t_i))
    t_s[cumsum(ls)] <- data$s[js]

    h <- 1:length(t_i)
    h[anc_1] <- 1

    t_i <- c(mcmc$seq[[i]][1], t_i)
    t_s <- c(data$s[i], t_s)
    h <- c(0, h)
  }else{
    t_i <- mcmc$seq[[i]][1]
    t_s <- data$s[i]
    h <- 0
  }

  ttree <- matrix(c(t_i, t_s, h), ncol = 3)

  if(nrow(ttree) == 0){
    stop("??")
  }

  # 1st param in NBin offspring distribution
  rho <- mcmc$R * mcmc$psi / (1 - mcmc$psi)

  if(!data$ongoing){
    out <- probTTree(
      ttree, rho, 1-mcmc$psi, mcmc$pi, mcmc$a_g, 1/mcmc$lambda_g, mcmc$a_s, 1/mcmc$lambda_s, Inf, wbar0 = 0, delta_t = 0.1 # wbar0 argument doesn't matter here; unused
    )
    return(out)
  }

  # Min time of infection
  tinfmin <- min(ttree[,1])
  # Length of wbar
  wbar_len <- round((-tinfmin)/0.1)

  # If we need to access values of wbar further back than the earliest entry of wbar0, and if wbar hasn't converged, reject the move (likelihood = -Inf)
  if((length(mcmc$wbar) < wbar_len) & (mcmc$wbar[1] != mcmc$wbar[2])){
    return(-Inf)
  }

  # Relevant entries of wbar
  if(wbar_len > length(mcmc$wbar)){
    wbar0 <- c(
      rep(mcmc$wbar[1], wbar_len - length(mcmc$wbar)),
      mcmc$wbar
    )
  }else{
    wbar0 <- mcmc$wbar[(length(mcmc$wbar) - wbar_len + 1):length(mcmc$wbar)]
  }

  out <- probTTree(
    ttree, rho, 1-mcmc$psi, mcmc$pi, mcmc$a_g, 1/mcmc$lambda_g, mcmc$a_s, 1/mcmc$lambda_s, 0, wbar0, delta_t = 0.1
  )

  return(out)

}

# Compute genomic log likelihood for each person
# This is now just the probability of the iSNVs layered on top
g_lik <- function(mcmc, data, i, js = NULL){

  #return(0)

  if(i > data$n_obs){
    return(0)
  }

  if(i %in% mcmc$external_roots){
    return(0)
  }

  if(!data$vcf_present[i]){
    return(0)
  }

  # Log likelihood
  out <- 0

  # iSNVs detected in i
  isnv_a1 <- data$snvs[[i]]$isnv$a1
  isnv_pos <- data$snvs[[i]]$isnv$pos
  isnv_a2 <- data$snvs[[i]]$isnv$a2

  # Proportion of the population that is a de novo iSNV
  freq <- data$snvs[[i]]$isnv$af1
  freq[mcmc$bot[[i]]] <- 1 - freq[mcmc$bot[[i]]]

  # If using split bottleneck feature: which of these bottlenecks are split?
  if(data$split_bottlenecks & i != 1){
    h <- mcmc$h[i]
    if(h <= data$n_obs){
      if(data$vcf_present[h]){

        # Get iSNVs in h
        h_a1 <- data$snvs[[h]]$isnv$a1
        h_pos <- data$snvs[[h]]$isnv$pos
        h_a2 <- data$snvs[[i]]$isnv$a2

        # Filter to iSNVs that are also in i
        shared <- which(h_pos %in% isnv_pos)
        h_a1 <- h_a1[shared]
        h_pos <- h_pos[shared]
        h_a2 <- h_a2[shared]

        # Mapping onto positions of iSNVs in i
        mapping <- match(h_pos, isnv_pos)

        # Indices of isnv_pos that are the same iSNVs in h
        matched <- mapping[which(h_a1 == isnv_a1[mapping] & h_a2 == isnv_a2[mapping])]

        # Get list of shared iSNVs
        shared_a1 <- isnv_a1[matched]
        shared_pos <- isnv_pos[matched]
        shared_a2 <- isnv_a2[matched]

        if(length(shared_pos) > 0){
          # A split bottleneck is only possible when all mutations along the branch from h to i appear as shared iSNVs
          # There can be no fixed mutations along said branch

          ## Is a split bottleneck valid?
          # We will say that there is indeed a split bottleneck (valid = TRUE) whenever imposing one decreases the total number of needed mutations
          # The number of mutations that get eliminated is length(shared_pos)
          # The number of mutations we pay in penalty is the mutations along the branch from h to i that are not involved in the split bottleneck
          # We count these here
          involved <- 0
          if(length(mcmc$subs$pos[[i]]) > 0){
            for (k in 1:length(mcmc$subs$pos[[i]])) {
              if(mcmc$subs$pos[[i]][k] %in% shared_pos){
                ind <- match(mcmc$subs$pos[[i]][k], shared_pos)

                if(all(c(mcmc$subs$from[[i]][k], mcmc$subs$to[[i]][k]) %in% c(shared_a1[ind], shared_a2[ind]))){
                  # If the two nucleotides involved in the substitution on the branch from h to i match the two alleles in the shared iSNV, this substitution is involved in the split bottleneck
                  involved <- involved + 1
                }
              }
            }
          }
          uninvolved <- length(mcmc$subs$pos[[i]]) - involved

          # Declare a boolean named "valid": TRUE if we're taking the bottleneck infecting i to be split
          valid <- length(shared_pos) > uninvolved

          if(valid){
            # Revise list of iSNVs in i to those that are de novo in i
            keep <- setdiff(1:length(isnv_pos), matched)
            isnv_a1 <- isnv_a1[keep]
            isnv_pos <- isnv_pos[keep]
            isnv_a2 <- isnv_a2[keep]
            freq <- freq[keep]

            # As an approximation, penalize by cost of additional branch from h to i with uninvolved mutations (integrating out their time)
            delta_t <- mcmc$seq[[i]][1] - mcmc$seq[[h]][1]
            out <- out +
              (data$n_bases - uninvolved) * dpois(0, mcmc$mu * delta_t, log = T) + # Jukes-Cantor probability of 0 mutations on a branch of length delta_t
              uninvolved * dpois(1, mcmc$mu * delta_t, log = T) # Jukes-Cantor probability of 1 mutation on a branch of length delta_t
          }
        }else{
          # If no shared iSNVs, split bottleneck invalid
          valid <- F
        }
      }
    }
  }




  # Children of i
  if(is.null(js)){
    js <- which(mcmc$h == i)
  }


  # Time of emergence of SNVs in global phylogeny that might get picked up as iSNVs
  trans_isnv_pos <- integer(0)
  trans_isnv_from <- character(0)
  trans_isnv_to <- character(0)
  trans_isnv_time <- numeric(0)

  for (j in js) {

    # Which mutations actually occur in host i?
    keep <- which(mcmc$tmu[[j]] < mcmc$seq[[j]][length(mcmc$seq[[j]])])

    # Which sites have detected iSNVs that get passed on?
    trans_isnv_pos <- c(trans_isnv_pos, mcmc$subs$pos[[j]][keep])
    # What did said sites mutate from?
    trans_isnv_from <- c(trans_isnv_from, mcmc$subs$from[[j]][keep])
    # And what did said sites mutate into?
    trans_isnv_to <- c(trans_isnv_to, mcmc$subs$to[[j]][keep])
    # And when?
    trans_isnv_time <- c(trans_isnv_time, mcmc$tmu[[j]][keep])
  }

  ## If accounting for split bottlenecks, delete mutations that are no longer necessary along outgoing branches
  # Also, subtract out the likelihood of said mutations from m_lik calculation
  if(data$split_bottlenecks & i != 1){
    if(h <= data$n_obs){
      if(data$vcf_present[h]){
        if(valid){
          if(length(trans_isnv_pos) > 0){
            delete <- integer(0)

            for (k in 1:length(trans_isnv_pos)) {
              if(trans_isnv_pos[k] %in% shared_pos){

                # Index in shared isnvs of this mutation on outgoing branch
                ind <- match(trans_isnv_pos[k], shared_pos)

                # Check if to/from are the same
                if(all(c(trans_isnv_from[k], trans_isnv_to[k]) %in% c(shared_a1[ind], shared_a2[ind]))){

                  # If this passes, no longer need outgoing mutation k
                  delete <- c(delete, k)

                  # Adjust log likelihood to remove one (explicit) mutation
                  # Gives it a huge boost
                  out <- out - log(mcmc$mu/3)
                }
              }
            }

            # Each outgoing branch from i has 1/2 probability of choosing correct ancestral virion
            #out <- out + length(js)*log(1/2)

            # Delete iSNVs on the global phylogeny that are accounted for in split bottleneck
            keep <- setdiff(1:length(trans_isnv_pos), delete)
            trans_isnv_pos <- trans_isnv_pos[keep]
            trans_isnv_from <- trans_isnv_from[keep]
            trans_isnv_to <- trans_isnv_to[keep]
            trans_isnv_time <- trans_isnv_time[keep]

          }
        }
      }
    }
  }

  # For each position p, index of minimum time at which there's a mutation at site p
  keep <- integer(0)
  for (p in trans_isnv_pos) {
    keep <- c(
      keep,
      which(trans_isnv_pos == p)[which.min(trans_isnv_time[which(trans_isnv_pos == p)])]
    )
  }

  trans_isnv_pos <- trans_isnv_pos[keep]
  trans_isnv_from <- trans_isnv_from[keep]
  trans_isnv_to <- trans_isnv_to[keep]
  trans_isnv_time <- trans_isnv_time[keep]

  # Size of viral population at trans_isnv_time
  trans_isnv_size <- exp(mcmc$N_eff * (trans_isnv_time - mcmc$seq[[i]][1]))

  # Which iSNVs detected in i are also on the global phylogeny
  local_in_global <- integer(0)
  # Which iSNVs on the global phylogeny are also detected in i?
  global_in_local <- integer(0)

  # Which iSNVs detected in i aren't also on the global phylogeny, but occur at the same site as one on the global phylogeny?
  local_same_site <- integer(0)
  # Which iSNVs on the global phylogeny aren't also detected in i, but occur at the same site as one detected in i?
  global_same_site <- integer(0)

  for (p in isnv_pos) {
    if(p %in% trans_isnv_pos){
      ind_local <- match(p, isnv_pos)
      ind_global <- match(p, trans_isnv_pos)

      if(all(
        c(isnv_a1[ind_local], isnv_a2[ind_local]) %in% c(trans_isnv_from[ind_global], trans_isnv_to[ind_global])
      )){
        local_in_global <- c(local_in_global, ind_local)
        global_in_local <- c(global_in_local, ind_global)
      }else{
        local_same_site <- c(local_same_site, ind_local)
        global_same_site <- c(global_same_site, ind_global)
      }
    }
  }

  if(length(isnv_pos) > 0){
    local_alone <- setdiff(1:length(isnv_pos), c(local_in_global, local_same_site))
  }else{
    local_alone <- integer(0)
  }

  if(length(trans_isnv_pos) > 0){
    global_alone <- setdiff(1:length(trans_isnv_pos), c(global_in_local, global_same_site))
  }else{
    global_alone <- integer(0)
  }

  ## REMEMBER: N_eff = lambda / rho

  ## Log likelihood contribution

  # For iSNVs observed in i that aren't accounted for on the global phylogeny, compute marginal probability of the denovo frequency
  out <- out + sum(dprop(freq[local_alone], mcmc$mu / mcmc$N_eff, log = T) + log(1/3)) # Choice of "to" nucleotide is 1/3


  # For iSNVs observed in i where the global phylogeny has a different mutation at the same site, must have arisen before the one on the phylo tree)
  out <- out + sum(dprop_bounded(freq[local_same_site], trans_isnv_size[global_same_site], mcmc$mu / mcmc$N_eff, log = T) + log(1/3)) # Choice of "to" nucleotide is 1/3


  # For iSNVs observed in i that ARE accounted for on the global phylogeny, condition on when the first denovo SNV occurs (before/after the one on the phylo tree)
  out <- out + sum(log(
    dprop_bounded(freq[local_in_global], trans_isnv_size[global_in_local], mcmc$mu / mcmc$N_eff, log = F)/3 + # When there's an earlier emergence of this iSNV
      dbeta(freq[local_in_global], 1, trans_isnv_size[global_in_local]) * pgeom(trans_isnv_size[global_in_local] - 1, mcmc$mu / mcmc$N_eff, lower.tail = F) # When there's not
  ))


  # For iSNVs UNobserved in i that ARE accounted for on the global phylogeny, again condition on when the first denovo SNV occurs
  #print(global_alone)
  out <- out + sum(log(
    pprop_bounded(data$filters$af, trans_isnv_size[global_alone], mcmc$mu / mcmc$N_eff, log = F) + # When there's an earlier emergence of this iSNV
      pbeta(data$filters$af, 1, trans_isnv_size[global_alone]) * pgeom(trans_isnv_size[global_alone] - 1, mcmc$mu / mcmc$N_eff, lower.tail = F) # When there's not
  ))

  # And finally, all other sites
  out <- out + (data$n_bases - length(unique(c(isnv_pos, trans_isnv_pos))) - length(data$snvs[[i]]$missing)) * pprop(data$filters$af, mcmc$mu / mcmc$N_eff, log = T)

  return(out)
}
