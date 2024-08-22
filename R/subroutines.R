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

### Helper functions

# Functions to convert between raw format and nucleotide letters:
base_to_raw <- function(b){
  out <- rep(as.raw(04), length(b))
  out[b == "A"] <- as.raw(136)
  out[b == "C"] <- as.raw(40)
  out[b == "G"] <- as.raw(72)
  out[b == "T"] <- as.raw(24)
  return(out)
}

raw_to_base <- function(r){
  out <- rep("N", length(r))
  out[r == as.raw(136)] <- "A"
  out[r == as.raw(40)] <- "C"
  out[r == as.raw(72)] <- "G"
  out[r == as.raw(24)] <- "T"
  return(out)
}

genetic_info <- function(seq1, seq2, filters, vcf = NULL){
  # List output: list of SNVs, iSNVs, and positions with no information
  out <- list()

  ## Get iSNVs from VCF, if provided
  if(!is.null(vcf)){

    out$isnv <- list()

    # Position
    pos <- vcf$V2
    # Allele on reference genome
    ref <- vcf$V4
    # Alternate allele
    alt <- vcf$V5

    # Which allele is present in the root?
    root_allele <- raw_to_base(seq1[pos])

    # For which positions does the root allele match the ALT allele?
    # We will need to swap the ref and alt at such positions
    root_alt <- which(root_allele == alt)

    # What's the ref allele at these positions
    ref_swap <- ref[root_alt]

    # What's the alt allele at these positions
    alt_swap <- alt[root_alt]

    # Swap
    ref[root_alt] <- alt_swap
    alt[root_alt] <- ref_swap

    # Final scenario: the root allele is neither ref nor alt
    neither <- which(root_allele != alt & root_allele != ref)

    # Here, simply set the ref category to the root
    ref[neither] <- root_allele[neither]

    # Info column
    info <- vcf$V8

    # Read depth
    dp <- gsub(".*DP=", "", info)
    dp <- sub(";.*", "", dp)
    dp <- as.numeric(dp)

    # Allele frequency
    af <- gsub(".*AF=", "", info)
    af <- sub(";.*", "", af)
    af <- as.numeric(af)

    # Strand bias
    sb <- gsub(".*;SB=", "", info)
    sb <- sub(";.*", "", sb)
    sb <- as.numeric(sb)

    # Which sites pass the filters?
    keep <- which(dp >= filters$dp & sb < filters$sb & af >= filters$af & af <= 1 - filters$af & !(pos %in% filters$common))
    pos <- pos[keep]
    ref <- ref[keep]
    alt <- alt[keep]
    af <- af[keep]

    out$isnv$call <- paste0(ref, pos, alt)
    out$isnv$pos <- pos
    out$isnv$af <- af

  }

  out$snv <- list()

  ## Get SNVs from FASTA
  snv_pos <- which(
    seq1 != seq2 &
      seq1 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24)) &
      seq2 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24))
  )

  ## Get missing sites from FASTA in seq2
  missing_pos <- which(
    !(seq2 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24)))
  )

  # Remove positions already accounted for in VCF, if provided
  if(!is.null(vcf)){
    snv_pos <- setdiff(snv_pos, pos)
    missing_pos <- setdiff(missing_pos, pos)
  }

  old <- raw_to_base(seq1[snv_pos])
  new <- raw_to_base(seq2[snv_pos])

  out$snv$call <- paste0(old, snv_pos, new)
  out$snv$pos <- snv_pos
  out$missing <- list()
  out$missing$pos <- missing_pos

  return(out)

}

# Number of total cases created over an interval of length delta_t
tot_cases <- function(mcmc, delta_t){
  if(mcmc$R == 1){
    (delta_t / (mcmc$a_g / mcmc$lambda_g)) + 1
  }else{
    (mcmc$R^((delta_t / (mcmc$a_g / mcmc$lambda_g)) + 1) - 1) / (mcmc$R - 1)
  }
}

# Probability not sampled in 0, 1, ..., g_max generations
alpha_gs <- function(mcmc, g_max){
  out <- 1 - mcmc$alpha
  for (i in 1:g_max) {
    if(is.infinite(mcmc$rho)){
      out <- c(
        out,
        (1 - mcmc$alpha) * exp(R * (out[i] - 1))
      )

    }else{
      out <- c(
        out,
        (1 - mcmc$alpha) * (1 + out[i] * (mcmc$psi - 1))^(-mcmc$rho) * mcmc$psi^mcmc$rho
      )
    }
  }
  out
}

# log probability none of the cases created over an interval of length delta_t are sampled
log_p_unsampled <- function(mcmc, delta_t){
  tot_cases(mcmc, delta_t) * log(1 - mcmc$alpha)
}

# Convert adjacency matrix to ancestry vector
adj_to_anc <- function(adj, i, h = NULL){
  if(is.null(h)){
    h <- rep(0, ncol(adj))
    h[i] <- NA
  }
  children <- which(adj[i,] == 1 & h == 0)
  h[children] <- i
  for (j in children) {
    h <- adj_to_anc(adj, j, h)
  }
  return(h)
}


# Approximate parsimony tree
split_cluster <- function(cluster, past_muts, id){

  # Mutations in cluster
  muts <- table(unlist(mcmc$m01[cluster]))

  # Get rid of those that are already accounted for
  muts <- muts[!(names(muts) %in% past_muts)]

  # If no more mutations, nothing to do!
  if(length(muts) == 0){
    return(
      list(
        list(
          cluster,
          past_muts,
          id
        )
      )
    )
  }else{
    # Best mutation at which to split
    mut <- names(which.max(muts))

    # Who in the cluster has the mutation?
    who <- cluster[which(sapply(cluster, function(i){mut %in% mcmc$m01[[i]]}))]

    # If it's everyone, run split_cluster again, now updating past_muts
    if(length(who) == length(cluster)){
      return(
        split_cluster(who, c(past_muts, mut), id)
      )
    }else{
      # Return a list of the two resulting clusters from the split
      return(
        c(
          list(
            list(
              cluster,
              past_muts,
              id
            )
          ),
          split_cluster(setdiff(cluster, who), past_muts, c(id, "A")),
          split_cluster(who, c(past_muts, mut), c(id, "B"))
        )
      )
    }
  }
}

### Initialize to (approximate) parsimony tree
# if(FALSE){
#   pars <- split_cluster(2:n, character(0), character(0))
#
#   # For each entry of "pars", who is its ancestor?
#   pars_anc <- NA
#   # Mutations added since ancestor
#   added <- list(character(0))
#   for (i in 2:length(pars)) {
#     id <- pars[[i]][[3]]
#     if(length(id) == 1){
#       anc_id <- character(0)
#     }else{
#       anc_id <- id[1:(length(id) - 1)]
#     }
#
#     for (j in i:1) { # Probably faster to reverse order of inner for loop
#       if(identical(anc_id, pars[[j]][[3]])){
#         pars_anc[i] <- j
#         break
#       }
#     }
#
#     added[[i]] <- setdiff(pars[[i]][[2]], pars[[pars_anc[i]]][[2]])
#   }
#
#   # Which nodes are terminal?
#   terminal <- which(!(1:length(pars) %in% pars_anc))
#
#   # "pars" is a list of length length(pars).
#   # For certain elements of this list, we need to create a new unobserved host, which will be assigned some id.
#   # We track that with this vector:
#   new_hosts <- rep(NA, length(pars))
#
#   # Loop through each entry in "pars" and update transmission tree
#   for (k in 1:length(pars)) {
#
#     if(!(k %in% terminal)){
#
#       # Increase n
#       mcmc$n <- mcmc$n + 1
#       i <- mcmc$n
#
#       new_hosts[k] <- i
#
#       if(k == 1){
#         mcmc$h[i] <- 1
#
#         mcmc$m01[[i]] <- setdiff(
#           added[[k]],
#           snvs[[1]]$isnv$call
#         )
#
#         mcmc$m10[[i]] <- character(0)
#
#         mcmc$m0y <- character(0)
#         mcmc$m1y[[i]] <- character(0)
#         mcmc$mx0[[i]] <- setdiff(
#           snvs[[1]]$isnv$call,
#           union(
#             snvs[[i]]$snv$call,
#             snvs[[i]]$isnv$call
#           )
#         )
#         mcmc$mx1[[i]] <- intersect(
#           snvs[[1]]$isnv$call,
#           snvs[[i]]$snv$call
#         )
#         if(i==1){
#           mcmc$mxy[[i]] <- character(0)
#         }else{
#           mcmc$mxy[[i]] <- intersect(
#             snvs[[1]]$isnv$call,
#             snvs[[i]]$isnv$call
#           )
#         }
#
#       }else{
#         mcmc$h[i] <- new_hosts[pars_anc[k]]
#
#       }
#     }
#   }
# }


# Get the ancestry of a single node, down to the root
ancestry <- function(h, i){
  if(is.na(h[i])){
    return(i)
  }else{
    return(c(ancestry(h, h[i]), i))
  }
}

# Get the generations of an ancestor vector
generations <- function(h, i){
  if(length(which(h %in% i)) == 0){
    return(list(i))
  }else{
    return(c(list(i), generations(h, which(h %in% i))))
  }
}

# Get approximate time of infection for tracked and untracked hosts leading into a node
get_ts <- function(mcmc, i){
  if(i == 1){
    return(mcmc$t[i])
  }else{
    delta_t <- mcmc$t[i] - mcmc$t[mcmc$h[i]]
    # Increment back in time
    inc <- delta_t / (mcmc$w[i] + 1)
    seq(mcmc$t[i], by = -inc, length.out = mcmc$w[i] + 1)
  }

}

# Get whether someone is observed, in unlisted form
get_obs <- function(mcmc, data, i){
  if(i <= data$n_obs){
    c(T, rep(F, mcmc$w[i]))
  }else{
    rep(F, mcmc$w[i] + 1)
  }
}

# JC evolution
evolveJC <- function(init, mu, delta_t){
  1/4 + (init - 1/4)*exp(-(4*mu/3) * delta_t)
}

# Get initial bottleneck, given composition in the form (fraction of absent particles, fraction of present particles)
init_bot <- function(init_absent, init_present, mu, b, delta_t){
  init <- c(init_absent, init_present)
  # Evolve composition forward
  if(delta_t > 0){
    init <- evolveJC(init, mu, delta_t)
  }

  return(c(
    (1-b) * init[1] + b * init[1]^2,
    b * 2 * init[1] * init[2],
    (1-b) * init[2] + b * init[2]^2
  ))
}

# Evolve P(absent, isnv, present) from bottleneck to bottleneck
evolve_bot <- function(bot, mu, b, delta_t){
  # If the SNV starts as absent, what fraction of particles are still absent, and what fraction are present?
  probs_start_absent <- evolveJC(c(1,0), mu, delta_t)

  # If the SNV starts as present, what fraction of particles have the absent type, and what fraction are present?
  probs_start_present <- evolveJC(c(0,1), mu, delta_t)

  # Limitation: not accounting for types other than present/absent going thru bottleneck; hence, probabilities won't sum to 1

  # Probability absent at next bottleneck
  p_absent <-
    # Absent to absent
    bot[1] * ((1-b) * probs_start_absent[1] + b * probs_start_absent[1]^2) +
    # iSNV to absent
    bot[2] * (((1-b) / 2) + (b / 3)) +
    # Present to absent
    bot[3] * ((1-b) * probs_start_present[1] + b * probs_start_present[1]^2)

  # Probability iSNV at next bottleneck
  p_isnv <-
    # Absent to iSNV
    bot[1] * b * 2 * probs_start_absent[1] * probs_start_absent[2] +
    # iSNV to iSNV
    bot[2] * b / 3 +
    # Present to iSNV
    bot[3] * b * 2 * probs_start_present[1] * probs_start_present[2]

  # Probability present at next bottleneck
  p_present <-
    # Absent to present
    bot[1] * ((1-b) * probs_start_absent[2] + b * probs_start_absent[2]^2) +
    # iSNV to present
    bot[2] * (((1-b) / 2) + (b / 3)) +
    # Present to present
    bot[3] * ((1-b) * probs_start_present[2] + b * probs_start_present[2]^2)

  return(c(p_absent, p_isnv, p_present))
}

# Evolve P(absent, isnv, present) from bottleneck to bottleneck, repeatedly for a transmission chain
evolve_bot_repeatedly <- function(bot, mu, b, deltas){
  if(length(deltas) == 0){
    return(bot)
  }

  for (delta_t in deltas) {
    bot <- evolve_bot(bot, mu, b, delta_t)
  }

  return(bot)
}


# Log probability of genetic information in i, given bottleneck probs infecting i
log_p_given_bot <- function(bot, freq, p, p_new_isnv, log_p_no_isnv){
  if(freq == 0){
    log(bot[1]) + log_p_no_isnv
  }else if(freq == 1){
    log(bot[3]) + log_p_no_isnv
  }else{
    log(
      bot[1] * p_new_isnv * denovo(freq, p) +
        bot[2] +
        bot[3] * p_new_isnv * denovo(1 - freq, p)
    )
  }
}

# Probabilities of successive split bottlenecks, starting with init (initial frequency) (can be a vector)
p_all_split <- function(b, w, init){
  (b^(w + 1) * 2 * init * (1 - init) / 3^w)
}


# Distribution of de novo iSNVs
denovo <- function(x, p, log = FALSE){
  k <- 1/sqrt(p)
  if(log){
    log(1-(1-x)^k * (1 + k*x)) - log(k) - 2*log(x)
  }else{
    (1-(1-x)^k * (1 + k*x)) / (k*x^2)
  }
}

# CDF of distribution of de novo iSNVs
denovo_cdf <- function(x, p){
  k <- 1/sqrt(p)
  ((1-x)^(k+1) + k*x + x - 1)/(k*x)
}

# Distribution of de novo iSNVs, normalized
denovo_normed <- function(x, p, filters, log = FALSE){
  k <- 1/sqrt(p)
  if(log){
    log(1-(1-x)^k * (1 + k*x)) - log(k) - 2*log(x) - log(1 - denovo_cdf(filters$af, p))
  }else{
    ((1-(1-x)^k * (1 + k*x)) / (k*x^2)) / (1 - denovo_cdf(filters$af, p))
  }
}

# Mean of denovo_normed



## Maximum time of infection for a host i
# If fix_child_seq = F, we allow max_t to go all the way up to the min of mcmc$seq[[j]][1] for j child of i
get_max_t <- function(mcmc, data, i, fix_child_seq = TRUE){
  js <- which(mcmc$h == i)
  ts <- c()
  for (j in js) {
    if(fix_child_seq){
      ts <- c(ts, mcmc$seq[[j]][length(mcmc$seq[[j]])])
    }else{
      ts <- c(ts, mcmc$seq[[j]][1])
    }
  }
  if(i <= data$n_obs){
    ts <- c(ts, data$s[i])
  }
  if(length(ts) == 0){
    return(Inf)
  }else{
    return(min(ts))
  }
}

## Sample number of hosts along an edge (i.e. length of seq)
# Minimum of 1
# Choose the best value with probability 0.95, otherwise random value
# Support is non-negative integers when fix_latest_host = TRUE
# Support is positive integers when fix_latest_host = FALSE
rw <- function(mcmc, min_t, max_t, fix_latest_host){
  if(runif(1) < 0.95){
    if(fix_latest_host){
      return(max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)) - 1), 0))
    }else{
      return(max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)) - 1), 1))
    }
  }else{
    if(fix_latest_host){
      return(rpois(1, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g))))
    }else{
      return(1 + rpois(1, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g))))
    }
  }
}

# LOG density of number of hosts along an edge
dw <- function(w, mcmc, min_t, max_t, fix_latest_host){
  if(fix_latest_host){
    ideal <- max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)) - 1), 0)
  }else{
    ideal <- max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)) - 1), 1)
  }

  if(w == ideal){
    if(fix_latest_host){
      return(
        log(
          0.95 + 0.05 * dpois(w, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)))
        )
      )
    }else{
      return(
        log(
          0.95 + 0.05 * dpois(w - 1, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)))
        )
      )
    }
  }else{
    if(fix_latest_host){
      return(
        log(0.05) + dpois(w, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)), log = TRUE)
      )
    }else{
      return(
        log(0.05) + dpois(w - 1, ((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)), log = TRUE)
      )
    }
  }
}

# Sample times of infection of hosts along an edge
rseq <- function(w, min_t, max_t, mcmc){

  if(w==0){
    return(numeric(0))
  }

  ## FOR NOW: just sorted uniform draw
  return(sort(runif(w, min_t, max_t), decreasing = T))

  # draw <- cumsum(rgamma(w + 1, shape = mcmc$a_g, rate = mcmc$lambda_g))
  # return(rev(draw / draw[length(draw)])[-1] * (max_t - min_t) + min_t)

}

# LOG Density of sampled times
# Note seq does not include the latest host if fix_latest_host = TRUE
dseq <- function(seq, w, min_t, max_t, mcmc){
  if(any(seq >= max_t) | any(seq <= min_t)){
    return(-Inf)
  }

  out <- 0
  if(length(seq) == 0){
    return(out)
  }

  ## FOR NOW: just sorted uniform density
  return(lfactorial(w) + w * log(1 / (max_t - min_t)))

  # # Record max_t - min_t, for change of density
  # delta_t <- max_t - min_t
  #
  # for (i in 1:w) {
  #   # Get the fraction of the remaining interval max_t - min_t taken up by seq[i], going back in time
  #   #print((max_t - seq[i]))
  #   #print(max_t - min_t)
  #   frac <- (max_t - seq[i]) / (max_t - min_t)
  #
  #   # This is the ratio of a single Gamma(a_g, lambda_g) draw divided by itself plus another w - i + 1 i.i.d. Gammas
  #   # Hence the ratio is Beta(a_g, a_g * (w - i + 1))
  #   # Final term is change of variables
  #   out <- out + dbeta(frac, mcmc$a_g, mcmc$a_g * (w - i + 1), log = T) + log(delta_t / (max_t - min_t))
  #
  #   # Finally we reset max_t to be seq[i], so that we calculate the remainder of the interval for next iteration in loop
  #   max_t <- seq[i]
  # }
  #
  # return(out)
}

# # Test this out with simple mcmc sampler
# w <- 3
# max_t <- 2
# samps <- list(max_t * seq(1 - 1/(w+1), 1/(w+1), -1/(w+1)))
# for (i in 2:10000) {
#   pr <- samps[[i-1]] + rnorm(w, 0, 0.1)
#   if(log(runif(1)) < dseq(pr, w, 0, max_t, mcmc) - dseq(samps[[i-1]], w, 0, max_t, mcmc)){
#     samps[[i]] <- pr
#   }else{
#     samps[[i]] <- samps[[i-1]]
#   }
# }
# hist(unlist(samps), breaks = seq(0,2,0.02))

# dseq(c(1.3, 0.9, 0.6), 3, 0, 2, mcmc)
#
# ds <- c()
# for (i in seq(10,15,0.05)) {
#   for (j in seq(10, 15, 0.05)) {
#     ds <- c(ds, dseq(c(i,j), 2, 10, 15, mcmc))
#   }
# }
# sum(exp(ds))

# dseq(c(0, 0.5), 2, 0, 1, mcmc)




## Softmax function, used for choosing arbitrary new ancestors
softmax <- function(v, tau){
  exp(v/tau) / sum(exp(v/tau))
}

## Score function: approximates the utility of attaching i to j in terms of parsimony
score <- function(mcmc, i, j){
  sum(mcmc$m01[[i]] %in% union(mcmc$m01[[j]], mcmc$m0x[[j]])) +
    sum(mcmc$m10[[i]] %in% union(mcmc$m10[[j]], mcmc$m1x[[j]])) +
    sum(union(mcmc$m0y[[i]], mcmc$m1y[[i]]) %in% union(mcmc$m0y[[j]], mcmc$m1y[[j]]))
}

## Path from i to j, going down then up
paths <- function(h, i, j){
  anc_i <- ancestry(h, i)
  anc_j <- ancestry(h, j)
  overlap <- length(intersect(anc_i, anc_j))
  return(list(
    rev(anc_i[overlap:length(anc_i)]),
    anc_j[overlap:length(anc_j)]
  ))
}


# Update genetics for the following topological move:
# From g -> i, g -> h
# To g -> h -> i
update_genetics_upstream <- function(mcmc, i, h){

  # Proposal
  prop <- mcmc

  # Everything that doesn't stay the same in i
  all_i <- unique(c(
    mcmc$m01[[i]],
    mcmc$m10[[i]],
    mcmc$m0y[[i]],
    mcmc$m1y[[i]],
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    mcmc$mxy[[i]]
  ))

  # Everything that doesn't stay the same in h
  all_h <- unique(c(
    mcmc$m01[[h]],
    mcmc$m10[[h]],
    mcmc$m0y[[h]],
    mcmc$m1y[[h]],
    mcmc$mx0[[h]],
    mcmc$mx1[[h]],
    mcmc$mxy[[h]]
  ))

  prop$m01[[i]] <- setdiff(mcmc$m01[[i]], all_h) # 00 in h, 01 in i
  prop$m01[[i]] <- union(prop$m01[[i]], intersect(mcmc$mx0[[h]], mcmc$mx1[[i]])) # x0 in h, x1 in i
  prop$m01[[i]] <- union(prop$m01[[i]], setdiff(mcmc$m10[[h]], all_i)) # 10 in h, 11 in i

  prop$m10[[i]] <- setdiff(mcmc$m10[[i]], all_h) # 11 in h, 10 in i
  prop$m10[[i]] <- union(prop$m10[[i]], intersect(mcmc$mx1[[h]], mcmc$mx0[[i]])) # x1 in h, x0 in i
  prop$m10[[i]] <- union(prop$m10[[i]], setdiff(mcmc$m01[[h]], all_i)) # 01 in h, 00 in i

  prop$m0y[[i]] <- setdiff(mcmc$m0y[[i]], all_h) # 00 in h, 0y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m10[[h]], mcmc$m1y[[i]])) # 10 in h, 1y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$mx0[[h]], mcmc$mxy[[i]])) # x0 in h, xy in i

  prop$m1y[[i]] <- setdiff(mcmc$m1y[[i]], all_h) # 11 in h, 1y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m01[[h]], mcmc$m0y[[i]])) # 01 in h, 0y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$mx1[[h]], mcmc$mxy[[i]])) # x1 in h, xy in i

  prop$mx0[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx0[[i]]) # xy in h, x0 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], setdiff(mcmc$m0y[[h]], all_i)) # 0y in h, 00 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], intersect(mcmc$m1y[[h]], mcmc$m10[[i]])) # 1y in h, 10 in i

  prop$mx1[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx1[[i]]) # xy in h, x1 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], setdiff(mcmc$m1y[[h]], all_i)) # 1y in h, 11 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], intersect(mcmc$m0y[[h]], mcmc$m01[[i]])) # 0y in h, 01 in i

  prop$mxy[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mxy[[i]]) # xy in h, xy in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$m1y[[h]], mcmc$m1y[[i]])) # 1y in h, 1y in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$m0y[[h]], mcmc$m0y[[i]])) # 0y in h, 0y in i

  # Whew.
  return(prop)
}

# Update genetics for the following topological move:
# From g -> h -> i
# To g -> i, g -> h
update_genetics_downstream <- function(mcmc, i, h){

  prop <- mcmc

  # Everything that doesn't stay the same in i
  all_i <- unique(c(
    mcmc$m01[[i]],
    mcmc$m10[[i]],
    mcmc$m0y[[i]],
    mcmc$m1y[[i]],
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    mcmc$mxy[[i]]
  ))

  # Everything that doesn't stay the same in h
  all_h <- unique(c(
    mcmc$m01[[h]],
    mcmc$m10[[h]],
    mcmc$m0y[[h]],
    mcmc$m1y[[h]],
    mcmc$mx0[[h]],
    mcmc$mx1[[h]],
    mcmc$mxy[[h]]
  ))

  prop$m01[[i]] <- setdiff(mcmc$m01[[h]], all_i) # 01 in h, 11 in i
  prop$m01[[i]] <- union(prop$m01[[i]], intersect(mcmc$m0y[[h]], mcmc$mx1[[i]])) # 0y in h, x1 in i
  prop$m01[[i]] <- union(prop$m01[[i]], setdiff(mcmc$m01[[i]], all_h)) # 00 in h, 01 in i

  prop$m10[[i]] <- setdiff(mcmc$m10[[h]], all_i) # 10 in h, 00 in i
  prop$m10[[i]] <- union(prop$m10[[i]], intersect(mcmc$m1y[[h]], mcmc$mx0[[i]])) # 1y in h, x0 in i
  prop$m10[[i]] <- union(prop$m10[[i]], setdiff(mcmc$m10[[i]], all_h)) # 11 in h, 10 in i

  prop$m0y[[i]] <- setdiff(mcmc$m0y[[i]], all_h) # 00 in h, 0y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m01[[h]], mcmc$m1y[[i]])) # 01 in h, 1y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m0y[[h]], mcmc$mxy[[i]])) # 0y in h, xy in i

  prop$m1y[[i]] <- setdiff(mcmc$m1y[[i]], all_h) # 11 in h, 1y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m10[[h]], mcmc$m0y[[i]])) # 10 in h, 0y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m1y[[h]], mcmc$mxy[[i]])) # 1y in h, xy in i

  prop$mx0[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx0[[i]]) # xy in h, x0 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], setdiff(mcmc$mx0[[h]], all_i)) # x0 in h, 00 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], intersect(mcmc$mx1[[h]], mcmc$m10[[i]])) # x1 in h, 10 in i

  prop$mx1[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx1[[i]]) # xy in h, x1 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], setdiff(mcmc$mx1[[h]], all_i)) # x1 in h, 11 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], intersect(mcmc$mx0[[h]], mcmc$m01[[i]])) # x0 in h, 01 in i

  prop$mxy[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mxy[[i]]) # xy in h, xy in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$mx1[[h]], mcmc$m1y[[i]])) # x1 in h, 1y in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$mx0[[h]], mcmc$m0y[[i]])) # x0 in h, 0y in i

  # Whew.
  return(prop)
}

# Wrap as a function: switch from
# h_old -> i, h_old -> h_new to
# h_old -> h_new -> i
shift_upstream <- function(mcmc, data, i, h_old, h_new){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  mcmc <- update_genetics_upstream(mcmc, i, h_new) # Update genetics. i is inheriting from h_new.
  return(mcmc)
}

# Wrap as a function: switch from
# h_new -> h_old -> i
# h_old -> i, h_new -> i
shift_downstream <- function(mcmc, data, i, h_old, h_new){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  mcmc <- update_genetics_downstream(mcmc, i, h_old) # Update genetics. i is inheriting from h_new, but compared to genetics of h_old
  return(mcmc)
}

# Flip the genotype for a SNV
flip_genotype <- function(mcmc, i, js, snv){
  prop <- mcmc
  ## Run through cases of updating genetic info in i
  if(snv %in% mcmc$m01[[i]]){

    # Delete from 01 in i
    prop$m01[[i]] <- setdiff(mcmc$m01[[i]], snv)

    # Note that we're changing from 1 to 0 in i
    add <- FALSE
  }else if(snv %in% mcmc$mx1[[i]]){

    # Delete from x1 in i
    prop$mx1[[i]] <- setdiff(mcmc$mx1[[i]], snv)

    # Union to x0 in i
    prop$mx0[[i]] <- union(mcmc$mx0[[i]], snv)

    # Note that we're changing from 1 to 0 in i
    add <- FALSE
  }else if(snv %in% mcmc$m10[[i]]){

    # Delete from 01 in i
    prop$m10[[i]] <- setdiff(mcmc$m10[[i]], snv)

    # Note that we're changing from 0 to 1 in i
    add <- TRUE
  }else if(snv %in% mcmc$mx0[[i]]){

    # Delete from x1 in i
    prop$mx0[[i]] <- setdiff(mcmc$mx0[[i]], snv)

    # Union to x0 in i
    prop$mx1[[i]] <- union(mcmc$mx1[[i]], snv)

    # Note that we're changing from 0 to 1 in i
    add <- TRUE
  }else if(snv %in% c(unlist(mcmc$m10[js]), unlist(mcmc$m1y[js]))){ # 11 in i

    # Union to 10 in i
    prop$m10[[i]] <- union(mcmc$m10[[i]], snv)

    add <- FALSE
  }else if(snv %in% c(unlist(mcmc$m01[js]), unlist(mcmc$m0y[js]))){ # 00 in i

    # Union to 01 in i
    prop$m01[[i]] <- union(mcmc$m01[[i]], snv)

    add <- TRUE
  }else{
    # In the final case, we need to search the ancestry of i to determine whether the SNV is present or absent
    h <- mcmc$h[i]
    add <- T
    while (h != 1) {
      if(snv %in% c(unlist(mcmc$m01[[h]]), unlist(mcmc$mx1[[h]]))){
        add <- F
        break
      }else{
        h <- mcmc$h[h]
      }
    }
    if(add){
      # Union to 01 in i
      prop$m01[[i]] <- union(mcmc$m01[[i]], snv)
    }else{
      # Union to 10 in i
      prop$m10[[i]] <- union(mcmc$m10[[i]], snv)
    }
  }

  ## Now update genetic info for j in js
  if(add){
    for (j in js) {
      if(snv %in% mcmc$m01[[j]]){
        # Delete from 01 in j
        prop$m01[[j]] <- setdiff(mcmc$m01[[j]], snv)
      }else if(snv %in% mcmc$m0y[[j]]){
        # Delete from 0y in j
        prop$m0y[[j]] <- setdiff(mcmc$m0y[[j]], snv)
        # Add to 1y in j
        prop$m1y[[j]] <- union(mcmc$m1y[[j]], snv)
      }else{
        # Otherwise, it was in 00 in j
        prop$m10[[j]] <- union(mcmc$m10[[j]], snv)
      }
    }
  }else{
    for (j in js) {
      if(snv %in% mcmc$m10[[j]]){
        # Delete from 10 in j
        prop$m10[[j]] <- setdiff(mcmc$m10[[j]], snv)
      }else if(snv %in% mcmc$m0y[[j]]){
        # Delete from 1y in j
        prop$m1y[[j]] <- setdiff(mcmc$m1y[[j]], snv)
        # Add to 0y in j
        prop$m0y[[j]] <- union(mcmc$m0y[[j]], snv)
      }else{
        # Otherwise, it was in 11 in j
        prop$m01[[j]] <- union(mcmc$m01[[j]], snv)
      }
    }
  }

  return(prop)
}

# Change the genotype
change_genotype <- function(mcmc, data, snv, from, to, i, js, is_observed){
  if(from != to){
    # Figure out whether the snv is absent, present, or isnv in h[i]
    if(snv %in% c(mcmc$m01[[i]], mcmc$m0y[[i]])){
      old_state_anc <- "0"
    }else if(snv %in% c(mcmc$mx0[[i]], mcmc$mxy[[i]], mcmc$mx1[[i]])){
      old_state_anc <- "x"
    }else if(snv %in% c(mcmc$m10[[i]], mcmc$m1y[[i]])){
      old_state_anc <- "1"
    }else if(from == "absent"){
      old_state_anc <- "0"
    }else if(from == "present"){
      old_state_anc <- "1"
    }else{
      stop("Error!!!")
    }

    if(length(js) > 0){
      # Figure out whether the snv is absent, present, or isnv in each j
      old_state_kids <- c()
      for (j in 1:length(js)) {
        if(snv %in% c(mcmc$mx0[[js[j]]], mcmc$m10[[js[j]]])){
          old_state_kids[j] <- "0"
        }else if(snv %in% c(mcmc$m0y[[js[j]]], mcmc$mxy[[js[j]]], mcmc$m1y[[js[j]]])){
          old_state_kids[j] <- "y"
        }else if(snv %in% c(mcmc$m01[[js[j]]], mcmc$mx1[[js[j]]])){
          old_state_kids[j] <- "1"
        }else if(from == "absent"){
          old_state_kids[j] <- "0"
        }else if(from == "present"){
          old_state_kids[j] <- "1"
        }else{
          stop("Error!!!!!!")
        }
      }
    }



    # From which list do we delete in i and js?
    if(from == "absent"){
      delete_i <- paste0("m", old_state_anc, "0")
      if(length(js) > 0){
        delete_js <- paste0("m", "0", old_state_kids)
      }
    }else if(from == "isnv"){
      delete_i <- paste0("m", old_state_anc, "y")
      if(length(js) > 0){
        delete_js <- paste0("m", "x", old_state_kids)
      }
    }else if(from == "present"){
      delete_i <- paste0("m", old_state_anc, "1")
      if(length(js) > 0){
        delete_js <- paste0("m", "1", old_state_kids)
      }
    }

    # To which list do we add in i and js?
    if(to == "absent"){
      add_i <- paste0("m", old_state_anc, "0")
      if(length(js) > 0){
        add_js <- paste0("m", "0", old_state_kids)
      }
    }else if(to == "isnv"){
      add_i <- paste0("m", old_state_anc, "y")
      if(length(js) > 0){
        add_js <- paste0("m", "x", old_state_kids)
      }
    }else if(to == "present"){
      add_i <- paste0("m", old_state_anc, "1")
      if(length(js) > 0){
        add_js <- paste0("m", "1", old_state_kids)
      }
    }

    # Modify the lists
    if(!is.null(mcmc[[delete_i]])){
      mcmc[[delete_i]][[i]] <- setdiff(mcmc[[delete_i]][[i]], snv)
    }
    if(!is.null(mcmc[[add_i]])){
      mcmc[[add_i]][[i]] <- union(mcmc[[add_i]][[i]], snv)
    }
    if(length(js) > 0){
      for (j in 1:length(js)) {
        if(!is.null(mcmc[[delete_js[j]]])){
          mcmc[[delete_js[j]]][[js[j]]] <- setdiff(mcmc[[delete_js[j]]][[js[j]]], snv)
        }
        if(!is.null(mcmc[[add_js[j]]])){
          mcmc[[add_js[j]]][[js[j]]] <- union(mcmc[[add_js[j]]][[js[j]]], snv)
        }
      }
    }
  }

  # If from == "isnv", need to delete this isnv from mcmc$isnv
  if(from == "isnv"){
    # Index of isnv on mcmc$isnv
    ind <- match(snv, mcmc$isnv$call[[i]])
    mcmc$isnv$call[[i]] <- mcmc$isnv$call[[i]][-ind]
    mcmc$isnv$af[[i]] <- mcmc$isnv$af[[i]][-ind]
  }

  # Similarly, if to == "isnv", need to add to the list, AND pick its frequency
  if(to == "isnv"){
    mcmc$isnv$call[[i]] <- c(mcmc$isnv$call[[i]], snv)

    if(is_observed){
      if(data$vcf_present[i]){
        if(from == "absent"){
          mcmc$isnv$af[[i]] <- c(mcmc$isnv$af[[i]], runif(1, 0, data$filters$af))
        }else if(from == "present"){
          mcmc$isnv$af[[i]] <- c(mcmc$isnv$af[[i]], runif(1, 1 - data$filters$af, 1))
        }
        log_p <- log(1 / data$filters$af)
      }else{
        if(from == "absent"){
          mcmc$isnv$af[[i]] <- c(mcmc$isnv$af[[i]], runif(1, 0, 1/2))
        }else if(from == "present"){
          mcmc$isnv$af[[i]] <- c(mcmc$isnv$af[[i]], runif(1, 1/2, 1))
        }
        log_p <- log(2) # log(1 / 0.5)
      }
    }else{
      mcmc$isnv$af[[i]] <- c(mcmc$isnv$af[[i]], runif(1))
      log_p <- 0 # log(1)
    }
  }else{
    log_p <- 0
  }

  return(list(mcmc, log_p))
}

## Is a snv absent, isnv, or present in a host i with children js?
snv_status <- function(mcmc, i, js, snv){
  if(snv %in% c(
    mcmc$m0y[[i]],
    mcmc$mxy[[i]],
    mcmc$m1y[[i]]
  )){
    from <- "isnv"
  }else if(snv %in% c(
    mcmc$mx0[[i]],
    mcmc$m10[[i]],
    unlist(mcmc$m0y[js]),
    unlist(mcmc$m01[js])
  )){
    from <- "absent"
  }else if(snv %in% c(
    mcmc$mx1[[i]],
    mcmc$m01[[i]],
    unlist(mcmc$m1y[js]),
    unlist(mcmc$m10[js])
  )){
    from <- "present"
  }else{
    # In this case, we need to trace the ancestry until we find the most recent m01 or m10
    from <- "absent"
    h <- mcmc$h[i]
    while (!is.na(h)) {
      if(snv %in% c(unlist(mcmc$m01[[h]]), unlist(mcmc$mx1[[h]]))){
        from <- "present"
        break
      }else{
        h <- mcmc$h[h]
      }
    }
  }
  return(from)
}



# Accept / reject
accept_or_reject <- function(prop, mcmc, data, update, hastings = 0, check_parsimony = integer(0), noisy = FALSE){

  if(length(check_parsimony) > 0){
    for (k in check_parsimony) {
      if(
        !genotype(prop, data, k, check_parsimony = T)
      ){
        return(mcmc)
      }
    }
  }

  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  # Accept / reject
  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
    if(noisy){
      print("Move accepted")
    }
    return(prop)
  }else{
    return(mcmc)
  }
}


## Get list of all nodes upstream from a given node (including indirectly)
# We can do this using recursion!
get_upstream <- function(h, i){
  out <- which(h == i)
  for (j in out) {
    out <- c(out, get_upstream(h, j))
  }
  return(out)
}

## Efficiently compute total number of upstream nodes for each node (including self)
total_degree <- function(h, d){
  n <- length(h)
  out <- rep(1, length(h))
  frontier <- which(d == 0)

  while (length(frontier) > 0 & !identical(frontier, 1)) {
    new_frontier <- c()
    for (i in frontier) {
      out[h[i]] <- out[h[i]] + out[i] # Back up degree into parent
      d[h[i]] <- d[h[i]] - 1 # Prune child
      if(d[h[i]] == 0){
        if(!is.na(h[i])){
          new_frontier <- c(new_frontier, h[i])
        }
      }
    }
    frontier <- new_frontier
    if(length(frontier) == 1){
      if(frontier == 1){
        frontier <- integer(0)
      }
    }
  }

  return(out)
}

# BFS traversal of tree
bfs <- function(i, h){
  out <- i
  frontier <- which(h == i)
  while (length(frontier) > 0) {
    out <- c(out, frontier)
    frontier <- which(h %in% frontier)
  }
  return(out)
}

# DFS traversal of tree, lowest number first
dfs <- function(h){
  stack <- 1
  explored <- c()
  while (length(stack) > 0) {
    who <- sort(which(h == stack[1]))
    explored <- c(explored, stack[1])
    stack <- stack[-1]
    stack <- c(who, stack)
  }
  return(explored)
}

# Get generation of each node
gen2 <- function(mcmc){
  ord <- bfs(1,mcmc$h)
  out <- rep(NA, mcmc$n)
  for (i in ord) {
    if(i == 1){
      out[i] <- 0
    }else{
      out[i] <- out[mcmc$h[i]] + mcmc$w[i] + 1
    }
  }
  return(out)
}

chop <- function(mcmc, data, old_roots){

  # Initial tree (will change)
  h <- mcmc$h

  # Node degrees
  d <- sapply(1:mcmc$n, function(i){sum(mcmc$h[2:mcmc$n] == i)})


  # Traverse the tree in reverse-BFS order
  ord <- rev(bfs(1, h))

  # Minimum number of nodes per subtree
  lambda <- mcmc$n / data$n_subtrees

  # Tree outputs (not including roots)
  trees <- list()

  # Root outputs
  roots <- c()

  # All upstream nodes, not including self
  w <- rep(0, mcmc$n)

  for (v in ord) {
    if(v == 1){

      sub <- bfs(v, h)

      trees <- c(trees, list(sort(sub[-1])))
      roots <- c(roots, v)
    }else{
      # Update number of upstream nodes of vertex v
      if(d[v] > 0){
        kids <- which(h == v)
        if(length(kids) > 0){
          w[v] <- w[v] + length(kids) + sum(w[kids])
        }
      }
      # If weight is large enough, and root is observed, and root isn't previous root, hack off a piece of the tree
      if(
        w[v] >= lambda &
        #v <= data$n_obs &
        !(v %in% old_roots)
      ){
        if(mcmc$n - length(unlist(trees)) - w[v] >= lambda){

          sub <- bfs(v, h)

          trees <- c(trees, list(sort(sub[-1])))
          roots <- c(roots, v)

          # Delete nodes from tree, except root
          ## CHECK kids is correct
          h[kids] <- NA

          # Reset upstream nodes of root to nothing
          w[v] <- 0
        }
      }
    }
  }

  return(list(roots, trees))
}






## For debugging: consensus changes from root
cc_from_root <- function(mcmc, i){
  if(i == 1){
    return(0)
  }else{
    anc <- ancestry(mcmc$h, i)
    count <- 0
    for (j in anc[2:length(anc)]) {
      count <- count + length(mcmc$m01[[j]]) + length(mcmc$mx1[[j]]) - length(mcmc$m10[[j]]) - length(mcmc$m1y[[j]])
    }
    return(count)
  }
}






