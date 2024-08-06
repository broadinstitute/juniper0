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
get_max_t <- function(mcmc, data, i){
  ts <- mcmc$t[which(mcmc$h == i)]
  if(i <= data$n_obs){
    ts <- c(ts, data$s[i])
  }
  if(length(ts) == 0){
    return(Inf)
  }else{
    return(min(ts))
  }
}

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
update_genetics_upstream <- function(prop, mcmc, i, h){
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
update_genetics_downstream <- function(prop, mcmc, i, h){
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
shift_upstream <- function(mcmc, data, i, h_old, h_new, resample_t = FALSE, resample_w = FALSE){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  if(resample_t){
    max_t <- get_max_t(mcmc, data, i)
    min_t <- mcmc$t[h_new]
    mcmc$t[i] <- runif(1, min_t, max_t)
  }
  if(resample_w){
    mcmc$w[i] <- rpois(1, (mcmc$t[i] - mcmc$t[h_new]) * mcmc$lambda_g / mcmc$a_g) # Biased sample, but hopefully good enough
  }else{
    mcmc$w[i] <- mcmc$w[i] - mcmc$w[h_new] - 1
  }
  mcmc <- update_genetics_upstream(mcmc, mcmc, i, h_new) # Update genetics. i is inheriting from h_new.
  mcmc$d[h_old] <- mcmc$d[h_old] - 1
  mcmc$d[h_new] <- mcmc$d[h_new] + 1
  return(mcmc)
}

# Wrap as a function: switch from
# h_new -> h_old -> i
# h_old -> i, h_new -> i
shift_downstream <- function(mcmc, data, i, h_old, h_new, resample_t = FALSE, resample_w = FALSE){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  if(resample_t){
    max_t <- get_max_t(mcmc, data, i)
    min_t <- mcmc$t[h_new]
    mcmc$t[i] <- runif(1, min_t, max_t)
  }
  if(resample_w){
    mcmc$w[i] <- rpois(1, (mcmc$t[i] - mcmc$t[h_new]) * mcmc$lambda_g / mcmc$a_g) # Biased sample, but hopefully good enough
  }else{
    mcmc$w[i] <- mcmc$w[i] + mcmc$w[h_old] + 1
  }
  mcmc <- update_genetics_downstream(mcmc, mcmc, i, h_old) # Update genetics. i is inheriting from h_new, but compared to genetics of h_old
  mcmc$d[h_old] <- mcmc$d[h_old] - 1
  mcmc$d[h_new] <- mcmc$d[h_new] + 1
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

## Resample the genotype for an unobserved host, or for an observed host with missing sites, based on (approximate) parsimony
genotype <- function(mcmc, data, i, js, comparison = F){
  # Get all SNVs
  snvs <- data$all_snv

  # If i is observed, the only positions that can change are those with missing data
  if(i <= data$n_obs){
    snvs <- intersect(snvs, data$snvs[[i]]$missing$call)
  }

  # Which ones go 0y or 1y in j, or x0 or x1 in i? (with multiplicity)
  isnvs <- c(
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m0y[js]),
    unlist(mcmc$m1y[js])
  )
  # Ensure these are valid mutations to flip
  if(i <= data$n_obs){
    isnvs <- isnvs[isnvs %in% snvs]
  }

  tab_isnv <- table(isnvs)
  ind_isnv <- match(names(tab_isnv), snvs) # Indices of these iSNVs in "snvs"

  # And the others?
  non_isnvs <- c(
    mcmc$m10[[i]],
    mcmc$m01[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m10[js])
  )
  # Ensure these are valid mutations to flip
  if(i <= data$n_obs){
    non_isnvs <- non_isnvs[non_isnvs %in% snvs]
  }
  tab_non_isnv <- table(non_isnvs)
  ind_non_isnv <- match(names(tab_non_isnv), snvs) # Indices of these SNVs in "snvs"

  # Probability of swapping
  probs <- rep(0, length(snvs))
  probs[ind_non_isnv] <- probs[ind_non_isnv] + unname(tab_non_isnv)
  probs[ind_isnv] <- probs[ind_isnv] + unname(tab_isnv)/2
  probs <- probs / (length(js) + 1)

  # By parsimony, round probability
  probs[probs < 0.5] <- 0
  probs[probs > 0.5] <- 1

  # Add random noise
  probs <- (data$eps/2) + (1 - data$eps)*probs

  # If our goal is to compute the probability that a new host i has this genotype...
  if(comparison){
    # The probability of creating a new genotype at i equal to that in MCMC is the probability we don't swap anything in i
    return(sum(log(1-probs)))
  }else{
    # Which ones get swapped?
    which_swap <- runif(length(snvs)) < probs
    swaps <- snvs[which_swap]

    # Log probability of picking this genotype
    log_p <- sum(log(probs[which_swap])) + sum(log(1 - probs[!which_swap]))

    for (snv in swaps) {
      mcmc <- flip_genotype(mcmc, i, js, snv)
    }

    return(list(mcmc, log_p))
  }
}

# Accept / reject
accept_or_reject <- function(prop, mcmc, data, update, hastings = 0){
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  # Accept / reject
  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
    #print("coo-ee")
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
  d <- mcmc$d

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






