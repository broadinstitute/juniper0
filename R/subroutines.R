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

    # 1st allele, alphabetically
    out$isnv$a1 <- character(0)
    # 2nd allele, alphabetically
    out$isnv$a2 <- character(0)
    # Mutated proportion of first allele
    out$isnv$af1 <- numeric(0)
    # Position
    out$isnv$pos <- integer(0)

    for (p in unique(pos)) {
      out$isnv$pos <- c(out$isnv$pos, p)
      props <- rep(0, 4)
      rows <- which(pos == p)
      afs <- af[rows]
      alts <- alt[rows]
      refs <- ref[rows]
      props[match(alts, c("A", "C", "G", "T"))] <- afs
      props[match(refs[1], c("A", "C", "G", "T"))] <- 1 - sum(afs)

      # Which props are nonzero
      nonzero <- which(props > 0)
      out$isnv$a1 <- c(out$isnv$a1, c("A", "C", "G", "T")[nonzero[1]])
      out$isnv$a2 <- c(out$isnv$a2, c("A", "C", "G", "T")[nonzero[2]])
      out$isnv$af1 <- c(out$isnv$af1, props[nonzero[1]])

      if(sum(props > 0) > 2){
        print("Triallelic")
        print(props)
      }
    }
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

  old <- raw_to_base(seq1[snv_pos])
  new <- raw_to_base(seq2[snv_pos])

  out$snv$from <- old
  out$snv$pos <- snv_pos
  out$snv$to <- new
  out$missing <- missing_pos

  return(out)

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

# Get the ancestry of a single node, down to the root
ancestry <- function(h, i){
  if(is.na(h[i])){
    return(i)
  }else{
    return(c(ancestry(h, h[i]), i))
  }
}

# Get the MRCA node of a bunch of nodes
mrca <- function(h, js){
  # Ancestry of each j
  ancs <- lapply(js, ancestry, h=h)

  # Index of the ancestor we will try
  ind <- 1

  # Does this work as the MRCA?
  works <- TRUE
  while (works) {
    if(all(sapply(ancs, function(v){ancs[[1]][ind + 1] %in% v}))){
      ind <- ind + 1
    }else{
      works <- FALSE
    }
  }
  return(ancs[[1]][ind])
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
get_ts <- function(mcmc, data, i){

  max_t <- get_max_t(mcmc, data, i)
  min_t <- mcmc$seq[[mcmc$h[i]]][1]

  # Number of hosts along the edge
  n_hosts <- max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g)) - 1), 1)

  delta_t <- mcmc$t[i] - mcmc$t[mcmc$h[i]]
  # Increment back in time
  inc <- (max_t - min_t) / (n_hosts + 1)
  seq(max_t - inc, by = -inc, length.out = n_hosts)

}

# Get total evolutionary time
evo_time <- function(i, mcmc){
  mcmc$seq[[i]][1] - mcmc$seq[[mcmc$h[i]]][1]
}

tot_evo_time <- function(mcmc){
  sum(sapply(2:mcmc$n, evo_time, mcmc = mcmc))
}

## Get optimal values of R and pi, approximately

# Objective function of R and pi
obj <- function(params, t, a_g, lambda_g){

  #print(params)
  #print(optimize_tmrca)

  R <- params[1]
  pi <- params[2]

  if(pi <= 0 | pi > 1 | R <= 0){
    return(Inf)
  }

  max_t <- max(t)

  # Generation interval
  g <- a_g / lambda_g

  # Number of generations
  G <- max_t / g

  # Get rid of last time to be infected, since we say the process stops there
  t <- t[-which.max(t)]

  if(length(t) > round((R^(G+1) - 1) / (R-1))){
    return(Inf)
  }

  out <- dbinom(length(t), round((R^(G+1) - 1) / (R-1)), pi, log = T)

  if(is.nan(out)){
    print(params)
  }

  # Population growth curve is N(t) = exp(t*log(R))
  # Normalizing by its integral to get a PDF, we have f(t) = (exp(t*log(R)) * log(R)) / (exp(max_t * log(R)) - 1)
  # Taking the log...
  #out <- 0
  out <- out + sum(log(
    R^(t/g) * log(R) / ((R^(max_t/g) - 1) * g)
  ))

  return(-out) # Because optimize function finds minimum
}

opt_R_pi <- function(s, a_g, lambda_g, a_s, lambda_s){

  # Delete NAs
  s <- s[!is.na(s)]

  # Shift date of sample collection back by mean sojourn interval to get approximate time of infection
  t <- s - (a_s / lambda_s)

  init_vals <- c(2, 0.5)
  vals <- optim(init_vals, obj, t = t, a_g = a_g, lambda_g = lambda_g)

  return(vals$par)

}




# Probability density function for the waiting time to a denovo iSNV at a given site
ddenovo <- function(x, mu, N_eff, log){
  if(log){
    mu - exp(N_eff*x)*mu + N_eff*x + log(N_eff) + log(mu)
  }else{
    exp(mu - exp(N_eff*x)*mu + N_eff*x) * N_eff * mu
  }
}

# Cumulative density function
pdenovo <- function(x, mu, N_eff, log){
  out <- 1 - exp(-mu * (exp(N_eff *x) - 1))
  if(log){
    return(log(out))
  }else{
    return(out)
  }
}

# Survival function of denovo frequencies
sdenovo <- function(x, mu, N_eff, log){
  if(log){
    -mu * (exp(N_eff *x) - 1)
  }else{
    exp(-mu * (exp(N_eff *x) - 1))
  }
}

# Marginal PDF of frequencies of denovo iSNVs
dprop <- function(x, mu, log){
  if(log){
    log(mu) - 2 * log(mu + x - mu*x)
  }else{
    mu/(mu + x - mu *x)^2
  }
}

# CDF
pprop <- function(x, mu, log){
  if(log){
    log(x) - log(mu + x - mu*x)
  }else{
    x/(mu + x - mu * x)
  }
}

# Survival
sprop <- function(x, mu, log){
  if(log){
    log(mu - mu * x) - log(mu + x - mu*x)
  }else{
    (mu - mu * x)/(mu + x - mu * x)
  }
}

# Probability density that the iSNV frequency is x AND the first substituion occurred by the time the viral population is N
dprop_bounded <- function(x, N, mu, log){
  if(log){
    log(
      (mu + (1 - mu)^N * mu * (1 - x)^N * (-1 + mu * N * (-1 + x) - N * x))
    ) -
      2 * log(
        (mu + x - mu * x)
      )
  }else{
    (mu + (1 - mu)^N * mu * (1 - x)^N * (-1 + mu * N * (-1 + x) - N * x))/(mu + x - mu * x)^2
  }
}

# Probability that the iSNV frequency doesn't exceed x AND the first substituion occurred by the time the viral population is N
pprop_bounded <- function(x, N, mu, log){
  if(log){
    log(
      (-x + (1 - mu)^N * (mu * (-1 + (1 - x)^N) * (-1 + x) + x))
    ) -
      log(
        (mu * (-1 + x) - x)
      )
  }else{
    (-x + (1 - mu)^N * (mu * (-1 + (1 - x)^N) * (-1 + x) + x))/(mu * (-1 + x) - x)
  }
}

# Draw Beta(1, (1-mu)/mu) (mean mu), then rescaled to [t_min, t_max]
rbeta1_rescaled <- function(n, mu, t_min, t_max){
  rbeta(n, 1, (1-mu)/mu) * (t_max - t_min) + t_min
}

# Density for the above
dbeta1_rescaled <- function(x, mu, t_min, t_max, log){
  out <- dbeta((x - t_min) / (t_max - t_min), 1, (1-mu)/mu, log = log)
  if(log){
    return(out - log(t_max - t_min))
  }else{
    return(out / (t_max - t_min))
  }
}

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



## Softmax function, used for choosing arbitrary new ancestors
softmax <- function(v, tau){
  exp(v/tau) / sum(exp(v/tau))
}

## Score function: approximates the utility of attaching i to j in terms of parsimony
score <- function(mcmc, i, j){
  length(intersect(mcmc$subs$pos[[i]], mcmc$subs$pos[[j]]))
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
# Upstream:
# From g -> i, g -> h
# To g -> h -> i

# Downstream:
# From g -> h -> i
# To g -> i, g -> h
update_genetics <- function(mcmc, i, h, upstream){

  if(upstream){
    # These are reversed, because when we first move i onto h, we have to delete everything that was added from g to h, then adjust based on what's added h to i
    from <- mcmc$subs$to[[h]]
    pos <- mcmc$subs$pos[[h]]
    to <- mcmc$subs$from[[h]]
  }else{
    from <- mcmc$subs$from[[h]]
    pos <- mcmc$subs$pos[[h]]
    to <- mcmc$subs$to[[h]]
  }

  if(length(mcmc$subs$pos[[i]]) > 0){
    for (j in 1:length(mcmc$subs$pos[[i]])) {
      if(mcmc$subs$pos[[i]][j] %in% pos){
        ind <- which(pos == mcmc$subs$pos[[i]][j])
        # Update the "to" allele
        to[ind] <- mcmc$subs$to[[i]][j]
      }else{
        from <- c(from, mcmc$subs$from[[i]][j])
        pos <- c(pos, mcmc$subs$pos[[i]][j])
        to <- c(to, mcmc$subs$to[[i]][j])
      }
    }
  }

  # Clear out SNVs where from == to
  keep <- from != to

  mcmc$subs$from[[i]] <- from[keep]
  mcmc$subs$pos[[i]] <- pos[keep]
  mcmc$subs$to[[i]] <- to[keep]

  return(mcmc)

}



# Wrap as a function: switch from
# Upstream = T
# h_old -> i, h_old -> h_new to
# h_old -> h_new -> i
# Upstream = F
# h_new -> h_old -> i
# h_old -> i, h_new -> i


shift <- function(mcmc, data, i, h_old, h_new, upstream){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  # For update genetics upstream, h represents the new ancestor
  if(upstream){
    mcmc <- update_genetics(mcmc, i, h_new, TRUE) # Update genetics. i is inheriting from h_new.
  }else{
    # For update genetics downstream, h represents the old ancestor
    mcmc <- update_genetics(mcmc, i, h_old, FALSE)
  }

  return(mcmc)
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
  if(is.infinite(hastings)){
    stop("hastings error")
  }

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






