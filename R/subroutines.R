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

is_date <- function(x) {
  !is.na(tryCatch(as.Date(x), error = function(e) NA))
}

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

genetic_info <- function(seq_name, seq1, seq2, filters, vcf = NULL){
  # List output: list of SNVs, iSNVs, and positions with no information
  out <- list()

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
    keep <- which(dp >= filters$dp & sb < filters$sb & af >= filters$af & af <= 1 - filters$af & !(pos %in% filters$common) & !(pos %in% missing_pos))
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

      if(length(nonzero) < 2){
        #print(i)
        print(vcf)
        print(p)
        stop("???")
      }

      # Multiallelic case: mask alleles that aren't in the top 2 frequencies
      if(length(nonzero) > 2){
        # Get the indices of props going from biggest to smallest
        allele_ord <- sort.int(props, decreasing = T, index.return = T)$ix
        nonzero <- allele_ord[1:2]
      }

      out$isnv$a1 <- c(out$isnv$a1, c("A", "C", "G", "T")[nonzero[1]])
      out$isnv$a2 <- c(out$isnv$a2, c("A", "C", "G", "T")[nonzero[2]])
      out$isnv$af1 <- c(out$isnv$af1, props[nonzero[1]])

      # if(sum(props > 0) > 2){
      #   print("Triallelic")
      #   print(props)
      # }
    }
  }



  if(any(out$isnv$pos %in% out$missing)){
    stop("There are sites reported as missing in the FASTA and non-missing in the corresponding VCF")
  }

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

  max_t <- mcmc$seq[[i]][1]
  min_t <- mcmc$seq[[mcmc$h[i]]][1]

  # Number of hosts along the edge
  n_hosts <- max(round(((max_t - min_t) / (mcmc$a_g / mcmc$lambda_g))), 1)

  # Increment back in time
  seq(max_t, min_t, length.out = n_hosts + 1)[1:n_hosts]

}

# Get total evolutionary time
evo_time <- function(i, mcmc, data){
  # If unrooted and i is the only child of the root and i is unobserved, no contribution
  if(!data$rooted){
    if(i > data$n_obs){
      if(mcmc$h[i] == 1){
        if(length(which(mcmc$h == 1)) == 1){
          return(0)
        }
      }
    }
  }

  return(
    (mcmc$seq[[i]][1] - mcmc$seq[[mcmc$h[i]]][1]) * (data$n_bases - length(mcmc$dropout[[i]]))
  )
}

tot_evo_time <- function(mcmc, data){
  sum(sapply(2:mcmc$n, evo_time, mcmc = mcmc, data = data))
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



  # Now, the initial value of R needs to be large enough such that the expected cumulative number of cases generated over the epidemic exceeds the number of samples
  # Max time
  max_t <- max(t)
  # Generation interval
  g <- a_g / lambda_g
  # Number of generations
  G <- max_t / g

  # Number of total cases as a function of R: (R^(G+1) - 1) / (R-1)
  # This needs to exceed length(t) - 1
  # Can show: for G > 0, (R^(G+1) - 1) / (R-1) > R^(G+1) / R = R^G
  # So suffices to pick R such that R^G > length(t) - 1
  init_R <- (length(t) - 1)^(1/G)

  init_vals <- c(init_R, 0.5)
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
  if(i <= data$n_obs & i != 1){
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
  #return(sort(runif(w, min_t, max_t), decreasing = T))

  draw <- extraDistr::rdirichlet(1, rep(mcmc$a_g / mcmc$lambda_g, w + 1))[1,]
  draw <- cumsum(draw[1:(length(draw) - 1)])
  draw <- draw * (max_t - min_t) + min_t

  return(rev(draw))

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
  #return(lfactorial(w) + w * log(1 / (max_t - min_t)))

  # Record max_t - min_t, for change of density
  delta_t <- max_t - min_t

  seq <- rev((seq - min_t) / (max_t - min_t))
  seq <- diff(c(0,seq))
  seq <- c(seq, 1 - sum(seq))

  return(
    extraDistr::ddirichlet(seq, rep(mcmc$a_g / mcmc$lambda_g, w + 1), log = T) + w * log(1 / (max_t - min_t))
  )

}

## Draw from density that's uniform above 0 to max_delta, exponential below 0 with conditional mean mu_delta, and such that the PDF is continuous on (-Inf, max_delta)
# Solve for probability p of being negative
# p / mu_delta = (1 - p) / max_delta
# p = mu_delta / (mu_delta + max_delta)
rdelta <- function(mu_delta, max_delta){
  p <- mu_delta / (mu_delta + max_delta)
  if(runif(1) < p){
    return(-rexp(1, 1 / mu_delta))
  }else{
    return(runif(1, 0, max_delta))
  }
}

ddelta <- function(delta, mu_delta, max_delta, log){
  p <- mu_delta / (mu_delta + max_delta)
  if(log){
    if(delta <= 0){
      log(p) + dexp(-delta, 1 / mu_delta, log = TRUE)
    }else{
      log(1-p) + log(1 / max_delta)
    }
  }else{
    if(delta <= 0){
      p * dexp(-delta, 1 / mu_delta)
    }else{
      (1-p) / max_delta
    }
  }
}



## Log softmax function, used for choosing arbitrary new ancestors
lsoftmax <- function(v, tau){
  biggest <- max(v/tau)
  # Denominator
  denom <- biggest + log(sum(exp(v/tau - biggest)))
  return(v/tau - denom)

  #exp(v/tau) / sum(exp(v/tau))
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
  keep <- (from != to)

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

    # Update dropout in order from tips towards root
    mcmc$dropout[[h_new]] <- get_dropout(mcmc, data, h_new)
    mcmc$dropout[[h_old]] <- get_dropout(mcmc, data, h_old)

    mcmc <- update_genetics(mcmc, i, h_new, TRUE) # Update genetics. i is inheriting from h_new.

  }else{
    # Update dropout in order from tips towards root
    mcmc$dropout[[h_old]] <- get_dropout(mcmc, data, h_old)
    mcmc$dropout[[h_new]] <- get_dropout(mcmc, data, h_new)

    # For update genetics downstream, h represents the old ancestor
    mcmc <- update_genetics(mcmc, i, h_old, FALSE)
  }

  return(mcmc)
}

get_dropout <- function(mcmc, data, i){

  # If rooted and i is the root, nothing to do here
  # UPDATE if we decide to allow missing positions in the root
  if(data$rooted & i == 1){
    return(integer(0))
  }

  # Children of i
  js <- which(mcmc$h == i)

  # If no children, only dropouts are missing positions if i observed, or all positions if not
  if(length(js) == 0){
    if(i <= data$n_obs){
      return(data$snvs[[i]]$missing)
    }else{
      return(1:data$n_bases)
    }
  }

  # Dropout at i is dropout at kids of i
  out <- Reduce(intersect, mcmc$dropout[js])

  # If i unobserved, nothing left to do here
  if(i > data$n_obs | (!data$rooted & i == 1)){
    return(out)
  }

  # Only dropout positions are missing ones
  out <- intersect(out, data$snvs[[i]]$missing)
  return(out)

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

# Convert ancestry list to kids per node
anc_to_js <- function(h, is){
  out <- rep(list(integer(0)), length(h))
  for (j in 2:length(h)) {
    out[[h[j]]] <- c(out[[h[j]]], j)
  }
  return(out[is])
}

# Accept / reject
accept_or_reject <- function(prop, mcmc, data, update_e, update_g, update_m, hastings = 0, check_parsimony = integer(0), noisy = FALSE){
  if(is.infinite(hastings)){
    stop("hastings error")
  }

  # for (i in 1:prop$n) {
  #   if(
  #     any(prop$dropout[[i]] != get_dropout(prop, data, i)) | (!all(prop$dropout[[prop$h[i]]] %in% prop$dropout[[i]]))
  #   ){
  #     stop("Dropout updated incorrectly")
  #   }
  # }





  if(length(check_parsimony) > 0){
    for (k in check_parsimony) {
      if(
        !genotype(prop, data, k, check_parsimony = T)
      ){
        return(mcmc)
      }
    }
  }

  ## Check that no SNVs are listed in "dropout"
  # for (i in 1:prop$n) {
  #   if(any(
  #     prop$subs$pos[[i]] %in% prop$dropout[[i]]
  #   )){
  #     print(i)
  #     stop("No mutations should be listed at positions that drop out")
  #   }
  # }

  if(length(update_e) > 0){

    js <- anc_to_js(prop$h, update_e)
    prop$e_lik[update_e] <- mapply(e_lik_personal, update_e, js, MoreArgs = list(mcmc = prop, data = data))

  }

  if(length(update_g) > 0){

    js <- anc_to_js(prop$h, update_g)
    prop$g_lik[update_g] <- mapply(g_lik, update_g, js, MoreArgs = list(mcmc = prop, data = data))

  }

  if(length(update_m) > 0){

    js <- anc_to_js(prop$h, update_m)
    prop$m_lik[update_m] <- mapply(m_lik, update_m, js, MoreArgs = list(mcmc = prop, data = data))

  }

  prop$prior <- prior(prop)

  # Accept / reject
  if(
    log(runif(1)) <
    sum(prop$e_lik) +
    sum(prop$g_lik) +
    sum(prop$m_lik) +
    prop$prior -
    sum(mcmc$e_lik) -
    sum(mcmc$g_lik) -
    sum(mcmc$m_lik) -
    mcmc$prior +
    hastings
  ){
    if(noisy){
      print("Move accepted")
    }
    return(prop)
  }else{
    return(mcmc)
  }
}

# BFS traversal of tree
bfs <- function(h){
  n <- length(h)
  kids <- rep(list(integer(0)), n)

  for (i in 2:n) {
    kids[[h[i]]] <- c(kids[[h[i]]], i)
  }

  out <- rep(1, n)
  frontier <- kids[[1]]
  counter <- 2
  while (length(frontier) > 0) {
    out[counter:(counter + length(frontier) - 1)] <- frontier
    counter <- counter + length(frontier)
    frontier <- unlist(kids[frontier])
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

  # Traverse the tree in reverse-BFS order
  ord <- rev(bfs(h))

  # Minimum number of nodes per subtree
  lambda <- mcmc$n / data$n_subtrees

  # Root outputs
  roots <- c()

  # Childen of a node, not including self
  kids <- rep(list(integer(0)), mcmc$n)

  # How many nodes have been chopped off so far?
  n_cut <- 0

  for (v in ord) {
    if(v == 1){

      #trees <- c(trees, kids[1])
      roots <- c(roots, 1)
    }else{

      if(
        length(kids[[v]]) >= lambda & # Enough kids to make a subtree
        !(v %in% old_roots) & # Wasn't a previous root
        mcmc$n - n_cut - length(kids[[v]]) >= lambda # Root has enough kids
      ){
        # Make sure global root has enough kids
        #trees <- c(trees, kids[v])
        roots <- c(roots, v)
        n_cut <- n_cut + length(kids[[v]])

        # Back up only the root into its ancestor, since it's ALSO part of its parent subtree
        kids[[h[v]]] <- c(kids[[h[v]]], v)


      }else{
        # Otherwise, back up kids of v onto kids of h[v]
        kids[[h[v]]] <- c(kids[[h[v]]], v, kids[[v]])
      }
    }
  }

  trees <- lapply(kids[roots], sort)

  return(list(roots, trees))
}

## For debugging: infer mu just based on isnv data
infer_mu_from_isnvs <- function(mcmc, data){
  mus <- mcmc$mu
  for (r in 2:1000) {
    prop <- mcmc
    prop$mu <- mcmc$mu + rnorm(1, 0, 2e-6)
    prop$g_lik <- sapply(1:prop$n, g_lik, mcmc = prop, data = data)
    if(log(runif(1)) < sum(prop$g_lik) - sum(mcmc$g_lik)){
      mcmc <- prop
    }
    mus <- c(mus, mcmc$mu)
  }
  return(mus)
}

# Extinction probability objective to minimize
obj_ext <- function(x, rho, psi){
  (log(x) - rho * (log(psi) - log(1 - (1-psi) * x)))^2
}

p_ext <- function(rho, psi){
  optim(1e-5, obj_ext, method = "Brent", rho=rho, psi=psi, lower = 0, upper = 1)$par
}
