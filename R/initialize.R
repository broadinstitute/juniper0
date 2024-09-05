#' Initialize Outbreak Reconstruction
#'
#' This function initializes the outbreak reconstruction algo.
#'
#' @param n_subtrees How many subtrees to parallelize over. Defaults to NULL, in which case it's the number of cores available to the machine, so long as each subtree has at least 25 nodes.
#' @param n_global Number of global iterations. Defaults to 100.
#' @param n_local Number of local iterations per global iteration. Defaults to 100.
#' @param sample_every Number of local iterations per one sample. Must divide n_local.
#' @param init_mst Should we initialize to a minimum spanning tree? (Set to FALSE for large datasets due to long runtime.)
#' @param init_ancestry If TRUE, the initial ancestry is specified in a one-column .csv file called ancestry.csv, where the entry in the ith row is the initial ancestor of host i.
#' @param rooted If TRUE, the sequence ref.fasta is treated as the root of the transmission network. If FALSE, no root is prespecified. In the latter case, it is recommended to specify a fixed mutation rate value via fixed_mu to ensure convergence.
#' @param record Parameters to be recorded in the MCMC output.
#' @param filters Filters for within-host variation data. List consisting of three named values: af, dp, and sb, meaning minor allele frequency threshhold, read depth threshhold, and strand bias threshhold, respectively. Defaults to NULL, in which case these filters are set to 0.03, 100, and 10, respectively.
#' @param indir Name of the directory in which we have aligned.fasta, ref.fasta, date.csv, vcf folder, and other optional inputs. Defaults to "input_data".
#' @param a_g Shape parameter, generation interval. Defaults to 5.
#' @param lambda_g Rate parameter, generation interval. Defaults to 1.
#' @param a_s Shape parameter, sojourn interval. Defaults to 5.
#' @param lambda_s Rate parameter, sojourn interval. Defaults to 1.
#' @param psi Second parameter of the negative-binomially distributed offspring distribution. Defaults to 0.5.
#' @param init_mu Initial value of the mutation rate, in substitutions/site/day. May be fixed using fixed_mu = TRUE, or inferred otherwise. Defaults to 1e-6.
#' @param fixed_mu If FALSE (the default), the mutation rate is estimated. If TRUE, the mutation rate is fixed at this value for the duration of the algorithm.
#' @param N_eff The growth rate of the within-host effective population size. Specifically, at time t after inoculation, a host has exp(N_eff * t) virions in their body. Defaults to log(100).
#' @return The initial configuration of the Markov Chain.
#' @export
initialize <- function(
    n_subtrees = NULL,
    n_global = 100, # Number of global moves
    n_local = 100, # Number of local moves per global move
    sample_every = 100, # Per how many local moves do we draw one sample? Should be a divisor of n_local
    init_mst = TRUE, # Should we initialize to a minimum spanning tree?
    init_ancestry = FALSE, # Specify the starting ancestry
    rooted = TRUE, # Is the root of the transmission network fixed at the ref sequence?
    N = NA, # Population size
    record = c("n", "h", "seq", "N_eff", "mu", "pi", "R"), # Which aspects of mcmc do we want to record
    filters = NULL,
    # If FALSE, all names must match exactly, with names of VCFs being the same as the names on the FASTA, plus the .vcf suffix
    indir = "input_data", # Name of the directory in which we have aligned.fasta, ref.fasta, date.csv, vcf folder, and other optional inputs
    a_g = 5, # Shape parameter, generation interval
    lambda_g = 1, # Rate parameter, generation interval
    a_s = 5, # Shape parameter, sojourn interval
    lambda_s = 1, # Rate parameter, sojourn interval
    psi = 0.5, # Second parameter in negative binomial offspring distribution. E[NBin(rho, psi)] = R => rho*(1-psi)/psi = R => rho = R*psi / (1-psi)
    init_mu = 2e-5,
    fixed_mu = F, # Should mutation rate be fixed? Defaults to FALSE.
    N_eff = log(100) # Effective population size, within host
){

  ## Filters
  if(is.null(filters)){
    filters <- list(
      af = 0.03,
      dp = 100,
      sb = 10
    )
    if(file.exists(paste0("./", indir, "/problematic.csv"))){
      filters$common <- read.csv(paste0("./", indir, "/problematic.csv"))[,1]
    }else{
      filters$common <- integer(0)
    }
  }


  ### Data Processing

  ## FASTAs and dates

  # Load the FASTA of sequences
  fasta <- ape::read.FASTA(paste0("./", indir, "/aligned.fasta"))

  # Times of collection for the non-reference sequence
  s_nonref <- as.Date(gsub(".*\\|", "", names(fasta)))

  # If rooted, reference sequence is provided. Otherwise, initialized to earliest case in dataset.
  # In the latter case, the genome will be updated regardless, so this initialization doesn't matter

  if(rooted){
    # Load the reference sequence
    ref_genome <- ape::read.FASTA(paste0("./", indir, "/ref.fasta"))
    s_ref <- as.Date(gsub(".*\\|", "", names(ref_genome)))
  }else{
    earliest <- which.min(s_nonref)
    ref_genome <- fasta[earliest]
    names(ref_genome) <- "ref_genome"

    # Going to initialize time of collection to 10 days before earliest sequence collection
    s_ref <- min(s_nonref) - 10
  }

  # Length of genome
  n_bases <- length(ref_genome[[1]])

  # The first genome is itself the ref genome
  fasta <- c(ref_genome, fasta)
  s <- c(as.Date(NA), s_nonref) # Ref genome always treated as unsampled

  # The earliest case to be collected must be offset by 1, since we've appended the reference genome in front
  if(!rooted){
    earliest <- earliest + 1
  }

  # Shift all dates such that the last date of sample collection is always 0
  s_max <- max(s, na.rm = T)
  s <- as.numeric(difftime(s, s_max, units = 'days'))

  # Time of most recent common ancestor, before s_max
  tmrca <- as.numeric(difftime(s_ref, s_max, units = 'days'))





  # Number of samples
  n <- length(fasta)

  if(is.null(n_subtrees)){
    n_subtrees <- max(min(parallel::detectCores(), floor(n / 25)), 1) # Minimum 25 nodes per subtree
  }

  # Prevent MST initialization when dataset is large
  if(init_mst & n > 1000){
    stop("Dataset is too large to use init_mst = TRUE. Specify a strating configuration by using init_ancestry, or initialize to a star-shaped configuration by setting init_mst = FALSE.")
  }

  # Warn when using unrooted tree, unfixed mutation rate
  if(!rooted & !fixed_mu){
    message("Warning: Convergence may fail when rooted = FALSE and fixed_mu = FALSE, especially if samples span a short date range. Consider setting one of these arguments to TRUE.")
  }

  # Names of sequences
  names <- gsub("\\|.*", "", names(fasta))

  # VCF files present
  vcfs <- list.files(paste0("./", indir, "/vcf"))





  ## List of SNVs present per sample
  message("Processing FASTA and VCF files...")
  snvs <- list()
  pb = txtProgressBar(min = 0, max = n, initial = 0)

  vcfs_prefix <- gsub(".vcf", "", vcfs)
  vcf_present <- c()
  for (i in 1:n) {
    # Locate the correct vcf file
    who <- which(vcfs_prefix == names[i])
    if(length(who) == 1){
      vcf <- read.delim(
        paste0("./", indir, "/vcf/", vcfs[who]),
        colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character")
      )
      colnames(vcf) <- paste0("V", 1:ncol(vcf))
      snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters, vcf = vcf)
      vcf_present[i] <- TRUE
    }else{
      snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters)
      vcf_present[i] <- FALSE
    }
    setTxtProgressBar(pb,i)
  }

  close(pb)


  # Compile vectors of all positions with SNVs (may have duplicates) and all SNVs (unique)
  message("Compiling list of mutations...")
  all_pos <- c()
  all_snv <- c()
  pb = txtProgressBar(min = 0, max = n, initial = 0)
  for (i in 1:n) {
    calls <- snvs[[i]]$snv$call
    new <- which(!(calls %in% all_snv))
    all_snv <- c(all_snv, calls[new])
    all_pos <- c(all_pos, snvs[[i]]$snv$pos[new])
    if(!is.null(snvs[[i]]$isnv)){
      calls <- snvs[[i]]$isnv$call
      new <- which(!(calls %in% all_snv))
      all_snv <- c(all_snv, calls[new])
      all_pos <- c(all_pos, snvs[[i]]$isnv$pos[new])
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)

  # Come back to missing data

  ### Initialize MCMC and data

  data <- list()
  data$s <- s
  data$N <- N #population size
  data$n_obs <- n # number of observed hosts, plus 1 (index case)
  data$n_bases <- n_bases
  data$snvs <- snvs
  data$tau = 0.2
  # Number of subtrees to chop into is n_cores, as long as each subtree has at least 100 people
  data$n_subtrees <- n_subtrees
  data$n_global <- n_global
  data$n_local <- n_local
  data$sample_every <- sample_every
  data$record <- record
  data$filters <- filters
  data$vcf_present <- vcf_present
  data$rooted <- rooted
  data$init_mu <- init_mu
  data$fixed_mu <- fixed_mu
  data$names <- names

  # Later: could move elements of "data" that aren't used in MCMC to a new "config" list

  mcmc <- list()
  mcmc$n <- n # number of tracked hosts
  mcmc$h <- rep(1, n)
  mcmc$h[1] <- NA

  # SNVs: from, pos, to
  mcmc$subs <- list()
  mcmc$subs$from <- list()
  mcmc$subs$pos <- list()
  mcmc$subs$to <- list()

  # For each iSNV, TRUE if allele a1 in the bottleneck, FALSE otherwise
  mcmc$bot <- list()

  for (i in 1:n) {
    # Positions at which the consensus genome changes
    pos <- data$snvs[[i]]$snv$pos
    cons_change <- logical(0)

    if(length(pos) > 0){
      # Positions in which the initial bottleneck genome should indeed report a consensus change
      for (p in 1:length(pos)) {
        cons_change[p] <- TRUE
        if(pos[p] %in% data$snvs[[i]]$isnv$pos){
          # Index of which isnv
          ind <- which(data$snvs[[i]]$isnv$pos == pos[p])
          ref_allele <- raw_to_base(ref_genome[[1]][pos[p]])
          # If either allele in the iSNV matches the ref genome (root), no consensus change in i
          if(data$snvs[[i]]$isnv$a1[ind] == ref_allele | data$snvs[[i]]$isnv$a2[ind] == ref_allele){
            cons_change[p] <- FALSE
          }
        }
      }
    }

    mcmc$subs$from[[i]] <- data$snvs[[i]]$snv$from[cons_change]
    mcmc$subs$pos[[i]] <- data$snvs[[i]]$snv$pos[cons_change]
    mcmc$subs$to[[i]] <- data$snvs[[i]]$snv$to[cons_change]

    isnv_pos <- data$snvs[[i]]$isnv$pos
    mcmc$bot[[i]] <- logical(0)

    if(length(isnv_pos) > 0){
      for (p in 1:length(isnv_pos)) {
        ref_allele <- raw_to_base(ref_genome[[1]][isnv_pos[p]])

        if(data$snvs[[i]]$isnv$a1[p] == ref_allele){
          # Does a1 match the ref allele? If so, a1 is in the bottleneck
          mcmc$bot[[i]][p] <- T
        }else if(data$snvs[[i]]$isnv$a2[p] == ref_allele){
          # Does a2 match the ref allele? If so, it's in the bottleneck
          mcmc$bot[[i]][p] <- F
        }else{
          # If neither, then whichever is the consensus allele is in the bottleneck
          mcmc$bot[[i]][p] <- data$snvs[[i]]$isnv$af1[p] > 0.5
        }
      }
    }
  }

  mcmc$a_g <- a_g # shape parameter of the generation interval
  mcmc$lambda_g <- lambda_g # rate parameter of the generation interval. FOR NOW: fixing at 1.
  mcmc$a_s <- a_s # shape parameter of the sojourn interval
  mcmc$lambda_s <- lambda_s # rate parameter of the sojourn interval. FOR NOW: fixing at 1.
  mcmc$mu <- init_mu
  mcmc$psi <- psi # second parameter, NBin offspring distribution (computed in terms of R0)

  # Optimize R and pi
  #vals <- opt_R_pi(data$s - tmrca, mcmc$a_g, mcmc$lambda_g, mcmc$a_s, mcmc$lambda_s)
  #print(vals)
  vals <- c(2, 0.5)

  mcmc$R <- vals[1] # Reproductive number
  mcmc$pi <- vals[2] # Probability of sampling
  mcmc$N_eff <- N_eff

  # Sequence of times at which the hosts along the edge leading into i were sampled
  mcmc$seq <- list()
  mcmc$seq[[1]] <- tmrca
  mcmc$seq[2:n] <- lapply(2:n, get_ts, mcmc = mcmc, data = data)

  # Times at which mutations occur
  mcmc$tmu <- list()

  for (i in 1:n) {
    n_subs <- length(mcmc$subs$from[[i]])
    if(n_subs > 0){
      mcmc$tmu[[i]] <- runif(n_subs, tmrca, mcmc$seq[[i]][1])
    }else{
      mcmc$tmu[[i]] <- numeric(0)
    }
  }

  for (i in 1:n) {
    data$snvs[[i]]$snv <- NULL
  }

  print(length(unlist(mcmc$subs$from)))

  if(!data$rooted){
    mcmc <- genotype(mcmc, data, 1)[[1]]
  }

  print(length(unlist(mcmc$subs$from)))

  # Also track the epidemiological and genomic likelihoods, and prior
  # The genomic likelihood we will store on a per-person basis, for efficiency purposes
  mcmc$e_lik <- e_lik(mcmc, data)
  mcmc$g_lik <- sapply(1:n, g_lik, mcmc = mcmc, data = data)
  mcmc$prior <- prior(mcmc)

  return(list(mcmc, data))
}
