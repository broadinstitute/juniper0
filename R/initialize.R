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
#' @param check_names If TRUE, checks whether all of the names in the FASTA match the names of the VCFs and dates.
#' @param virus being studied (affects various presets). So far, options are "SARS-CoV-2", "H5N1", and "HepA". More to come soon!
#' @param a_g Shape parameter, generation interval. Defaults to 5.
#' @param lambda_g Rate parameter, generation interval. Defaults to 1.
#' @param a_s Shape parameter, sojourn interval. Defaults to 5.
#' @param lambda_s Rate parameter, sojourn interval. Defaults to 1.
#' @param rho Overdispersion parameter. Defaults to Inf, indicating the offspring distribution is Poisson.
#' @param R Reproductive number (average over entire outbreak). Defaults to 1. Exactly one of R and growth may be specified; the other must be set to NULL.
#' @param growth Exponential growth rate of cases. Defaults to NULL. Exactly one of R and growth may be specified; the other must be set to NULL.
#' @param fixed_mu If NA (the default), the mutation rate is estimated. If a numeric value (in substitutions/site/day), the mutation rate is fixed at this value for the duration of the algorithm.
#' @return The initial configuration of the Markov Chain.
#' @export
initialize <- function(
    n_subtrees = NULL,
    n_global = 100, # Number of global moves
    n_local = 100, # Number of local moves per global move
    sample_every = 100, # Per how many local moves do we draw one sample? Should be a divisor of n_local
    init_mst = FALSE, # Should we initialize to a minimum spanning tree? Bad idea if dataset is large.
    init_ancestry = FALSE, # Specify the starting ancestry
    rooted = TRUE, # Is the root of the transmission network fixed at the ref sequence?
    N = NA, # Population size
    record = c("n", "h", "w", "t", "b", "mu", "p"), # Which aspects of mcmc do we want to record
    filters = NULL,
    check_names = TRUE, # Should we check to make sure all of the names in the FASTA match the names of the VCFs and dates?
    # If FALSE, all names must match exactly, with names of VCFs being the same as the names on the FASTA, plus the .vcf suffix
    virus = "SARS-CoV-2", # Pathogen being studied
    a_g = 5, # Shape parameter, generation interval
    lambda_g = 1, # Rate parameter, generation interval
    a_s = 5, # Shape parameter, sojourn interval
    lambda_s = 1, # Rate parameter, sojourn interval
    rho = Inf, # Overdispersion parameter (Inf indicates Poisson distribution)
    R = 1, # Reproductive number (average over entire outbreak)
    growth = NULL, # Exponential growth rate of cases, i.e. # of cases at time t is exp(growth * t)
    fixed_mu = NA # Should mutation rate be fixed? If a number is given, fixed at that. NA means we estimate mu
){

  # n_subtrees = NULL
  # n_global = 100 # Number of global moves
  # n_local = 100 # Number of local moves per global move
  # sample_every = 100 # Per how many local moves do we draw one sample? Should be a divisor of n_local
  # init_mst = FALSE # Should we initialize to a minimum spanning tree? Bad idea if dataset is large.
  # init_ancestry = FALSE # Specify the starting ancestry
  # rooted = TRUE # Is the root of the transmission network fixed at the ref sequence?
  # N = NA # Population size
  # record = c("n", "h", "w", "t", "b", "mu", "p") # Which aspects of mcmc do we want to record
  # filters = NULL
  # check_names = TRUE # Should we check to make sure all of the names in the FASTA match the names of the VCFs and dates?
  # # If FALSE, all names must match exactly, with names of VCFs being the same as the names on the FASTA, plus the .vcf suffix
  # virus = "SARS-CoV-2" # Pathogen being studied
  # a_g = 5 # Shape parameter, generation interval
  # lambda_g = 1 # Rate parameter, generation interval
  # a_s = 5 # Shape parameter, sojourn interval
  # lambda_s = 1 # Rate parameter, sojourn interval
  # rho = Inf # Overdispersion parameter (Inf indicates Poisson distribution)
  # R = 1 # Reproductive number (average over entire outbreak)
  # growth = NULL # Exponential growth rate of cases, i.e. # of cases at time t is exp(growth * t)
  # fixed_mu = NA # Should mutation rate be fixed? If a number is given, fixed at that. NA means we estimate mu

  # if(!rooted & is.na(N)){
  #   stop("If the root of the tree is not fixed, N must be specified.")
  # }

  if((!is.null(R) & !is.null(growth)) | (is.null(R) & is.null(growth))){
    stop("Exactly one of the inputs R and growth may be specified. The other must be set to NULL.")
  }else if(!is.null(R)){
    growth <- log(R) / (a_g / lambda_g)
  }else if(!is.null(growth)){
    R <- exp(growth * a_g / lambda_g)
  }

  ## Filters
  if(is.null(filters)){
    filters <- list(
      af = 0.03,
      dp = 100,
      sb = 10
    )
    if(file.exists("./input_data/problematic.csv")){
      filters$common <- read.csv("./input_data/problematic.csv")[,1]
    }else{
      filters$common <- integer(0)
    }
  }


  ## Data Processing

  # Load the reference sequence
  ref_genome <- ape::read.FASTA("./input_data/ref.fasta")

  # Length of genome
  n_bases <- length(ref_genome[[1]])

  # Load the FASTA of sequences
  fasta <- ape::read.FASTA("./input_data/aligned.fasta")

  # The first genome is itself the ref genome
  fasta <- c(ref_genome, fasta)

  # Number of samples
  n <- length(fasta)

  if(is.null(n_subtrees)){
    n_subtrees <- max(min(parallel::detectCores(), floor(n / 25)), 1) # Minimum 25 nodes per subtree
  }

  # Prevent MST initialization when dataset is large
  if(init_mst & n > 1000){
    stop("Dataset is too large to use init_mst = TRUE. Specify a strating configuration by using init_ancestry, or initialize to a star-shaped configuration by setting init_mst = FALSE.")
  }

  # Names of sequences
  names <- names(fasta)

  # VCF files present
  vcfs <- list.files("./input_data/vcf/")

  # Date
  date <- read.csv("input_data/date.csv")

  if(check_names){
    s <- c()
    for (i in 1:n) {
      # Check if we have a date for the sample
      included <- grepl(names[i], date[,1])
      if(sum(included) >= 2){
        stop(paste("Multiple sample collection dates found for sequence", names[i]))
      }else if(sum(included) == 1){
        s[i] <- date[,2][included]
      }else{
        s[i] <- NA
      }
    }
    s[1] <- 0
  }else{
    s <- date$date[match(names, date$cases)]
    s[1] <- 0
  }


  ## List of SNVs present per sample
  message("Processing FASTA and VCF files...")
  snvs <- list()
  pb = txtProgressBar(min = 0, max = n, initial = 0)
  if(check_names){
    vcf_present <- c()
    for (i in 1:n) {
      # Check if we have a VCF file
      included <- grepl(names[i], vcfs)
      if(sum(included) >= 2){
        stop(paste("Multiple VCF files found for sequence", names[i]))
      }else if(sum(included) == 1){
        vcf <- read.table(paste0("./input_data/vcf/", vcfs[included]), row.names = NULL)
        snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters, vcf = vcf)
        vcf_present[i] <- TRUE
      }else{
        snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters)
        vcf_present[i] <- FALSE
      }
      setTxtProgressBar(pb,i)
    }
  }else{
    vcfs_prefix <- gsub(".vcf", "", vcfs)
    vcf_present <- c()
    for (i in 1:n) {
      # Locate the correct vcf file
      who <- which(vcfs_prefix == names[i])
      if(length(who) == 1){
        vcf <- read.table(paste0("./input_data/vcf/", vcfs[who]), row.names = NULL)
        snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters, vcf = vcf)
        vcf_present[i] <- TRUE
      }else{
        snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters)
        vcf_present[i] <- FALSE
      }
      setTxtProgressBar(pb,i)
    }
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

  # Remove irrelevant missing site info; add SNVs with missing data
  message("Identifying missing data in samples...")
  pb = txtProgressBar(min = 0, max = n, initial = 0)
  for (i in 1:n) {
    # Which positions in all_pos are detected in the missing sites in the sample?
    present <- which(all_pos %in% snvs[[i]]$missing$pos)
    # Change missing$pos to be a vector of these sites (may have duplicates)
    snvs[[i]]$missing$pos <- all_pos[present]
    # Add a new vector of the SNVs for which there's no information
    snvs[[i]]$missing$call <- all_snv[present]
    setTxtProgressBar(pb,i)
  }
  close(pb)

  ### Initialize MCMC and data

  ## For MCMC initialization: minimum spanning tree
  if(init_mst){
    snv_dist <- ape::dist.dna(fasta, "N")
    tree <- ape::mst(snv_dist)
    init_h <- adj_to_anc(tree, 1)
  }

  if(init_ancestry){
    init_h <- read.table("./input_data/ancestry.csv")[,1]
  }

  data <- list()
  data$s <- s
  data$N <- N #population size
  data$n_obs <- n # number of observed hosts, plus 1 (index case)
  data$n_bases <- n_bases
  data$snvs <- snvs
  data$eps <- 0.005 # Explore/exploit tradeoff for genotypes of new nodes
  data$p_move <- 0.6
  data$tau = 0.2
  #data$n_cores <-
  # Number of subtrees to chop into is n_cores, as long as each subtree has at least 100 people
  data$n_subtrees <- n_subtrees
  data$n_global <- n_global
  data$n_local <- n_local
  data$sample_every <- sample_every
  data$record <- record
  #data$n_subtrees <- 3
  data$filters <- filters
  data$virus <- virus
  data$growth <- growth
  data$R <- R
  data$vcf_present <- vcf_present
  data$rooted <- rooted
  if(is.na(fixed_mu)){
    data$fixed_mu <- F
  }else{
    data$fixed_mu <- T
  }

  # Old feature from previous version
  data$pooled_coalescent = T
  data$disjoint_coalescent = F

  mcmc <- list()
  mcmc$n <- n # number of tracked hosts
  mcmc$h <- rep(1, n) # ancestors; initialized to index case
  # mcmc$w <- rep(0, n) # edge weights; initialized to 0
  # mcmc$w[1] <- 0 # For convenience
  mcmc$h[1] <- NA
  mcmc$t <- s - 5 # time of contracting
  mcmc$t[2:n] <- pmax(mcmc$t[2:n], -3.5)
  mcmc$m01 <- list() # fixed mutations added in each transmission link
  mcmc$m10 <- list() # fixed mutations deleted in each transmission link
  mcmc$m0y <- list() # 0% -> y%, 0 < y < 100
  mcmc$m1y <- list() # 100% -> y%, 0 < y < 100
  mcmc$mx0 <- list() # x% -> 0%, 0 < x < 100
  mcmc$mx1 <- list() # x% -> 100%, 0 < x < 100
  mcmc$mxy <- list() # x% -> y%, 0 < x < 100, 0 < y < 100
  for (i in 1:n) {
    mcmc$m01[[i]] <- snvs[[i]]$snv$call
    mcmc$m10[[i]] <- character(0)
    if(is.null(snvs[[i]]$isnv$call)){
      mcmc$m0y[[i]] <- character(0)
    }else{
      mcmc$m0y[[i]] <- snvs[[i]]$isnv$call
    }
    mcmc$m1y[[i]] <- character(0)
    mcmc$mx0[[i]] <- character(0)
    mcmc$mx1[[i]] <- character(0)
    mcmc$mxy[[i]] <- character(0)
  }

  if(init_mst | init_ancestry){
    gens <- generations(init_h, 1)
    max_t <- min(s[2:n] - 5)

    message("Initializing transmission network...")
    progress <- 0
    pb = txtProgressBar(min = 0, max = n, initial = 0)
    for (g in 2:length(gens)) {
      for (i in gens[[g]]) {

        if(g >= 3){
          anc <- ancestry(init_h, i)
          for (j in 2:(length(anc) - 1)) {
            mcmc <- update_genetics_upstream(mcmc, mcmc, i, anc[j])
            mcmc$m01[[j]] <- setdiff(mcmc$m01[[j]], snvs[[j]]$missing$call) # Remove calls for missing positions
            mcmc$m10[[j]] <- setdiff(mcmc$m10[[j]], snvs[[j]]$missing$call)
          }
        }



        progress <- progress + 1
        setTxtProgressBar(pb,progress)
      }

    }
    close(pb)
    mcmc$h <- init_h

    # Initialize time of infection
    ord <- rev(bfs(1, mcmc$h))
    mcmc$t <- rep(NA, n)
    mcmc$t[1] <- -5
    for (i in ord) {
      mcmc$t[i] <- min(c(data$s[i] - 5, mcmc$t[which(mcmc$h == i)] - 5))
    }
  }

  # If not rooted, initialize root to earliest case
  if(!rooted){
    # Earliest case, besides root
    earliest <- which.min(s[2:n]) + 1

    # Set infection time one generation earlier, to ensure infectious at time of transmission
    mcmc$t[earliest] <- mcmc$t[earliest] - (a_g/lambda_g)

    # Time of infection of ref is arbitrarily early
    mcmc$t[1] <- -Inf
    data$s[1] <- -Inf

    # Update ancestry
    mcmc$h[earliest] <- 1
    mcmc$h[setdiff(2:n, earliest)] <- earliest
    mcmc$w[earliest] <- Inf

    # Update genetics
    for (i in setdiff(2:n, earliest)) {
      mcmc <- update_genetics_upstream(mcmc, mcmc, i, earliest)
    }

  }

  mcmc$b <- 0.05 # Probability bottleneck has size >1
  mcmc$a_g <- a_g # shape parameter of the generation interval
  mcmc$lambda_g <- lambda_g # rate parameter of the generation interval. FOR NOW: fixing at 1.
  mcmc$a_s <- a_s # shape parameter of the sojourn interval
  mcmc$lambda_s <- lambda_s # rate parameter of the sojourn interval. FOR NOW: fixing at 1.
  if(is.na(fixed_mu)){
    if(virus == "SARS-CoV-2"){
      mcmc$mu <- 2.7e-6 # mutation rate, sites/day
      mcmc$p <- 1e-6 # mutation rate, sites/cycle
    }else if(virus == "H5N1"){
      mcmc$mu <- 1e-5 # mutation rate, sites/day
      mcmc$p <- 5e-6 # mutation rate, sites/cycle
    }else if(virus == "HepA"){
      mcmc$mu <- 2e-6
      mcmc$p <- 1e-6
    }
  }else{
    mcmc$mu <- fixed_mu
    mcmc$p <- fixed_mu / 2
  }

  mcmc$v <- 1000 # burst size
  mcmc$lambda <- 1 # expo growth rate of bursts
  mcmc$rho <- rho # first parameter, NBin offspring distribution (overdispersion param)

  # Mean number of cases at time t is R^(t/g) = exp(t * log(R) / g), g mean generation interval
  # growth = log(R) / g => R = exp(g * growth)
  mcmc$psi <- mcmc$rho / (data$R + mcmc$rho) # second parameter, NBin offspring distribution (computed in terms of R0)

  mcmc$w <- round(((mcmc$t - mcmc$t[mcmc$h]) / (mcmc$a_g / mcmc$lambda_g)))
  mcmc$w[1] <- 0

  # Functions of MCMC params
  mcmc$d <- sapply(1:n, function(x){sum(mcmc$h[2:n] == x)}) # Node degrees

  # Also track the epidemiological and genomic likelihoods, and prior
  # The genomic likelihood we will store on a per-person basis, for efficiency purposes
  mcmc$e_lik <- e_lik(mcmc, data)
  mcmc$g_lik <- c(NA, sapply(2:n, g_lik, mcmc = mcmc, data = data))
  mcmc$prior <- prior(mcmc)

  return(list(mcmc, data))
}
