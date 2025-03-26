#' Initialize Outbreak Reconstruction
#'
#' This function initializes the outbreak reconstruction algo.
#'
#' @param n_subtrees How many subtrees to parallelize over. Defaults to NULL, in which case it's the number of cores available to the machine, so long as each subtree has at least 25 nodes.
#' @param n_global Number of global iterations. Defaults to 100.
#' @param n_local Number of local iterations per global iteration. Defaults to 100.
#' @param sample_every Number of local iterations per one sample. Must divide n_local.
#' @param root If not NA, the sequence name given is treated as the root of the transmission network. Defaults to NA
#' @param record Parameters to be recorded in the MCMC output.
#' @param filters Filters for within-host variation data. List consisting of three named values: af, dp, and sb, meaning minor allele frequency threshhold, read depth threshhold, and strand bias threshhold, respectively. Defaults to NULL, in which case these filters are set to 0.03, 100, and 10, respectively.
#' @param indir Name of the directory in which we have fasta(s), metadata table, and (optionally) vcf files. Defaults to "input_data."
#' @param a_g Shape parameter, generation interval. Defaults to 5.
#' @param lambda_g Rate parameter, generation interval. Defaults to 1.
#' @param a_s Shape parameter, sojourn interval. Defaults to 5.
#' @param lambda_s Rate parameter, sojourn interval. Defaults to 1.
#' @param psi Second parameter of the negative-binomially distributed offspring distribution. Defaults to 0.5.
#' @param init_mu Initial value of the mutation rate, in substitutions/site/day. May be fixed using fixed_mu = TRUE, or inferred otherwise. Defaults to 1e-6.
#' @param fixed_mu If FALSE (the default), the mutation rate is estimated. If TRUE, the mutation rate is fixed at its initial value for the duration of the algorithm.
#' @param upper_mu Upper bound on the mutation rate, in substitutions/site/day. Defaults to 1.
#' @param init_N_eff The growth rate of the within-host effective population size. Specifically, at time t after inoculation, a host has exp(N_eff * t) virions in their body. Defaults to log(100).
#' @param fixed_N_eff If FALSE (the default), the within-host population growth rate is estimated. If TRUE, the within-host population growth rate is fixed at its initial value for the duration of the algorithm.
#' @param upper_N_eff Upper bound on the growth rate of the within-host effective population size. Defaults to 100.
#' @param init_R Initial value of the reproductive number. Defaults to 1.
#' @param fixed_R If FALSE (the default), the reproductive number is estimated. If TRUE, the reproductive number is fixed at its initial value for the duration of the algorithm.
#' @param init_pi Initial value of the sampling rate. Defaults to 0.5.
#' @param fixed_pi If FALSE (the default), the sampling rate is estimated. If TRUE, the sampling rate is fixed at its initial value for the duration of the algorithm.
#' @param ongoing TRUE if the outbreak is ongoing, or FALSE if has it concluded (i.e. no more samples will ever be collected from this outbreak). Defaults to TRUE.
#' @param safety Either NA or a non-negative value. If NA, safety mode (checking that the likelihood is correct after each global iteration) is disabled. If numeric, the maximum tolerance for a difference in the re-computed versus stored likelihood before throwing an error. Defaults to NA.
#' @param split_bottlenecks If TRUE, likelihood accounts for split bottlenecks (new beta feature). Defaults to FALSE.
#' @return The initial configuration of the Markov Chain.
#' @export
initialize <- function(
    n_subtrees = NULL,
    n_global = 100, # Number of global moves
    n_local = 100, # Number of local moves per global move
    sample_every = 100, # Per how many local moves do we draw one sample? Should be a divisor of n_local
    root = NA, # Is the root of the transmission network fixed at the specified sequence?
    N = NA, # Population size
    record = c("n", "h", "seq", "N_eff", "mu", "pi", "R"), # Which aspects of mcmc do we want to record
    filters = NULL,
    indir = "input_data", # Name of the directory in which we have fasta(s), metadata table, and (optionally) vcf files.
    a_g = 5, # Shape parameter, generation interval
    lambda_g = 1, # Rate parameter, generation interval
    a_s = 5, # Shape parameter, sojourn interval
    lambda_s = 1, # Rate parameter, sojourn interval
    psi = 0.5, # Second parameter in negative binomial offspring distribution. E[NBin(rho, psi)] = R => rho*(1-psi)/psi = R => rho = R*psi / (1-psi)
    init_mu = 2e-5,
    fixed_mu = FALSE, # Should mutation rate be fixed? Defaults to FALSE.
    upper_mu = 1,
    init_N_eff = log(100), # Effective population size, within host
    fixed_N_eff = FALSE,
    upper_N_eff = 100,
    init_R = 1,
    fixed_R = FALSE,
    init_pi = 0.5,
    fixed_pi = FALSE,
    ongoing = TRUE,
    safety = NA,
    split_bottlenecks = FALSE
){

  # Does the cluster have a prespecified root?
  rooted <- !is.na(root)

  if(!rooted & init_pi == 1 & !ongoing){
    stop("If treating all cases as sampled, root must be provided.")
  }

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

  # Read in metadata file
  meta_names <- list.files(paste0("./", indir), pattern = "metadata")

  if(length(meta_names) == 0){
    stop("No metadata file found.")
  }

  if(length(meta_names) > 1){
    stop(paste0("Multiple possible metadata files found: ", paste(meta_names, collapse = ", "), "."))
  }

  # Get the suffix of metadata file
  meta_extension <- sub(".*\\.", "", meta_names)

  if(meta_extension == "tsv"){
    meta <- read.csv(paste0("./", indir, "/", meta_names), sep = "\t")
  }else if(meta_extension == "csv"){
    meta <- read.csv(paste0("./", indir, "/", meta_names))
  }else{
    stop("Metadata file must be in .tsv or .csv format.")
  }

  # Check metadata file has appropriate columns
  if(!("sample" %in% colnames(meta))){
    stop("Metadata file must contain a column named 'sample' with the names of the sequences.")
  }
  if(!("date" %in% colnames(meta))){
    stop("Metadata file must contain a column named 'date' with the sequence collection dates.")
  }

  if(!all(is_date(meta$date))){
    stop("The 'date' column in the metadata file is not in a standard date format.")
  }

  if(any(duplicated(meta$sample))){
    stop(paste0("The following sequence names are duplicated in the metadata table: ", paste(meta$sample[duplicated(meta$sample)], collapse = ", "), "."))
  }

  if(rooted){
    if(!(root %in% meta$sample)){
      stop("The specified root does not match any of the sequence names in the metadata table.")
    }
  }

  meta$date <- as.Date(meta$date)




  ## If user hasn't provided single aligned.fasta file, and has instead provided individual fastas...
  if(!file.exists(paste0("./", indir, "/aligned.fasta"))){
    fasta_names <- list.files(paste0("./", indir), pattern = ".fasta")
    all_fasta_names <- fasta_names

    # Keep only the fastas that are found in the metadata list
    keep <- integer(0)
    for (i in 1:nrow(meta)) {
      # Line of metadata table that matches the fasta name
      ind <- which(grepl(meta$sample[i], fasta_names))

      if(length(ind) == 0){
        stop(paste0("No FASTA file for sequence name ", meta$sample[i], "."))
      }

      if(length(ind) > 1){
        stop(paste0("Multiple FASTA files for sequence name ", meta$sample[i], ": ", paste(fasta_names[ind], collapse = " "), "."))
      }

      keep <- c(keep, ind)
    }

    # Subset to fasta files kept
    removed <- setdiff(1:length(fasta_names), keep)
    if(length(removed) > 0){
      message(paste0(
        "The following FASTA files could not be matched to the sequence names in the metadata table, and will not be included in the analysis: ",
        paste(fasta_names[removed], collapse = ", ")
      ))
    }


    fasta_names <- fasta_names[keep]

    if(length(fasta_names) == 0){
      stop("No FASTA files that match the metadata table sequence names were found.")
    }

    # Read in each fasta
    fasta <- ape::read.FASTA(paste0("./", indir, "/", fasta_names[1]))
    if(length(fasta_names) > 1){
      for (i in 2:length(fasta_names)) {
        fasta <- c(fasta, ape::read.FASTA(paste0("./", indir, "/", fasta_names[i])))
      }
    }

    # Get length of each fasta
    lengths <- sapply(fasta, length)

    if(!all(lengths == lengths[1])){
      stop("FASTA files are not aligned.")
    }

    names(fasta) <- paste0(meta$sample, "|", meta$date)

    # Bookkeeping: move all fastas to separate folder
    if(!dir.exists(paste0("./", indir, "/fasta"))){
      dir.create(paste0("./", indir, "/fasta"))
    }

    for (i in 1:length(all_fasta_names)) {
      file.rename(paste0("./", indir, "/", fasta_names[i]), paste0("./", indir, "/fasta/", fasta_names[i]))
    }

    # Write aligned fasta
    ape::write.FASTA(fasta, file = paste0("./", indir, "/aligned.fasta"))

  }else{
    fasta <- ape::read.FASTA(paste0("./", indir, "/aligned.fasta"))

    # Keep only the fastas that are found in the metadata list
    keep <- integer(0)
    for (i in 1:nrow(meta)) {
      # Line of metadata table that matches the fasta name
      ind <- which(grepl(meta$sample[i], names(fasta)))

      if(length(ind) == 0){
        stop(paste0("No sequence in aligned.fasta found for sequence name ", meta$sample[i], "."))
      }

      if(length(ind) > 1){
        stop(paste0("Multiple sequences in aligned.fasta found for sequence name ", meta$sample[i], ": ", paste(names(fasta)[ind], collapse = " "), "."))
      }

      keep <- c(keep, ind)
    }

    # Subset to fasta files kept
    removed <- setdiff(1:length(fasta), keep)
    if(length(removed) > 0){
      message(paste0(
        "The following entries of the FASTA could not be matched to the sequence names in the metadata table, and will not be included in the analysis: ",
        paste(names(fasta)[removed], collapse = ", ")
      ))
    }


    fasta <- fasta[keep]

    if(length(fasta) == 0){
      stop("No sequences names in the FASTA match the sequence names in the metadata table.")
    }

    names(fasta) <- paste0(meta$sample, "|", meta$date)

  }

  ## FASTAs and dates

  # If cluster is rooted, remove root from fasta
  if(rooted){
    root_ind <- which(meta$sample == root)
    root_genome <- fasta[root_ind]
    fasta <- fasta[-root_ind]
    s_root <- as.Date(gsub(".*\\|", "", names(root_genome)))
  }

  # Times of collection for the non-root sequence
  s_nonroot <- as.Date(gsub(".*\\|", "", names(fasta)))

  # If rooted, root sequence is provided. Otherwise, initialized to earliest case in dataset.
  # In the latter case, the genome will be updated regardless, so this initialization doesn't matter

  if(!rooted){
    earliest <- which.min(s_nonroot)
    root_genome <- fasta[earliest]
    names(root_genome) <- "ref_genome"

    # Going to initialize time of collection to 1 generation interval and 1 sojourn interval before earliest sequence collection
    s_root <- min(s_nonroot) - (a_g / lambda_g) - (a_s / lambda_s)
  }

  # Fill the root_genome's missing sites with "A". This is simply a placeholder, since genotype() will update all of them
  root_missing <- which(!(root_genome[[1]] %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24))))
  root_genome[[1]][!(root_genome[[1]] %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24)))] <- base_to_raw("A")

  # Length of genome
  n_bases <- length(root_genome[[1]])

  # The first genome is itself the root genome
  fasta <- c(root_genome, fasta)
  if(rooted){
    s <- c(s_root, s_nonroot)
  }else{
    s <- c(as.Date(NA), s_nonroot)
  }


  # The earliest case to be collected must be offset by 1, since we've appended the root genome in front
  if(!rooted){
    earliest <- earliest + 1
  }

  # Shift all dates such that the last date of sample collection is always 0
  s_max <- max(s, na.rm = T)
  s <- as.numeric(difftime(s, s_max, units = 'days'))

  # Time of most recent common ancestor, before s_max
  tmrca <- as.numeric(difftime(s_root, s_max, units = 'days'))





  # Number of samples
  n <- length(fasta)

  if(is.null(n_subtrees)){
    n_subtrees <- max(min(parallel::detectCores(), floor(n / 25)), 1) # Minimum 25 nodes per subtree
  }

  if(split_bottlenecks & n_subtrees > 1){
    stop("Split bottlenecks are not yet available in parallel. Set split_bottlenecks = FALSE or n_subtrees = 1.")
  }

  # Names of sequences
  names <- gsub("\\|.*", "", names(fasta))

  # Bookkeeping: move all VCF files into a vcf folder
  vcfs <- list.files(paste0("./", indir), pattern = ".vcf")

  if(length(vcfs) > 0){

    if(!dir.exists(paste0("./", indir, "/vcf"))){
      dir.create(paste0("./", indir, "/vcf"))
    }

    for (i in 1:length(vcfs)) {
      file.rename(paste0("./", indir, "/", vcfs[i]), paste0("./", indir, "/vcf/", vcfs[i]))
    }

  }

  # VCF files present in vcf folder
  vcfs <- list.files(paste0("./", indir, "/vcf"))

  ## List of SNVs present per sample
  message("Processing FASTA and VCF files...")
  snvs <- list()
  pb = txtProgressBar(min = 0, max = n, initial = 0)

  vcf_present <- c()
  messages <- character(0)
  for (i in 1:n) {
    # Locate the correct vcf file
    who <- which(grepl(names[i], vcfs))

    if(length(who) > 1){
      stop(paste0("Multiple VCF files found for sequence ", names[i], ": ", paste(vcfs[who], collapse =  ", "), "."))
    }
    if(length(who) == 1){
      #print(paste0("./", indir, "/vcf/", vcfs[who]))
      # First see if there are any lines to read
      test <- readLines(
        paste0("./", indir, "/vcf/", vcfs[who])
      )
      first_char <- substring(test, 1, 1)
      if(all(first_char == "#")){
        # We use the 8-column LoFreq VCF format, and provide an empty file
        vcf <- data.frame(
          V1 = character(0),
          V2 = integer(0),
          V3 = character(0),
          V4 = character(0),
          V5 = character(0),
          V6 = character(0),
          V7 = character(0),
          V8 = character(0)
        )
      }else{
        # Figure out whether we're using LoFreq VCF or standard format
        col_names_index <- 1
        while (first_char[col_names_index] == "#") {
          col_names_index <- col_names_index + 1
        }
        col_names_index <- col_names_index - 1 # Col names row starts with #CHROM

        col_names <- test[col_names_index]
        col_names <- unlist(strsplit(col_names, "\t"))

        if(identical(
          col_names,
          c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
        )){
          #message("Reading VCF files in LoFreq format...")

          vcf <- read.delim(
            paste0("./", indir, "/vcf/", vcfs[who]),
            colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character"),
            comment.char = "#",
            header = F
          )
          colnames(vcf) <- paste0("V", 1:ncol(vcf))

        }else if(identical(
          col_names,
          c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")

        )){
          #message("Reading VCF files in standard format...")

          vcf <- read.delim(
            paste0("./", indir, "/vcf/", vcfs[who]),
            colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character", "character", "character"),
            comment.char = "#",
            header = F
          )

          # Assemble info column per LoFreq format
          af <- rep("0", nrow(vcf))
          sb <- rep("0", nrow(vcf))
          dp <- rep("0", nrow(vcf))

          for (j in 1:nrow(vcf)) {
            metrics <- strsplit(vcf$V9[j], "\\:")[[1]]
            measurements <- strsplit(vcf$V10[j], "\\:")[[1]]

            for (k in 1:length(metrics)) {
              if(metrics[k] == "AF"){
                af[j] <- measurements[k]
              }
              if(metrics[k] == "SB"){
                sb[j] <- measurements[k]
              }
              if(metrics[k] == "DP"){
                dp[j] <- measurements[k]
              }
            }
          }

          vcf$V8 = paste0("DP=", dp, ";AF=", af, ";SB=", sb)
          vcf <- vcf[, 1:8, drop = F]
          colnames(vcf) <- paste0("V", 1:ncol(vcf))



        }else{
          stop("VCF files are in an unknown format. The supported formats are standard VCF or LoFreq VCF. A standard VCF has column names #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and SAMPLE. A LoFreq VCF has the same column names, except without FORMAT and SAMPLE. Please check VCF files and try again.")

        }


      }

      ## Filter out duplicated iSNV calls
      isnvs <- paste0(vcf$V4, vcf$V2, vcf$V5)
      dup <- duplicated(isnvs)


      if(any(dup)){
        # Get rid of duplicate iSNV calls
        vcf <- vcf[!dup, ]

        # Save a message

        messages <- c(messages, paste0(
          "In the VCF file for sequence ",
          names[i],
          " multiple entries were found for the following iSNV(s): ",
          paste(isnvs[dup], collapse = ", "),
          "."
        ))

      }

      snvs[[i]] <- genetic_info(names[i], root_genome[[1]], fasta[[i]], filters = filters, vcf = vcf)

      # Even if the VCF has no lines, it's present, i.e. we observed no iSNVs
      vcf_present[i] <- TRUE
    }else{
      snvs[[i]] <- genetic_info(names[i], root_genome[[1]], fasta[[i]], filters = filters)
      vcf_present[i] <- FALSE
    }
    setTxtProgressBar(pb,i)

  }

  close(pb)

  if(length(messages) > 0){
    for (m in messages) {
      message(m)
    }
    message("Duplicated iSNV calls will be masked.")
  }

  # Update missing positions in 1st sequence
  if(rooted){
    snvs[[1]]$missing <- root_missing
  }

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
  data$upper_mu <- upper_mu
  data$init_N_eff <- init_N_eff
  data$fixed_N_eff <- fixed_N_eff
  data$upper_N_eff <- upper_N_eff
  data$init_R <- init_R
  data$fixed_R <- fixed_R
  data$init_pi <- init_pi
  data$fixed_pi <- fixed_pi
  data$ongoing <- ongoing
  data$names <- names
  data$s_max <- s_max
  data$safety <- safety
  data$split_bottlenecks <- split_bottlenecks

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
          root_allele <- raw_to_base(root_genome[[1]][pos[p]])
          # If either allele in the iSNV matches the root genome (root), no consensus change in i
          if(data$snvs[[i]]$isnv$a1[ind] == root_allele | data$snvs[[i]]$isnv$a2[ind] == root_allele){
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
        root_allele <- raw_to_base(root_genome[[1]][isnv_pos[p]])

        if(data$snvs[[i]]$isnv$a1[p] == root_allele){
          # Does a1 match the root allele? If so, a1 is in the bottleneck
          mcmc$bot[[i]][p] <- T
        }else if(data$snvs[[i]]$isnv$a2[p] == root_allele){
          # Does a2 match the root allele? If so, it's in the bottleneck
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

  mcmc$R <- init_R
  mcmc$pi <- init_pi

  #mcmc$R <- vals[1] # Reproductive number
  #mcmc$pi <- vals[2] # Probability of sampling
  mcmc$N_eff <- init_N_eff

  # Sequence of times at which the hosts along the edge leading into i were sampled
  mcmc$seq <- list()
  if(rooted){
    mcmc$seq[[1]] <- data$s[1] - (mcmc$a_s / mcmc$lambda_s) + rnorm(1, 0, 0.01)
  }else{
    mcmc$seq[[1]] <- tmrca
  }

  # Time of first infection is one average sojourn interval pre-test
  # Data jittered for intialization of transmission network
  # Cannot be earlier than tmrca!
  mcmc$seq[2:n] <- as.list(
    pmax(
      data$s[2:n] - (mcmc$a_s / mcmc$lambda_s) + rnorm(n-1, 0, 0.01),
      tmrca + 0.01
    )
  )

  # Times at which mutations occur
  mcmc$tmu <- list()

  # Which positions are missing in i and all (direct and indeirect) children of i?
  mcmc$dropout <- list()

  # Which positions are missing everywhere?
  dropout_everywhere <- 1:n_bases

  for (i in 1:n) {
    n_subs <- length(mcmc$subs$from[[i]])
    if(n_subs > 0){
      mcmc$tmu[[i]] <- runif(n_subs, tmrca, mcmc$seq[[i]][1])
    }else{
      mcmc$tmu[[i]] <- numeric(0)
    }

    # Also note which positions drop out at i
    mcmc$dropout[[i]] <- data$snvs[[i]]$missing
    dropout_everywhere <- intersect(dropout_everywhere, mcmc$dropout[[i]])
  }
  # Whichever positions dropout everywhere also dropout in the root
  if(!rooted){
    mcmc$dropout[[1]] <- dropout_everywhere
  }else{
    mcmc$dropout[[1]] <- integer(0)
  }


  # No longer any need to store SNVs relative to root
  for (i in 1:n) {
    data$snvs[[i]]$snv <- NULL
  }

  mcmc <- genotype(mcmc, data, 1)[[1]]


  ## Re-initialize h[i] to nearest genetic neighbor to i infected before i
  ## Can use shift_upstream

  # Order in which to shift nodes
  ord <- sort.int(unlist(mcmc$seq), decreasing = T, index.return = T)$ix

  mut_names <- list()
  for (i in 1:n) {
    mut_names[[i]] <- paste0(mcmc$subs$from[[i]], mcmc$subs$pos[[i]], mcmc$subs$to[[i]])
  }


  for (i in 1:(n-1)) {

    # Which person are we updating the ancestor of?
    who <- ord[i]

    # Distances, genetically, to each other person
    dists <- rep(0, n-i)

    for (j in (i+1):n) {

      # Who is the proposed ancestor?
      anc <- ord[j]

      # Distance from i to j
      # First compute number of shared mutations relative to root
      n_shared <- length(intersect(
        mut_names[[who]], mut_names[[anc]]
      ))

      # Distance is number of mutations in who + number in anc - 2*n_shared
      dists[j-i] <- length(mut_names[[who]]) + length(mut_names[[anc]]) - 2*n_shared
    }

    # Get the indices of "dists" going from lowest/best to highest/worst
    # We will try each one until we find a move that doesn't violate local parsimony
    # Usually the first move will be fine, but occasionally not

    # Choose the ancestor that minimizes the distance
    anc <- (ord[(i+1):n])[which.min(dists)]

    # Update ancestor
    mcmc <- shift(mcmc, data, who, 1, anc, upstream = T)

    # Clear out SNVs that dropout in "who"
    keep <- which(!(mcmc$subs$pos[[who]] %in% mcmc$dropout[[who]]))
    mcmc$subs$from[[who]] <- mcmc$subs$from[[who]][keep]
    mcmc$subs$pos[[who]] <- mcmc$subs$pos[[who]][keep]
    mcmc$subs$to[[who]] <- mcmc$subs$to[[who]][keep]

    # Update genotype at who
    # This will also update mutation times leading into "who"
    mcmc <- genotype(mcmc, data, who)[[1]]

  }

  ## Resample seq for each i, and at the same time, check parsimony at each i

  # Get non-parsimonious hosts
  nonparsimonious <- integer(0)

  for (i in ord) {

    if(i != 1){
      # Only infer unsampled hosts along edge when pi < 1, or when outbreak is ongoing
      if(mcmc$pi < 1 | data$ongoing){
        mcmc$seq[[i]] <- get_ts(mcmc, data, i)
      }
    }

    if(!identical(get_dropout(mcmc, data, i), mcmc$dropout[[i]])){
      print(i)
      print("dropout error 2")
    }

    # Check parsimony while we're at it
    if(!genotype(mcmc, data, i, check_parsimony = T)){
      nonparsimonious <- c(nonparsimonious, i)
    }
  }

  # Correct parsimony errors
  nonparsimonious <- list(sort(nonparsimonious))

  while (length(nonparsimonious[[length(nonparsimonious)]]) > 0) {
    # New not parsimonious
    new <- integer(0)
    for(i in nonparsimonious[[length(nonparsimonious)]]){
      mcmc <- genotype(mcmc, data, i)[[1]]
      neighbors <- which(mcmc$h == i)
      if(i != 1){
        neighbors <- c(mcmc$h[i], neighbors)
      }

      for (j in neighbors) {
        if(!genotype(mcmc, data, j, check_parsimony = T)){
          new <- c(new, j)
        }
      }
    }

    if(any(sapply(nonparsimonious, function(x){identical(x, new)}))){
      print(nonparsimonious)
      print(new)
      stop("could not find locally parsimonious initial configuration.")
    }

    nonparsimonious <- c(nonparsimonious, list(sort(unique(new))))
  }

  # Double-check we've indeed achieved local parsimony
  for (i in 1:n) {
    if(!genotype(mcmc, data, i, check_parsimony = T)){
      stop("Failed local parsimony")
    }
  }


  ## Check that no SNVs are listed in "dropout"
  for (i in 1:n) {
    if(any(
      mcmc$subs$pos[[i]] %in% mcmc$dropout[[i]]
    )){
      print(i)
      print("No mutations should be listed at positions that drop out")
    }
  }

  if(data$rooted){
    data$t_min <- min(data$s, na.rm = T) - 20 * (mcmc$a_s / mcmc$lambda_s)
  }else{
    data$t_min <- (min(data$s, na.rm = T) - 10) * 10 # Set minimum time of anything happening to 10 times earlier than 10 less than the min sampling time
  }

  if(data$ongoing){
    # Can consider updating this by forcing it to converge?
    mcmc$wbar <- wbar(data$t_min, 0, mcmc$R * mcmc$psi / (1 - mcmc$psi), 1 - mcmc$psi, mcmc$pi, mcmc$a_g, 1 / mcmc$lambda_g, mcmc$a_s, 1 / mcmc$lambda_s, 0.1)
  }else{
    mcmc$wbar <- 0 # Placeholder
  }

  # Also track the epidemiological and genomic likelihoods, and prior
  # The genomic likelihood we will store on a per-person basis, for efficiency purposes
  mcmc$e_lik <- sapply(1:n, e_lik_personal, mcmc = mcmc, data = data)
  mcmc$g_lik <- sapply(1:n, g_lik, mcmc = mcmc, data = data)
  mcmc$m_lik <- sapply(1:n, m_lik, mcmc = mcmc, data = data)
  mcmc$prior <- prior(mcmc)

  return(list(mcmc, data))
}
