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

### Simulate outbreak

set.seed(211)

library(reshape2)
library(stringr)
library(ape)

source("main.R")

### Generate transmission network

# Epi params
a_g <- 5
lambda_g <- 1
a_s <- 5
lambda_s <- 1
mu <- 3e-6
p <- 1.5e-6
v <- 1000 # virions produced per cycle
b <- 0.05 # probability bottleneck has size 2
k <- round(1/sqrt(p)) # number of virions at end of exponential growth phase
N_bases <- 29903 # length of viral genome
rho <- 0.5
psi <- 0.5 / (1.5 + 0.5)
N <- 1e9

# Probability of a mutation in exponential growth phase
p_growth_mut <- 1 - (1-p)^k

# Number of generations to simulate
G <- 7

get_adj <- function(run, n_obs){
  out <- matrix(0, nrow = n_obs, ncol = n_obs)
  len <- length(run)
  for (i in (len*0.1 + 1):len) {
    h <- run[[i]]$h
    present <- which(!is.na(h) & h <= n_obs)
    present <- present[present <= n_obs]
    out[cbind(h[present], present)] <- out[cbind(h[present], present)] + 1
  }
  return(out / (0.9*len))
}

run_sim <- function(index, misspecify = 0){

  n <- 0
  while (n < 10 | n > 25) {
    h <- NA
    t <- -5
    s <- 0

    # Proportions of particles that are mutated
    props <- list(rep(0, N_bases))

    for (g in 1:G) {
      # Who's in the current generation?
      if(g == 1){
        genr <- 1
      }else{
        genr <- which(h %in% genr)
      }

      # Loop over people in current generation
      for (i in genr) {
        # Inherited genotype
        if(i != 1){
          # Proportions of inherited mutations
          anc_props <- props[[h[i]]]

          # Evolve this per JC
          delta_t <- t[i] - (t[h[i]] + log(1/sqrt(p)) / (mu / p) / log(v))

          # Evolved anc_props
          #anc_props[anc_props == 0] <- 1/4 - (1/4)*exp((-4/3)*mu*delta_t)
          #anc_props[anc_props == 1] <- 1 - (1/4 - (1/4)*exp((-4/3)*mu*delta_t))
          anc_props <- 1/4 + (anc_props - 1/4)*exp((-4/3)*mu*delta_t)

          props[[i]] <- rep(0, N_bases)

          if(runif(1) < b){
            # Which sites have a split bottleneck?
            split <- which(runif(N_bases) < 2*anc_props*(1 - anc_props))
            props[[i]][split] <- runif(length(split))

            # Which sites inherit the alternate allele?
            alt <- which(runif(N_bases) < anc_props^2 / (anc_props^2 + (1 - anc_props)^2))
            alt <- setdiff(alt, split)
          }else{
            split <- integer(0)
            alt <- which(runif(N_bases) < anc_props)
          }

          # Which sites pick up an iSNV?
          isnv <- which(runif(N_bases) < p_growth_mut)
          isnv <- setdiff(isnv, split)

          isnv_props <- rbeta(length(isnv), 1, sample(1:k, length(isnv), replace = T))

          props[[i]][setdiff(isnv, alt)] <- isnv_props[!(isnv %in% alt)]
          props[[i]][intersect(isnv, alt)] <- 1 - isnv_props[isnv %in% alt]
          props[[i]][setdiff(alt, isnv)] <- 1



        }

        if(g != G){
          # Generate kids
          #n_kids <- ifelse(runif(1) < 0.1, 2, 1)
          n_kids <- rnbinom(1, rho, psi)
          #n_kids <- 1

          if(n_kids > 0){
            # What are their indices?
            who <- (length(h) + 1):(length(h) + n_kids)

            h[who] <- i
            t[who] <- t[i] + pmax(rgamma(length(who), a_g, lambda_g), log(1/sqrt(p)) / (mu / p) / log(v) + 1/2) # To ensure positive evolutionary time
            s[who] <- t[who] + rgamma(length(who), a_s, lambda_s)

          }
        }
      }
    }

    # Number of people
    n <- length(s)
    print("L")
  }



  names <- paste(1:n)
  names <- str_pad(names, 3, pad = "0")
  names <- paste0("sim_", names)

  # For writing files, move into new directory
  setwd("./simulated-outbreaks/")
  dirname <- paste0("outbreak_", index)
  dir.create(paste0("./", dirname), showWarnings = F)
  setwd(paste0("./", dirname))
  dir.create("./input_data", showWarnings = F)
  dir.create("./input_data/vcf", showWarnings = F)

  # Write fastas
  cons <- list()
  for (i in 2:n) {
    newseq <- rep("A", N_bases)
    newseq[props[[i]] > 1/2] <- "C"
    cons <- c(cons, list(newseq))
  }

  names(cons) <- names[2:n]
  write.dna(cons, file = "./input_data/aligned.fasta", format = "fasta")

  # Write ref genome
  ref <- list(rep("A", N_bases))
  names(ref) <- "ref"
  write.dna(ref, file = "./input_data/ref.fasta", format = "fasta")

  # Write test date table
  dates <- cbind(names(cons), round(s[2:n]))
  write.csv(dates, "./input_data/date.csv", row.names = F, quote = F)

  # First, remove existing vcfs, just in case
  do.call(file.remove, list(list.files("./input_data/vcf/", full.names = TRUE)))

  # Write VCFs
  for (i in 2:n) {
    # Which positions have iSNVs in read data?
    pos <- which(props[[i]] > 0 & props[[i]] < 1)
    af <- props[[i]][pos]
    alt_reads <- round(10000 * af)
    info <- paste0(
      "DP=",
      10000,
      ";AF=",
      format(round(alt_reads / 10000, 6), scientific = F),
      ";SB=0;DP4=0,", # strand bias wouldn't actually be 0 here, but assuming no strand bias in synthetic data
      10000 - alt_reads,
      ",0,",
      alt_reads
    )
    vcf <- data.frame(
      CHROM = "ref",
      POS = pos,
      ID = ".",
      REF = "A",
      ALT = "C",
      QUAL = 5000,
      FILTER = "PASS",
      INFO = info
    )

    write.table(vcf, file = paste0("./input_data/vcf/", names[i], ".vcf"), col.names = F, row.names = F)
  }

  ### Reconstruction


  init <- initialize(
    pooled_coalescent = F,
    disjoint_coalescent = T
  )
  mcmc <- init[[1]]
  data <- init[[2]]

  mcmc$rho <- 0.5
  mcmc$psi <- 0.5/(0.5 + 1.5)

  if(misspecify == 1){
    mcmc$psi = 0.5/(0.5 + 2.5)
  }

  if(misspecify == 2){
    mcmc$a_g <- 6
  }

  if(misspecify == 3){
    mcmc$a_s = 6
  }

  data$n_local = 10
  data$sample_every = 10
  data$n_global = 100

  output <- run_mcmc(mcmc, data, noisy = T)

  ## Analysis

  adj <- get_adj(output[[2]], data$n_obs)

  # Overall accuracy
  acc <- mean(adj[cbind(h[2:n], 2:n)])

  # What proportion of ancestral choices with a probability >50% are correct?
  hits <- c()
  for (i in 1:n) {
    if(any(adj[,i] > 0.5)){
      who <- which(adj[,i] > 0.5)
      if(h[i] == who){
        hits <- c(hits, 1)
      }else{
        hits <- c(hits, 0)
      }
    }
  }
  hits <- ifelse(length(hits) == 0, NA, mean(hits))

  # Does a (minimum-cardinality) 50% posterior interval contain the correct ancestor?
  cover <- c()
  for (i in 1:n) {
    ord <- sort.int(adj[,i], index.return = T, decreasing = T)$ix
    csum <- cumsum(adj[ord,i])
    if(max(csum) > 0.5){
      first <- min(which(csum > 0.5))
      who <- ord[1:first]
      if(h[i] %in% who){
        cover <- c(cover, 1)
      }else{
        cover <- c(cover, 0)
      }
    }
  }
  cover <- ifelse(length(cover) == 0, NA, mean(cover))

  # for (i in 1:(n-1)) {
  #   print(sum(cons[[i]] == "C"))
  # }

  ## Go back to root directory
  setwd('..')
  setwd('..')

  # Return accuracy
  return(
    c(
      acc,
      hits,
      cover
    )
  )

}

plots <- list()

for (i in 0:3) {
  sim_stats <- mclapply(1:100, run_sim, misspecify = i, mc.cores = 12)
  sim_stats <- as.data.frame(matrix(unlist(sim_stats), ncol = 3, byrow = T))
  colnames(sim_stats) <- c("Accuracy", "Hit Rate", "Coverage")
  sim_stats <- melt(sim_stats)
  colnames(sim_stats)[1] <- "Metric"

  plots[[i+1]] <- ggplot(sim_stats, aes(x = Metric)) +
    geom_boxplot(aes(y=value, color = Metric, fill = Metric)) +
    scale_fill_manual(values = c("pink", "lightgreen", "lightblue")) +
    scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
    xlab("Metric") +
    ylab("Value") +
    theme_minimal() +
    theme(legend.position = 'none')
}

combined <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], labels = LETTERS)

ggsave(combined, file = "./figs/simulations.pdf", width = 6.5, height = 6.5)
ggsave(combined, file = "./figs/simulations.png", width = 6.5, height = 6.5)
