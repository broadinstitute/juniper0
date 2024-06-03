# Demonstration of outbreak reconstruction algorithm
set.seed(1)
source("main.R")

## Unhash these for regular run

init <- initialize(
  init_mst = T,
  pooled_coalescent = T,
  disjoint_coalescent = F,
  n_global = 100,
  n_subtrees = 1,
  virus = "H5N1"
)

mcmc <- init[[1]]
data <- init[[2]]

out <- run_mcmc(mcmc, data, noisy = T)
#save(out, file = "output.RData")



n_reps <- length(out[[1]])
plot(out[[1]])

adj <- matrix(0, ncol = data$n_obs, nrow = data$n_obs)
direct <- matrix(0, ncol = data$n_obs, nrow = data$n_obs)
burnin <- n_reps * 0.2 + 1

# Generations per transmission
bs <- c()
mus <- c()
ps <- c()
for (i in burnin:n_reps) {
  bs <- c(bs, out[[2]][[i]]$b)
  mus <- c(mus, out[[2]][[i]]$mu)
  ps <- c(ps, out[[2]][[i]]$p)

  trans <- cbind(out[[2]][[i]]$h, 1:out[[2]][[i]]$n)
  trans[1,1] <- 0
  direct_trans <- trans[trans[,1] <= data$n_obs & trans[,1] >= 1 & trans[,2] <= data$n_obs & trans[,2] >= 1 & out[[2]][[i]]$w == 0, ]
  trans <- trans[trans[,1] <= data$n_obs & trans[,1] >= 1 & trans[,2] <= data$n_obs & trans[,2] >= 1, ]

  direct[direct_trans] <- direct[direct_trans] + 1
  adj[trans] <- adj[trans] + 1
}
adj <- adj / (n_reps - burnin + 1)
direct <- direct / (n_reps - burnin + 1)

# Load the reference sequence
ref_genome <- read.FASTA("./input_data/ref.fasta")

# Load the FASTA of sequences
fasta <- read.FASTA("./input_data/aligned.fasta")

# The first genome is itself the ref genome
fasta <- c(ref_genome, fasta)

View(as.matrix(dist.dna(fasta, "N", pairwise.deletion = T)))

rownames(adj) <- names(fasta)
colnames(adj) <- names(fasta)
rownames(direct) <- names(fasta)
colnames(direct) <- names(fasta)

write.csv(direct, "direct_transmissions.csv")
