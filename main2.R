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

### Execute large-scale outbreak reconstruction algorithm, in parallel over many cores
set.seed(220)
## Libraries
library(ape)
library(Rcpp)
#library(igraph)
#library(ggraph)
library(parallel)

source("likelihood.R")
source("moves.R")
source("prior.R")
source("subroutines.R")
source("initialize.R")
source("global_mcmc.R")
source("local_mcmc.R")

load("data_init.RData")
load("mcmc_init.RData")

### M-H algo
output <- list()

liks <- c()

for (r in 1:data$n_global) {

  # For reproducible results
  set.seed(r)

  # Make global moves
  mcmc <- global_mcmc(mcmc, data)

  # Chop up the tree into pieces
  breakdowns <- breakdown(mcmc, data)
  mcmcs <- breakdowns[[1]]
  datas <- breakdowns[[2]]

  message(paste("Parallelizing over", length(mcmcs), "cores..."))

  ## Write each mcmcs[[i]] and datas[[i]] to its own directory
  dir.create("./state/substate/")
  cmd <- mclapply(1:length(mcmcs), function(i){
    mcmc_tmp <- mcmcs[[i]]
    data_tmp <- datas[[i]]
    dir.create(paste0("./state/substate/tree_", i))
    save(mcmc_tmp, file = paste0("./state/substate/tree_", i, "/mcmc.RData"))
    save(data_tmp, file = paste0("./state/substate/tree_", i, "/data.RData"))
    return(paste0("Rscript move_local.R ", i))
  },
  mc.cores = data$n_subtrees)


  cmd <- paste(unlist(cmd), collapse = " & ")
  cmd <- paste(cmd, "& wait")

  # Run on subdirectories in parallel
  system(cmd)

  ## Compile results. Parallelize

  all_res <- mclapply(1:length(mcmcs), function(i){
    load(paste0("./state/substate/tree_", i, "/res.RData"))
    return(res)
  },
  mc.cores = data$n_subtrees)

  # (or in series....)
  # all_res <- list()
  # for (i in 1:length(mcmcs)) {
  #   load(paste0("./state/substate/tree_", i, "/res.RData"))
  #   all_res[[i]] <- res
  # }

  # Amalgamate results of parallel MCMC run
  amalgam <- amalgamate(all_res, mcmcs, datas, mcmc, data)

  # For cleanup: delete substate directory
  unlink("./state/substate/", recursive=TRUE)

  # Record amalgamated results, filtering to parameters of interest
  for (i in 1:length(amalgam)) {
    output <- c(output, list(
      amalgam[[i]][data$record]
    ))
  }

  # "mcmc" is now the most recent result
  mcmc <- amalgam[[length(amalgam)]]

  #print(r)

  liks <- c(liks, mcmc$e_lik + sum(mcmc$g_lik[2:mcmc$n]) + mcmc$prior)


  message(paste(r, "global iterations complete. Log-likelihood =", round(liks[r], 2)))
  #print(plot_current(mcmc$h, data$n_obs))
  #print(mcmc$w)

  # if(r == 10){
  #   data$n_subtrees <- 3
  # }

}

save(output, file = "output.RData")
