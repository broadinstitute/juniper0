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

### Execute large-scale outbreak reconstruction algorithm
set.seed(232)
## Libraries
library(ape)
library(Rcpp)
library(ggplot2)
library(igraph)
library(ggraph)
library(cowplot)
library(parallel)

source("likelihood.R")
source("moves.R")
source("prior.R")
source("subroutines.R")
sourceCpp("cpp_subroutines.cpp")
source("initialize.R")
source("global_mcmc.R")
source("local_mcmc.R")

### M-H algo
run_mcmc <- function(mcmc, data, noisy = F){
  output <- list()
  liks <- c()

  for (r in 1:data$n_global) {

    # For reproducible results
    #set.seed(r)

    # Make global moves
    mcmc <- global_mcmc(mcmc, data)

    # Chop up the tree into pieces
    breakdowns <- breakdown(mcmc, data)
    mcmcs <- breakdowns[[1]]
    datas <- breakdowns[[2]]

    if(noisy){
      message(paste("Parallelizing over", length(mcmcs), "cores..."))
    }

    all_res <- parallel::mclapply(
      1:length(mcmcs),
      function(i, mcmcs, datas){
        local_mcmc(mcmcs[[i]], datas[[i]])
      },
      mcmcs = mcmcs,
      datas = datas,
      mc.set.seed = F,
      mc.cores = length(mcmcs)
    )
    #...or run in series
    # all_res <- list()
    # for (j in 1:length(mcmcs)) {
    #   all_res[[j]] <- local_mcmc(mcmcs[[j]], datas[[j]])
    # }

    # Amalgamate results of parallel MCMC run
    amalgam <- amalgamate(all_res, mcmcs, datas, mcmc, data)

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


    if(noisy){
      message(paste(r, "global iterations complete. Log-likelihood =", round(liks[r], 2)))
      print(plot_current(mcmc$h, data$n_obs))
      print(mcmc$w)
       print(mcmc$mu)
       print(mcmc$p)
       print(mcmc$a_g)
       print(mcmc$rho * (1-mcmc$psi) / mcmc$psi)
       print(mcmc$prior)
       #print(mcmc$rho * (1 - mcmc$psi) / mcmc$psi)
      # print(length(unlist(mcmc$m01)) + length(unlist(mcmc$m10)))
      # print(length(unlist(mcmc$mx1)))
      #print(data$s - mcmc$t[1:data$n_obs])
      #print(mcmc$lambda)
      #print(mcmc$h)
      # print(mcmc$a_g)
    }


    # if(r == 10){
    #   data$n_subtrees <- 3
    # }

  }
  return(list(
    liks, output
  ))
}



