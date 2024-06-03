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

# Update global parameters

global_mcmc <- function(mcmc, data){

  # Move 1
  mcmc <- moves$b(mcmc, data)

  mcmc <- moves$a_g(mcmc, data)

  # Move 6
  mcmc <- moves$mu(mcmc, data)

  # Move 7
  mcmc <- moves$p(mcmc, data)

  # Move 8
  #mcmc <- moves$v(mcmc, data)
  #mcmc <- moves$lambda(mcmc, data)

  # mcmc <- moves$rho(mcmc, data)
  # mcmc <- moves$psi(mcmc, data)

  # We are fixing parameters associated with moves 2-5, 9-10
  return(mcmc)
}

## After making global moves, chop up the tree, and make one copy of "mcmc" for each subtree
breakdown <- function(mcmc, data){

  subtrees <- chop(mcmc, data)
  n_subtrees <- length(subtrees[[1]])
  mcmcs <- list()
  datas <- list()
  for (i in 1:n_subtrees) {
    mcmcs[[i]] <- mcmc
    mcmcs[[i]]$root <- subtrees[[1]][i]
    mcmcs[[i]]$cluster <- subtrees[[2]][[i]]
    if(mcmcs[[i]]$root == 1){
      mcmcs[[i]]$anc_cluster <- NA
    }else{
      for (j in 1:n_subtrees) {
        if(mcmcs[[i]]$root %in% c(subtrees[[1]][j], subtrees[[2]][[j]])){
          mcmcs[[i]]$anc_cluster <- j
        }
      }
    }


    joined <- c(mcmcs[[i]]$root, mcmcs[[i]]$cluster)

    # To save memory: extract only necessary components of MCMC and data
    mcmcs[[i]]$n <- length(joined)
    mcmcs[[i]]$h <- mcmcs[[i]]$h[joined]
    mcmcs[[i]]$h <- match(mcmcs[[i]]$h, joined)
    mcmcs[[i]]$w <- mcmcs[[i]]$w[joined]
    mcmcs[[i]]$t <- mcmcs[[i]]$t[joined]
    mcmcs[[i]]$m01 <- mcmcs[[i]]$m01[joined]
    mcmcs[[i]]$m10 <- mcmcs[[i]]$m10[joined]
    mcmcs[[i]]$m0y <- mcmcs[[i]]$m0y[joined]
    mcmcs[[i]]$m1y <- mcmcs[[i]]$m1y[joined]
    mcmcs[[i]]$mx0 <- mcmcs[[i]]$mx0[joined]
    mcmcs[[i]]$mx1 <- mcmcs[[i]]$mx1[joined]
    mcmcs[[i]]$mxy <- mcmcs[[i]]$mxy[joined]
    mcmcs[[i]]$d <- mcmcs[[i]]$d[joined]
    mcmcs[[i]]$g_lik <- mcmcs[[i]]$g_lik[joined]

    datas[[i]] <- data
    datas[[i]]$s <- datas[[i]]$s[joined]
    datas[[i]]$n_obs <- sum(joined <= data$n_obs)
    datas[[i]]$snvs <- datas[[i]]$snvs[joined]
    datas[[i]]$frozen <- setdiff(which(joined %in% subtrees[[1]]), 1)

    # Degree of roots of other trees needs to be set to 0
    mcmcs[[i]]$d[datas[[i]]$frozen] <- 0

    # e_lik needs to be re-computed based on the smaller set
    mcmcs[[i]]$e_lik <- e_lik(mcmcs[[i]], datas[[i]])

  }

  return(list(mcmcs, datas))

}
