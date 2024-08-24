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

  mcmc <- move_b(mcmc, data)

  if(!data$fixed_mu){
    mcmc <- move_mu(mcmc, data)
  }

  mcmc <- move_p(mcmc, data)
  mcmc <- move_pi(mcmc, data)
  mcmc <- move_R(mcmc, data)

  return(mcmc)
}

## After making global moves, chop up the tree, and make one copy of "mcmc" for each subtree
breakdown <- function(mcmc, data, old_roots){

  subtrees <- chop(mcmc, data, old_roots)
  n_subtrees <- length(subtrees[[1]]) # number of subtrees, indexed 1, 2, ..., n_subtrees
  mcmcs <- list()
  datas <- list()

  for (i in 1:n_subtrees){
    # Initialize mcmc for each subtree
    mcmcs[[i]] <- mcmc
    # Initialize a vector of nodes in subtree i that are roots of OTHER subtrees
    mcmcs[[i]]$external_roots <- integer()
    # Initialize a vector of WHICH OTHER subtree each one is a root of
    mcmcs[[i]]$external_subtrees <- integer()
  }
  for (i in 1:n_subtrees) {

    # Root of subtree i
    root <- subtrees[[1]][i]

    # People in subtree i, excluding root
    cluster <- subtrees[[2]][[i]]

    if(root == 1){
      # Record which subtree is ancestral to subtree i
      mcmcs[[i]]$anc_cluster <- NA
    }else{
      for (j in 1:n_subtrees) {
        if(root %in% c(subtrees[[1]][j], subtrees[[2]][[j]]) & j != i){

          # Record which subtree is ancestral to subtree i
          mcmcs[[i]]$anc_cluster <- j

          # Record which node in subtree j is a root, and of which subtree it's a root
          mcmcs[[j]]$external_roots <- c(mcmcs[[j]]$external_roots, root)
          mcmcs[[j]]$external_subtrees <- c(mcmcs[[j]]$external_subtrees, i)
        }
      }
    }
  }

  for (i in 1:n_subtrees) {

    # Root of subtree i
    root <- subtrees[[1]][i]

    # People in subtree i, excluding root
    cluster <- subtrees[[2]][[i]]

    joined <- c(root, cluster)

    # To save memory: extract only necessary components of MCMC and data
    mcmcs[[i]]$n <- length(joined)
    mcmcs[[i]]$h <- mcmcs[[i]]$h[joined]
    mcmcs[[i]]$h <- match(mcmcs[[i]]$h, joined)
    mcmcs[[i]]$seq <- mcmcs[[i]]$seq[joined]
    mcmcs[[i]]$m01 <- mcmcs[[i]]$m01[joined]
    mcmcs[[i]]$m10 <- mcmcs[[i]]$m10[joined]
    mcmcs[[i]]$m0y <- mcmcs[[i]]$m0y[joined]
    mcmcs[[i]]$m1y <- mcmcs[[i]]$m1y[joined]
    mcmcs[[i]]$mx0 <- mcmcs[[i]]$mx0[joined]
    mcmcs[[i]]$mx1 <- mcmcs[[i]]$mx1[joined]
    mcmcs[[i]]$mxy <- mcmcs[[i]]$mxy[joined]
    mcmcs[[i]]$d <- mcmcs[[i]]$d[joined]
    mcmcs[[i]]$g_lik <- mcmcs[[i]]$g_lik[joined]

    # The roots of external clusters have been relabeled
    mcmcs[[i]]$external_roots <- match(mcmcs[[i]]$external_roots, joined)

    datas[[i]] <- data
    datas[[i]]$s <- datas[[i]]$s[joined]
    datas[[i]]$n_obs <- sum(joined <= data$n_obs)
    # If the root is unobserved, we treat it as observed anyway, because it cannot be added / deleted
    if(root > data$n_obs){
      datas[[i]]$n_obs <- datas[[i]]$n_obs + 1
    }

    datas[[i]]$snvs <- datas[[i]]$snvs[joined]
    datas[[i]]$frozen <- setdiff(which(joined %in% subtrees[[1]]), 1)

    datas[[i]]$observed_root <- (root < data$n_obs)

    # Root always treated as fixed within subtree, UNLESS it's the root subtree
    if(root != 1){
      datas[[i]]$rooted <- T
    }

    # e_lik needs to be re-computed based on the smaller set
    mcmcs[[i]]$e_lik <- e_lik(mcmcs[[i]], datas[[i]])

    # Record only observed hosts in each cluster
    mcmcs[[i]]$cluster <- joined[joined <= data$n_obs]

  }

  return(list(mcmcs, datas, subtrees[[1]])) # Last entry is the vector of roots, used later

}
