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

# Schedule of moves within a subtree

j <- as.numeric(commandArgs(TRUE))
load(paste0("tmp/mcmcs_", j, ".RData"))
load(paste0("tmp/datas_", j, ".RData"))

res <- list()

for (r in 1:data$n_local) {

  # Move 4
  mcmc <- move_seq(mcmc, data, also_resample_tmu = F)

  # Move 5
  mcmc <- move_seq(mcmc, data, also_resample_tmu = T)

  # Move 6
  mcmc <- move_w_t(mcmc, data)

  # Move 7
  mcmc <- move_w_t(mcmc, data, recursive = T)

  # Move 8
  mcmc <- move_genotype(mcmc, data)

  # Move 9
  if(runif(1) < 1/2){
    mcmc <- move_h_step(mcmc, data)
  }else{
    mcmc <- move_h_step(mcmc, data, upstream = F)
  }

  # Move 10
  mcmc <- move_h_global(mcmc, data, biassed = F)

  # Move 11
  mcmc <- move_h_global(mcmc, data)

  # Move 12
  mcmc <- move_swap(mcmc, data)

  # Move 13
  mcmc <- move_swap(mcmc, data, exchange_children = T)

  # Move 14
  if(runif(1) < 1/2){
    mcmc <- move_create(mcmc, data)
  }else{
    mcmc <- move_delete(mcmc, data)
  }

  # Move 15
  if(runif(1) < 1/2){
    mcmc <- move_create(mcmc, data, upstream = F)
  }else{
    mcmc <- move_delete(mcmc, data, upstream = F)
  }

  # Move 16
  if(runif(1) < 1/2){
    mcmc <- move_create(mcmc, data, upstream = T, biassed = T)
  }else{
    mcmc <- move_delete(mcmc, data, upstream = T, biassed = T)
  }

  # Append new results
  if(r %% data$sample_every == 0){
    res <- c(res, list(mcmc))
  }
}

save(res, file = paste0("tmp/res_", j, ".RData"))









