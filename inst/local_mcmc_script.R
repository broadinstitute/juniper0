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

for (r in 1:d0$n_local) {

  # Move 4
  m0 <- move_seq(m0, d0, also_resample_tmu = F)

  # Move 5
  m0 <- move_seq(m0, d0, also_resample_tmu = T)

  # Move 6
  m0 <- move_w_t(m0, d0)

  # Move 7
  m0 <- move_w_t(m0, d0, recursive = T)

  # Move 8
  m0 <- move_genotype(m0, d0)

  # Move 9
  if(runif(1) < 1/2){
    m0 <- move_h_step(m0, d0)
  }else{
    m0 <- move_h_step(m0, d0, upstream = F)
  }

  # Move 10
  m0 <- move_h_global(m0, d0, biassed = F)

  # Move 11
  m0 <- move_h_global(m0, d0)

  # Move 12
  m0 <- move_swap(m0, d0)

  # Move 13
  m0 <- move_swap(m0, d0, exchange_children = T)

  # Move 14
  if(runif(1) < 1/2){
    m0 <- move_create(m0, d0)
  }else{
    m0 <- move_delete(m0, d0)
  }

  # Move 15
  if(runif(1) < 1/2){
    m0 <- move_create(m0, d0, upstream = F)
  }else{
    m0 <- move_delete(m0, d0, upstream = F)
  }

  # Move 16
  if(runif(1) < 1/2){
    m0 <- move_create(m0, d0, upstream = T, biassed = T)
  }else{
    m0 <- move_delete(m0, d0, upstream = T, biassed = T)
  }

  # Append new results
  if(r %% d0$sample_every == 0){
    res <- c(res, list(m0))
  }
}

save(res, file = paste0("tmp/res_", j, ".RData"))









