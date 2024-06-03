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

### Prior function and defaults

prior <- function(mcmc, priors = NULL){
  if(is.null(priors)){
    priors <- list()
  }
  if(is.null(priors$b)){
    priors$b <- function(x){
      # Unif(0,1)
      if(0 <= x & x <= 1){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$a_g)){
    priors$a_g <- function(x){
      # Unif(0,100)
      # if(0 <= x & x <= 100){
      #   return(-log(100))
      # }else{
      #   return(-Inf)
      # }
      return(dgamma(x, 5, 1, log = T))
    }
  }

  if(is.null(priors$lambda_g)){
    priors$lambda_g <- function(x){
      # Unif(0,5) on variance
      x <- sqrt(5/x)
      if(0 <= x & x <= 5){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$a_s)){
    priors$a_s <- function(x){
      # Unif(0,20)
      if(0 <= x & x <= 20){
        return(-log(20))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$lambda_s)){
    priors$lambda_s <- function(x){
      # Unif(0,10) on variance
      x <- sqrt(5/x)
      if(0 <= x & x <= 10){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$mu)){
    priors$mu <- function(x){
      # Unif(0,0.01)
      if(0 <= x & x <= 0.01){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$p)){
    priors$p <- function(x){
      # Unif(0,0.01)
      if(0 <= x & x <= 0.01){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$v)){
    priors$v <- function(x){
      # Unif(1,1e4)
      if(0 <= x & x <= 1e4){
        return(-log(1e4))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$lambda)){
    priors$lambda <- function(x){
      dnorm(x, 2, 1, log = T)
    }
  }

  if(is.null(priors$rho)){
    priors$rho <- function(x){
      x <- x*(1-0.1)/0.1
      # Unif(0,5)
      if(0 <= x & x <= 5){
        return(-log(5))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$psi)){
    priors$psi <- function(x){
      # Unif(0,100)
      if(0 <= x & x <= 100){
        return(-log(100))
      }else{
        return(-Inf)
      }
    }
  }

  priors$b(mcmc$b) +
    priors$a_g(mcmc$a_g) +
    priors$a_s(mcmc$a_s) +
    priors$mu(mcmc$mu) +
    priors$p(mcmc$p) +
    priors$v(mcmc$v) +
    priors$lambda(mcmc$lambda) +
    priors$rho(mcmc$rho) +
    priors$psi(mcmc$psi)

}






