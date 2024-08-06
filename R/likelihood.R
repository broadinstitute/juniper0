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

# Compute epidemiological log likelihood
e_lik <- function(mcmc, data){

  if(
    any(mcmc$w < 0) |
      (!data$rooted & mcmc$d[1] > 1)

    # mcmc$a_g < 0 |
    # mcmc$lambda_g < 0 |
    # mcmc$a_s < 0 |
    # mcmc$lambda_s < 0 |
    # mcmc$rho < 0 |
    # mcmc$psi < 0
  ){
    return(-Inf)
  }else{
    # Which people do we compute the epi prior for?
    if(data$rooted){
      who <- 2:mcmc$n
    }else{
      who <- which(mcmc$h != 1)

    }

    if(data$experimental){

      if(any(mcmc$t > max(data$s, na.rm = T))){
        return(-Inf)
      }

      t_max <- max(data$s, na.rm = T)

      ids <- c()
      hs <- c()
      counter <- mcmc$n
      for (i in 1:mcmc$n) {
        if(mcmc$w[i] == 0){
          ids <- c(ids, i)
          hs <- c(hs, mcmc$h[i])
        }else{
          ids <- c(ids, i, (counter+1):(counter+mcmc$w[i]))
          hs <- c(hs, (counter+1):(counter+mcmc$w[i]), mcmc$h[i])
          counter <- counter + mcmc$w[i]
        }
      }
      t_inf <- unlist(mcmc$seq)

      hs[ids] <- hs
      t_inf[ids] <- t_inf
      #ids[ids] <- ids
      t_samp <- c(data$s[1:data$n_obs], rep(NA, mcmc$n - data$n_obs + sum(mcmc$w)))

      ttree <- matrix(c(t_inf, t_samp, hs), ncol = 3, byrow = F)
      ttree[1, 3] <- 0 # Ancestor of person 1 designated as 0

      for (i in 2:nrow(ttree)) {
        print(ttree[i, 1] - ttree[ttree[i,3], 1])
      }

      return(
        probTTree(
          ttree, mcmc$rho, 1-mcmc$psi, mcmc$pi, mcmc$a_g, 1/mcmc$lambda_g, mcmc$a_s, 1/mcmc$a_s, t_max, delta_t = 0.1
        )
      )

    }else{
      return(
        # Generation intervals
        sum(dgamma(mcmc$t[who] - mcmc$t[mcmc$h[who]], shape = (mcmc$w[who] + 1) * mcmc$a_g, rate = mcmc$lambda_g, log = T)) +

          # Sojourn intervals
          sum(dgamma(data$s[2:data$n_obs] - mcmc$t[2:data$n_obs], shape = mcmc$a_s, rate = mcmc$lambda_s, log = T)) +

          # xi-coalescent
          ifelse(
            mcmc$rho != Inf,
            sum(lfactorial(mcmc$d + mcmc$rho - 1)) - mcmc$n * lfactorial(mcmc$rho - 1) + sum(mcmc$d) * log((1-mcmc$psi) / mcmc$psi) + sum(mcmc$w[who]) * log((mcmc$rho * (1 - mcmc$psi) / mcmc$psi)),
            (sum(mcmc$d) + sum(mcmc$w[who])) * log(data$R)
          ) #+
        #ifelse(data$rooted, 0, correction)


      )
    }




  }
}

# Compute genomic log likelihood for each person

g_lik <- function(mcmc, data, i){

  if(mcmc$v < 0 | mcmc$mu < 0 | mcmc$p < 0 | mcmc$b < 0 | mcmc$lambda < 0 | mcmc$b > 1 | mcmc$w[i] < 0){
    return(-Inf)
  }else{
    # Time of end of expo growth phase for ancestor of i, which is when we reach k = 1/sqrt(p) virions
    #g <- mcmc$t[mcmc$h[i]] - (mcmc$p * log(mcmc$p) * (-1 + 0.577215664901533 + log(mcmc$v)))/(2*mcmc$mu * log(mcmc$v))

    #g <- mcmc$t[mcmc$h[i]] - log(mcmc$p) / (2*log(mcmc$mu / mcmc$p))
    #g <- mcmc$t[mcmc$h[i]] + (mcmc$p *(-(log(mcmc$p)/2) + log(mcmc$v)))/(mcmc$mu* log(mcmc$v))
    #mcmc$p <- mcmc$mu/mcmc$lambda
    h <- mcmc$h[i]

    # Time of end of exponential growth phase in h
    t_g <- mcmc$t[h] + log(1/sqrt(mcmc$p)) / (mcmc$mu / mcmc$p) / log(mcmc$v)
    #g <- mcmc$t[h] + 1

    # Evolutionary time
    delta_t <- mcmc$t[i] - mcmc$t[h]

    # If i unobserved, or has no iSNV info provided, simply evolve from end of growth phase in h[i] to end of growth phase in i
    isnv_info <- TRUE
    if(i > data$n_obs){
      isnv_info <- FALSE
    }
    if(i <= data$n_obs){
      if(!data$vcf_present[i]){
        isnv_info <- FALSE
      }
    }

    # if(!isnv_info){
    #   delta_t <- mcmc$t[i] - mcmc$t[h]
    # }

    # Evolutionary time from first downstream host of h to infection time of i, approx
    if(is.infinite(delta_t)){
      delta_t_prime <- Inf
    }else{
      if(data$experimental){
        # Time of first transmission on the chain from h to i
        t_1st_trans <- mcmc$seq[[i]][length(mcmc$seq[[i]])]
        delta_t_prime <- mcmc$t[i] - t_1st_trans

        if(t_1st_trans < mcmc$t[h]){
          return(-Inf)
        }
      }else{
        delta_t_prime <- delta_t * mcmc$w[i] / (mcmc$w[i] + 1)
      }
    }




    if(delta_t <= 0){
      -Inf
    }else{
      # log probability of SPECIFIC iSNV in expo growth phase
      log_p_isnv <- log(
        (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * (1 - denovo_cdf(data$filters$af, mcmc$p)) / 3
      )

      # log probability of no iSNV in expo growth phase
      log_p_no_isnv <- log(
        (1 - mcmc$p)^(1/sqrt(mcmc$p)) +
          (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * denovo_cdf(data$filters$af, mcmc$p)
      )

      # Number of sites without a mutation
      no_mut <- data$n_bases -
        length(mcmc$m01[[i]]) -
        length(mcmc$m10[[i]]) -
        length(mcmc$m0y[[i]]) -
        length(mcmc$m1y[[i]]) -
        length(mcmc$mx0[[i]]) -
        length(mcmc$mx1[[i]]) -
        length(mcmc$mxy[[i]]) #-
        #length(data$filters$common)

      # print(1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t))
      # print(evolveJC(1, mcmc$mu, delta_t))

      # Likelihood from x = 0, y = 0 or x = 1, y = 1
      out <- no_mut * (log(evolveJC(1, mcmc$mu, delta_t)) + ifelse(isnv_info, log_p_no_isnv, 0)) +

        # Likelihood from x = 0, y = 1 and x = 1, y = 0
        (length(mcmc$m01[[i]]) + length(mcmc$m10[[i]])) * (log(evolveJC(0, mcmc$mu, delta_t)) + ifelse(isnv_info, log_p_no_isnv, 0))

      # If there actually are any reported iSNVs...
      if(i <= data$n_obs){
        if(length(data$snvs[[i]]$isnv$call) > 0){

          # Frequencies of added iSNVs
          freq_0y <- data$snvs[[i]]$isnv$af[match(mcmc$m0y[[i]], data$snvs[[i]]$isnv$call)]

          # Frequencies of deleted iSNVs
          freq_1y <- data$snvs[[i]]$isnv$af[match(mcmc$m1y[[i]], data$snvs[[i]]$isnv$call)]

          out <- out +
            # Likelihood from x = 0, 0 < y < 1
            length(mcmc$m0y[[i]]) * log_p_isnv +
            sum(log(
              (evolveJC(1, mcmc$mu, delta_t)) * denovo_normed(freq_0y, mcmc$p, data$filters) +
                (evolveJC(0, mcmc$mu, delta_t)) * denovo_normed(1 - freq_0y, mcmc$p, data$filters)
            )) +

            # Likelihood from x = 1, 0 < y < 1
            length(mcmc$m1y[[i]]) * log_p_isnv +
            sum(log(
              (evolveJC(1, mcmc$mu, delta_t))*denovo_normed(1 - freq_1y, mcmc$p, data$filters) +
                (evolveJC(0, mcmc$mu, delta_t))*denovo_normed(freq_1y, mcmc$p, data$filters)
            ))

          if(h <= data$n_obs){
            if(length(data$snvs[[h]]$isnv$call) > 0){

              # Frequency of shared iSNV in ancestor of case i
              freq_xy_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[h]]$isnv$call)]

              # If transmission occurred before end of exponential growth phase, need to back-mutate these frequencies.
              if(t_1st_trans < t_g){
                # To do this, we need to figure out what frequency they started at at t[h], on average
                start_freqs_xy <- rep(0, length(mcmc$mxy[[i]]))

                # Which of mcmc$mx0 are 1y in h?
                start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$m1y[[h]])] <- 1

                # If xy in h: start freq is 1/2 (approximation; can try improving this)
                start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$mxy[[h]])] <- 1/2

                # Proportion of exponential growth phase completed
                prop_exp_complete <- (t_1st_trans - mcmc$t[h]) / (t_g - mcmc$t[h])

                # Linearly interpolate frequencies
                freq_xy_anc <- freq_xy_anc * prop_exp_complete + start_freqs_xy * (1 - prop_exp_complete)
              }

              # Frequency of shared iSNV in case i
              freq_xy <- data$snvs[[i]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[i]]$isnv$call)]


              out <- out +
                # Likelihood from 0 < x < 1, 0 < y < 1
                sum(log(
                  (evolveJC(1, mcmc$mu, delta_t_prime)*(1 - freq_xy_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*freq_xy_anc) * exp(log_p_isnv) * denovo_normed(freq_xy, mcmc$p, data$filters) * (1 - p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)) +
                    (evolveJC(1, mcmc$mu, delta_t_prime)*freq_xy_anc + evolveJC(0, mcmc$mu, delta_t_prime)*(1 - freq_xy_anc)) * exp(log_p_isnv) * denovo_normed(1 - freq_xy, mcmc$p, data$filters) * (1 - p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)) +
                    p_all_split(mcmc$b, mcmc$w[i], freq_xy_anc)
                ))
            }
          }
        }
      }

      if(h <= data$n_obs){
        if(length(data$snvs[[h]]$isnv$call) > 0){

          # Frequency of iSNV in ancestor of case i
          freq_x0_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mx0[[i]], data$snvs[[h]]$isnv$call)]
          freq_x1_anc <- data$snvs[[h]]$isnv$af[match(mcmc$mx1[[i]], data$snvs[[h]]$isnv$call)]

          # If transmission occurred before end of exponential growth phase, need to back-mutate these frequencies.
          if(t_1st_trans < t_g){
            # To do this, we need to figure out what frequency they started at at t[h], on average
            start_freqs_x0 <- rep(0, length(mcmc$mx0[[i]]))

            # Which of mcmc$mx0 are 1y in h?
            start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$m1y[[h]])] <- 1

            # If xy in h: start freq is 1/2 (approximation; can try improving this)
            start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$mxy[[h]])] <- 1/2

            # Repeat for mcmc$mx1
            start_freqs_x1 <- rep(0, length(mcmc$mx1[[i]]))

            # Which of mcmc$mx0 are 1y in h?
            start_freqs_x1[which(mcmc$mx1[[i]] %in% mcmc$m1y[[h]])] <- 1

            # If xy in h: start freq is 1/2 (approximation; can try improving this)
            start_freqs_x1[which(mcmc$mx1[[i]] %in% mcmc$mxy[[h]])] <- 1/2

            # Proportion of exponential growth phase completed
            prop_exp_complete <- (t_1st_trans - mcmc$t[h]) / (t_g - mcmc$t[h])

            # Linearly interpolate frequencies
            freq_x0_anc <- freq_x0_anc * prop_exp_complete + start_freqs_x0 * (1 - prop_exp_complete)
            freq_x1_anc <- freq_x1_anc * prop_exp_complete + start_freqs_x1 * (1 - prop_exp_complete)

            if(any(freq_x1_anc < 0) | any(freq_x0_anc < 0)){
              print("warning2")
              print(mcmc$seq[[i]])
              print(t_1st_trans)
              print(mcmc$t[h])
              print(freq_x0_anc)
              print(freq_x1_anc)
            }
          }



          out <- out +
            # Likelihood from 0 < x < 1, y = 0
            length(mcmc$mx0[[i]]) * ifelse(isnv_info, log_p_no_isnv, 0) +
            sum(log(
              evolveJC(1, mcmc$mu, delta_t_prime)*(1 - freq_x0_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*freq_x0_anc
            )) + sum(log(1 - p_all_split(mcmc$b, mcmc$w[i], freq_x0_anc))) + # probability we don't transmit successive split bottlenecks

            # Likelihood from 0 < x < 1, y = 1
            length(mcmc$mx1[[i]]) * ifelse(isnv_info, log_p_no_isnv, 0) +
            sum(log(
              evolveJC(1, mcmc$mu, delta_t_prime)*(freq_x1_anc) + evolveJC(0, mcmc$mu, delta_t_prime)*(1 - freq_x1_anc)
            )) + sum(log(1 - p_all_split(mcmc$b, mcmc$w[i], freq_x1_anc))) # probability we don't transmit successive split bottlenecks

        }
      }

      if(delta_t_prime < 0){
        print("warning")
        print(delta_t)
        print(delta_t_prime)
        print(mcmc$seq[[i]])
        print(out)
      }



      return(out)
    }
  }
}









