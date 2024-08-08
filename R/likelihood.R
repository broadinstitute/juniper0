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

    # Ancestor of host i
    h <- mcmc$h[i]

    # Time of end of exponential growth phase in h
    t_g <- mcmc$t[h] + log(1/sqrt(mcmc$p)) / (mcmc$mu / mcmc$p) / log(mcmc$v)

    # Evolutionary time intervals for each transmission on the chain from h to i, bottleneck to bottleneck
    # Doesn't include h -> last element of seq[[i]]
    # Oldest to newest, hence rev()
    deltas <- diff(rev(mcmc$seq[[i]]))
    t_1st_trans <- mcmc$seq[[i]][length(mcmc$seq[[i]])]

    # if(t_1st_trans < mcmc$t[h]){
    #   return(-Inf)
    # }

    delta_init <- t_1st_trans - t_g # Initial end of expo phase to first bottleneck; may be negative

    # probability of SPECIFIC iSNV in expo growth phase
    p_new_isnv <- (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) / 3
    log_p_new_isnv <- log(p_new_isnv)

    # log probability of no iSNV in expo growth phase
    log_p_no_isnv <- (1/sqrt(mcmc$p)) * log(
      (1 - mcmc$p)
    )

    # Number of sites without a mutation
    n_no_mut <- data$n_bases -
      length(mcmc$m01[[i]]) -
      length(mcmc$m10[[i]]) -
      length(mcmc$m0y[[i]]) -
      length(mcmc$m1y[[i]]) -
      length(mcmc$mx0[[i]]) -
      length(mcmc$mx1[[i]]) -
      length(mcmc$mxy[[i]])

    # Bottleneck infecting first host on the chain. Format: (P(absent), P(iSNV), P(present))
    # Since it's unclear whether we're going 0 to 0 or 1 to 1 for the non-mutant sites, assume 0 -> 0

    # bot_absent is the probability of each bottleneck infecting i, starting with the SNV absent in h
    # bot_present is the probability of each bottleneck infecting i, starting with the SNV present in h
    bot_absent <- c(1, 0, 0)
    bot_present <- c(0, 0, 1)
    # Evolve over first bottleneck
    if(delta_init > 0){
      bot_absent <- evolve_bot(bot_absent, mcmc$mu, mcmc$b, delta_init)
      bot_present <- evolve_bot(bot_present, mcmc$mu, mcmc$b, delta_init)
    }
    # Evolve over remaining bottlenecks
    bot_absent <- evolve_bot_repeatedly(bot_absent, mcmc$mu, mcmc$b, deltas)
    bot_present <- evolve_bot_repeatedly(bot_present, mcmc$mu, mcmc$b, deltas)

    out <- 0

    # Likelihood from x = 0, y = 0 or x = 1, y = 1
    if(n_no_mut > 0){
      out <- out + n_no_mut * (log(bot_absent[1]) + log_p_no_isnv)
    }

    # Likelihood from x = 0, y = 1
    if(length(mcmc$m01[[i]]) > 0){
      out <- out + length(mcmc$m01[[i]]) * (log(bot_absent[3]) + log_p_no_isnv)
    }

    # Likelihood from x = 1, y = 0
    if(length(mcmc$m10[[i]]) > 0){
      out <- out + length(mcmc$m10[[i]]) * (log(bot_present[1]) + log_p_no_isnv)
    }

    # Frequencies of added iSNVs
    freq_0y <- mcmc$isnv$af[[i]][match(mcmc$m0y[[i]], mcmc$isnv$call[[i]])]

    # Frequency of shared iSNV in case i
    freq_xy <- mcmc$isnv$af[[i]][match(mcmc$mxy[[i]], mcmc$isnv$call[[i]])]

    # Frequencies of deleted iSNVs
    freq_1y <- mcmc$isnv$af[[i]][match(mcmc$m1y[[i]], mcmc$isnv$call[[i]])]

    out <- out +
      # Likelihood from x = 0, 0 < y < 1
      sum(
        log(
          bot_absent[1] * p_new_isnv * denovo(freq_0y, mcmc$p) +
            bot_absent[2] +
            bot_absent[3] * p_new_isnv * denovo(1 - freq_0y, mcmc$p)
        )
      ) +

      # Likelihood from x = 1, 0 < y < 1
      sum(
        log(
          bot_present[1] * p_new_isnv * denovo(freq_1y, mcmc$p) +
            bot_present[2] +
            bot_present[3] * p_new_isnv * denovo(1 - freq_1y, mcmc$p)
        )
      )

    # Frequency of iSNV in ancestor of case i
    freq_x0_anc <- mcmc$isnv$af[[h]][match(mcmc$mx0[[i]], mcmc$isnv$call[[h]])]
    freq_xy_anc <- mcmc$isnv$af[[h]][match(mcmc$mxy[[i]], mcmc$isnv$call[[h]])]
    freq_x1_anc <- mcmc$isnv$af[[h]][match(mcmc$mx1[[i]], mcmc$isnv$call[[h]])]

    #print(freq_x1_anc)
    #print(t_1st_trans)
    #print(mcmc$t[h])

    # If transmission occurred before end of exponential growth phase, need to back-mutate these frequencies.
    if(delta_init < 0){
      # To do this, we need to figure out what frequency they started at at t[h], on average
      start_freqs_x0 <- rep(0, length(mcmc$mx0[[i]]))

      # Which of mcmc$mx0 are 1y in h?
      start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$m1y[[h]])] <- 1

      # If xy in h: start freq is 1/2 (approximation; can try improving this)
      start_freqs_x0[which(mcmc$mx0[[i]] %in% mcmc$mxy[[h]])] <- 1/2


      # Repeat for mcmc$mxy
      start_freqs_xy <- rep(0, length(mcmc$mxy[[i]]))

      # Which of mcmc$mx0 are 1y in h?
      start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$m1y[[h]])] <- 1

      # If xy in h: start freq is 1/2 (approximation; can try improving this)
      start_freqs_xy[which(mcmc$mxy[[i]] %in% mcmc$mxy[[h]])] <- 1/2


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
      freq_xy_anc <- freq_xy_anc * prop_exp_complete + start_freqs_xy * (1 - prop_exp_complete)
      freq_x1_anc <- freq_x1_anc * prop_exp_complete + start_freqs_x1 * (1 - prop_exp_complete)
    }

    # Bottlenecks for each of the three cases
    bot_x0 <- mapply(init_bot, 1 - freq_x0_anc, freq_x0_anc, MoreArgs = list(mu = mcmc$mu, b = mcmc$b, delta_t = delta_init), SIMPLIFY = F)
    bot_xy <- mapply(init_bot, 1 - freq_xy_anc, freq_xy_anc, MoreArgs = list(mu = mcmc$mu, b = mcmc$b, delta_t = delta_init), SIMPLIFY = F)
    bot_x1 <- mapply(init_bot, 1 - freq_x1_anc, freq_x1_anc, MoreArgs = list(mu = mcmc$mu, b = mcmc$b, delta_t = delta_init), SIMPLIFY = F)

    #print(bot_x0)
    #print(bot_xy)
    #print(freq_x1_anc)
    #print(bot_x1)

    bot_x0 <- lapply(bot_x0, evolve_bot_repeatedly, mu = mcmc$mu, b = mcmc$b, deltas = deltas)
    bot_xy <- lapply(bot_xy, evolve_bot_repeatedly, mu = mcmc$mu, b = mcmc$b, deltas = deltas)
    bot_x1 <- lapply(bot_x1, evolve_bot_repeatedly, mu = mcmc$mu, b = mcmc$b, deltas = deltas)



    # Log probabilities for each observation in each case
    p_x0 <- unlist(lapply(bot_x0, log_p_given_bot, freq = 0, p = mcmc$p, p_new_isnv = p_new_isnv, log_p_no_isnv = log_p_no_isnv))
    p_xy <- unlist(mapply(log_p_given_bot, bot_xy, freq_xy, MoreArgs = list(p = mcmc$p, p_new_isnv = p_new_isnv, log_p_no_isnv = log_p_no_isnv), SIMPLIFY = F))
    p_x1 <- unlist(lapply(bot_x1, log_p_given_bot, freq = 1, p = mcmc$p, p_new_isnv = p_new_isnv, log_p_no_isnv = log_p_no_isnv))

    out <- out + sum(p_x0) + sum(p_xy) + sum(p_x1)

    return(out)

  }
}









