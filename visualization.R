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

### Visualize MCMC output

library(ggplot2)
library(cowplot)
library(ggraph)
library(igraph)
library(ape)
library(parallel)

source("subroutines.R")

# Load in results
load("output.RData")

mus <- c()
ps <- c()
bs <- c()
vs <- c()
for (i in 1001:10000) {
  mus[i] <- output[[i]]$mu
  ps[i] <- output[[i]]$p
  bs[i] <- output[[i]]$b
  vs[i] <- output[[i]]$v
}
mus <- mus[!is.na(mus)]
ps <- ps[!is.na(ps)]
bs <- bs[!is.na(bs)]

h <- output[[10000]]$h
t <- output[[10000]]$t
w <- output[[10000]]$w
n <- length(h)


## Lineage by case
# Read in sequence names and metadata
cons <- read.FASTA("~/Desktop/input_data_huge/aligned.fasta")
names <- names(cons)
meta <- read.csv("~/Desktop/input_data_huge/metadata.csv")
lin <- meta$Sequence.Nextclade.Clade[match(names, meta$case)]
lineage <- rep("Unsampled", n)
lineage[2:(length(lin) + 1)] <- lin

n_obs <- length(cons) + 1

# Data cleaning
lineage[lineage == "20C"] <- "20C (Epsilon)"
lineage[lineage == "21L"] <- "21L (Omicron)"
lineage[lineage == "22A"] <- "22A (Omicron)"
lineage[lineage == "22B"] <- "22B (Omicron)"
lineage[lineage == "22D"] <- "22D (Omicron)"
lineage[lineage == "22E"] <- "22E (Omicron)"
lineage[lineage == "22F"] <- "22F (Omicron)"
lineage[lineage == "23A"] <- "23A (Omicron)"


# Reverse BFS order
bfs2 <- function(h){
  queue <- 1
  explored <- c()
  while (length(queue) > 0) {
    queue <- c(queue, which(h == queue[1]))
    explored <- c(explored, queue[1])
    queue <- queue[-1]
    #print(length(explored))
  }
  return(explored)
}

# Impute lineages of unobserved nodes
for (i in rev(bfs(1, h))) {
  if(lineage[i] == "Unsampled"){
    kids <- which(h == i)
    lineage[i] <- names(which.max(table(lineage[kids])))
  }
}

# BFS order, expanding nodes by lineage
bfs_lineage <- function(h){
  queue <- 1
  explored <- c()
  while (length(queue) > 0) {
    who <- which(h == queue[1])
    who <- who[sort.int(lineage[who], index.return = T)$ix]
    queue <- c(queue, which(h == queue[1]))
    explored <- c(explored, queue[1])
    queue <- queue[-1]
    #print(length(explored))
  }
  return(explored)
}

# DFS order, expanding nodes by lineage
dfs_lineage <- function(h){
  stack <- 1
  explored <- c()
  while (length(stack) > 0) {
    who <- which(h == stack[1])
    who <- who[sort.int(lineage[who], index.return = T)$ix]
    explored <- c(explored, stack[1])
    stack <- stack[-1]
    stack <- c(who, stack)

    print(length(explored))
  }
  return(explored)
}



ord <- rev(dfs_lineage(h))

# Leaves
leaves <- which(!(1:n %in% h))
n_leaves <- length(leaves)

# Lineage of each leaf
#leaf_lineage <- lineage[leaves]




# Angle of each node
thetas <- c()
#thetas[leaves[sort.int(leaf_lineage, index.return = T)$ix]] <- (1:n_leaves) * 2* pi / n_leaves

leaf_count <- 0

for (i in ord) {
  if(!(i %in% leaves)){
    kids <- which(h == i)
    thetas[i] <- mean(thetas[kids])
  }else{
    thetas[i] <- leaf_count * 2 * pi / n_leaves
    leaf_count <- leaf_count + 1
  }
}

# Reset unsampled lineages to unsampled
lineage[(length(lin) + 2):length(lineage)] <- "Unsampled"
lineage[1] <- "Unsampled"
lineage[lineage == "recombinant"] <- "Recombinant"

offset <- 550
rotation <- 0.68 #radians
df <- data.frame(x = ((t-offset)*cos(thetas + rotation))[t>offset], y = ((t-offset)*sin(thetas + rotation))[t>offset], lineage = lineage[t>offset])
df_standard <- data.frame(x = t, y = thetas, lineage = lineage)

colors <- c(
  '#482a43',
  '#baa0b5',
  '#37495e',
  '#2b24a9',
  '#41481a',
  '#94558a',
  '#442b33',
  '#5a6a0a',
  '#ff7566',
  '#640308',
  '#ab7971',
  '#00607f',
  '#404b89',
  '#ddb4e1',
  '#000000',
  '#cccccc'
)

ggplot(df, aes(x=x,y=y,color=lineage)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = colors, name = NULL) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.direction = "horizontal",
    legend.position=c(0.45, -0.05),
    legend.text = element_text(size=6),
    legend.key.width= unit(0.1, 'cm'),
    legend.key.height= unit(0.3, 'cm'),
    plot.margin = margin(0,-1.4,2,0, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1)))

ggsave("./figs/radial.pdf", width=6.5, height = 8)
ggsave("./figs/radial.png", width=6.5, height = 8)



# vertical segments
xs <- c()
ystart <- c()
yend <- c()
for (i in 1:n) {
  kids <- which(h == i)
  if(length(kids) > 0){
    xs <- c(xs, t[i])
    ystart <- c(ystart, min(thetas[kids]))
    yend <- c(yend, max(thetas[kids]))
  }
}

big <- ggplot() +
  geom_segment(mapping = aes(x = t[h], xend = t, y = thetas, yend = thetas), linewidth = 0.05) +
  geom_segment(mapping = aes(x = xs, xend = xs, y = ystart, yend = yend), linewidth = 0.05) +
  geom_point(mapping= aes(x = df_standard$x, y = df_standard$y, color = df_standard$lineage), alpha = 0.5, size = 0.01) +
  scale_color_manual(values = colors, name = "Lineage") +
  xlab("Evolutionary Time (days)") +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    legend.position=c(0.85, 0.65)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1)))

big

ggsave("./figs/mass_network.pdf", width = 6.5, height = 8, limitsize = F)
ggsave("./figs/mass_network.png", width = 6.5, height = 8, limitsize = F)


mu_hist <- ggplot(data.frame(x = mus), aes(x = x)) +
  geom_histogram(aes(y=after_stat(density)), bins = 10, color = "white", fill = "grey") +
  xlab(paste("Value of mu")) +
  ylab("Probability Density") +
  #xlim(2.1e-6, 2.4e-6) +
  scale_x_continuous(breaks = seq(min(mus), max(mus), by = (max(mus) - min(mus)) /2), labels = function(x) format(x, scientific = T, digits = 3)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0.2, r = 0.3, b = 0, l = 0, unit = "in"))
p_hist <- ggplot(data.frame(x = ps), aes(x = x)) +
  geom_histogram(aes(y=after_stat(density)), bins = 10, color = "white", fill = "grey") +
  xlab(paste("Value of p")) +
  ylab("Probability Density") +
  scale_y_continuous(labels = function(x) format(x, scientific = T, digits = 2)) +
  scale_x_continuous(breaks = seq(min(ps), max(ps), by = (max(ps) - min(ps)) / 2), labels = function(x) format(x, scientific = T, digits = 3)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0.2, r = 0.3, b = 0, l = 0, unit = "in"))
b_hist <- ggplot(data.frame(x = bs), aes(x = x)) +
  geom_histogram(aes(y=after_stat(density)), breaks = seq(0,0.01, 0.001), color = "white", fill = "grey") +
  xlab(paste("Value of b")) +
  ylab("Probability Density") +
  scale_x_continuous(breaks = seq(0, 0.01, by = 0.01/2), labels = function(x) format(x, scientific = T, digits = 2)) +
  theme_minimal() +
  theme(plot.margin = margin(t = 0.2, r = 0.3, b = 0, l = 0, unit = "in"))

hists <- plot_grid(
  mu_hist, p_hist, b_hist,
  ncol = 3,
  labels = "AUTO"
)

ggsave("./figs/mass_hists.pdf", width = 8, height = 3, limitsize = F)
ggsave("./figs/mass_hists.png", width = 8, height = 3, limitsize = F)

## Next figure: genetic distance between case and ancestor
dists <- c()
times <- c()
gens <- c()
for (i in 2:n) {
  if(h[i] <= n_obs & i <= n_obs){
    dists[i] <- dist.dna(cons[c(i-1, h[i] - 1)], "N")[1]
  }
  times[i] <- t[i] - t[h[i]]
  gens[i] <- w[i] + 1

  if(i %% 100 == 0){
    print(i)
  }
}

hist(dists[w == 0])

# Is facility of ancestor same?
anc_fac <- c()
for (i in 2:n) {
  if(h[i] <= n_obs & i <= n_obs){
    anc_fac[i] <- meta$Practice.Name[match(names[i-1], meta$case)] == meta$Practice.Name[match(names[h[i]-1], meta$case)]
  }
}

mean(anc_fac, na.rm = T)

# Okay, and what if we randomly assigned each person someone of the same lineage who was sampled 0-10 days before?
# Facility by case
facs <- meta$Practice.Name[match(names, meta$case)]
lins <- meta$pango_lineage_full[match(names, meta$case)]
dats <- as.numeric(difftime(
  as.Date(meta$Collection.Date[match(names, meta$case)]),
  as.Date("2020-01-01"),
  units = "days"
))

get_new_h <- function(xx){
  new_h <- c()
  for (i in 1:length(facs)) {
    choices <- which(dats <= dats[i] & dats >= dats[i] - 10 & lins == lins[i])
    choices <- setdiff(choices, i)
    if(length(choices) == 1){
      new_h[i] <-choices
    }else if(length(choices) > 1){
      new_h[i] <- sample(choices, 1)
    }else{
      new_h[i] <- NA
    }
  }
  return(mean(facs == facs[new_h], na.rm = T))
}


rand <- mclapply(1:10, get_new_h, mc.cores = 10)
mean(unlist(rand))
range(unlist(rand))







# table(lineage[h[which(lineage == "22E")]])
# weird <- which(lineage == "22E" & lineage[h] == "22B (Omicron)")
#
# i <- h[weird[3]]
# j <- weird[3]
#
# mini_cons <- cons[c(i-1, j-1)]
# dist.dna(mini_cons, "N")
#
# mcmc <- output[[10000]]

