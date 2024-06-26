% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_network}
\alias{plot_network}
\title{Plot Current Network Configuration}
\usage{
plot_network(init, type = "standard")
}
\arguments{
\item{init}{State of the Markov chain, in the same format as generated by initialize().}

\item{type}{Either "standard" (default), which returns a standard phylogenetic tree layout, or "radial", which returns a radial tree layout.}
}
\value{
The network plot.
}
\description{
This function plots the current network topology as the Markov chain runs.
}
