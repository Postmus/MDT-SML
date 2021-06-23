library("hitandrun")
library("DirichletReg")
library("plyr")
library("abind")
library("maxLik")

### Functions for constructing the simulated log-likelihood function for Dirichlet models without covariates ###

getHarSamplesSubject <- function(constr, n.har=1e3) {
  samples <- hitandrun(constr, n.har)
  n.crit <- ncol(constr$constr)
  array(samples, dim=c(1, n.har, n.crit))
}

getSMLinputs <- function(constrSample, n.har=1e3) {
  har.samples <- c()
  for (constr in constrSample) {
    har.samples <- abind(har.samples, getHarSamplesSubject(constr, n.har), along=1)
  }
  har.samples
}

genLogLikelihood <- function(har.samples) {
  n.patients <- dim(har.samples)[1]
  function(log.alpha) {
    alpha <- exp(log.alpha)
    logLik <- 0
    for (i in 1:n.patients) {
      logLik <- logLik + log(mean(ddirichlet(har.samples[i,,], alpha)))
    }
    logLik
  }
}

performSimulatedMLE <- function(constrSample, n.har=1e3) {
  har.samples <- getSMLinputs(constrSample, n.har)
  logLikelihood <- genLogLikelihood(har.samples)
  centroids <- aaply(har.samples, 1, colMeans)
  starting.values <- log(MCMCprecision::fit_dirichlet(centroids)$alpha)
  list(centroidsMLE=starting.values, simulatedMLE=maxLik(logLik=logLikelihood, start=starting.values, method="BFGS"))
}

