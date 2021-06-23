source("Simulated ML procedure.R")
source("MDT algorithm.R")

euclideanDistance <- function(x, y) {
  stopifnot(length(x) == length(y))
  sqrt(sum((x-y)^2))
}

maxAboluteDifference <- function(x, y) {
  max(abs(x - y))
}

runSingleExperiment <- function(alpha, weights.sample, n.har=1e2, n.steps=2) {
  
  n.subjects <- nrow(weights.sample)
  
  # Perform MDT
  constrSample <- vector("list",n.subjects)
  for (subject in 1:n.subjects) {
    constrSample[[subject]] <- performMDT(as.vector(weights.sample[subject,]), n.steps.after.ranking=n.steps)
  }
  
  # Perform simulated ML estimation
  results <- performSimulatedMLE(constrSample, n.har) 
  
  # Calculate estimated population mean and precision parameters
  true.precision <- sum(alpha)
  true.weights <- alpha/true.precision
  
  precision.SMLE <- sum(exp(results$simulatedMLE$estimate))
  weights.SMLE <- exp(results$simulatedMLE$estimate) / precision.SMLE
  
  precision.centroids <- sum(exp(results$centroidsMLE))
  weights.centroids <- exp(results$centroidsMLE) / precision.centroids
  
  alpha.Dirichlet.MLE <- MCMCprecision::fit_dirichlet(weights.sample)$alpha
  precision.Dirichlet.MLE <- sum(alpha.Dirichlet.MLE)
  weights.Dirichlet.MLE <- alpha.Dirichlet.MLE / precision.Dirichlet.MLE
  
  # Assess accuracy
  eucl.dist <- c(euclideanDistance(true.weights, weights.Dirichlet.MLE), euclideanDistance(true.weights, weights.SMLE), euclideanDistance(true.weights, weights.centroids))
  max.abs.diff <- c(maxAboluteDifference(true.weights, weights.Dirichlet.MLE), maxAboluteDifference(true.weights, weights.SMLE), maxAboluteDifference(true.weights, weights.centroids))
  precision.diff <- (c(precision.Dirichlet.MLE, precision.SMLE, precision.centroids) - true.precision)/true.precision
  
  # Return results
  data.frame(n.samples=rep(n.subjects, 3), method=c("exact", "SML", "centroids"), eucl.dist=eucl.dist, max.abs.diff=max.abs.diff, precision.diff=precision.diff)
  
}


PerformExperiments <- function(n.crit, n.subjects, n.questions, n.experiments=100, n.har=1e2) {
  results <- c()
  for (experiment in 1:n.experiments) {
    print(experiment)
    
    # Randomly generate parameters of the Dirichlet distribution; reject alpha if ratio highest to lowest weight is > 10
    alpha <- as.vector((runif(1)*5 + 5)*rdirichlet(n=1, alpha=rep(1, n.crit)))  
    while(max(alpha/sum(alpha))/min(alpha/sum(alpha))>10) {
      alpha <- as.vector((runif(1)*5 + 5)*rdirichlet(n=1, alpha=rep(1, n.crit)))  
    }

    for (cur.n.subjects in n.subjects) {
      weights.sample <- rdirichlet(n=cur.n.subjects, alpha=alpha) 
      for (cur.n.questions in n.questions) {
        cur.results <- runSingleExperiment(alpha, weights.sample, n.har, cur.n.questions/(n.crit-1))
        cur.results$experiment.id <- paste(n.crit, cur.n.subjects, cur.n.questions, sep="_")
        cur.results$n.subjects <- rep(cur.n.subjects, nrow(cur.results))
        cur.results$n.questions <- rep(cur.n.questions, nrow(cur.results))
        cur.results$run.id <- rep(experiment, nrow(cur.results))
        results <- rbind(results, cur.results)
      }
    }
  }
  results
}


### Conduct computational experiments and store results into an RData object ###

set.seed(1234) 
results.crit4 <- PerformExperiments(n.crit=4, n.subjects=c(100, 300), n.questions=c(3, 6, 9), n.experiments=100, n.har=1e2)
save(results.crit4, file="./compExperiments/compExp_4_100.RData")

set.seed(14321)
results.crit6 <- PerformExperiments(n.crit=6, n.subjects=c(100, 300), n.questions=c(5, 10, 15), n.experiments=100, n.har=1e2)
save(results.crit6, file="./compExperiments/compExp_6_100.RData")
