library("hitandrun")

performMDT <- function(weights, n.steps.after.ranking=2) {
  n.crit <- length(weights)
  constr <- simplexConstraints(n.crit)
  ranking <- order(weights, decreasing=T)
  for (index in 1:(n.crit-1)) {
    constr <- mergeConstraints(constr, performBisection(weights, ranking[index], ranking[index+1], n.steps.after.ranking))
  }
  constr
}

performBisection <- function(weights, index.i, index.j, n.steps) {
  n.crit <- length(weights)
  cutoffs <- 1/seq(0, 1, 1/(2^n.steps))
  weight.ratio <- weights[index.i] / weights[index.j]
  upper.bound <- cutoffs[max(which(cutoffs>=weight.ratio))]
  lower.bound <- cutoffs[min(which(cutoffs<=weight.ratio))]
  if (upper.bound!=Inf) {
    constr <- mergeConstraints(upperRatioConstraint(n.crit, index.i, index.j, upper.bound), lowerRatioConstraint(n.crit, index.i, index.j, lower.bound))
  } else {
    constr <- lowerRatioConstraint(n.crit, index.i, index.j, lower.bound)
  }
  constr
}
