library("ggplot2")
library("gridExtra")
library("ggpubr")
library("plyr")

generateBoxPlot <- function(results, measure="eucl.dist", plot.title=NULL) {
  results$n.questions <- as.factor(results$n.questions)
  results$group <- interaction(results$method, results$n.questions)
  results$method <- as.factor(results$method)
  lvls <- levels(results$method)
  ## redo levels
  lvls = c('Centroid', 'Exact', 'SML')
  levels(results$method) <- lvls


  if (measure=="eucl.dist") {
      p <- ggplot(data=results,
                  mapping=aes(x=reorder(method, eucl.dist, FUN=function(x){-mean(x)}, order=TRUE),
                              y=eucl.dist, fill=n.questions, group=group))
  } 
  if (measure=="max.abs.dist") {
      p <- ggplot(data=results,
                  mapping=aes(x=reorder(method, max.abs.diff, FUN=function(x){-mean(x)}, order=TRUE),
                              y=max.abs.diff, fill=n.questions, group=group))
  }
  if (measure=="precision") {
    p <- ggplot(data=results,
                mapping=aes(x=reorder(method, precision.diff, FUN=function(x){-mean(x)}, order=TRUE),
                            y= precision.diff, fill=n.questions, group=group))
  }

  p <- p + geom_boxplot()
  p <- p + labs(x=NULL, fill="Number of questions")
  if (!is.null(plot.title)) {
    p <- p + labs(title=plot.title)
  }
  if (measure=="eucl.dist") {
    p <- p + labs(y="Euclidian distance")
  } 
  if (measure=="max.abs.dist") {
    p <- p + labs(y="Maximum absolute difference")
  }
  if (measure=="precision") {
    p <- p + labs(y="Relative difference")
  }
  p <- p + theme(plot.title = element_text(hjust=0.5), legend.position="bottom")
  p
}

### Load simulation resutls and produce plots ###

load(file="./compExperiments/compExp_4_100.RData")
load(file="./compExperiments/compExp_6_100.RData")

# Euclidean distance
p100.4 <- generateBoxPlot(results.crit4[results.crit4$n.samples==100,], measure="eucl.dist",  plot.title="4 attributes, n=100")
p300.4 <- generateBoxPlot(results.crit4[results.crit4$n.samples==300,], measure="eucl.dist",  plot.title="4 attributes, n=300")
p100.6 <- generateBoxPlot(results.crit6[results.crit6$n.samples==100,], measure="eucl.dist",  plot.title="6 attributes, n=100")
p300.6 <- generateBoxPlot(results.crit6[results.crit6$n.samples==300,], measure="eucl.dist",  plot.title="6 attributes, n=300")

pdf("./Figures/results_euclDist.pdf", width=12, height=8, onefile=FALSE)
print(ggarrange(p100.4, p300.4, p100.6, p300.6, nrow=2, ncol=2, legend='bottom'))
dev.off()

# Precision
pre100.4 <- generateBoxPlot(results.crit4[results.crit4$n.samples==100,], measure="precision",  plot.title="4 attributes, n=100")
pre300.4 <- generateBoxPlot(results.crit4[results.crit4$n.samples==300,], measure="precision",  plot.title="4 attributes, n=300")
pre100.6 <- generateBoxPlot(results.crit6[results.crit6$n.samples==100,], measure="precision",  plot.title="6 attributes, n=100")
pre300.6 <- generateBoxPlot(results.crit6[results.crit6$n.samples==300,], measure="precision",  plot.title="6 attributes, n=300")

pdf("./Figures/results_precision.pdf", width=12, height=8, onefile=FALSE)
print(ggarrange(pre100.4, pre300.4, pre100.6, pre300.6, nrow=2, ncol=2, legend='bottom'))
dev.off()
