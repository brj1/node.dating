library(ape)
library(data.table)

tree.like <- function(tree, node.dates, mu) {
	if (mu < 0)
		return(-Inf)
	
	scale.lik <- sum(-lgamma(tree$edge.length+1)+(tree$edge.length+1)*log(mu))
		
	calc.lik <- function(ch.node, edge, par.node) {
		tim <- ch.node - par.node
				
		edge*log(tim)-mu*tim
	}
			
	sum(calc.lik(node.dates[tree$edge[,2]], tree$edge.length, node.dates[tree$edge[,1]])) + scale.lik
}

cat("generating trees...\n")
trees <- lapply(unlist(lapply(4^seq(1, 5), rep, 10)), rtree)

sds <- .1^seq(1, 5)

rmse <- function(x, y) sqrt(sum((x-y)^2)/length(x))

clock <- function(tree, rate, noise) rnorm(tree$edge.length, mean=tree$edge.length * rate, sd=noise)

write.csv(data.frame(tips=c(), noise=c(), steps=c(), likelihood=c(), rmse=c(), time=c()), file="trial.csv")

supress <- lapply(sds, function(s) {
		cat(paste0("\nnoise: ", s,"\n"))

		tab <- rbindlist(lapply(trees, function(tree) {
				n <- length(tree$tip.label)
				
				cat(paste0(n, " "))
				
				s.tree <- tree
				s.tree$edge.length <- clock(tree, 1, s)
				tree$time <- node.depth.edgelength(s.tree)
				tree$mu <- estimate.mu(tree, tree$time[1:n])
				
				times <- rep(0, 11)
				liks <- rep(0, 11)
				rmses <- rep(0, 11)
				waste.time <- 0
				
				start.time <- Sys.time()
				dates <- estimate.dates(tree, tree$time[1:n], mu=tree$mu, opt.tol=1e-8, lik.tol=0, nsteps=0, is.binary=T)
				times[1] <- Sys.time() - start.time
				
				rmses[1] <- rmse(dates, tree$time)
				liks[1] <- tree.like(tree, dates, tree$mu)				
												
				for (i in 1:10) {
					waste.time <- (Sys.time() - start.time) - times[i]
					dates <- estimate.dates(tree, dates, mu=tree$mu, opt.tol=1e-8, lik.tol=0, nsteps=100, is.binary=T)
					times[i + 1] <- (Sys.time() - start.time) - waste.time
					rmses[i + 1] <- rmse(dates, tree$time)
					liks[i + 1] <- tree.like(tree, dates, tree$mu)
				}
				
				data.frame(tips=rep(n, 11), noise=rep(s, 11), steps=100*(0:10), likelihood=liks, rmse=rmses, time=times)
			}))
			
		write.csv(tab, file="trial.csv", col.names=F, row.names=F, append=T)
	})