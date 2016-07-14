library(ape)
library(data.table)

cat("generating trees...")
trees <- lapply(unlist(lapply(4^seq(1, 6), rep, 10)), rtree)

sds <- .1^seq(1, 5)

rmse <- function(x, y) sqrt(sum((x-y)^2)/length(x))

clock <- function(tree, rate, noise) rnorm(tree$edge.length, mean=tree$edge.length * rate, sd=noise)

tab <- rbindlist(lapply(sds, function(s) {
		cat(paste0("noise: ", s,"\n"))

		rbindlist(lapply(trees, function(tree) {
				n <- length(tree$tip.label)
				s.tree <- tree
				s.tree$edge.length <- clock(tree, 1, s)
				tree$time <- node.depth.edgelength(s.tree)
				tree$mu <- estimate.mu(tree, tree$time[1:n])
				start.time <- Sys.time()
				times <- rep(0, 11)
				dates <- estimate.dates(tree, tree$time[1:n], mu=tree$mu, opt.tol=1e-8, lik.tol=0, nsteps=0)
				rmses <- c(rmse(dates, tree$time), rep(0, 10))
				times[1] <- Sys.time() - start.time
								
				for (i in 1:10) {
					dates <- estimate.dates(tree, dates, mu=tree$mu, opt.tol=1e-8, lik.tol=1000, nsteps=0)
					rmses[i + 1] <- rmse(dates, tree$time)
					times[i + 1] <- Sys.time() - start.time
				}
				
				data.frame(tips=rep(n, 11), noise=rep(s, 11), steps=1000*(0:10), rmse=rmses, time=times)
			}))
	}))

cat("write...")
write.csv(tab, file="trial.csv")