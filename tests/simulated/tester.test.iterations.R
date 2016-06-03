source("../../src/node.dating.R")

extract_dates <- function(x) as.numeric(gsub("(.+)_([0123456789.]+)", "\\2", x, perl=T))
n.trees = 50

trees <- c(lapply(grep("HIV_[123456789]_rooted.tre", dir("."), value=T, perl=T), read.tree), lapply(grep("HIV_[123456789][1234567890]_rooted.tre", dir("."), value=T, perl=T), read.tree))

trees <- lapply(trees, function(t) {t$tip.date <- extract_dates(t$tip.label); t})
trees <- lapply(trees, function(t) {t$mu <- estimate.mu(t, t$tip.date); t})
trees <- lapply(trees, function(t) {mrcas <- mrca(t); o <- order(row.names(mrcas)); t$mrcas <- mrcas[o, o]; t})

mrcas <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca.csv"), header=T)[,2:101])
dates <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates.csv"), header=T)[,2:101])

o <- lapply(1:n.trees, function(n) order(names(mrcas[[n]])))
mrcas <- lapply(1:n.trees, function(n) mrcas[[n]][o[[n]], o[[n]]])
dates <- lapply(1:n.trees, function(n) dates[[n]][o[[n]], o[[n]]])

occ <- lapply(1:n.trees, function(n) table(unlist(mrcas[[n]])))

rmse <- function(mrca.test, dates.test, n, f) {
	occ.test <- table(unlist(mrca.test))
	
	w <- 0
	d <- 0
	
	for (i in 1:99) {
		for (j in (i+1):100) {
			weight <- f(occ[[n]][mrcas[[n]][i, j]], occ.test[mrca.test[i, j]])
			w <- w + weight
			d <- d + weight * (dates[[n]][i, j] - dates.test[mrca.test[i, j]])^2
		}
	}	
			
	as.numeric(sqrt(d / w))
}

sqrt.weight <- function(x, y) 1 / sqrt(x * y)

rmses <- matrix(0, nrow=n.trees, ncol=102)

cat("Init\n")

trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$tip.date, t$mu, opt.tol=1e-8, nsteps=0, lik.tol=0, is.binary=T); t})

rmses[, 1] <- unlist(lapply(1:n.trees, function(i) rmse(trees[[i]]$mrcas, trees[[i]]$node.date, i, sqrt.weight)))

cat("First\n")

node.date.1 <- lapply(trees, function(t) estimate.dates(t, t$node.date, t$mu, opt.tol=1e-8, nsteps=0, lik.tol=0, is.binary=T))

rmses[, 2] <- unlist(lapply(1:n.trees, function(i) rmse(trees[[i]]$mrcas, node.date.1[[i]], i, sqrt.weight)))

for (i in 1:100) {
	cat(as.character(i * 100))
	cat("\n")
	
	trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$node.date, t$mu, opt.tol=1e-8, nsteps=100, lik.tol=0, is.binary=T); t})
	
	rmses[, i + 2] <- unlist(lapply(1:n.trees, function(i) rmse(trees[[i]]$mrcas, trees[[i]]$node.date, i, sqrt.weight)))
}

