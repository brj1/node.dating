source("../../src/node.dating.R")
source("../../src/extra.R")

extract_dates <- function(x) as.numeric(gsub("(.+)_([0123456789.]+)", "\\2", x, perl=T))
n.trees = 50

cat("reading trees\n")

trees <- c(lapply(grep("HIV_[123456789]_rooted.tre", dir("."), value=T, perl=T), read.tree), lapply(grep("HIV_[123456789][1234567890]_rooted.tre", dir("."), value=T, perl=T), read.tree))

trees <- lapply(trees, function(t) {t$tip.date <- extract_dates(t$tip.label); t})
trees <- lapply(trees, function(t) {t$mu <- estimate.mu(t, t$tip.date); t})
trees <- lapply(trees, function(t) {mrcas <- mrca(t); o <- order(row.names(mrcas)); t$mrcas <- mrcas[o, o]; t})
trees <- lapply(trees, function(t) {t$occ <- table(unlist(t$mrcas)); t})

cat("reading true data\n")

mrcas <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca.csv"), header=T)[,2:101])
dates <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates.csv"), header=T)[,2:101])

o <- lapply(1:n.trees, function(n) order(names(mrcas[[n]])))
mrcas <- lapply(1:n.trees, function(n) mrcas[[n]][o[[n]], o[[n]]])
dates <- lapply(1:n.trees, function(n) dates[[n]][o[[n]], o[[n]]])

occ <- lapply(1:n.trees, function(n) table(unlist(mrcas[[n]])))

cat("computing rmse weights\n")

sqrt.weight <- function(x, y) 1 / sqrt(x * y)

rmse.weight <- lapply(1:n.trees, function(n) matrix(unlist(lapply(1:100, function(i) lapply(1:100, function(j) if (i == j) 0 else sqrt.weight(occ[[n]][mrcas[[n]][i, j]], trees[[n]]$occ[trees[[n]]$mrcas[i, j]])))), nrow=100, ncol=100))

rmse <- function(n) {
	w <- 0
	d <- 0
	 	
	d <- rmse.weight[[n]] * (dates[[n]] - trees[[n]]$node.date[trees[[n]]$mrcas])^2
				
	as.numeric(sqrt(sum(d) / sum(rmse.weight[[n]])))
}

#rmses <- matrix(0, nrow=n.trees, ncol=102)
#lik <- matrix(0, nrow=n.trees, ncol=102)
rmses <- matrix(0, nrow=n.trees, ncol=1001)
lik <- matrix(0, nrow=n.trees, ncol=1001)

cat("init\n")

trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$tip.date, t$mu, opt.tol=1e-8, nsteps=0, lik.tol=0, is.binary=T); t})

rmses[, 1] <- unlist(lapply(1:n.trees, rmse))
lik[, 1] <- unlist(lapply(trees, function(t) tree.like(t, t$node.date, t$mu)))

#cat("First\n")

#node.date.1 <- lapply(trees, function(t) estimate.dates(t, t$node.date, t$mu, opt.tol=1e-8, nsteps=0, lik.tol=0, is.binary=T))

#rmses[, 2] <- unlist(lapply(1:n.trees, function(i) rmse(trees[[i]]$mrcas, node.date.1[[i]], i, sqrt.weight)))
#lik[, 2] <- unlist(lapply(1:n.trees, function(i) tree.like(trees[[i]], node.date.1[[i]], trees[[i]]$mu)))

#for (i in 1:100) {
for (i in 1:1000) {
#	cat(as.character(i * 100))
	cat(as.character(i))
	cat("\n")
	
#	trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$node.date, t$mu, opt.tol=1e-8, nsteps=100, lik.tol=0, is.binary=T); t})
	trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$node.date, t$mu, opt.tol=1e-8, nsteps=0, lik.tol=0, is.binary=T); t})
	
#	rmses[, i + 2] <- unlist(lapply(1:n.trees, function(j) rmse(trees[[j]]$mrcas, trees[[j]]$node.date, j, sqrt.weight)))
#	lik[, i + 2] <- unlist(lapply(trees, function(t) tree.like(t, t$node.date, t$mu)))
	rmses[, i + 1] <- unlist(lapply(1:n.trees, rmse))
	lik[, i + 1] <- unlist(lapply(trees, function(t) tree.like(t, t$node.date, t$mu)))
}

pdf("iter.pdf", width=8, height=8)
boxplot(rmses[,seq(0, 100) + 1], xlab="Steps", ylab="RMSE (days)", main="RMSE with different numbers of iterations", xaxt='n')
axis(1, seq(0, 10)*10 + 1, seq(0, 10)*10)
dev.off()