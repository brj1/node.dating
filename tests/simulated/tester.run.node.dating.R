options(digits.secs = 4)
st <- Sys.time()

source("../../src/node.dating.R")

args <- commandArgs(trailing=T)

nsteps <- as.numeric(args[1])

extract_dates <- function(x) as.numeric(gsub("(.+)_([0123456789.]+)", "\\2", x, perl=T))

n.trees <- 50

trees <- lapply(1:n.trees, function(i) read.tree(paste0("HIV_", i, "/HIV_", i,"_rooted.tre")))

trees <- lapply(trees, function(t) {t$tip.date <- extract_dates(t$tip.label); t})
trees <- lapply(trees, function(t) {t$mu <- estimate.mu(t, t$tip.date); t})
trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$tip.date, t$mu, opt.tol=1e-8, nsteps=nsteps, show.steps=0, lik.tol=0, is.binary=T); t})

Sys.time() - st

for (i in 1:n.trees) {
	m <- mrca(trees[[i]])
	dates <- trees[[i]]$node.date
	
	write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_mrca_node.dating.", nsteps, ".csv"))
	
	m <- apply(m, c(1, 2), function(n) dates[n])
	
	write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_dates_node.dating.", nsteps, ".csv"))
}

Sys.time() - st