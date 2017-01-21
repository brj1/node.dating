library(ape)

n.trees <- 50
n.tips <- 100

mrcas <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "/HIV_", n, "_mrca.csv"), header=T)[,2:101])
dates <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "/HIV_", n, "_dates.csv"), header=T)[,2:101])
mrcas.node.dating <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "/HIV_", n, "_mrca_node.dating.csv"), header=T)[,2:101])

o <- lapply(1:n.trees, function(n) match(names(mrcas[[n]]), gsub("(.+)_.+", "\\1", names(mrcas.node.dating[[n]]))))
mrcas.node.dating <- lapply(1:n.trees, function(n) mrcas.node.dating[[n]][o[[n]], o[[n]]])

occ <- lapply(1:n.trees, function(n) table(unlist(mrcas[[n]])))
occ.node.dating <- lapply(1:n.trees, function(n) table(unlist(mrcas.node.dating[[n]])))

sqrt.weight <- function(x, y) 1 / sqrt(x * y)

rmse.weight <- lapply(1:n.trees, function(n) {
		m <- matrix(sqrt.weight(occ[[n]][unlist(mrcas[[n]])],occ.node.dating[[n]][unlist(mrcas.node.dating[[n]])]), nrow=100, ncol=100)
		diag(m) <- 0
		m
	})

get.min <- function(node, n) {
	mask <- unlist(mrcas.node.dating[[n]]) == node
	sum(unlist(rmse.weight[[n]])[mask] * unlist(dates[[n]])[mask]) / sum(unlist(rmse.weight[[n]])[mask])
}

dates.nodes <- lapply(1:n.trees, function(n) unlist(lapply(1:(n.tips - 1) + n.tips, get.min, n=n)))

dates.node.dating <- lapply(1:n.trees, function(n) apply(mrcas.node.dating[[n]], c(1, 2), function(m) if (m <= n.tips) 0 else dates.nodes[[n]][m - n.tips]))

s <- lapply(1:n.trees, function(n) write.csv(dates.node.dating[[n]], paste0("HIV_", n, "/HIV_", n, "_dates_best.csv")))
s <- lapply(1:n.trees, function(n) write.csv(mrcas.node.dating[[n]], paste0("HIV_", n, "/HIV_", n, "_mrca_best.csv")))