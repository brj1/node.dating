args <- commandArgs(trailing=T)

program <- args[1]

n.trees <- 50

mrcas <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca.csv"), header=T)[,2:101])
dates <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates.csv"), header=T)[,2:101])
mrcas.node.dating <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca_", program, ".csv"), header=T)[,2:101])
dates.node.dating <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates_", program, ".csv"), header=T)[,2:101])

o <- lapply(1:n.trees, function(n) order(names(mrcas[[n]])))
mrcas <- lapply(1:n.trees, function(n) mrcas[[n]][o[[n]], o[[n]]])
dates <- lapply(1:n.trees, function(n) dates[[n]][o[[n]], o[[n]]])

o <- lapply(1:n.trees, function(n) order(names(mrcas.node.dating[[n]])))
mrcas.node.dating <- lapply(1:n.trees, function(n) mrcas.node.dating[[n]][o[[n]], o[[n]]])
dates.node.dating <- lapply(1:n.trees, function(n) dates.node.dating[[n]][o[[n]], o[[n]]])

occ <- lapply(1:n.trees, function(n) table(unlist(mrcas[[n]])))
occ.node.dating <- lapply(1:n.trees, function(n) table(unlist(mrcas.node.dating[[n]])))

sqrt.weight <- function(x, y) 1 / sqrt(x * y)

rmse.weight <- lapply(1:n.trees, function(n) {
		m <- matrix(sqrt.weight(occ[[n]][unlist(mrcas[[n]])],occ.node.dating[[n]][unlist(mrcas.node.dating[[n]])]), nrow=100, ncol=100)
		diag(m) <- 0
		m
	})

rmse <- function(n) {	 	
	d <- rmse.weight[[n]] * (dates[[n]] - dates.node.dating[[n]])^2
				
	sqrt(sum(d) / sum(rmse.weight[[n]]))
}

rmses <- lapply(1:n.trees, function(i) rmse(i))

mean(unlist(rmses))