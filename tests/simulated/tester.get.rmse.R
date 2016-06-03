n.trees = 50

mrcas <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca.csv"), header=T)[,2:101])
dates <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates.csv"), header=T)[,2:101])
mrcas.node.dating <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_mrca_node.dating.csv"), header=T)[,2:101])
dates.node.dating <- lapply(1:n.trees, function(n) read.csv(paste0("HIV_", n, "_dates_node.dating.csv"), header=T)[,2:101])

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

rmse <- apply(1:n.trees, function(i) rmse(mrca.test, dates.test, i, sqrt.weight))