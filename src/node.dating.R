## node.dating.R (2016-06-21)
## This file is part of the R-package `ape'.
## See the file COPYING in the package ape available at cran.r-project.org for licensing issues.

# Copyright (c) 2016, Bradley R. Jones, BC Centre for Excellence in HIV/AIDS
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the BC Centre for Excellence in HIV/AIDS nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL The BC Centre for Excellence in HIV/AIDS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Estimate the mutation rate and node dates based on tip dates.
#
# Felsenstein, Joseph. "Evolutionary trees from DNA sequences: A maximum
# likelihood approach." Journal of Molecular Evolution 17 (1981):368-376.
#
# Rambaut, Andrew. "Estimating the rate of molecular evolution: incorporating 
# non-contemporaneous sequences into maximum likelihood phylogenies." 
# Bioinformatics 16.4 (2000): 395-399.

library(ape)

# Estimate the mutation rate of a phylogenetic tree from the tip dates using 
# linear regression. This model assumes that the tree follows a molecular 
# clock.
#
# t: rooted tree with edge lengths equal to genetic distance
#
# tip.dates: vector of dates for the tips, in the same order as t$tip.label.
#            Tip dates can be censored with NA values
#
# returns a list containing the tree, the date of the root, the mutation rate,
# the log likelihood of the linear regression and the log likelihood of the
# null model (mu=0)
estimate.mu <- function(t, node.dates) {
	# fit linear model
	g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~ node.dates, na.action=na.omit)
	null.g <- glm(node.depth.edgelength(t)[1:length(node.dates)] ~ 1, na.action=na.omit)
	
	list(tree=t, root.date=coef(g)[[1]], mu=ceof(g)[[2]], logLik=logLik(g), null.logLik(null.g))
}

# Estimate the dates of the internal nodes of a phylogenetic tree.
#
# t: rooted tree with edge lengths equal to genetic distance
#
# node.dates: either a vector of dates for the tips, in the same order as 
#             t$tip.label; or a vector of dates to initalize each node
#
# mu: mutation rate, either a vector of size one for a strict molecular clock
#     or a vector with a local molecular clock along each edge
#
# min.date: the minimum date that a node can have (needed for optimize()). The 
#           default is -.Machine$double.xmax
#
# show.steps: set to print the log likelihood every show.steps. Set to 0 to 
#             supress output
#
# opt.tol: tolerance for optimization precision. By default, the optimize()
#          function uses a tolerance of .Machine$double.eps^0.25 (see ?optimize)
#
# lik.tol: tolerance for likelihood comparison. estimate.dates will stop when
#          the log likelihood between successive trees is less than like.tol. If
#          0 will stop after nsteps steps.
#
# nsteps: the maximum number of steps to run. If 0 will run until the log 
#         likelihood between successive runs is less than lik.tol. The default 
#         is 1000.
#
# is.binary: if the phylogentic tree is binary, setting is.binary to TRUE, will 
#            run a optimization method
#
# If lik.tol and nsteps are both 0 then estimate.dates will only run the inital 
# step.
#
# returns a list containing the tree, a vector of the estimated dates of the 
# tips and internal nodes, the mutation rate and the log likelihood of the tree
estimate.dates <- function(t, node.dates, mu = estimate.mu(t, node.dates)$mu, min.date = -.Machine$double.xmax, show.steps = 0, opt.tol = 1e-8, nsteps = 1000, lik.tol = if (nsteps <= 0) opt.tol else 0, is.binary = is.binary.tree(t)) {
	if (any(mu < 0))
		stop(paste("mu (", mu, ") less than 0", sep=""))
		
	mu <- if (length(mu) == 1) rep(mu, length(t$edge.length)) else mu

	# init vars
	n.tips <- length(t$tip.label)
	nodes <- unique(reorder(t)$edge[,1])
	dates <- if (length(node.dates) == n.tips) {
			c(node.dates, rep(NA, t$Nnode))
		} else {
			node.dates
		}
		
	lik.sens <- if (lik.tol == 0) opt.tol else lik.tol
	
	# Don't count initial step if all values are seeded
	iter.step <-  if (any(is.na(dates))) 0 else 1
	
	children <- lapply(1:t$Nnode,
		function(x) {
			which(t$edge[,1] == x + n.tips)
		})
	parent <- lapply(1:t$Nnode,
		function(x) {
			which(t$edge[,2] == x + n.tips)
		})
		
	# to process children before parents
	nodes <- c(1)
	i <- 1
	
	while (i <= length(nodes)) {                                                                                                                                                                                                                                                                                                                                                                            
		nodes <- c(nodes, t$edge[t$edge[,1] == nodes[i] + n.tips, 2] - n.tips)
		
		i <- i + 1
	}
	
	nodes <- nodes[nodes > 0]
	nodes <- rev(nodes)
		
	# calculate likelihood
	scale.lik <- sum(-lgamma(t$edge.length+1)+(t$edge.length+1)*log(mu))
		
	calc.Like <- function(ch.node, ch.edge, x) {
		tim <- ch.node - x
						
		t$edge.length[ch.edge]*log(tim)-mu[ch.edge]*tim
	}
		
	opt.fun <- function(x, ch, p, ch.edge, p.edge, use.parent=T) {	
		sum(if (!use.parent || length(dates[p]) == 0 || is.na(dates[p])) {		
				calc.Like(dates[ch], ch.edge, x)
			} else {
				calc.Like(c(dates[ch], x), c(ch.edge, p.edge), c(rep(x, length(dates[ch])), dates[p]))
			})
	}
	
	solve.bin <- function(bounds, ch.times, ch.edge) {
		ch.edge.length <- t$edge.length[ch.edge]
		
		a <- sum(mu[ch.edge])
		b <- ch.edge.length[1] + ch.edge.length[2] - a * (ch.times[1] + ch.times[2])
		c.0 <- a*ch.times[1] * ch.times[2] - ch.times[1] * ch.edge.length[2] - ch.times[2] * ch.edge.length[1]
						
		if (b ^ 2 - 4 * a * c.0 < 0) {		
			return(bounds[1 + (sum(calc.Like(ch.times, ch.edge, bounds[2] - opt.tol)) > sum(calc.Like(ch.times, ch.edge, bounds[1] + opt.tol)))])
		}
		else {
			x.1 <- (-b + sqrt(b ^ 2 - 4 * a * c.0)) / (2 * a)
			x.2 <- (-b - sqrt(b ^ 2 - 4 * a * c.0)) / (2 * a)
			
			x <- c(bounds[1] + opt.tol, bounds[2] - opt.tol)
			if (bounds[1] < x.1 && x.1 < bounds[2])
				x <- c(x, x.1)
			if (bounds[1] < x.2 && x.2 < bounds[2])
				x <- c(x, x.2)
		 			 	
			return(x[which.max(unlist(lapply(x, function(y) sum(calc.Like(ch.times, ch.edge, y)))))])
		}
	}
	
	solve.cube <- function(bounds, ch.times, ch.edge, par.time, par.edge) {
		ch.edge.length <- t$edge.length[ch.edge]
		par.edge.length <- t$edge.length[par.edge]
	
		a <- sum(mu[ch.edge]) - mu[par.edge]
		b <- sum(ch.edge.length) + par.edge.length - a * (sum(ch.times) + par.time)
		c.0 <- a * (ch.times[1] * ch.times[2] + ch.times[1] * par.time + ch.times[2] * par.time) - (ch.times[1] + ch.times[2]) * par.edge.length - (ch.times[1] + par.time) * ch.edge.length[2] - (ch.times[2] + par.time) * ch.edge.length[1]
		d <- ch.edge.length[1] * ch.times[2] * par.time + ch.edge.length[2] * ch.times[1] * par.time + par.edge.length * ch.times[1] * ch.times[2] - a * prod(ch.times) * par.time
				
		delta.0 <- complex(real=b^2 - 3 * a * c.0)
		delta.1 <- complex(real=2 * b^3 - 9 * a * b * c.0 + 27 * a^2 * d)
		C <- ((delta.1 + sqrt(delta.1^2 - 4 * delta.0^3)) / 2)^(1/3)
		
		x.1 <- Re(-1 / (3 * a) * (b + complex(real=1) * C + delta.0 / (complex(real=1) * C)))
		x.2 <- Re(-1 / (3 * a) * (b + complex(real=-1/2, imaginary=sqrt(3)/2) * C + delta.0 / (complex(real=-1/2, imaginary=sqrt(3)/2) * C)))
		x.3 <- Re(-1 / (3 * a) * (b + complex(real=-1/2, imaginary=-sqrt(3)/2) * C + delta.0 / (complex(real=-1/2, imaginary=-sqrt(3)/2) * C)))
		
		x <- c(bounds[1] + opt.tol, bounds[2] - opt.tol)
		if (x.1 && bounds[1] < x.1 && x.1 < bounds[2])
			x <- c(x, x.1)
		if (bounds[1] < x.2 && x.2 < bounds[2])
			x <- c(x, x.2)
		if (bounds[1] < x.3 && x.3 < bounds[2])
			x <- c(x, x.3)
		
		return(x[which.max(unlist(lapply(x, function(y) sum(calc.Like(c(ch.times, y), c(ch.edge, par.edge), c(y, y, par.time))))))])
	}
	
	estimate <- function(node) {
		ch.edge <- children[[node]]
		ch <- t$edge[ch.edge, 2]
		
		p.edge <- parent[[node]]
		p <- t$edge[p.edge, 1]
		
		m <- if (length(p) == 0 || is.na(dates[p])) {
				min.date
			} else {
				dates[p]
			}
			
		if (is.binary) {
			if (length(dates[p]) == 0 || is.na(dates[p]))
				solve.bin(c(m, min(dates[ch])), dates[ch], ch.edge)
			else 
				solve.cube(c(m, min(dates[ch])), dates[ch], ch.edge, dates[p], p.edge)
		}
		else {				
			res <- optimize(opt.fun, c(m, min(dates[ch])), ch, p, ch.edge, p.edge, maximum=T)
		
			res$maximum
		}
	}
	
	# iterate to estimate dates
	lik <- NA
	
	repeat
	{
		for (n in nodes) {		
			dates[n + n.tips] <- estimate(n)
		}
				
		new.lik <- sum(calc.Like(dates[tree$edge[,2]], 1:length(tree$edge.length), dates[tree$edge[,1]])) + scale.lik
		
		if (show.steps > 0 && ((iter.step %% show.steps) == 0)) {
			cat(paste("Step: ", iter.step, ", Likelihood: ", new.lik, "\n", sep=""))
		}
		
		if ((lik.tol > 0 && (!is.na(lik) && (is.infinite(lik) || is.infinite(new.lik) || new.lik - lik < lik.tol))) || (nsteps > 0 && iter.step >= nsteps) || (lik.tol <= 0 && nsteps <= 0)) {
			if (is.infinite(lik) || is.infinite(new.lik)) {
				warning("Likelihood infinite")
			}
			else if (!is.na(lik) && new.lik + lik.sens < lik) {			
				warning("Likelihood less than previous estimate")
			}
					
			break
		} else {
			lik <- new.lik
		}
		
		iter.step <- iter.step + 1
	}
	
	if (show.steps > 0) {
		cat(paste("Step: ", iter.step, ", Likelihood: ", new.lik, "\n", sep=""))
	}
	
	list(tree=t, node.date=dates, mu=mu, logLik=new.lik)
}


