# Copyright (c) 2016, Bradley R. Jones, BC Centre for Excellence in HIV/AIDS
# All rights reserved.
#
# Extra methods for node.dating

# Calculate the tree log likelihood 
#
# t: rooted tree with edge lengths equal to genetic distance
#
# tip.dates: vector of dates for the tips, in the same order as t$tip.label
#
# mu: mutation rate
#
# returns the log likelihood as a double
tree.like <- function(tree, node.dates, mu) {
	if (any(mu < 0))
		return(-Inf)
	
	scale.lik <- sum(-lgamma(tree$edge.length+1)+(tree$edge.length+1)*log(mu))
		
	calc.lik <- function(ch.node, edge, par.node) {
		tim <- ch.node - par.node
						
		edge*log(tim)-mu*tim
	}
			
	sum(calc.lik(node.dates[tree$edge[,2]], tree$edge.length, node.dates[tree$edge[,1]])) + scale.lik
}

# Auxiliary method for plot.tree
col.nodes <- function(t, col.tip, col.node) {
	tip.cols <- if (length(col.tip) == length(t$tip.label)) col.tip else rep(col.tip[[1]], length(t$tip.label))
	node.cols <- if (length(col.node) == t$Nnode) col.node else rep(col.node[[1]], t$Nnode)
			
	c(tip.cols, node.cols)
}

# Plots a tree's node in a distance versus time
plot.time.tree <- function(t, col.tip="#00aa6666", col.node="#aa660066", pch.tip=16, pch.node=5, ...) {
	plot(t$node.date, node.depth.edgelength(t), col=col.nodes(t, col.tip, col.node), pch=col.nodes(t, pch.tip, pch.node), ...)
	apply(t$edge, 1, function(e) {
			arrows(t$node.date[e[1]], node.depth.edgelength(t)[e[1]], t$node.date[e[2]], node.depth.edgelength(t)[e[2]], col="#00000040", length=.1)
		})
	abline(-min(t$node.date)*t$mu, t$mu, lty=3)
}
