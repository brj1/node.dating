source("../../src/node.dating.R")
source("../../src/extra.R")

extract_dates <- function(t) as.numeric(gsub("(.+)_([0123456789\\.]+)$", "\\2", t$tip.label, perl=T))
extract_types <- function(t) gsub("(.+)_((PBMC|PLASMA))_(.+)$", "\\2", t$tip.label, perl=T)

mark.dna <- function(t, col.rna, col.dna) {
	x <- rep(col.rna, length(t$tip.label))
	x[t$tip.type == "PBMC"] <- col.dna
	x
}

tree <- read.tree("patient_16617.tre")
tree$tip.date <- extract_dates(tree)

print(tree)

tree <- rtt(tree, tree$tip.date, opt.tol=1e-8)

tree$tip.date <- extract_dates(tree)
tree$tip.type <- extract_types(tree)

tree$mu <- estimate.mu(tree, tree$tip.date)
tree$node.date <- estimate.dates(tree, tree$tip.date, tree$mu, nsteps=1000, show.step=100, is.binary=T, opt.tol=1e-8)

pdf("patient_16617.pdf", height=8.5, width=11)

plot.time.tree(tree, col.tip=mark.dna(tree, "#6600aa66", "#00aa6666"), pch.tip=mark.dna(tree, 16, 15), main="Genetic distance versus time of Patient 16617 from LANL", xlab="Time since sero-conversion (days)", ylab="Genetic distance from the root (substitions per base)")

legend(2170, .02, c("RNA samples", "DNA samples", "Internal nodes"), col=c("#00aa66ff", "#6600aaff", "#aa6600ff"), pch=c(16, 15, 5))