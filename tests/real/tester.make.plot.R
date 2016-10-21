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

par(mar=c(4.5,5.1,2.2,1))

plot.time.tree(tree, col.tip=mark.dna(tree, "#00aa6666", "#00aa6666"), pch.tip=mark.dna(tree, 15, 15), xlab="Time since first sample (days)", ylab="Genetic distance from root (subs. per base)", cex.axis=2, cex.lab=2.5, cex=1.2)

legend(-240, .165, c("HIV samples", "Internal nodes"), col=c("#00aa66ff", "#aa6600ff"), pch=c(15, 5), cex=2.5)