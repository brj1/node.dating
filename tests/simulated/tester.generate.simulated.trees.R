library(ape)
library(TreeSim)
library(NELSI)

source("node.dating.R")

extract_dates <- function(x) as.numeric(gsub("(.+)_([0123456789.]+)", "\\2", x, perl=T))

date.branches <- function(s.tree) {
	tree <- s.tree$phylogram
	subsr <- s.tree$tree.data.matrix[,6]
	times <- s.tree$tree.data.matrix[,7]
	
	# Calculate the cumulative times
	tree$edge.length <- times
	times <- node.depth.edgelength(tree)
	
	tree$edge.length <- subsr

	dates <- unlist(Map(toString, times))[1:n.tips]
	tree$tip.label <- paste(tree$tip.label, dates, sep='_')
	tree
}

raxml <- function(seq, parsimony.seed=NULL, bootstrap.seed=NULL, executable='raxmlHPC8-pthreads', threads=4, N=10, name='raxml', clear=T, model="GTRGAMMA") {
	if(class(seq) != "DNAbin") stop('seq should be of class \'DNAbin\'')
	if(is.null(parsimony.seed)) {
		parsimony.seed <- as.integer(sample(2**31,1))
		warning(sprintf('parsimony.seed should be fixed for debugging!\nSetting to %d!', parsimony.seed))
	}

	if(is.null(bootstrap.seed)) {
		bootstrap.seed <- as.integer(sample(2**31,1))
		warning(sprintf('bootstrap.seed should be fixed for debugging!\nSetting to %d!', bootstrap.seed))
	}

	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}

	wd <- path.expand(getwd())
	executable <- path.expand(executable)
	r.name <- paste0(name,"_R")
	b.name <- paste0(name,"_B")

	dnafile <- tempfile()
	write.dna(seq, dnafile)

	cmd <- sprintf('%s -T %d -b %d -m %s -p %d -N %d -s %s -n %s -w %s -O', 
					executable, threads, bootstrap.seed, model, parsimony.seed, N, dnafile, r.name, wd)

	cmd2 <- sprintf('%s -T %d -f d -m %s  -s %s -N %d -n %s -w %s -p %d -O', 
					executable, threads, model,  dnafile, N, b.name, wd, parsimony.seed)
	
	cmd3 <- sprintf('%s -T %d -f b -n %s -m %s -t %s/RAxML_bestTree.%s -z %s/RAxML_bootstrap.%s -s %s -w %s -O', 
					executable, threads, name, model, wd, b.name, wd, r.name, dnafile, wd)

	system(cmd)
	system(cmd2)
	system(cmd3)

	tr <- read.tree(sprintf('RAxML_bipartitions.%s', name))
	if(clear) {
		cwd <- dir()
		unlink(cwd[grep('^RAxML_', cwd)])
	}
	unlink(dnafile)
	
	tr
}

n.trees <- 50
n.partitions <- n.trees
n.replicates <- 1
n.tips <- 100

clock.rate <- c(0.0001964) # 0.0028
noise.rate <- c(0.00001417)
sampprob <-c(0.005237) # 2.78e-4 / ( 0.88e-4 + 2.78e-4)
lambda <- c(0.05116) # 8.23e-4
mu <- c(0.05006) # 0.88e-4 + 2.78e-4
times<-c(0)

trees <- apply(matrix(rep(n.tips,n.trees)), 1, sim.bdsky.stt, lambdasky=lambda, deathsky=mu, timesky=times, sampprobsky=sampprob, rho=0, timestop=0)
trees <- lapply(trees, function(x) {unroot(x[[1]])})

for (i in 1:n.trees) {
	m <- mrca(trees[[i]])
	dates <- node.depth.edgelength(trees[[i]])
		
	write.csv(m, file=paste0("HIV_", i, "_mrca.csv"))
		
	m <- apply(m, c(1, 2), function(n) dates[n])
		
	write.csv(m, file=paste0("HIV_", i, "_dates.csv"))
}

sim.trees <- lapply(trees, simulate.clock, params=list(rate=clock.rate, noise=noise.rate))
trees <- lapply(sim.trees, date.branches)

indel_control <- sprintf(
"
[TYPE] NUCLEOTIDE 2

[SETTINGS]
  [output]                   FASTA 
  [randomseed]               %s

[MODEL]    HKY_HIV
  [submodel] HKY 8.5               
  [statefreq] 0.42 0.15 0.15 0.28
"
, 1989)


for(i in 1:n.trees) {
	tree_dat <- write.tree(trees[[i]])
#	print(tree_dat)
	indel_control <- paste0(indel_control, sprintf("[TREE] tree_%d %s \n", i, tree_dat))
}

#index <- sample(1:n.trees, n.partitions, replace = (n.partitions > n.trees))
index <- 1:n.trees
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("[PARTITIONS] pHKY_%d [tree_%d HKY_HIV 700] \n", i, index[i]))
}

indel_control <- paste0(indel_control, "[EVOLVE] \n")
for(i in 1:n.partitions) {
	indel_control <- paste0(indel_control, sprintf("    pHKY_%d %d HIV_%d \n", i, n.replicates, i))
}

write(indel_control, 'control.unfixed.txt')

system("python fix_control.py control.unfixed.txt control.txt")
system("indelible")

trees <- lapply(seq(1, n.trees), function(i) {
		dna <- read.FASTA(paste0("HIV_", i, "_TRUE.fas"))
		tree <- raxml(dna)
		write.tree(tree, file=paste0("HIV_", i , "_unrooted.tre"))
		tree
	})

#source("node.dating.R")
#n.trees <- 50
#trees <- c(lapply(grep("HIV_[123456789]_unrooted.tre", dir("."), value=T, perl=T), read.tree), lapply(grep("HIV_[123456789][1234567890]_unrooted.tre", dir("."), value=T, perl=T), read.tree))

trees <- lapply(trees, function(t) {t$tip.date <- extract_dates(t$tip.label); t})
trees <- lapply(trees, function(t) {rtt(t, t$tip.date, opt.tol=1e-8)})
lapply(seq(1, n.trees), function(i) write.tree(trees[[i]], file=paste0("HIV_", i , "_rooted.tre")))

trees <- lapply(trees, function(t) {t$tip.date <- extract_dates(t$tip.label); t})
trees <- lapply(trees, function(t) {t$mu <- estimate.mu(t, t$tip.date); t})
trees <- lapply(trees, function(t) {t$node.date <- estimate.dates(t, t$tip.date, t$mu, opt.tol=1e-8, lik.tol=1e-10, nsteps=1000, show.steps=100); t})

for (i in 1:n.trees) {
	m <- mrca(trees[[i]])
	dates <- trees[[i]]$node.date
	
	write.csv(m, file=paste0("HIV_", i, "_mrca_node.dating.csv"))
	
	m <- apply(m, c(1, 2), function(n) dates[n])
	
	write.csv(m, file=paste0("HIV_", i, "_dates_node.dating.csv"))
}

for (i in 1:n.trees) {
	system(paste0("java -classpath ~/random/tempest dr.app.tempest.Hacker HIV_", i, "_rooted.tre HIV_", i, "_rooted.nex"))
		
	tree <- read.nexus(file=paste0("HIV_", i, "_rooted.nex"))
	
	m <- mrca(tree)
	dates <- node.depth.edgelength(tree)
	shift <- extract_dates(tree$tip.label)[1] - dates[1]
	dates <- dates + shift
	
	write.csv(m, file=paste0("HIV_", i, "_mrca_tempest.csv"))
	
	m <- apply(m, c(1, 2), function(n) {dates[n]})
	
	write.csv(m, file=paste0("HIV_", i, "_dates_tempest.csv"))
}

for (i in 1:n.trees) {
	write.table(100, file=paste0("HIV_", i, "_dates.txt"), row.names=F, col.names=F)
	write.table(data.frame(x=trees[[i]]$tip.label, y=trees[[i]]$tip.date), file=paste0("HIV_", i, "_dates.txt"), row.names=F, col.names=F, sep="\t", append=T, quote=F)

	system(paste0("~/Downloads/lsd-0.2/bin2/src/lsd -i HIV_", i, "_rooted.tre -o HIV_", i, "_rooted_lsd -d HIV_", i, "_dates.txt -c -n 1"))
		
	tree <- read.tree(file=paste0("HIV_", i, "_rooted_lsd.date.newick"))
	
	m <- mrca(tree)
	dates <- node.depth.edgelength(tree)
	shift <- extract_dates(tree$tip.label)[1] - dates[1]
	dates <- dates + shift
	
	write.csv(m, file=paste0("HIV_", i, "_mrca_lsd.csv"))
	
	m <- apply(m, c(1, 2), function(n) {dates[n]})
	
	write.csv(m, file=paste0("HIV_", i, "_dates_lsd.csv"))
}

for (i in 1:n.trees) {		
	tree <- read.nexus(file=paste0("~/beast_ancre/node.dating/HIV_", i, "/HIV_", i, ".nex"))
	
	m <- mrca(tree)
	dates <- node.depth.edgelength(tree)
	shift <- extract_dates(tree$tip.label)[1] - dates[1]
	dates <- dates + shift
	
	write.csv(m, file=paste0("HIV_", i, "_mrca_beast.csv"))
	
	m <- apply(m, c(1, 2), function(n) {dates[n]})
	
	write.csv(m, file=paste0("HIV_", i, "_dates_beast.csv"))
}