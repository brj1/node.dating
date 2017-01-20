library(ape)

extract_dates <- function(x) as.numeric(gsub("(.+)_([0123456789.]+)", "\\2", x, perl=T))
n.trees <- 50

for (i in 1:n.trees) {
      tree <- read.nexus(file=paste0("HIV_", i, "/HIV_", i, "_rooted.nex"))
       
      m <- mrca(tree)
      dates <- node.depth.edgelength(tree)
      shift <- extract_dates(tree$tip.label)[1] - dates[1]
      dates <- dates + shift
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_mrca_tempest.csv"))
      
      m <- apply(m, c(1, 2), function(n) {dates[n]})
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_dates_tempest.csv"))
                    
      tree <- read.tree(file=paste0("HIV_", i, "/HIV_", i, "_rooted_lsd.date.newick"))
      
      m <- mrca(tree)
      dates <- node.depth.edgelength(tree)
      shift <- extract_dates(tree$tip.label)[1] - dates[1]
      dates <- dates + shift
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_mrca_lsd.csv"))
      
      m <- apply(m, c(1, 2), function(n) {dates[n]})
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_dates_lsd.csv"))
      
      tree <- read.nexus(file=paste0("HIV_", i, "/HIV_", i, ".beast.10e4.nex"))
      
      m <- mrca(tree)
      dates <- node.depth.edgelength(tree)
      shift <- extract_dates(tree$tip.label)[1] - dates[1]
      dates <- dates + shift
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_mrca_beast.10e4.csv"))
      
      m <- apply(m, c(1, 2), function(n) {dates[n]})
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_dates_beast.10e4.csv"))
      
      tree <- read.nexus(file=paste0("HIV_", i, "/HIV_", i, ".beast.10e6.nex"))
      
      m <- mrca(tree)
      dates <- node.depth.edgelength(tree)
      shift <- extract_dates(tree$tip.label)[1] - dates[1]
      dates <- dates + shift
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_mrca_beast.10e6.csv"))
      
      m <- apply(m, c(1, 2), function(n) {dates[n]})
      
      write.csv(m, file=paste0("HIV_", i, "/HIV_", i, "_dates_beast.10e6.csv"))
}