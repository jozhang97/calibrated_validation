# This file has a few important outputs
#  1. results/diversitree_ancestral_states.pdf - tree visualization with character state likelyhoods
#  2. sim_data/classe.tree
#  3. sim_data/tips.csv - list of node to character mappings

library(diversitree)
library(ape)
set.seed(1234)

# simulate tree and character under CLaSSE

# used similar parameters as the BiSSE script 
# (difference is low unqiue cladogenetic speciation rates)
lambda_111 = 0.2
lambda_112 = 0.000001
lambda_122 = 0.000001
lambda_211 = 0.000001
lambda_212 = 0.000001
lambda_222 = 0.4

mu_1 = 0.01
mu_2 = 0.1

q_12 = 0.1
q_21 = 0.4

# Use diversitree:::default.argnames.classe(3) to see arg names
## pars_test = c(runif(27)) # classe(3) has 27 pars, for a separate test
## phy_test = tree.classe(pars_test, max.t=2, max.taxa=Inf, 
##                  include.extinct=FALSE, x0=NA)  
## print(phy_test)

pars = c(lambda_111, lambda_112, lambda_122,
         lambda_211, lambda_212, lambda_222,
         mu_1, mu_2, q_12, q_21)
pars_bisse = c(lambda_111, lambda_222, mu_1, mu_2, q_12, q_21)

# simulates and returns a phylogeny given parameters and 
# initial state for max.t amount of time 
# No cap on the number of taxa
# Do not need information on extinct taxa
# Sample initial state from equilibrium distribution
phy = tree.classe(pars, max.t=22, max.taxa=Inf, 
                  include.extinct=FALSE, x0=NA)  
print("CLaSSE Tree")
print(phy)
states = phy$tip.state

phy_bisse = tree.bisse(pars_bisse, max.t=22, max.taxa=Inf, 
                  include.extinct=FALSE, x0=NA)  
print("BiSSE Tree")
print(phy_bisse)
# ------------------------------------

num_iter = 25
num_total_nodes_bisse = 0
num_total_tips_bisse = 0
num_total_nodes = 0
num_total_tips = 0
for (v in c(runif(num_iter))){
    phy_bisse = tree.bisse(pars_bisse, max.t=22, max.taxa=Inf, 
                          include.extinct=FALSE, x0=NA)  
    num_nodes_bisse = phy_bisse$Nnode
    num_tips_bisse = length(phy_bisse$tip.state)
    if (!is.null(num_nodes_bisse)) {
        num_total_nodes_bisse = num_total_nodes_bisse + num_nodes_bisse
    }
    num_total_tips_bisse = num_total_tips_bisse + num_tips_bisse

    phy = tree.classe(pars, max.t=22, max.taxa=Inf, 
                      include.extinct=FALSE, x0=NA)  
    num_nodes = phy$Nnode
    num_tips = length(phy$tip.state)
    if (!is.null(num_nodes)) {
        num_total_nodes = num_total_nodes + num_nodes
    }
    num_total_tips = num_total_tips + num_tips
}
print(num_total_nodes_bisse)
print(num_total_tips_bisse)
print(num_total_nodes)
print(num_total_tips)

print(states)
write.tree(phy, file="sim_data/classe.tree")
write.tree(phy_bisse, file="sim_data/bisse.tree")
q()

# calculate likelihood
# outputs likelyhood function which takes in parameters and outputs the likelyhood of the parameters
lik = make.bisse(phy, states, strict=FALSE)
rate = pars[1]
num_taxa = length(states)
lnl = lik(pars,root=ROOT.FLAT, condition.surv=!TRUE, intermediates=TRUE) - log( rate ) + num_taxa * log(2) - sum(log(1:num_taxa)) - log(2)
#cat("diversitree lnl =", lnl, "\n", file="results/likelihoods.txt", append=TRUE)

# infer marginal ancestral state reconstructions
# output is N x V where 
#   N number of character states (e.g. 2)
#   V number of nodes in the tree
# This is where diversitree comes up with the likelyhoods of the character states per node 
anc_states = asr.marginal(lik, pars)

# plot the inferred ancestral states
pdf("results/diversitree_ancestral_states.pdf")
plot(phy, cex=.5, label.offset=0.2)
col = c("#004165", "#eaab00")
nodelabels(pie=t(anc_states), piecol=col, cex=.5)
dev.off() # This line gives me an Null Device thing

# now write the tree to a newick file
write.tree(phy, file="sim_data/bisse.tree")

# now write tip data file
a = FALSE
for (i in 1:num_taxa) {
    cat(c(names(states[i]), states[i]), file="sim_data/tips.csv", append=a, sep=",")
    a = TRUE
    cat("\n", append=a, file="sim_data/tips.csv")
}
