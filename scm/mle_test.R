library(ape)
library(diversitree)
set.seed(1234)

## HELPERS
# Compares the states in the truth and the reconstruction for given states/indices of interest
compare.states = function(truth, pred, indices) {
    num.correct = 0
    for (i in indices) {
        if (truth[i] == pred[i]) {
            num.correct  = num.correct + 1
        }
    }
    num.samples = length(indices)
    return (num.correct / num.samples)
}


# Set parameters
lambda = 1/20
mu = 1/80
q = 1/20
sim.time = 50
pars = c(lambda, lambda, mu, mu, q, q)


# Get a ground truth tree
phy = tree.bisse(pars, sim.time, max.taxa=1000, include.extinct=FALSE, x0=NA)
if (is.null(phy)) {
    print("bad tree")
    q()
}
tips = phy$tip.state
ntips = length(tips)
node.truth = phy$node.state


# Calculate MLE on tree
lik = make.bisse(phy, phy$tip.state)
asr.marginal = asr.marginal(lik, pars)


# Compare truth tips and MLE tips for accuracy
node.marginal = vector("list", phy$Nnode)
for (i in 1:phy$Nnode) {
    if (asr.marginal[1,i] > asr.marginal[2,i]) {
        node.marginal[i] = 1 - 1
    } else {
        node.marginal[i] = 2 - 1
    }
}

acc = compare.states(node.truth, node.marginal, 1:phy$Nnode)
cat("Total accuracy: ", acc, "\n")


# Compare accuracy on the more recent and ancient nodes
all.node.depth = node.depth.edgelength(phy)  # Note this includes tips too... time from beginning of time 
node.sorted = order(all.node.depth)[1 : phy$Nnode] # gives us the order of internal nodes from most ancient to least
node.ancient = node.sorted[1 : as.integer(phy$Nnode/4)] - ntips # quartile of most ancient 
node.recent = node.sorted[as.integer(phy$Nnode * 3/4) : phy$Nnode] - ntips # quartile of least ancient

acc.ancient = compare.states(node.truth, node.marginal, node.ancient)
cat("Ancient accuracy: ", acc.ancient, "\n")
acc.recent = compare.states(node.truth, node.marginal, node.recent)
cat("Recent accuracy: ", acc.recent, "\n")



