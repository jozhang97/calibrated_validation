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
# pars = c(lambda, lambda, mu, mu, q, q)

pars = c(0.2, 0.4, 0.001, 0.1, 0.1, 0.4)  # from RevBayes exp


# Get a ground truth tree
phy = tree.bisse(pars, max.taxa=22, x0=0)
if (is.null(phy)) {
    print("bad tree")
    q()
}
tips = phy$tip.state
print("Tree tips")
print(tips)
print("Tree in Newick format")
write.tree(phy)
ntips = length(tips)
node.truth = phy$node.state


# Calculate likelyhood on tree
sampling.f = c(1,1)
lik = make.bisse(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE)
tree.lik = lik(pars=pars, root.p=NULL, intermediates=TRUE, condition.surv=FALSE) 
# prior on root is the weighted average of D0 and D1, i.e., ROOT.OBS = D = D0 * (D0/(D0+D1)) + D1 * (D1/(D0+D1))
print("Tree log likeyhood")
print(tree.lik[1])

# Calculate ancestral MLE
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
print("Compare asr to ground truth")
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


# Write the results out
print("asr marginal likelyhoods")
print(asr.marginal)
write.csv(asr.marginal, file = "diversitree_anc_states.csv")


