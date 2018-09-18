library(ape)
library(diversitree)
set.seed(1234)

lambda = 1/20
mu = 1/80
q = 1/20
sim.time = 50
pars = c(lambda, lambda, mu, mu, q, q)

# Get a ground truth tree
phy = tree.bisse(pars, sim.time, max.taxa=10000000, include.extinct=FALSE, x0=NA)
if (is.null(phy)) {
    print("bad tree")
    q()
}
tips = phy$tip.state
node.truth = phy$node.state
num.nodes = length(node.truth)

# Calculate MLE on tree
lik = make.bisse(phy, phy$tip.state)
# fit = find.mle(lik, pars)
# history = history.from.sim.discrete(phy, 0:1)
asr.marginal = asr.marginal(lik, pars)

# Compare truth tips and MLE tips for accuracy
node.marginal = vector("list", num.nodes)
for (i in 1:num.nodes) {
    if (asr.marginal[1,i] > asr.marginal[2,i]) {
        node.marginal[i] = 1 - 1
    } else {
        node.marginal[i] = 2 - 1
    }
}

num.correct = 0
for (i in 1:num.nodes) {
    if (node.marginal[i] == node.truth[i]) {
        num.correct  = num.correct + 1
    }
}
print(num.correct/num.nodes)
