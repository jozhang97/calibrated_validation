library(ggplot2)
library(gridExtra)
library(sjPlot)
library(diversitree)
library(HDInterval)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)
output.dir <- args[1]
prefix <- args[2]
n.sim <- args[3]
param.names <- args[4]
prior.dists <- args[5]
prior.params <- args[6]
sim.time <- as.numeric(args[7])

param.names.vec <- unlist(strsplit(param.names, split=","))
prior.dists.vec <- unlist(strsplit(prior.dists, split=","))
prior.params.vec <- unlist(strsplit(prior.params, split=";"))

# --- START: Prior sampling stuff --- #
get.95 <- function(a.vector) {
    res = hdi(a.vector, cresMass=.95)
    return(c(res[[1]], res[[2]], mean(a.vector)))
}

sample.parameter <- function(output.dir, prefix, n.sim, param.name, prior.dist, prior.params, plot.flag) {
## sample.parameter <- function(output.dir, prefix, n.sim, mu, std, param.name) {
    # Setting up arbitrary simulation values, and making objects
    ignored.flag = FALSE
    s = rep(0, n.sim)
    ps = rep(0, 1000)
    xaxis.min = 0
    xaxis.max = 0
    if (prior.dist == "lnorm") {
        mu.std = unlist(strsplit(prior.params, split=","))
        mu = as.numeric(mu.std[1]); std = as.numeric(mu.std[2])
        s = rlnorm(n.sim, meanlog=mu, sdlog=std)
        ps = rlnorm(1000, meanlog=mu, sdlog=std)
        xaxis.min = min(ps) - 0.05
        xaxis.max = max(ps) + 0.05
    }
    else if (prior.dist == "exp") {
        r = as.numeric(prior.params)
        s = rexp(n.sim, rate=r)
        ps = rexp(1000, rate=r)
        xaxis.max = max(ps)
    }
    else if (prior.dist == "NA") {
        ignored.flag = TRUE # parameters (from triplets) that are necessary to simulate under CLaSSE, but which we do not allow for
        s = rep(0, 1000)
    } 
    df <- data.frame(s); names(df) <- "value" # samples df

    # only make plot of parameter (triplet) that we allowed for 
    if (!ignored.flag & plot.flag) {
        cat(paste0("Drawing param ", param.name, " from ", prior.dist, " prior (moments: ", prior.params, ")"))
        cat(paste0(" \\ Prior mean: ", mean(ps), "\n"))
    
    ## after ggtitle
    ## stat_function(fun=dlnorm,
    ##           args=list(r[1]:r[2], meanlog=mean(log(p.df$value)), sdlog=sd(log(p.df$value)))
    ##           ) +

        prior.samples.plot <- ggplot(df, aes(x=value)) +
            geom_histogram(aes(y=stat(density)), alpha=.4, bins=100) +
            scale_x_continuous(limits=c(xaxis.min, xaxis.max)) +
            xlab(paste0(param.name, " value")) + ylab("Density") + 
            theme(
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.background = element_blank(),
                plot.title = element_text(hjust=0.5),
                axis.line = element_line(),
                axis.ticks = element_line(color="black"),
                axis.text.x = element_text(color="black", size=10),
                axis.text.y = element_text(color="black", size=10),
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12)
            )

        return(list(prior.samples.plot, df))
    }
    ## ggsave(file=paste0(output.dir, param.name, ".pdf"), plot=prior.samples.plot, width=3, height=3)
    ## prior.samples.plot # uncomment and run to see graph
    else {
        ## cat(paste0("Did not plot ", param.name, ".\n"))
        return(list(NA, df))
    }
}

sample.parameters <- function(output.dir, prefix, n.sim, param.names.vec, prior.dists.vec, prior.params.vec, plot.flag) {
    num_params = length(param.names.vec)
    params = vector("list", length = num_params)
    plots = vector("list", length = num_params)

    for (i in 1 : num_params) {
        param.name = param.names.vec[i]
        ## mu = mus.vec[i]
        ## std = stds.vec[i]
        prior.dist = prior.dists.vec[i]
        prior.params = prior.params.vec[i]
        res = sample.parameter(output.dir, prefix, n.sim, param.name, prior.dist, prior.params, plot.flag)
        plots[[i]] = res[[1]]
        params[[i]] = res[[2]]
    }

    params.df = as.data.frame(params)
    names(params.df) = param.names.vec

    # printing all panels in single graph (png() won't work)
    if (plot.flag) {
        cat("Plotting prior_samples.pdf\n")
        pdf(paste0(output.dir, "prior_samples.pdf"), width=6, height=7)
        plot_grid(plots[!is.na(plots)])
        dev.off()
    }

    return(params.df)
}

# ----------- Test Run --------------------------
# output.dir <- "~/Desktop/"
# prefix <- "test"
# param.name <- "parameter1"
# n.sim <- 10000 # number of simulations
# sample.parameters(output.dir, prefix, n.sim, param.name)

# ----------- True Run --------------------------
if (length(args) < 6) {
  print("Not enough args. Output directory, Prefix, Num simulations, mu, sigma, parameter names")
}

params.df <- sample.parameters(output.dir, prefix, n.sim, param.names.vec, prior.dists.vec, prior.params.vec, TRUE) # table with all true parameters for all simulations
init.params.df <- sample.parameters(output.dir, prefix, min(100, n.sim), param.names.vec, prior.dists.vec, prior.params.vec, FALSE) # table with all initialization (for .xml) parameters for all simulations
write.table(init.params.df, file=paste0(output.dir, "data_param_inits.csv"), row.names=FALSE, quote=FALSE, sep="|")
# --- END: Prior sampling stuff --- #

# --- START: Simulations --- #
plot.ntips <- function(my.df) {
    cat("\nPlotting ntips.pdf\n")
    res <- ggplot(my.df, aes(x=ntips)) +
        geom_histogram(aes(y=stat(density)), alpha=.4, bins=100) +
        ggtitle(paste("Median number of tips =", median(my.df$ntips), sep=" ")) +
        xlab("Number of tips") + ylab("Density") + 
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(hjust=0.5),
            axis.line = element_line(),
            axis.ticks = element_line(color="black"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12)
        ) + xlim(0,1000)

    ggsave(paste0(output.dir, "ntips.pdf"), plot=res, width=6, height=3)
}

# BiSSE
## initializing cols in df
params.df$tree <- NA
params.df$ntips <- 0
params.df$tipstates <- NA

## getting prior HDIs
many.samples.from.priors <- vector("list", length(param.names.vec))
for (i in 1:length(prior.params.vec)) {
    ignored = FALSE
    prior.dist = prior.dists.vec[i]

    if (prior.dist == "lnorm") {
        mu.std = unlist(strsplit(prior.params.vec[i], split=","))
        mu = as.numeric(mu.std[1]); std = as.numeric(mu.std[2])
        many.samples.from.priors[[i]] = rlnorm(20000, meanlog=as.numeric(mu), sdlog=as.numeric(std))
    }

    else if (prior.dist == "exp") {
        rate = prior.params.vec[i]
        many.samples.from.priors[[i]] = rexp(20000, rate=as.numeric(rate))
    }

    else if (prior.dist == "NA") { ignored = TRUE } 
    
    if (!ignored) {
        prior.hdi = get.95(many.samples.from.priors[[i]])
        lower.colname = paste0("hdilower.", param.names.vec[i])
        upper.colname = paste0("hdiupper.", param.names.vec[i])
        params.df = cbind(params.df, prior.hdi[[1]])
        params.df = cbind(params.df, prior.hdi[[2]])
        names(params.df)[c(ncol(params.df)-1,ncol(params.df))] = c(lower.colname, upper.colname)
        ## cat(paste0(param.names.vec[i], " hdilower ", prior.hdi[[1]], " hdiupper ", prior.hdi[[2]], " mean ", prior.hdi[[3]], "\n")) # verify they match samples
    }
}

## simulating, storing and printing
too.large <- 0
too.small <- 0
simulated.trees <- 0
save(params.df, file="test_sim_classe.RData")
for (i in 1:nrow(params.df)) {
    pars = unlist(params.df[i,1:length(param.names.vec)], use.names=FALSE)
    ## cat(names(params.df)[1:length(param.names.vec)]); cat("\n")
    ## print(pars)
    ## phy = tree.bisse(pars, sim.time, max.taxa=10000, include.extinct=FALSE, x0=NA)
    
    if (length(param.names.vec) == 6) {
        phy = tree.bisse(pars, sim.time, max.taxa=10000, include.extinct=FALSE, x0=NA)
    } else {
        phy = tryCatch(tree.classe(pars, sim.time, max.taxa=10000, include.extinct=FALSE, x0=NA),
                       error = function(e) {
                           cat("tree.classe() bombed\n."); return(NULL)
                       }
                       )
    }

    if (!is.null(phy)) {
        ## cat(paste0("N tips=", length(phy$tip.state), "\n"))
        if (length(phy$tip.state) > 1000) { too.large = too.large + 1; next }
        if (length(phy$tip.state) < 30) { too.small = too.small + 1; next }

        simulated.trees = simulated.trees + 1
        cat(paste0("Simulated ",simulated.trees," trees.\r"))
        params.df[i,"tree"] = write.tree(phy)
        params.df[i,"ntips"] = length(phy$tip.state)
        params.df[i,"tipstates"] = paste(
            paste(names(phy$tip.state), phy$tip.state + 1, sep="="), # need to index by 1 
            collapse=",")

        if (simulated.trees == 100) {
            plot.ntips(params.df[!is.na(params.df$"tree"),]) # plotting ntips.pdf
            break
        }
    }
}
cat(paste0("Threw away ", too.large, " too large simulations.\n"))
cat(paste0("Threw away ", too.small, " too small simulations.\n"))
write.table(params.df[!is.na(params.df$"tree"),], file=paste0(output.dir, "data_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
write.table(params.df[is.na(params.df$"tree"),], file=paste0(output.dir, "discarded_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
