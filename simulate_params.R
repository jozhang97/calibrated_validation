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
mus <- args[4]
stds <- args[5]
param.names <- args[6] 
sim.time <- as.numeric(args[7])

mus.vec <- as.numeric(unlist(strsplit(mus, split=",")))
stds.vec <- as.numeric(unlist(strsplit(stds, split=",")))
param.names.vec <- unlist(strsplit(param.names, split=","))
print(mus.vec)
print(stds.vec)
print(param.names.vec)

# --- START: Prior sampling stuff --- #
get.95 <- function(a.vector) {
    res = hdi(a.vector, cresMass=.95)
    return(c(res[[1]], res[[2]]))
}

sample.parameter <- function(output.dir, prefix, n.sim, mu, std, param.name) {
    # Setting up arbitrary simulation values, and making objects
    ps <- rlnorm(1000, meanlog=mu, sdlog=std) # plotting samples: default mean and sd in log scale, and =0 and 1, respectively
    s <- rlnorm(n.sim, meanlog=mu, sdlog=std) # actual samples we take from the prior
    p.df <- data.frame(ps); names(p.df) <- "value" # df for curve
    df <- data.frame(s); names(df) <- "value" # samples df
    r <- range(ps) # min and max from samples, for plotting
    xaxis.max <- 3 # for plotting

    print(paste0("Drawing param ", param.name, " with mean ", mu, " and stdev ", std))
    print(mean(ps))
    ## after ggtitle
    ## stat_function(fun=dlnorm,
    ##           args=list(r[1]:r[2], meanlog=mean(log(p.df$value)), sdlog=sd(log(p.df$value)))
    ##           ) +

    
    prior.samples.plot <- ggplot(df, aes(x=value)) +
        geom_histogram(aes(y=stat(density)), alpha=.4, bins=100) +
        scale_x_continuous(breaks=seq(0, xaxis.max, by=1), limits=c(0, xaxis.max)) +
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

    ## ggsave(file=paste0(output.dir, param.name, ".pdf"), plot=prior.samples.plot, width=3, height=3)
    ## prior.samples.plot # uncomment and run to see graph
    return(list(prior.samples.plot, df))
}

sample.parameters <- function(output.dir, prefix, n.sim, mus.vec, stds.vec, param.names.vec, plot.flag) {
    num_params = length(param.names.vec)
    params = vector("list", length = num_params)
    plots = vector("list", length = num_params)

    for (i in 1 : num_params) {
        param.name = param.names.vec[i]
        mu = mus.vec[i]
        std = stds.vec[i]
        res = sample.parameter(output.dir, prefix, n.sim, mu, std, param.name)
        plots[[i]] = res[[1]]
        params[[i]] = res[[2]]
    }

    params.df = as.data.frame(params)
    names(params.df) = param.names.vec

    # printing all panels in single graph (png() won't work)
    if (plot.flag) {
        cat("Plotting prior_samples.pdf\n")
        pdf(paste0(output.dir, "prior_samples.pdf"), width=6, height=7)
        plot_grid(plots)
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

params.df <- sample.parameters(output.dir, prefix, n.sim, mus.vec, stds.vec, param.names.vec, TRUE) # table with all true parameters for all simulations
init.params.df <- sample.parameters(output.dir, prefix, min(100, n.sim), mus.vec, stds.vec, param.names.vec, FALSE) # table with all initialization (for .xml) parameters for all simulations
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
        ) + xlim(0,750)

    ggsave(paste0(output.dir, "ntips.pdf"), plot=res, width=6, height=3)
}

# BiSSE
## initializing cols in df
params.df$tree <- NA
params.df$ntips <- 0
params.df$tipstates <- NA

## getting prior HDIs
many.samples.from.priors <- vector("list", length(param.names.vec))
for (i in 1:length(param.names.vec)) {
    many.samples.from.priors[[i]] <- rlnorm(20000, meanlog=mus.vec[i], sdlog=stds.vec[i])
    prior.hdi = get.95(many.samples.from.priors[[i]])
    lower.colname = paste0("hdilower.", param.names.vec[i])
    upper.colname = paste0("hdiupper.", param.names.vec[i])
    params.df = cbind(params.df, prior.hdi[[1]])
    params.df = cbind(params.df, prior.hdi[[2]])
    names(params.df)[c(ncol(params.df)-1,ncol(params.df))] = c(lower.colname, upper.colname)
}

## simulating, storing and printing
too.large <- 0
simulated.trees <- 0
for (i in 1:nrow(params.df)) {
    pars = unlist(params.df[i,1:6], use.names=FALSE)
    phy = tree.bisse(pars, sim.time, max.taxa=10000, include.extinct=FALSE, x0=NA)

    if (!is.null(phy)) {
        if (length(phy$tip.state) > 750) { too.large = too.large + 1; next }
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
    # else { cat("died\n") }
}
cat(paste0("Threw away ", too.large, " simulations.\n"))
write.table(params.df[!is.na(params.df$"tree"),], file=paste0(output.dir, "data_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
write.table(params.df[is.na(params.df$"tree"),], file=paste0(output.dir, "discarded_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
