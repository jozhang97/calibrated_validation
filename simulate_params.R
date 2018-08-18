library(ggplot2)
library(gridExtra)
library(sjPlot)
library(diversitree)
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)
output.dir <- args[1]
prefix <- args[2]
n.sim <- args[3]
mu <- as.numeric(args[4])
std <- as.numeric(args[5])
sim.time <- as.numeric(args[6])
param.names <- args[7:length(args)]

# --- START: Prior sampling stuff --- #
sample.parameter <- function(output.dir, prefix, n.sim, mu, std, param.name) {
    # Setting up arbitrary simulation values, and making objects
    ps <- rlnorm(20000, meanlog=mu, sdlog=std) # plotting samples: default mean and sd in log scale, and =0 and 1, respectively
    s <- rlnorm(n.sim, meanlog=mu, sdlog=std) # actual samples we take from the prior
    p.df <- data.frame(ps); names(p.df) <- "value" # df for curve
    df <- data.frame(s); names(df) <- "value" # samples df
    r <- range(ps) # min and max from samples, for plotting
    xaxis.max <- 3 # for plotting

    prior.samples.plot <- ggplot(df, aes(x=value)) +
        geom_histogram(aes(y=stat(density)), alpha=.4, bins=100) +
        ggtitle(param.name) +
        stat_function(fun=dlnorm,
                      args=list(r[1]:r[2], meanlog=mean(log(p.df$value)), sdlog=sd(log(p.df$value)))
                      ) +
        scale_x_continuous(breaks=seq(0, xaxis.max, by=1), limits=c(0, xaxis.max)) +
        xlab("Parameter value") + ylab("Density") + 
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

    ## prior.samples.plot # uncomment and run to see graph
    return(list(prior.samples.plot, df))
}

sample.parameters <- function(output.dir, prefix, n.sim, mu, std, param.names, plot.flag) {
    num_params = length(param.names)
    params = vector("list", length = num_params)
    plots = vector("list", length = num_params)

    for (i in 1 : num_params) {
        param.name = param.names[i]
        res = sample.parameter(output.dir, prefix, n.sim, mu, std, param.name)
        plots[[i]] = res[[1]]
        params[[i]] = res[[2]]
    }

    params.df = as.data.frame(params)
    names(params.df) = param.names
    save.path = paste(paste0(output.dir, prefix), "all_params.csv", sep="_")

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

params.df <- sample.parameters(output.dir, prefix, n.sim, mu, std, param.names, TRUE) # table with all true parameters for all simulations
init.params.df <- sample.parameters(output.dir, prefix, min(100, n.sim), mu, std, param.names, FALSE) # table with all initialization (for .xml) parameters for all simulations
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
        )

    ggsave(paste0(output.dir, "ntips.pdf"), plot=res, width=6, height=3)
}

# BiSSE
params.df$tree <- NA
params.df$ntips <- 0
params.df$tipstates <- NA
simulated.trees <- 0
for (i in 1:nrow(params.df)) {
    pars = unlist(params.df[i,1:6], use.names=FALSE)
    phy = tree.bisse(pars, sim.time, max.taxa=1000, include.extinct=FALSE, x0=NA)

    if (!is.null(phy)) {        
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
write.table(params.df[!is.na(params.df$"tree"),], file=paste0(output.dir, "data_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
write.table(params.df[is.na(params.df$"tree"),], file=paste0(output.dir, "discarded_param_tree.csv"),
          row.names=FALSE, quote=FALSE, sep="|")
