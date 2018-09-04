## args <- commandArgs(trailingOnly=TRUE)
## beast.output.folder <- args[1]
## burnin.n <- args[2]
## csv.folder <- args[3]

library(ggplot2)
library(sjPlot)

true.names <- unlist(strsplit("l0,l1,m0,m1,q01,q10", split=","))
beast.names <- unlist(strsplit("Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix2,FlatQMatrix3", split=","))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names), paste0("lower",beast.names), paste0("upper",beast.names)))
names(name.df) <- c("true.names", "beast.names", "beast.lower", "beast.upper")
csv.dir <- "/home/fkur465/Documents/uoa/calibrated_validation/csvs_plots/"

make.regression.plot <- function(true.param.name, beast.param.mean, beast.param.lower, beast.param.upper, true.df, beast.df, hpd=FALSE) {
    true.name = as.character(true.param.name)
    beast.name = as.character(beast.param.mean)
    beast.lower = as.character(beast.param.lower)
    beast.upper = as.character(beast.param.upper)
    x = as.numeric(true.df[,names(true.df)==true.name])
    min.x = 0
    ## min.x = min(x)
    max.x = max(x)
    n = length(x)
    ## max.x = sort(x, partial=n-1)[n-2] # throwing out one outlier
    y = as.numeric(beast.df[,names(beast.df)==beast.name])
    lower = as.numeric(beast.df[,names(beast.df)==beast.lower])
    upper = as.numeric(beast.df[,names(beast.df)==beast.upper])
    ## min.y = min(lower)
    ## max.y = max(upper)
    ## min.y = min(y)
    min.y = 0
    max.y = max(y)
    reg.df = data.frame(cbind(x,y,lower,upper))
    print(reg.df)

    plot = ggplot() + geom_point(data=reg.df, mapping=aes(x=x, y=y), shape=20) + xlim(min.x,max.y) + coord_cartesian(ylim=c(min.y, max.y)) +
        xlab(paste0("Simulated ",true.name)) + ylab("Posterior mean") + geom_abline(slope=1, linetype="dotted") +
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

    if (hpd) { plot = plot + geom_linerange(data=reg.df, mapping=aes(x=x, ymax=upper, ymin=lower), color="lightgray", alpha=.4) }
    return(plot)
}

load(paste0(csv.dir, "hpds.RData"))
true.df <- read.table(paste0(csv.dir, "data_param_tree.csv"), sep="|", head=TRUE)

large.idxs <- true.df[,"ntips"]>=150
true.df <- true.df[large.idxs,]
df <- df[large.idxs,]

small.idxs <- true.df[,"ntips"]<=10
true.df <- true.df[small.idxs,]
df <- df[small.idxs,]

all.plots <- vector("list", nrow(name.df))
for (r in 1:nrow(name.df)) {
    ## print(head(true.df[,c(-7,-8,-9)]))
    ## print(head(df))
    all.plots[[r]] = make.regression.plot(name.df$true.names[r], name.df$beast.names[r], name.df$beast.lower[r], name.df$beast.upper[r], true.df, df, hpd=TRUE)
}
plot_grid(all.plots)

pdf(paste0(csv.dir, "true_vs_postmean.pdf"), width=6, height=7)
pdf(paste0(csv.dir, "true_vs_postmean_large.pdf"), width=6, height=7)
pdf(paste0(csv.dir, "true_vs_postmean_small.pdf"), width=6, height=7)
plot_grid(all.plots)
dev.off()

