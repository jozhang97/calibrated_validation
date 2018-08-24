## args <- commandArgs(trailingOnly=TRUE)
## beast.output.folder <- args[1]
## burnin.n <- args[2]
## csv.folder <- args[3]

library(ggplot2)
library(sjPlot)

true.names <- unlist(strsplit("l0,l1,m0,m1,q01,q10", split=","))
beast.names <- unlist(strsplit("Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix2,FlatQMatrix3", split=","))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names)))
names(name.df) <- c("true.names", "beast.names")
csv.dir <- "/home/fkur465/Documents/uoa/calibrated_validation/csvs_plots/"

make.regression.plot <- function(true.param.name, beast.param.name, true.df, beast.df) {
    true.name = as.character(true.param.name)
    beast.name = as.character(beast.param.name)
    x = true.df[,names(true.df)==true.name]
    min.x = min(x)
    max.x = max(x)
    y = as.numeric(beast.df[,names(beast.df)==beast.name])
    reg.df = data.frame(cbind(x,y))

    plot = ggplot(reg.df, aes(x=x, y=y)) + geom_point(shape=20) + geom_smooth(method=lm, se=FALSE) + xlim(min.x,max.x) + ylim(min.x,max.x) +
    xlab(paste0("Simulated ",true.name)) + ylab("Posterior mean") +
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

    return(plot)
}

load(paste0(csv.dir, "hpds.RData"))
true.df <- read.table(paste0(csv.dir, "data_param_tree.csv"), sep="|", head=TRUE)

## large.idxs <- true.df[,"ntips"]>150
## true.df <- true.df[large.idxs,]
## df <- df[large.idxs,]
# segments(x0=0,y0=0,x1=45,y1=45)

all.plots <- vector("list", nrow(name.df))
for (r in 1:nrow(name.df)) {
    ## print(head(true.df[,c(-7,-8,-9)]))
    ## print(head(df))
    all.plots[[r]] = make.regression.plot(name.df$true.names[r], name.df$beast.names[r], true.df, df)
}
plot_grid(all.plots)

pdf(paste0(csv.dir, "true_vs_postmean.pdf"), width=6, height=7)
plot_grid(all.plots)
dev.off()

