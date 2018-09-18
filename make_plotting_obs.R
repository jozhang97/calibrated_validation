## args <- commandArgs(trailingOnly=TRUE)
## beast.output.folder <- args[1]
## burnin.n <- args[2]
## csv.folder <- args[3]

library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(stats)
library(gridExtra)

# BiSSE
true.names <- unlist(strsplit("l0,l1,m0,m1,q01,q10", split=","))
beast.names <- unlist(strsplit("Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix1,FlatQMatrix2", split=","))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names), paste0("lower",beast.names), paste0("upper",beast.names)))
names(name.df) <- c("true.names", "beast.names", "beast.lower", "beast.upper")
csv.dir <- "/home/fkur465/Documents/uoa/calibrated_validation/csvs_plots/"

# ClaSSE
true.names <- unlist(strsplit("l_111,l_313,l_312,m1,m2,m3,q01,q02,q10,q12,q20,q21", split=","))
beast.names <- unlist(strsplit("SympatricRate,SubsympatricRate,VicariantRate,Mu1,Mu2,Mu3,FlatQMatrix1,FlatQMatrix2,FlatQMatrix3,FlatQMatrix4,FlatQMatrix5,FlatQMatrix6", split=","))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names), paste0("lower",beast.names), paste0("upper",beast.names)))
names(name.df) <- c("true.names", "beast.names", "beast.lower", "beast.upper")
csv.dir <- "/home/fkur465/Documents/uoa/calibrated_validation/csvs_plots_classe/"

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

## large.idxs <- true.df[,"ntips"]>=150
## true.df <- true.df[large.idxs,]
## df <- df[large.idxs,]

## small.idxs <- true.df[,"ntips"]<=10
## true.df <- true.df[small.idxs,]
## df <- df[small.idxs,]

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

write.table(file=paste0(csv.dir, "one_arbitrary_simulation.txt"), true.df[true.df$ntips==250,c("tree","ntips", "tipstates")], quote=FALSE, row.names=FALSE, col.names=FALSE)

## CLaSSE + HKY
seqgen.dir <- "/home/fkur465/Documents/uoa/calibrated_validation/seqgen_seqs/"
log.df <- read.table(paste0(seqgen.dir, "seq0.log"), header=TRUE)

brewer.pal(8, "Set1")

lighten <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col*factor
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
    col
}

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

make.post.plot <- function(a.df, var.col, param.name, x.min.max.vec, a.color.fill, a.color.line, a.truth) {
    post.samples.plot = ggplot(a.df, aes(x=a.df[,var.col])) +
        geom_histogram(aes(y=stat(density)), color=a.color.line, fill=a.color.fill, alpha=.01, bins=100, size=2) +
        geom_histogram(aes(y=stat(density)), fill=a.color.fill, bins=100) +
        geom_vline(xintercept=a.truth, color="red") +
        scale_x_continuous(limits=x.min.max.vec) +
        xlab(param.name) + ylab("Density") + 
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
    return(post.samples.plot)
}

pal <- brewer.pal(8, "Set2")
pal <- colorRampPalette(pal)(9)

truth <- true.df[true.df$"ntips"==250,]$"l_111"*35
s.min.max <- quantile(log.df$SympatricRate, c(.005,.9999))
fill.color = lighten(pal[1], 1.3)
s.plot <- make.post.plot(log.df, 6, "Sympatric Rate", s.min.max, fill.color, pal[1], truth)
s.plot

truth <- true.df[true.df$"ntips"==250,]$"l_313"*35
ss.min.max <- quantile(log.df$SubsympatricRate, c(.001,.9999))
fill.color = lighten(pal[2], 1.4)
ss.plot <- make.post.plot(log.df, 7, "Subsympatric Rate", ss.min.max, fill.color, pal[2], truth)
ss.plot

truth <- true.df[true.df$"ntips"==250,]$"l_312"*35
v.min.max <- quantile(log.df$VicariantRate, c(.00001,.999))
fill.color = lighten(pal[3], 1.3)
v.plot <- make.post.plot(log.df, 8, "Vicariant Rate", v.min.max, fill.color, pal[3], truth)
v.plot

truth <- true.df[true.df$"ntips"==250,]$"q01"
q01.min.max <- quantile(log.df$FlatQMatrix1, c(.00001,.999))
fill.color = lighten(pal[4], 1.3)
q01.plot <- make.post.plot(log.df, 12, "q12", q01.min.max, fill.color, pal[4], truth)
q01.plot

truth <- true.df[true.df$"ntips"==250,]$"q02"
q02.min.max <- quantile(log.df$FlatQMatrix2, c(.00001,.999))
fill.color = lighten(pal[5], 1.3)
q02.plot <- make.post.plot(log.df, 13, "q13", q02.min.max, fill.color, pal[5], truth)
q02.plot

truth <- true.df[true.df$"ntips"==250,]$"q10"
q10.min.max <- quantile(log.df$FlatQMatrix3, c(.00001,.999))
fill.color = lighten(pal[6], 1.3)
q10.plot <- make.post.plot(log.df, 14, "q21", q10.min.max, fill.color, pal[6], truth)
q10.plot

truth <- true.df[true.df$"ntips"==250,]$"q12"
q12.min.max <- quantile(log.df$FlatQMatrix4, c(.00001,.999))
fill.color = lighten(pal[7], 1.3)
q12.plot <- make.post.plot(log.df, 15, "q13", q12.min.max, fill.color, pal[7], truth)
q12.plot

truth <- true.df[true.df$"ntips"==250,]$"q20"
q20.min.max <- quantile(log.df$FlatQMatrix5, c(.00001,.99999))
fill.color = lighten(pal[8], 1.3)
q20.plot <- make.post.plot(log.df, 16, "q31", q20.min.max, fill.color, pal[8], truth)
q20.plot

truth <- true.df[true.df$"ntips"==250,]$"q21"
q21.min.max <- quantile(log.df$FlatQMatrix6, c(.00001,.99999))
fill.color = lighten(pal[9], 1.3)
q21.plot <- make.post.plot(log.df, 17, "q32", q21.min.max, fill.color, pal[9], truth)
q21.plot 

grid.arrange(s.plot, ss.plot, v.plot,
             q01.plot, q02.plot,
             q10.plot, q12.plot,
             q20.plot, q21.plot)
