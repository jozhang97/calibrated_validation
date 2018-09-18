library(HDInterval)
library(naturalsort)
library(gtools)

args <- commandArgs(trailingOnly=TRUE)
beast.output.folder <- args[1]
burnin.n <- args[2]
csv.folder <- args[3]

get.95 <- function(a.vector) {
    res = hdi(a.vector, cresMass=.95)
    return(c(res[[1]], res[[2]], mean(a.vector)))
}

parse.log <- function(beast.output.folder, burnin.n, csv.folder) {
    header = ""
    burnin.n = as.numeric(burnin.n)
    all.dfs = vector("list", 100)
    summaries = vector("list", 100)
    all.files = list.files(beast.output.folder)
    ## file.order = naturalorder(all.files)
    count <- 1
    for (i in mixedsort(all.files)) {
        if (endsWith(i, ".log")) {
            cat(paste0("Reading table ", i, "\n"))
            df = read.table(paste0(beast.output.folder,i), header=TRUE)
	    w.o.burnin = (burnin.n/5000+1):nrow(df)
	    summaries[[count]] = lapply(df[w.o.burnin,5:ncol(df)], get.95) # gets lower and upper of 95 HPD, and mean posterior
	    if (count == 1) { header = c(as.vector(outer(c("lower", "upper", "mean"),names(summaries[[1]]), paste0)),"file") }
	    summaries[[count]]$file = gsub("_sim.log", "", i)
	    all.dfs[[count]] = df
	    count <- count+1
        }
    }

    df = data.frame(matrix(unlist(summaries), nrow=100, byrow=T), stringsAsFactors=FALSE)
    names(df) = header
    or = naturalorder(df$file)
    df = df[or,]
    ## print(df[or,]) # final result
    write.csv(df, file=paste0(csv.folder, "hpds.csv"), quote=FALSE, row.names=FALSE)
    save(df, file=paste0(csv.folder, "hpds.RData"))
}

parse.log(beast.output.folder, burnin.n, csv.folder)
