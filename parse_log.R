library(HDInterval)
library(naturalsort)

args <- commandArgs(trailingOnly=TRUE)
beast.output.folder <- args[1]
burnin.n <- args[2]
csv.folder <- args[3]

get.95 <- function(a.vector) {
    res = hdi(a.vector, cresMass=.95)
    return(c(res[[1]], res[[2]]))
}

parse_log <- function(beast.output.folder, burnin.n, csv.folder) {
    header <- ""
    burnin.n <- 500000
    all.dfs <- vector("list", 100)
    summaries <- vector("list", 100)

    count <- 1
    for (i in list.files(beast.output.folder)) {
        if (endsWith(i, ".log")) {
            df = read.table(paste0(beast.output.folder,i), header=TRUE)
	    w.o.burnin = (burnin.n/5000+1):nrow(df)
	    summaries[[count]] = lapply(df[w.o.burnin,5:ncol(df)], get.95)
	    if (count == 1) { header = c(as.vector(outer(c("lower", "upper"),names(summaries[[1]]), paste0)),"file") }
	    summaries[[count]]$file = gsub("_sim.log", "", i)
	    # all.dfs[[count]] = df
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

parse_log(beast.output.folder, burnin.n, csv.folder)
