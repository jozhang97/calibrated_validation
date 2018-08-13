library(ape)
library(ggplot2)

tip_analysis <- function (path) {
    files <- list.files(path=path, pattern="*tips.csv", full.names=TRUE, recursive=FALSE)

    rows <- lapply(files, function(x) {
        t <- read.csv(x, header=FALSE) # load file
        NROW(t)
    })

    tips = as.numeric(as.vector(rows))

    # TODO ask fabio to help with ggplot
    hist(tips, nclass = 10, xlim = c(0, 1000), prob = TRUE)
    hist(tips, nclass = 10, xlim = c(0, 100), prob = TRUE)

    hist_path = paste0(path, "tips_hist_HI.png")
    print(hist_path)
    png(hist_path)
    dev.off()

    num.inside = sum(5 <= tips & tips < 1000)
    num.outside = sum(tips < 5 | tips > 1000)
    num.total = length(tips)
    print("Inside/Outside")
    print(1.0 * num.inside / num.total)
    print(1.0 * num.outside / num.total)
    print("Mean/Median")
    print(mean(tips))
    print(median(tips))
}

args = commandArgs(trailingOnly=TRUE)

dir_path = args[1]
tip_analysis(dir_path)
