library(ape)
files <- list.files(path="/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/sim_data", pattern="*tips.csv", 
                    full.names=TRUE, recursive=FALSE)


rows <- lapply(files, function(x) {
        t <- read.csv(x, header=FALSE) # load file
        NROW(t)
})


tips = as.numeric(as.vector(rows))
tips.freq = table(tips) / sum(tips)
print(tips.freq)
hist(tips, nclass = 100, xlim = c(0, 100), prob = TRUE)
hist(tips, nclass = 100, xlim = c(0, 1000), prob = TRUE)
print(mean(tips))
print(median(tips))