args <- commandArgs(trailingOnly = TRUE)
ANALYSIS_DIR <- args[1]
WINDOWSIZE <- args[2]
BUFFER <- args[3]
PANEL_NAME <- args[4]

for(package in c("data.table", "rjson")) {
    if (!require(package, character.only = TRUE)) {
        print(paste0("installing ", package))
        install.packages(package, repos='http://cran.us.r-project.org')
        library(package)
    }
}
## install.packages("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.20.tar.gz")
library("parallel")


#ANALYSIS_DIR <- "/well/band/users/rbx225/GGVP/results/imputation/"
#WINDOWSIZE <- 5000000
#BUFFER <- 1000000
#CHRLIST <- 1:22

data <- mclapply(
    CHRLIST,
    mc.cores = 4,
    function(chr) {
        file <- file.path(ANALYSIS_DIR, "refs", paste0(PANEL_NAME, ".chr", chr, ".legend.gz"))
        data <- fread(
            cmd = paste0("gunzip -c ", file), ##  | grep -v '#' | cut -f2"
            data.table = FALSE
        )

        ## This alread checks for them in there
        range <- unique(floor(data[, "position"] / WINDOWSIZE))

        startV <- pmax((range) * WINDOWSIZE + 1 - BUFFER, data[1, "position"])
        endV <- pmin((range + 1) * WINDOWSIZE + BUFFER, data[nrow(data), "position"])

        list(
            start = startV,
            end = endV
        )

    }
)

print(data)

names(data) <- CHRLIST


write(toJSON(data), file.path(ANALYSIS_DIR, "regions.json"))
