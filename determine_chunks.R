for(package in c("data.table", "rjson")) {
    if (!require(package, character.only = TRUE)) {
        print(paste0("installing ", package))
        install.packages(package, repos='http://cran.us.r-project.org')
        library(package)
    }
}
## install.packages("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.20.tar.gz")
library("parallel")


ANALYSIS_DIR <- Sys.getenv("ANALYSIS_DIR")
WINDOWSIZE <- as.integer(Sys.getenv("WINDOWSIZE"))
BUFFER <- as.integer(Sys.getenv("BUFFER"))

CHRLIST <- 1:22

if (1 == 0) {

    ANALYSIS_DIR <- "/well/ansari/shared/lcwgs/results/2022_11_02/"
    WINDOWSIZE <- 5000000
    BUFFER <- 500000

}


data <- mclapply(
    CHRLIST,
    mc.cores = 4,
    function(chr) {

        
        ## data <- fread(
        ##     cmd = paste0("gunzip -c ", file, " | grep -v '#' | cut -f2"),
        ##     data.table = FALSE
        ## )

        ##file <- file.path(ANALYSIS_DIR, "refs", paste0("oneKG.chr", chr, ".vcf.gz"))
        ##start <- as.integer(system(paste0("bcftools query -f '%POS\n' ", file, " | head -n 1"), intern = TRUE))
        ##end <- as.integer(system(paste0("bcftools query -f '%POS\n' ", file, " | tail -n 1"), intern = TRUE))

        file <- file.path(ANALYSIS_DIR, "refs", paste0("oneKG.chr", chr, ".legend.gz"))
        data <- fread(
            cmd = paste0("gunzip -c ", file), ##  | grep -v '#' | cut -f2"
            data.table = FALSE
        )

        ##start <- as.integer(strsplit(system(paste0("gunzip -c ", file, " | head -n2 | tail -n 1"), intern = TRUE), " ")[[1]][2])
        ## end <- as.integer(strsplit(system(paste0("gunzip -c ", file, " | tail -n 1"), intern = TRUE), " ")[[1]][2])

        ##startI <- floor(start / WINDOWSIZE)
        ##endI <- ceiling(end / WINDOWSIZE)

        ## This alread checks for them in there
        range <- unique(floor(data[, "position"] / WINDOWSIZE))

        startV <- pmax((range) * WINDOWSIZE + 1 - BUFFER, data[1, "position"])
        endV <- pmin((range + 1) * WINDOWSIZE + BUFFER, data[nrow(data), "position"])

        ## check they have SNPs
        ## t <- table(floor(data[, "position"] / WINDOWSIZE))
        ## keep <- !is.na(match(min(range):(max(range) - 1), as.integer(names(t))))
        ## stopifnot(length(keep) == length(startV))
        ## startV <- startV[keep]
        ## endV <- endV[keep]
        
        list(
            start = startV,
            end = endV
        )
            
    }
)

print(data)

names(data) <- CHRLIST


write(toJSON(data), file.path(ANALYSIS_DIR, "regions.json"))

