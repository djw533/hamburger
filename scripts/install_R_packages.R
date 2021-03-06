

#########================= preliminary - read arguments =================#######

#### filter input - feed the whole tree to be used:
# args = commandArgs(trailingOnly=TRUE)
#
#
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(args) > 1) {
#   node_to_root <- args[2]
# } else if (length(args) > 2) {
#   stop("Too many arguments supplied", call.=FALSE)
# }
#
# hamburger_base_directory <- args[1]

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.r-project.org", dependencies = TRUE)
BiocManager::install()


ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, repos="http://cran.r-project.org", dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

bio_ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        BiocManager::install("ggtree")
    sapply(pkg, require, character.only = TRUE)
}



# usage
packages <- c("ggplot2", "dplyr", "gggenes", "RColorBrewer", "castor", "ape","glue")
ipak(packages)

bioconductor_packages <- c("ggtree")
bio_ipak(bioconductor_packages)
