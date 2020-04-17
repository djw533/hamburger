

dependencies_check = data.frame(matrix(nrow = 0, ncol = 1))
colnames(dependencies_check) <- c("package")


# ipak <- function(pkg){
#     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#     if (length(new.pkg)) {
#         print(paste0(new.pkg," not installed. Please install before running this script.",sep=""))
#         check_temp <- data.frame(c(new.pkg))
#         colnames(check_temp) <- colnames(dependencies_check)
#         dependencies_check <- rbind(dependencies_check,check_temp)
#       }
# }


# usage
packages <- c("ggplot2", "dplyr", "gggenes", "RColorBrewer", "castor", "ape", "ggtree","glue")
# ipak(packages)

for (pkg in packages) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
      #print(paste0(new.pkg," not installed. Please install before running this script.",sep=""))
      check_temp <- data.frame(c(new.pkg))
      colnames(check_temp) <- colnames(dependencies_check)
      dependencies_check <- rbind(dependencies_check,check_temp)
    }
}


write.csv(file = "dependencies_check.csv", dependencies_check, row.names = F, quote = F)
