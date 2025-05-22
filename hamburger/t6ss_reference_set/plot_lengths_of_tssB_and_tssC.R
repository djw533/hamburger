
library(ggplot2)

tssB_lengths <- read.table("tssB_lengths.txt",
                           sep="\t",
                           header = F)

tssC_lengths <- read.table("tssC_lengths.txt",
                           sep="\t",
                           header = F)


ggplot(tssB_lengths, aes(V2)) + geom_histogram(binwidth = 1) + 
  theme_bw() +
  xlab("length") +
  ggtitle("TssB lengths") +
  ggsave("tssB_lengths.png")


ggplot(tssC_lengths, aes(V2)) + geom_histogram(binwidth = 1) + 
  theme_bw() +
  xlab("length") +
  ggtitle("TssC lengths")+
  ggsave("tssC_lengths.png")
