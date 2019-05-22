library(ggplot2)
library(ggtree)
library(gggenes)
library(cowplot)


args <- commandArgs(TRUE)
tree_file <- args[1]


tree <- read.tree(tree_file)  ## input


genes <- read.csv(file="same_direction_gggenes.csv", header = T, comment.char = "", quote = "")

new_data <- genes[,c(2,1,3:7)]



t1 <- ggtree(tree, ladderize = T, right = T) + geom_rootedge(0.2) +  geom_tiplab(size=4) + xlim_tree(1.5) #+ xlim(0,2) #
t1 <- flip(t1, 4, 3) %>% flip(7,8)


l1 <- ggplot2::ggplot(genes, ggplot2::aes(xmin = start, xmax = end,y = operon, fill = gene)) +
  geom_gene_arrow() 
legend_1 <- cowplot::get_legend(l1)

ggdraw() +
  draw_plot(legend_1, x = 0.4, y = 0.45, width = .5, height = .5) +
  ggsave("legends.png",dpi=300, width=10,height=8, units=c("in"))


facet_plot(t1, panel='Operons',
           mapping = aes(xmin = start, xmax = end, fill = gene, x=start, forward = direction, y=y),
           data=new_data, geom=geom_gene_arrow, arrow_body_height = grid::unit(2, "mm"), arrowhead_height = grid::unit(3, "mm"), arrowhead_width = grid::unit(3,"mm")) +
  #geom_treescale(width = 2, offset = -1.2) +
  ggsave("blastburger_output.png", dpi=300, width=8,height=5, units=c("in"))
