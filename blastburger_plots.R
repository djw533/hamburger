library(ggplot2)
library(ggtree)
library(gggenes)
library(cowplot)


# args <- commandArgs(TRUE)
# tree_file <- args[1]


tree <- read.tree("../genus_tree.newick")  ## input


genes <- read.csv(file="same_direction_gggenes_input.csv", header = T, comment.char = "", quote = "")

new_data <- genes[,c(2,1,3:7)]


taxa_to_color<-c("26968_7#114","26968_7#200","26968_7#89","26968_7#170","26968_7#265","26968_7#264","26968_7#51","26968_7#4","26968_7#110","26968_7#117","26968_7#34","26968_7#211","26968_7#186","26968_7#181","26968_7#188","26968_7#199","26968_7#95","26968_7#73","26968_7#48","26968_7#47")
tree<-groupOTU(tree, taxa_to_color)


t1 <- ggtree(tree, ladderize = T, right = T) + 
  #geom_text2(aes(label=label, subset=isTip, color=group)) +
  geom_tiplab(size=2, aes(label=label, color=group)) + 
  scale_color_manual(values = c("black", "red"))
  geom_rootedge(0.2) +  
  xlim_tree(1.5) 
#+ xlim(0,2) #
#t1 <- flip(t1, 4, 3) %>% flip(7,8)


l1 <- ggplot2::ggplot(genes, ggplot2::aes(xmin = start, xmax = end,y = operon, fill = gene)) +
  geom_gene_arrow() +
  ggsave("gene_plots.svg")
legend_1 <- cowplot::get_legend(l1)

ggdraw() +
  draw_plot(legend_1, x = 0.4, y = 0.45, width = .5, height = .5) +
  ggsave("legends.png",dpi=300, width=10,height=8, units=c("in"))


facet_plot(t1, panel='Operons',
           mapping = aes(xmin = start, xmax = end, fill = gene, x=start, forward = direction, y=y),
           data=new_data, geom=geom_gene_arrow, arrow_body_height = grid::unit(1.8, "mm"), arrowhead_height = grid::unit(2.5, "mm"), arrowhead_width = grid::unit(2.5,"mm")) +
  theme(panel.grid.minor = element_line(colour="grey", size=0.1)) + 
  scale_y_continuous(minor_breaks = seq(1, 550, 2)) +
  #geom_treescale(width = 2, offset = -1.2) +
  ggsave("blastburger_output.pdf", dpi=300, width=30, height=50, units=c("in"), limitsize = F)
