library(ggplot2)
library(gggenes)






##colours:
t6ss_cols <- c("TssA" = "#3cb44b",
               "TssB" = "#ffe119",
               "TssC" = "#e6194b",
               "TssD" = "#4363d8",
               "TssE" = "#ff1493",
               "TssF" = "#911eb4",
               "TssG" = "#46f0f0",
               "TssH" = "#f032e6",
               "TssI" = "#bcf60c",
               "TssJ" = "#fabebe",
               "TssK" = "#008080",
               "TssL" = "#e6beff",
               "TssM" = "#9a6324",
               "PAAR_motif" = "#808000",
               "X-unknown" = "#ffffff")



genes_data <- read.csv("master_GGgenes.csv", header = T, stringsAsFactors = F, quote = "")



genes <- genes_data[,c(2,1,3:7)]


#plot first:
ggplot2::ggplot(genes, ggplot2::aes(xmin = start, xmax = end,y = operon, fill = gene, forward = direction)) +
  geom_gene_arrow(arrowhead_height = grid::unit(3, "mm"),
                  arrow_body_height = grid::unit(2,"mm"),
                  arrowhead_width = grid::unit(2, "mm")) +
 scale_fill_manual(values=c(t6ss_cols)) +
 geom_text(aes(x = (start+end)/2 , label = gene), fontface = "italic", angle = 45, hjust = 0, size = 4, nudge_y = 0.05) +
  theme_genes()
