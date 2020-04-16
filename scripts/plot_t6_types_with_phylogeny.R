#!/usr/bin/env Rscript

library(ggplot2)
#library(ggtree)
library(phytools)
library(dplyr)

####  colours for T6SS subtypes:
cols <- c("i1" = "#3cb44b",
          "i2" = "#ffe119",
          "i3" = "#e6194b",
          "i4a" = "#4363d8",
          "i4b" = "#ff1493",
          "i5" = "#911eb4",
          "0" = "black",
          "Undetermined" = "grey",
          "Genes num < threshold" = "grey",
          "Broken over contig" = "grey")



#########================= 1 - read arguments =================#######

#### filter input - feed the whole tree to be used:
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args) > 1) {
  node_to_root <- args[2]
} else if (length(args) > 2) {
  stop("Too many arguments supplied", call.=FALSE)
}

#########================= 2 - read in tree =================#######


input_tree <- args[1]

#input_tree <- "~/Documents/202002-Feb/20200214-BAPS_clustering_panaroo_serratia_alignment/panaroo_serratia.newick"

strains_tree <- read.tree(input_tree)

if (is.null(strains_tree)) {
  stop("Please provide a treefile for plotting", call.=FALSE)
}

if (length(args)== 2) {
  rooted_strains_tree <- root(strains_tree, outgroup = gsub("#","_",args[2]))
} else if (length(args) == 1) {
  rooted_strains_tree <- strains_trees
}

g1 <- ggtree(rooted_strains_tree, ladderize = T, right = T)


#########================= 3 - read in T6SS subtype data  =================#######

t6_subtypes <- read.csv("T6SS_cluster_types.csv",
                        header = T,
                        quote = "",
                        stringsAsFactors = F,
                        comment.char = "")

## fix hashes in the strain:
t6_subtypes$strain <- gsub("#","_",t6_subtypes$strain)



### test plotting first:
rownames(t6_subtypes) <- t6_subtypes$strain
#remove strain column:
t6_subtypes <- t6_subtypes %>% select(-strain)


g1 %>% gheatmap(t6_subtypes, color = NULL,
                colnames_position = "top",
                colnames_offset_y = 20,
                width = 0.18,
                colnames_angle = 90,
                offset = 0,
                hjust = 0) +
  scale_fill_gradient(low = "white", high = "blue") +
  ylim(0,length(strains_tree$tip.label)+100)

#### first take all of the T6SSs:

all_t6s <- t6_subtypes %>% select(predicted_T6SSs)
colnames(all_t6s) <- c("All T6SSs")

###sort colours :
max_high <- max(all_t6s)
high_cols <- colorRampPalette(c("white","#0000ff"))
highcols<- high_cols(max_high+1)
highn <- seq(0,max_high,1)
highn <- paste("h",highn,sep='')
names(highcols) <- highn
cols <- append(cols,highcols)


###sort out the number s again so that they are factors with h infront of them:
all_t6s$`All T6SSs` <- paste0("h",all_t6s$`All T6SSs`,sep="")


# plot high
t1 <- g1 %>% gheatmap(all_t6s, color = NULL,
                colnames_position = "top",
                colnames_offset_y = 20,
                width = 0.18,
                colnames_angle = 0,
                offset = 0,
                hjust = 0.5) +
  #scale_fill_gradient(low="white", high="#ff8000") +
  scale_fill_manual(values = c(cols),name="") +
  theme(legend.position = "right")


### now get the rest - minus the ones at the end:

specific_t6_subtypes <- t6_subtypes %>% select(-predicted_T6SSs,-Undetermined,-rejected_T6SSs,-contig_broken_T6SS)

####
for (name in colnames(specific_t6_subtypes)) {
  max_col <- max(specific_t6_subtypes[[name]])
  cols_ramp <- colorRampPalette(c("white",cols[name]))
  new_cols<- cols_ramp(max_col+1)
  col_seq <- seq(0,max_col,1)
  col_seq <- paste(name,"_",col_seq,sep='')
  names(new_cols) <- col_seq
  cols <- append(cols,new_cols)

  ### now put the prefix in the gheatmap....:

  specific_t6_subtypes[[name]] <- paste0(name,"_",specific_t6_subtypes[[name]],sep="")
}

### plot this
t2 <- t1 %>% gheatmap(specific_t6_subtypes, color = NULL,
                colnames_position = "top",
                colnames_offset_y = 20,
                width = 0.18*ncol(specific_t6_subtypes),
                colnames_angle = 0,
                offset = 0.15,
                hjust = 0.5) +
  #scale_fill_gradient(low="white", high="#ff8000") +
  scale_fill_manual(values = c(cols),name="") +
  theme(legend.position = "right")

### now do the undetermined ones:
scraps <- t6_subtypes %>% select(Undetermined,rejected_T6SSs,contig_broken_T6SS)
colnames(scraps) <- c("Undetermined","Genes num < threshold","Broken over contig")

for (name in colnames(scraps)) {
  max_col <- max(scraps[[name]])
  cols_ramp <- colorRampPalette(c("white",cols[name]))
  new_cols<- cols_ramp(max_col+1)
  col_seq <- seq(0,max_col,1)
  col_seq <- paste("grey_",col_seq,sep='')
  names(new_cols) <- col_seq
  cols <- append(cols,new_cols)

  ### now put the prefix in the gheatmap....:

  scraps[[name]] <- paste0("grey_",scraps[[name]],sep="")
}

t2 %>% gheatmap(scraps, color = NULL,
                      colnames_position = "top",
                      colnames_offset_y = 20,
                      width = 0.18*ncol(scraps),
                      colnames_angle = 90,
                      offset = 0.15 * ncol(specific_t6_subtypes),
                      hjust = 0) +
  #scale_fill_gradient(low="white", high="#ff8000") +
  scale_fill_manual(values = c(cols),name="") +
  theme(legend.position = "right") +
  ggsave("tree_w_T6SS_types.png",dpi = 300, height = 10, width = 20)
