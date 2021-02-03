#!/usr/bin/env Rscript


library(ggtree)
library(ape)
library(ggplot2)
library(castor)
library(gggenes)
library(dplyr)
library(RColorBrewer)

#### add in colours
cols <- c("i1" = "#3cb44b",
          "i2" = "#ffe119",
          "i3" = "#e6194b",
          "i4a" = "#4363d8",
          "i4b" = "#ff1493",
          "i5" = "#911eb4",
           "0" = "black")

tree_cols <- c("TssA" = "#3cb44b",
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
               "Non-model" = "#ffffff",
               "DUF4150" = "#808000",
               "FHA" = "red",
               "ImpE" = "blue",
               "PAAR_motif" = "pink",
               "Pkinase" = "green",
               "PP2C" = "purple",
               "TagF_N" = "yellow")


 #########================= preliminary - read arguments =================#######

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

hamburger_base_directory <- args[1]

###########=============== 1 - read in tssBC tree and associate reference set subtypes =================##############

tssBC_tree <- read.tree("tssBC_alignment.fasta.treefile")
tssBC_tree <- root(tssBC_tree, outgroup = "YP_211657.1_iii")


### load in the subtypes:
t6ss_subtypes <- read.csv(paste0(hamburger_base_directory,"/t6ss_reference_set/subtypes.csv",sep=""), # need to set this so that it works for anyone - or within this script:
                          header = T,
                          quote = "",
                          stringsAsFactors = F)

### group t6Ss subtypes for tree:
groups <- list()
i <- 0
###################### Loop to set species groups in tree
for (x in as.data.frame(table(unlist(t6ss_subtypes$subtype)))$Var1){
  i <- i + 1
  new_list <- c(subset.data.frame(t6ss_subtypes, subtype == x)$T6SS_example)
  print(x)
  groups[[x]] <- new_list
}

## add groups with OTU
groups_tree <- groupOTU(tssBC_tree, groups, overlap='abandon',connect = T)


###########============== 2 - get MRCA of the subtypes and descendents including hamburger observed T6SSs ===============##############


## get mrca and sort out the t6ss types
clades_colours_names <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(clades_colours_names) <- c("node_num","name","color_hex")

groups_mrca <- list()
i <- 0
for (i in 1:length(groups)){
  # print(i)
  x <- names(groups)[[i]]
  #print(x)
  appended <- append(groups_tree$tip.label, groups[[i]])
  new_list <- c(appended[duplicated(appended)])
  y <- new_list
  groups_mrca[[x]] <- new_list
  common_node <- getMRCA(groups_tree, new_list)
  #print(common_node)
  ### get data for the clade labeling:
  df_to_append <- data.frame(common_node,x,cols[x])
  colnames(df_to_append) <- c("node_num","name","color_hex")
  clades_colours_names <- rbind(clades_colours_names,df_to_append)
}

#### now plot the number of each t6ss subtype that they have....

#create vector list for the node numbers with each different t6ss subtype:

nodes_and_names <- clades_colours_names$node_num
names(nodes_and_names) <- clades_colours_names$name


### get descendents of the mrca for each t6ss subtype
i1_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i1"]-length(groups_tree$tip.label))$subtree)$tip.label#
i2_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i2"]-length(groups_tree$tip.label))$subtree)$tip.label#
i3_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i3"]-length(groups_tree$tip.label))$subtree)$tip.label#
i4a_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i4a"]-length(groups_tree$tip.label))$subtree)$tip.label#
i4b_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i4b"]-length(groups_tree$tip.label))$subtree)$tip.label#
i5_clusters <- (get_subtree_at_node(groups_tree, nodes_and_names["i5"]-length(groups_tree$tip.label))$subtree)$tip.label#

##### append and print these all out:
different_types <- list(i1_clusters,i2_clusters,i3_clusters,i4a_clusters,i4b_clusters,i5_clusters)
names(different_types) <- c("i1","i2","i3","i4a","i4b","i5")


all_types <- data.frame(matrix(nrow=0, ncol=2))
colnames(all_types) <- c("operon","type")

for (i in 1:length(different_types)) {
  df_to_append <- data.frame(different_types[[i]])
  t6_type_for_df <- names(different_types)[[i]]
  colnames(df_to_append) <- c("cluster")
  df_to_append$type <- t6_type_for_df

  all_types <- rbind(all_types,df_to_append)

}

### read in the statistics for merging:
### get in the gheatmap of presence / absence :
cluster_stats <- read.csv("cluster_stats.csv",
                                           comment.char = "",
                                           header = T,
                                           stringsAsFactors = F)


# get strain and cluster:
strain_and_cluster <- subset.data.frame(cluster_stats, select = c("gene_cluster","strain"))
strain_and_cluster$gene_cluster_backup <- strain_and_cluster$gene_cluster
#replace the hashes with underscores:
#strain_and_cluster$gene_cluster <- gsub("#","_",strain_and_cluster$gene_cluster)
### merge in the t6ss subtype:
strain_cluster_and_subtype <- merge(strain_and_cluster, all_types, by.x = "gene_cluster",by.y = "cluster", all.x = T)
## if NA for substype - call undetermined:
strain_cluster_and_subtype$type[is.na(strain_cluster_and_subtype$type)] <- "Undetermined"

## just take strain and subtype:
strain_and_subtype <- subset.data.frame(strain_cluster_and_subtype, select = c("strain", "type"))

## change to a table:
strain_and_subtype <- as.data.frame.matrix(table(strain_and_subtype))
strain_and_subtype <- data.frame(strain = rownames(strain_and_subtype), strain_and_subtype)

### now also add in total observed T6SSs, number rejected, and number rejected due to contig breaks
strain_stats <- read.csv("strain_statistics.csv",
                         comment.char = "",
                         header = T,
                         stringsAsFactors = F)

colnames(strain_stats) <- c("strain","predicted_T6SSs","rejected_T6SSs","contig_broken_T6SS")


strain_stats_ordered <- merge(strain_stats, strain_and_subtype, by.x = "strain", by.y = "strain", all.x = T)
strain_stats_ordered[is.na(strain_stats_ordered)] <- 0

### move some to the end:
strain_stats_ordered <- strain_stats_ordered%>%select(-rejected_T6SSs,rejected_T6SSs)
strain_stats_ordered <- strain_stats_ordered%>%select(-contig_broken_T6SS,contig_broken_T6SS)


### add the cluster type to the cluster_stats :
write.csv(file = "T6SS_cluster_types.csv",
          strain_stats_ordered,
          row.names = F,
          quote = F)


###### now add the subtype information back to the specific cluster information:
cluster_and_subtype <- subset.data.frame(strain_cluster_and_subtype, select=c("gene_cluster_backup","type"))
colnames(cluster_and_subtype) <- c("gene_cluster","T6SS_subtype")
cluster_stats <- merge(cluster_stats, cluster_and_subtype, by.x = "gene_cluster",by.y ="gene_cluster")

###write this back out and overwrite:
write.csv(file = "cluster_stats_with_subtype.csv",
          cluster_stats,
          row.names = F,
          quote = F)



###########============== 3 - plot tree of tssBC with the subtypes  ===============##############

### plot tree with labeled clades
g1 <- ggtree(groups_tree, ladderize = T, right = T, size=1, aes(color=group)) +
  scale_color_manual(values = c(cols)) +
  geom_treescale(x=0, y=0, offset = 1, width = 0.1, linesize = 0.5) +
  geom_cladelabel(node=clades_colours_names[1,]$node_num, label=clades_colours_names[1,]$name,
                  color=clades_colours_names[1,]$color_hex, offset=0,align = T) +
  geom_cladelabel(node=clades_colours_names[2,]$node_num, label=clades_colours_names[2,]$name,
                  color=clades_colours_names[2,]$color_hex, offset=0,align = T) +
  geom_cladelabel(node=clades_colours_names[3,]$node_num, label=clades_colours_names[3,]$name,
                  color=clades_colours_names[3,]$color_hex, offset=0,align = T) +
  geom_cladelabel(node=clades_colours_names[4,]$node_num, label=clades_colours_names[4,]$name,
                  color=clades_colours_names[4,]$color_hex, offset=0,align = T) +
  geom_cladelabel(node=clades_colours_names[5,]$node_num, label=clades_colours_names[5,]$name,
                  color=clades_colours_names[5,]$color_hex, offset=0,align = T) +
  geom_cladelabel(node=clades_colours_names[6,]$node_num, label=clades_colours_names[6,]$name,
                  color=clades_colours_names[6,]$color_hex, offset=0,align = T) +
  geom_rootedge(rootedge = 0.05) +
  ggsave("tssBC_tree_with_types.png",dpi = 300)


###########============== 4 - Draw tssBC tree with coloured gene clusters ===========##############

###draw the tree with the gggenes input

ref_set_genes <- read.csv(paste0(hamburger_base_directory,"/t6ss_reference_set/gggenes_input.csv",sep=""))
observed_set_genes <- read.csv("gggenes_input.csv")
## remove all hashes from the operon name:
observed_set_genes$operon <- gsub("#","_",observed_set_genes$operon)

##rbind these into a single set:
combined_genes <- rbind(ref_set_genes, observed_set_genes)

# now create a df with all genes in the same orientation and aligned to tssE:
##### prepare to plot all of these:
centred_complete_data <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(centred_complete_data) <- c("operon","number","start","end","gene","strand","direction")
average_position <- as.integer(mean(subset.data.frame(combined_genes, gene == "TssE")$start)) # get average starting position for "Ketoacyl-synt_C"


### add operons:
for (operon_num in unique(combined_genes$operon)) {

  temp_df <- subset.data.frame(combined_genes, operon == operon_num, select = c("operon","number","start","end","gene","strand","direction"))
  tssE_direction <- subset.data.frame(temp_df, gene == "TssE")$direction
  #print(tssE_direction)
  #
  if (nrow(subset.data.frame(temp_df, gene == "TssE")) > 0 ) {
    if (tssE_direction == "-1") {
      ##### need to reverse the numbers and the directions
      max_length <- max(temp_df$end)
      ## minus the length from the start and stop
      temp_df$start <- max_length - temp_df$start
      temp_df$end <- max_length - temp_df$end
      ## change colnames to switch the start and the stop over
      colnames(temp_df) <- c("operon","number","end","start","gene","strand","direction")
      ## change direction:
      temp_df$direction <- temp_df$direction * -1
    }

    # ##### correct the position of each of these:
    difference <- average_position + subset.data.frame(temp_df, gene == "TssE")$start
    temp_df$start <- temp_df$start - difference
    temp_df$end <- temp_df$end - difference

  }
  ## concatenate to new df
  centred_complete_data <- rbind(centred_complete_data ,temp_df)
}


###### plot the operons to the right of the tree:
max_length <- ((max(centred_complete_data$end) - min(centred_complete_data$start) ) / 2) / 1000


facet1 <- facet_plot(g1, panel='Operons',
           mapping = aes(xmin = start, xmax = end, fill = gene, x=start, forward = direction, y=y),
           data=centred_complete_data, color = "black", geom=geom_gene_arrow, size = 1,#, 0.0025 * length(unique(centred_complete_data$operon)) , #use this if want to scale the line thickness on the gene arrows according to the number of tips in the tssBC tree
           arrow_body_height = grid::unit(0.20, "cm"), arrowhead_height = grid::unit(0.25, "cm"), arrowhead_width = grid::unit(0.25,"cm")) +
  scale_fill_manual(values=c(tree_cols)) +
  # geom_treescale(width = 2, offset = -1.2) +
  ggsave("tssBC_tree_with_T6SS_operons.png",dpi = 300, units = c("cm"), width = max_length * 2, height = 0.3 * length(unique(centred_complete_data$operon)), limitsize = F) +
  ggsave("tssBC_tree_with_T6SS_operons.pdf",dpi = 300, units = c("cm"), width = max_length * 2, height = 0.3 * length(unique(centred_complete_data$operon)), limitsize = F)
