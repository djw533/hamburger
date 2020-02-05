library(ggplot2)
library(ggtree)
library(phytools)
library(castor)

#### add in colours
cols <- c("i1" = "#3cb44b",
          "i2" = "#ffe119",
          "i3" = "#e6194b",
          "i4a" = "#4363d8",
          "i4b" = "#ff1493",
          "i5" = "#911eb4",
           "0" = "black")

###load in tree
tssBC_tree <- read.tree("tssBC_alignment.fasta.treefile")
tssBC_tree <- midpoint.root(tssBC_tree)
ggtree(tssBC_tree, ladderize = T, right= T)

### load in the subtypes:
t6ss_subtypes <- read.csv("~/github/hamburger/t6ss_reference_set/subtypes.csv",
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

### plot tree with groups
ggtree(groups_tree, ladderize = T, right = T, size=0.5, aes(color=group)) +
  scale_color_manual(values = c(cols)) +
  geom_treescale(x=0, y=0, offset = 1, width = 0.1, linesize = 0.5)


## get mrca and sort out the t6ss types 
clades_colours_names <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(clades_colours_names) <- c("node_num","name","color_hex")

groups_mrca <- list()
i <- 0
###################### 
for (i in 1:length(groups)){
  # print(i)
  x <- names(groups)[[i]]
  print(x)
  appended <- append(groups_tree$tip.label, groups[[i]])
  new_list <- c(appended[duplicated(appended)])
  y <- new_list
  groups_mrca[[x]] <- new_list
  common_node <- getMRCA(groups_tree, new_list)
  print(common_node)
  ### get data for the clade labeling:
  df_to_append <- data.frame(common_node,x,cols[x])
  colnames(df_to_append) <- c("node_num","name","color_hex")
  clades_colours_names <- rbind(clades_colours_names,df_to_append)
}

### plot tree with labeled clades
g1 <- ggtree(groups_tree, ladderize = T, right = T, size=0.5, aes(color=group)) +
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
  geom_rootedge(rootedge = 0.05)



### now put these on the provided phylogeny

strains_tree <- read.tree("~/Documents/thesis/chapter_1/5_All_serratia_tr/serratia_tree_rooted.newick")

t1 <- ggtree(strains_tree, ladderize = T, right = T, size=0.5)


### get in the gheatmap of presence / absence :
hamburger_cluster_strain_stats <- read.csv("strain_statistics.csv",
                        comment.char = "",
                        header = T,
                        stringsAsFactors = F) 


hamburger_cluster_strain_stats$strain <- gsub("#","_", hamburger_cluster_strain_stats$strain) # remove hashes and replace with underscores like iqtree does

## remove the strain part and number_rejected_clusters
rownames(hamburger_cluster_strain_stats) <- hamburger_cluster_strain_stats$strain
hamburger_cluster_strain_stats <- subset.data.frame(hamburger_cluster_strain_stats, select=c("number_of_gene_clusters"))

###sort colours :
max_high <- max(hamburger_cluster_strain_stats) 
high_cols <- colorRampPalette(c("white","#0000ff"))
highcols<- high_cols(max_high+1)
highn <- seq(0,max_high,1)
highn <- paste("h",highn,sep='')
names(highcols) <- highn
cols <- append(cols,highcols)

###sort out the number s again so that they are factors with h infront of them:
hamburger_cluster_strain_stats$number_of_gene_clusters <- paste0("h",hamburger_cluster_strain_stats$number_of_gene_clusters,sep="")


# plot high
t1 %>% gheatmap(hamburger_cluster_strain_stats, color = NULL,
                                    colnames_position = "top",
                                    colnames_offset_y = 20,
                                    width = 0.18,
                                    colnames_angle = 0,
                                    offset = 0.14,
                                    hjust = 0.5) +
  #scale_fill_gradient(low="white", high="#ff8000") +
  scale_fill_manual(values = c(cols),name="") +
  theme(legend.position = "none") +
  ylim(0,750) 



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
names(different_types) <- c("i1","i2","i3","i4b","i4a","i5")


all_types <- data.frame(matrix(nrow=0, ncol=2))
colnames(all_types) <- c("operon","type")

for (i in 1:length(different_types)) {
  df_to_append <- data.frame(different_types[[i]])
  t6_type_for_df <- names(different_types)[[i]]
  colnames(df_to_append) <- c("cluster")
  df_to_append$type <- t6_type_for_df

  all_types <- rbind(all_types,df_to_append)
    
}

write.csv(file = "T6SS_cluster_types.csv",
          all_types,
          row.names = F,
          quote = F)

#### get this into an acceptable data frame

all_types$cluster <- sapply(strsplit(as.character(all_types$cluster),'_cluster_'), "[", 1)
colnames(all_types) <- c("strain","subtype")

####create df of the tips in the tree so can merge out the reference set types:
cluster_types_df <- data.frame(matrix(nrow = length(strains_tree$tip.label), ncol = 1))
colnames(cluster_types_df) <- c("strain")
cluster_types_df$strain <- strains_tree$tip.label

#merge
all_types <- merge(all_types, cluster_types_df, by.x = "strain", by.y = "strain")

###now separate
new_all_types <- as.data.frame.matrix(table(all_types))



#### Set the colours for these different subtypes:

for (name in colnames(new_all_types)) {
  max_col <- max(new_all_types[[name]]) 
  cols_ramp <- colorRampPalette(c("white",cols[name]))
  new_cols<- cols_ramp(max_col+1)
  col_seq <- seq(0,max_col,1)
  col_seq <- paste(name,"_",col_seq,sep='')
  names(new_cols) <- col_seq
  cols <- append(cols,new_cols)
  
  ### now put the prefix in the gheatmap....:
  
  new_all_types[[name]] <- paste0(name,"_",new_all_types[[name]],sep="")
}

### plot this
t1 %>% gheatmap(new_all_types, color = NULL,
                      colnames_position = "top",
                      colnames_offset_y = 10,
                      width = 0.18,
                      colnames_angle = 90,
                      offset = 0,
                      hjust = 0) +
  #scale_fill_gradient(low="white", high="#ff8000") +
  scale_fill_manual(values = c(cols),name="") +
  theme(legend.position = "none") +
  ylim(0,750) #+ ggsave("panel_d.svg")
