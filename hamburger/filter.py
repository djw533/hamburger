def clustering_func(numbers_list,max_gap,groupsize_min):
    ### Setting maximum gap between genes in the gene cluster

    maxgap = int(max_gap)

    ### clustering CDSs based on position and maxgap value

    numbers_list.sort()
    #print numbers_list
    groups = [[numbers_list[0]]]
    for x in numbers_list[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])

    ### filter groups based on the min_genes value

    filtered_groups = []


    for group in groups:
        if len(group) >= int(groupsize_min):
            filtered_groups.append(group)
    return filtered_groups


def same_contigs_check(gene_numbers, strain, genes_and_contig):

    contigs = []

    for num in gene_numbers:
        gene = "{strain}_{num}".format(strain = strain, num = str(num)) # covert the number back into the format in the annotation

        ### get the contig:
        contig = genes_and_contig[gene]

        contigs.append(contig)

    ### check if all the contigs are the same:
    if all(x == contigs[0] for x in contigs) == True:
        return(True, contig)

    else:
        return(False, "No contig")


def gene_names_in_cluster(gene_num_in_cluster, hmmer_output, strain):

    ###first get a dictionary with the gene number and the hmmer gene that was hit:

    ### also think - what if two hmmer queries match the same protein at some level?

    names_to_genes = {}

    for line in hmmer_output:
        gene_num = int(line.split()[0].split(strain+"_")[1]) # take the number after the strain, e.g. {strain}_0234 or {strain}_05422
        hmmer_query = line.split()[2]
        #put in dict:
        names_to_genes[gene_num] = hmmer_query


    ## cycle through the gene_numbers in the gene cluster and count number of each hmmer query type

    gene_names = {}


    for gene in gene_num_in_cluster:
        hmmer_query_name = names_to_genes[gene]

        ## add to the names_to_genes, or create key if not in the dict:
        if hmmer_query_name not in gene_names:
            gene_names[hmmer_query_name] = []

        gene_names[hmmer_query_name].append(gene)


    ### now report the number of unique gene hits found in the cluster:
    return(gene_names)
