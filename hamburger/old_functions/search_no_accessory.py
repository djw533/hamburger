def search_single_genome_no_accessory(mandatory,accessory,min_genes,genes_gap,upstream,downstream,cutoff,t6ss,output,gff_file):

    #args = parseArgs()

    ##### repeat the name changing here....
    ### if t6ss flag requested then set these using the t6ss default setting
    # if t6ss == True:
    #     min_genes_num = 8
    #     genes_gap_num = 12
    #     mandatory_models = T6SS_core#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_core.hmm"
    # elif t6ss == False:
    #     if mandatory == False:
    #         print("Need to provide input mandatory hmm profile(s): --mandatory")
    #         sys.exit()
    #     if min_genes == False:
    #         print("Need to set the minimum number of hits: --min_genes)")
    #         sys.exit()
    #     if genes_gap == False:
    #         print("Need to set the maximum genes gap: --genes_gap)")
    #         sys.exit()
        # now change the variable names from args.xxx
    min_genes_num = min_genes
    genes_gap_num = genes_gap
    mandatory_models = mandatory



    mandatory_names = hmmer_functions.read_hmms(mandatory_models)


    number_clusters_in_strain = 0
    number_rejected_clusters_in_strain = 0
    number_contig_break_clusters_in_strain = 0


    output_dir = output

    #print(gff_file)

    ### take the strain name (last field from "_" separator and remove the .gff from the end:)

    strain = gff_file.split('/')[-1].replace(".gff","")

    ### now made directory for this strain:
    strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
    os.makedirs(strain_dir)

#3 - extract protein sequences from the gff file

### extract the protein sequences
    prot_seqs = "{output_dir}/{strain}/{strain}.faa".format(output_dir = output_dir, strain = strain)
    gff2faa_output = gff_parsing.gff2faa(gff_file,prot_seqs,strain,output_dir)
    #altered_gene_names = gff2faa_output[0] # a "yes" or "no" - is only no because this is currently set - look at later??
    list_of_genes = gff2faa_output[1] # list of the new names to the names in the gff file
    genes_and_contig = gff2faa_output[2]

#4 - create TEMPORARY singlefasta files, fasta of the whole sequence, and separate annotation files

    gff_split_output = gff_parsing.gff_splitter(gff_file, strain_dir)

    annotation = gff_split_output[0]
    GC_genome = gff_split_output[1]

    ##set singlefasta directory
    singlefasta_dir = "{strain_dir}/singlefastas".format(strain_dir = strain_dir)
    fasta_parsing.multifasta_to_singlefasta("contigs.fna", strain_dir, singlefasta_dir)

#5 - Run HMMER and filter according to the cutoff
#- for both mandatory and accessory

    mandatory_hmmer_tuple = hmmer_functions.run_hmmer(hmmsearch, mandatory_models, prot_seqs, cutoff, strain_dir, "mandatory")
    mandatory_hmmer_output = mandatory_hmmer_tuple[0]
    mandatory_hmmer_genes = mandatory_hmmer_tuple[1]


    ###### then combine the accessory and mandatory results if required.... (think about adding any genes that are not allowed in gene clusters?)

    total_hmmer_output = mandatory_hmmer_output
    total_hmmer_genes = mandatory_hmmer_genes

#5 - Cluster the HMMER hits according to the parameters

    ### introduce break if there are no hmmer hits at all!:

    if len(mandatory_hmmer_genes) == 0:
        ##no hits
        with open(output_dir+"/log_file.txt", "a") as output:
            output.write("No hmmer hits in {strain}".format(strain =strain))
        ### and write to strain_statistics.csv
        with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
            output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = "0",rejected_clusters = "0",contig_break_clusters = "0"))
        return

    filtered_groups = filter.clustering_func(total_hmmer_genes,genes_gap_num,min_genes_num)

#6 - Check that clusters are on the the same contigs


    # print(annotation) # lines of the annotation from the gff file
    # print(filtered_groups) # list of lists for each gene cluster from the (simple) algorithm
    # print(strain) # name of the strain
    # print(list_of_genes) # dictionary with the "new" gene number/identifier as the key and the "old" gene number/identifier as the value
    # print(genes_and_contig) #dictionary with the "new" gene number/identifier as the key and the contig as the value

    ##set counter for the number of the gene cluster
    counter = 0

    #### check to see if there are no gene clusters that are reported:
    #if len(filtered_groups) == 0:
        #print("No gene clusters found")
        #continue

    for group in filtered_groups: # now getting into working with each gene cluster/gene cluster
        #print(group)
        ## make name

        contig_check =  filter.same_contigs_check(group, strain, genes_and_contig)
        #print(contig_check)

        if contig_check[0] == False:
            ### the hits are spread across more than one contigs

            with open(output_dir+"/log_file.txt", "a") as output:
                output.write("Contig break over gene cluster in {strain} \n\n".format(strain = strain))
            ### need to output this in some manner

            continue

        elif contig_check[0] == True: # if true, hits are on a single contig


            contig = contig_check[1] # store the name of the contig
            #update counter

            ### now check for the presence / absence of mandatory / accessory genes:




            query_types = filter.gene_names_in_cluster(group, mandatory_hmmer_output, strain)

            ## do the same again - are there enough genes:
            #print(number_query_types)
            if len(query_types) < int(min_genes_num):
            #    print("Not enough unique mandatory genes in the cluster found, not reporting gene cluster")

                continue


#7 - Extract the gff subsequence - with and without the extra sequences

        #### now set a name for the gene cluster:
        counter += 1 # update counter now that the identified gene cluster has passed filtering
        number_clusters_in_strain += 1 # and same for counting clusters / strain

        cluster_name = "{strain}_cluster_{counter}".format(strain=strain, counter = str(counter))

        #### find numbers for the beginning and the end of the gene cluster:
        min_gene = list_of_genes[strain+"_"+str(min(group))]
        max_gene = list_of_genes[strain+"_"+str(max(group))]


        ## get orientation, or if single gene:
        if min_gene == max_gene:
            orientation = "same"
        elif min_gene < max_gene:
            orientation = "forward"
        elif min_gene > max_gene:
            orientation = "reverse"



        extraction_details = extract.extract_a2b(min_gene, max_gene, annotation, strain, str(counter), upstream, downstream, contig, output_dir,orientation)

        GC_cluster = extraction_details[0]
        cluster_start = extraction_details[1]
        cluster_stop = extraction_details[2]
        extended_start = extraction_details[3]
        extended_stop = extraction_details[4]
        cluster_annotation = extraction_details[5]


        length_of_operon = cluster_start - cluster_stop
        if length_of_operon < 0:
            length_of_operon = length_of_operon * -1  ### change if the direction isn't left to write (i.s. on the reverse strand if you were to look at it in artemis)



#8 Write out the stats


        ### get dictionary for gene names:
        mandatory_genes_in_cluster = filter.gene_names_in_cluster(group, mandatory_hmmer_output, strain)


        #### get lists of the numbers of genes for each hmmer query, for the genes that are present in the operon
        list_of_mandatory  = []
        for gene in mandatory_names:
            if gene not in mandatory_genes_in_cluster:
                list_of_mandatory.append("0")
            else:
                list_of_mandatory.append(str(len(mandatory_genes_in_cluster[gene])))

        #### start creating the output:

        with open(output_dir+"/"+strain+"/cluster_stats.csv", "a") as output:
            output.write("""{strain}_cluster_{counter},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_mandatory_genes},{found_number_of_mandatory_genes},{percent_of_mandatory_genes_in_query},{hmm_genes},{GC_cluster},{GC_genome},{GCoperonbyGCgenome}\n""".format(strain=strain,counter=str(counter),contig=contig,start_on_contig=str(cluster_start),stop_on_contig=str(cluster_stop),length=abs(int(length_of_operon)),number_of_mandatory_genes=len(mandatory_names),found_number_of_mandatory_genes=len(mandatory_genes_in_cluster),percent_of_mandatory_genes_in_query=(float(len(mandatory_genes_in_cluster))/float(len(mandatory_names)))*100,hmm_genes=','.join(list_of_mandatory),GC_cluster=GC_cluster,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_cluster/GC_genome)))


        ####now work out and write out the input for gggenes:

        ### first, if using the accessory flag, need to create a new dictionary:


        all_genes_in_cluster = mandatory_genes_in_cluster

        ### now run a function to extract the gggenes information

        gggenes_input = extract.extract_gggenes_info(cluster_annotation, all_genes_in_cluster,list_of_genes, extended_start, extended_stop, strain, output_dir, str(counter))



        with open(output_dir+"/"+strain+"/gggenes_input.csv","a") as output:
            for line in gggenes_input:
                output.write(line+"\n")




    ####### record number of clusters / strain:
    with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
        output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = number_clusters_in_strain,rejected_clusters = number_rejected_clusters_in_strain,contig_break_clusters = str(number_contig_break_clusters_in_strain)))


    ### remove singlefastas folder:
    os.system("rm -r {output_dir}/{strain}/singlefastas".format(output_dir=output_dir,strain=strain))
    os.system("rm  {output_dir}/{strain}/contigs.fna".format(output_dir=output_dir,strain=strain))
