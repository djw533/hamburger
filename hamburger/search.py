


import os
import sys

#import functions within this directory
from hamburger import extract
from hamburger import fasta_parsing
from hamburger import filter
from hamburger import hmmer_functions
from hamburger import tool_check
from hamburger import gff_parsing
from hamburger import run_prodigal
from hamburger import clean

hmmsearch='hmmsearch'



def search_single_genome(mandatory_models,accessory_models,min_genes_num,genes_gap_num,upstream,downstream,cutoff,t6ss,output_dir,keep_files,gff_file,fasta_file):


    # # if no gff file supplied then run prodigal and create the gff file
    if gff_file == None and fasta_file != None:
        strain = fasta_file.split('/')[-1].replace(".fasta","")

        ### now made directory for this strain:
        strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
        os.makedirs(strain_dir)

        gff_file = "{strain_dir}/{strain}.gff".format(strain_dir = strain_dir, strain = strain) # set output gff file to be put into the overall output folder for posterity

        #try:
        run_prodigal.run_prodigal(fasta_file,gff_file) # now run prodigal, add fasta sequence write out as gff

        #except:
        #sys.exit("Problem with running prodigal")

    # if both gff and fasta file are supplied, then concatenate them together:
    elif gff_file != None and fasta_file != None:

        strain = gff_file.split('/')[-1].replace(".gff","")

        ### now made directory for this strain:
        strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
        os.makedirs(strain_dir)

        output_gff_file = "{strain_dir}/{strain}.gff".format(strain_dir = strain_dir, strain = strain) # set output gff file to be put into the overall output folder for posterity

        gff_parsing.concat_gff_and_fasta(fasta_file,gff_file,output_gff_file)

        gff_file  = output_gff_file

    else:
        # continue using the gff file that is passed across
        ### take the strain name (last field from "_" separator and remove the .gff from the end:)
        print(gff_file)
        strain = gff_file.split('/')[-1].replace(".gff","")

        ### now made directory for this strain:
        strain_dir = "{output_dir}/{strain}".format(output_dir = output_dir, strain = strain)
        os.makedirs(strain_dir)


    # read in the mandatory names
    mandatory_names = hmmer_functions.read_hmms(mandatory_models)

    #parse accessory (if they exist)
    if accessory_models is not None or t6ss == True:
        accessory_names = hmmer_functions.read_hmms(accessory_models)

    number_clusters_in_strain = 0
    number_rejected_clusters_in_strain = 0
    number_contig_break_clusters_in_strain = 0



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
    if accessory_models:
        accessory_hmmer_tuple = hmmer_functions.run_hmmer(hmmsearch, accessory_models, prot_seqs, cutoff, strain_dir, "accessory")
        accessory_hmmer_output = accessory_hmmer_tuple[0]
        accessory_hmmer_genes = accessory_hmmer_tuple[1]

        ### add accessory to the mandatory
        total_hmmer_output = mandatory_hmmer_output + accessory_hmmer_output
        total_hmmer_genes = mandatory_hmmer_genes + accessory_hmmer_genes

    else:
        total_hmmer_output = mandatory_hmmer_output
        total_hmmer_genes = mandatory_hmmer_genes

#5 - Cluster the HMMER hits according to the parameters

    #print("here {strain}".format(strain=gff_file))


    if len(mandatory_hmmer_genes) == 0:
        ##no hits
        with open(output_dir+"/log_file.txt", "a") as output:
            output.write("No hmmer hits in {strain}".format(strain =strain))
        ### and write to strain_statistics.csv
        with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
            output.write("{strain},{clusters},{rejected_clusters},{number_contig_break_clusters_in_strain}\n".format(strain=strain,clusters = "0",rejected_clusters = "0",number_contig_break_clusters_in_strain = number_contig_break_clusters_in_strain))
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
                number_contig_break_clusters_in_strain += 1

            ### need to output this in some manner

            continue

        elif contig_check[0] == True: # if true, hits are on a single contig


            contig = contig_check[1] # store the name of the contig
            #update counter

            ### now check for the presence / absence of mandatory / accessory genes:

            if accessory_models: # if accessory was used - remove these hits (if present) from the file:

                mandatory_in_cluster =  [x for x in group if x not in accessory_hmmer_genes]

                ### now associate each hit with the hmmer output to get the number of queries that were satisfied from the mandatory, concatenated hmm model
                mandatory_genes_in_cluster = filter.gene_names_in_cluster(mandatory_in_cluster, mandatory_hmmer_output, strain)


                #get the names of the accessory genes for later:
                accessory_in_cluster =  [x for x in group if x not in mandatory_hmmer_genes]

                ### now associate each hit with the hmmer output to get the number of queries that were satisfied from the mandatory, concatenated hmm model
                accessory_genes_in_cluster = filter.gene_names_in_cluster(accessory_in_cluster, accessory_hmmer_output, strain)

                ### cancel the search if the number of unique genes is not less than or equal to the min number set in the search:
                if len(mandatory_genes_in_cluster) < int(min_genes_num):
                #    print("Not enough unique mandatory genes in the cluster found, not reporting gene cluster")
                    number_rejected_clusters_in_strain += 1

                    continue

                    #### could also still extract the region, but report that t didn't pass threshold??

            else: ### if accessory models weren't used - but still need to query whether enough unique genes were found..

                mandatory_genes_in_cluster = filter.gene_names_in_cluster(group, mandatory_hmmer_output, strain)

                ## do the same again - are there enough genes:
                if len(mandatory_genes_in_cluster) < int(min_genes_num):
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


        #### get lists of the numbers of genes for each hmmer query, for the genes that are present in the operon
        list_of_mandatory  = []
        for gene in mandatory_names:
            if gene not in mandatory_genes_in_cluster:
                list_of_mandatory.append("0")
            else:
                list_of_mandatory.append(str(len(mandatory_genes_in_cluster[gene])))

        #### and now for accessory
        if accessory_models:
            list_of_accessory  = []
            for gene in accessory_names:
                if gene not in accessory_genes_in_cluster:
                    list_of_accessory.append("0")
                else:
                    list_of_accessory.append(str(len(accessory_genes_in_cluster[gene])))

        #### start creating the output:
        if accessory_models:
            with open(output_dir+"/"+strain+"/cluster_stats.csv", "a") as output:
                output.write("""{strain}_cluster_{counter},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_mandatory_genes},{found_number_of_mandatory_genes},{percent_of_mandatory_genes_in_query},{number_of_accessory_genes},{found_number_of_accessory_genes},{percent_of_accessory_genes_in_query},{hmm_genes},{GC_cluster},{GC_genome},{GCoperonbyGCgenome}\n""".format(strain=strain,counter=str(counter),contig=contig,start_on_contig=str(cluster_start),stop_on_contig=str(cluster_stop),length=abs(int(length_of_operon)),number_of_mandatory_genes=len(mandatory_names),found_number_of_mandatory_genes=len(mandatory_genes_in_cluster),percent_of_mandatory_genes_in_query=(float(len(mandatory_genes_in_cluster))/float(len(mandatory_names)))*100,number_of_accessory_genes=len(accessory_names),found_number_of_accessory_genes=len(accessory_genes_in_cluster),percent_of_accessory_genes_in_query=(float(len(accessory_genes_in_cluster))/float(len(accessory_names)))*100,hmm_genes=','.join(list_of_mandatory+list_of_accessory),GC_cluster=GC_cluster,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_cluster/GC_genome)))

        else:
            with open(output_dir+"/"+strain+"/cluster_stats.csv", "a") as output:
                output.write("""{strain}_cluster_{counter},{strain},{contig},{start_on_contig},{stop_on_contig},{length},{number_of_mandatory_genes},{found_number_of_mandatory_genes},{percent_of_mandatory_genes_in_query},{hmm_genes},{GC_cluster},{GC_genome},{GCoperonbyGCgenome}\n""".format(strain=strain,counter=str(counter),contig=contig,start_on_contig=str(cluster_start),stop_on_contig=str(cluster_stop),length=abs(int(length_of_operon)),number_of_mandatory_genes=len(mandatory_names),found_number_of_mandatory_genes=len(mandatory_genes_in_cluster),percent_of_mandatory_genes_in_query=(float(len(mandatory_genes_in_cluster))/float(len(mandatory_names)))*100,hmm_genes=','.join(list_of_mandatory),GC_cluster=GC_cluster,GC_genome=GC_genome,GCoperonbyGCgenome=str(GC_cluster/GC_genome)))


        ####now work out and write out the input for gggenes:

        ### first, if using the accessory flag, need to create a new dictionary:

        if accessory_models:
            all_genes_in_cluster = dict(mandatory_genes_in_cluster)
            all_genes_in_cluster.update(accessory_genes_in_cluster)

        else:
            all_genes_in_cluster = mandatory_genes_in_cluster

        ### now run a function to extract the gggenes information

        gggenes_input = extract.extract_gggenes_info(cluster_annotation, all_genes_in_cluster,list_of_genes, extended_start, extended_stop, strain, output_dir, str(counter))



        with open(output_dir+"/"+strain+"/gggenes_input.csv","a") as output:
            for line in gggenes_input:
                output.write(line+"\n")




    ####### record number of clusters / strain:
    with open(output_dir + "/"+strain+"/strain_statistics.csv","a") as output:
        output.write("{strain},{clusters},{rejected_clusters},{contig_break_clusters}\n".format(strain=strain,clusters = number_clusters_in_strain,rejected_clusters = number_rejected_clusters_in_strain,contig_break_clusters = number_contig_break_clusters_in_strain))

    ### remove singlefastas folder:
    #print("{output_dir}/{strain}/singlefastas".format(output_dir=output_dir,strain=strain))

    os.system("rm -r {output_dir}/{strain}/singlefastas".format(output_dir=output_dir,strain=strain))
    os.system("rm  {output_dir}/{strain}/contigs.fna".format(output_dir=output_dir,strain=strain))


    if keep_files is False:
        clean.clean_up_files("{output_dir}/{strain}".format(output_dir=output_dir,strain=strain), strain)
