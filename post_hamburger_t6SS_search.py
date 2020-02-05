#!/usr/bin/env python3

import os, time, multiprocessing

###### parse arguments

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="")
        parser.add_argument('-i',
                '--input_dir',
                action='store',
                required=True,
                help='Original output directory of hamburger')
    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



##### functions

def __extract_tssBC_sequences__(directory): # looks for tssBC, finds the closet copies, makes sure there is one sequence each and the genes are an acceptable length
    #### get the strain name from the directory:

    strain = directory.split('/')[-1]
    #print(strain)

    ## check to see if any gene clusters were extracted:

    strain_statistics_data = open("{directory}/strain_statistics.csv".format(directory = directory))
    strain_statistics = strain_statistics_data.readlines()
    strain_statistics_data.close()

    number_of_gene_clusters = int(strain_statistics[0].split(',')[1])

    if number_of_gene_clusters == 0:
        ### no gene clusters were found in this genome - so exit from the function:
        return


    ### open the protein sequence file now - so don't have to open it more than once:
    ### open file with protein sequences:
    protein_sequences_data = open("{directory}/{strain}.faa".format(directory = directory, strain = strain))
    protein_sequences_lines = protein_sequences_data.readlines()
    protein_sequences_data.close()
    entries_list = (''.join(protein_sequences_lines)).split('>')

    #### create files to write in the protein sequences:
    open("{directory}/{strain}_tssB.faa".format(directory = directory, strain = strain),"w")
    open("{directory}/{strain}_tssC.faa".format(directory = directory, strain = strain),"w")


    #### find tssB and tssC sequences in the directory - looking in the gggenes input csv file:
    genes_list_data = open("{directory}/gggenes_input.csv".format(directory = directory))
    genes_list = genes_list_data.readlines()
    genes_list_data.close()

    # make dictionary to put the different clusters in
    clusters = {}

    for line in genes_list:
        toks = line.split(',')

        cluster = toks[0]
        gene = toks[4]
        id = toks[9]

        if gene == "TssB" or gene == "TssC": # check to see if the gene is either tssB or tssC:
            if cluster not in clusters:
                clusters[cluster] = {} # create new dict in dict for this cluster to then add gene and gene numbers:
            ## add gene and the gene number
            if gene not in clusters[cluster]:
                clusters[cluster][gene] = [] # if gene not already in there - add a key for the gene, with the value being an empty list to put the CDS_id's in for that gene
            ##add the gene to this dictionary value
            clusters[cluster][gene].append(id)


    for cluster, gene in clusters.items():
        requires_sorting  = False
        if len(gene) < 2: # if there are not both genes for TssB and TssC - then can't analyse this cluster
        #    print("Not the genes required in {cluster}".format(cluster = cluster))
            continue
        for gene_name, ids in gene.items():

            if len(ids) > 1: # i.e. if there is more than one tssB or tssC gene in the cluster - need to pick the two genes that want in the alignment (want them to be adjacent)
                # need to find the two closest genes: - so mark as a cluster that needs to be changed:
                requires_sorting = True

        if requires_sorting == True:
            ### get the gene prefix:
            prefix = "_".join(id.strip("\n").split("_")[:-1])
            tssB_gene_nums = []
            tssC_gene_nums = []
            ### now get the integers of the gene numbers:
            for id in gene["TssB"]:
                tssB_gene_nums.append(int(id.strip("\n").split("_")[-1])) # take the gene number, convert to integer, and add to a new list:

            ### do the same for tssC
            for id in gene["TssC"]:
                tssC_gene_nums.append(int(id.strip("\n").split("_")[-1])) # take the gene number, convert to integer, and add to a new list:

            # now look for the two closest numbers: (i.e. tssB and tssC genes that are adjacent)
            closest_genes = None
            chosen_TssB = None
            chosen_TssC = None
            for tssB_num in tssB_gene_nums:
                for tssC_num in tssC_gene_nums:
                    distance = abs(tssB_num - tssC_num)
                    if closest_genes == None:
                        closest_genes = distance
                        chosen_TssB = tssB_num
                        chosen_TssC = tssC_num
                    elif distance < closest_genes:
                        closest_genes = distance
                        chosen_TssB = tssB_num
                        chosen_TssC = tssC_num


            chosen_TssB = "{prefix}_{num}".format(prefix = prefix, num = chosen_TssB)
            chosen_TssC = "{prefix}_{num}".format(prefix = prefix, num = chosen_TssC)

        if requires_sorting == False:
            chosen_TssB = gene["TssB"][0].strip("\n")
            chosen_TssC = gene["TssC"][0].strip("\n")
            ### now we have selected the tssB and tssC genes that we want - now we should check whether they are the right length:

        #### check the length of the genes within this cluster.

        for line in genes_list:
            toks = line.split(',')

            id = toks[9]
            start = int(toks[2])
            stop = int(toks[3])

        ### boundaries to set for TssB and TssC:
        ###  150aa < TssB < 200aa
        ###  400aa < TssC < 520aa

        # get length of tssB
            if chosen_TssB == id.strip("\n"):
                tssB_length = int(abs(stop-(start - 1)) / 3)
        # get length of tssC
            if chosen_TssC == id.strip("\n"):
                tssC_length = int(abs(stop-(start - 1)) / 3)

        if tssC_length < 400 or tssC_length > 520 or tssB_length < 150 or tssB_length > 200:
            ## not the right length:
            #print("not the right length in {cluster}".format(cluster=cluster))
            continue


        ### extract the protein sequences if passed from all of these...

        tssB_extracted = False
        tssC_extracted = False

        for entry in entries_list:
            entry_name = entry.split("\n")[0]
            #print(entry_name)
            if entry_name == chosen_TssB:
                # got the tssB sequence - now change the name and write out..
                # the first value if split by newline character will just be the sequence header, and will just replace this , but save the names they correspond to!
                sequence = entry.split("\n")[1:]
                # concatenate new name and sequence, then write out:
                tssB_sequence_to_write = ">{cluster_name}\n{sequence}".format(cluster_name = cluster, sequence = "\n".join(sequence))
                with open("{directory}/{strain}_tssB.faa".format(directory = directory, strain = strain),"a") as output:
                    output.write(tssB_sequence_to_write)
                tssB_extracted = True

            ## now do tssC
            if entry_name == chosen_TssC:
                # got the tssC sequence - now change the name and write out..
                # the first value if split by newline character will just be the sequence header, and will just replace this , but save the names they correspond to!
                sequence = entry.split("\n")[1:]
                # concatenate new name and sequence, then write out:
                tssC_sequence_to_write = ">{cluster_name}\n{sequence}".format(cluster_name = cluster, sequence = "\n".join(sequence))
                with open("{directory}/{strain}_tssC.faa".format(directory = directory, strain = strain),"a") as output:
                    output.write(tssC_sequence_to_write)
                tssC_extracted = True

            ### stop the loop if both tssC and tssB have been found:

            if tssB_extracted == True and tssC_extracted == True:
                break # done, can now move to the next cluster


def __concatenate_alignments__(aln_1, aln_2, concatenated_output_aln): # must have the same names for this to work
    alignment_1 = open(aln_1).readlines()

    alignment_2 = open(aln_2).readlines()



    sequences = {}


    for line in alignment_1:
        new_sequence = False

        if line.startswith(">"):
            new_sequence = True

        if new_sequence == True:
            identifier = line.strip().split('>')[1]
            if identifier not in sequences:
                sequences[identifier] = []

        if new_sequence == False:
            sequences[identifier].append(line)

    ## now do same with second file:

    for line in alignment_2:
        new_sequence = False

        if line.startswith(">"):
            new_sequence = True

        if new_sequence == True:
            identifier = line.strip().split('>')[1]
            if identifier not in sequences:
                sequences[identifier] = []

        if new_sequence == False:
            sequences[identifier].append(line)


    ### now write out:sys.argv[3],"w"

    open(concatenated_output_aln,"w")

    for identifier, sequence in sequences.items():
        with open(concatenated_output_aln,"a") as output:
            #write header
            output.write('>'+identifier+'\n')
            #write sequence:
            output.write(''.join(sequence))


### main script
def main():

    start = time.time()


    args = parseArgs()

    list_of_dirs = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(args.input_dir):
        if r.endswith("/"):
            continue # don't want to include the root directory name
        list_of_dirs.append(r)

    #print(list_of_dirs)


    # for dir in list_of_dirs:
    #
    #     __extract_tssBC_sequences__(dir)


    #or multiprocess - go through each dir of the hamburger output, and extract a tssB and tssC sequence from each cluster, if they are there


    pool = multiprocessing.Pool()

    result = pool.map(__extract_tssBC_sequences__, list_of_dirs)

    #### once extracted from each - as will multiprocess these stages - concatenate all of these, align, concatenate again, draw tree and then pass to R:


    ### adding in the reference set as well ---- need to sort out the absolute links here when putting onto github!!
    os.system("cat {output_dir}/*/*_tssC.faa ~/github/hamburger/t6ss_reference_set/tssC.fasta > {output_dir}/all_observed_tssC.faa".format(output_dir=args.input_dir))
    os.system("cat {output_dir}/*/*_tssB.faa ~/github/hamburger/t6ss_reference_set/tssB.fasta > {output_dir}/all_observed_tssB.faa".format(output_dir=args.input_dir))

    ## concatenate the reference set as well:


    ### align:
    os.system("muscle -in {output_dir}/all_observed_tssC.faa -out {output_dir}/all_observed_tssC_aligned.fasta".format(output_dir=args.input_dir))
    os.system("muscle -in {output_dir}/all_observed_tssB.faa -out {output_dir}/all_observed_tssB_aligned.fasta".format(output_dir=args.input_dir))


    #### concatenate alignments:


    __concatenate_alignments__("{output_dir}/all_observed_tssB_aligned.fasta".format(output_dir=args.input_dir), "{output_dir}/all_observed_tssC_aligned.fasta".format(output_dir=args.input_dir), "{output_dir}/tssBC_alignment.fasta".format(output_dir=args.input_dir))

    ### draw tree:
    os.system("iqtree -s {output_dir}/tssBC_alignment.fasta -bb 1000".format(output_dir=args.input_dir))

    ## now - pass to R to group together the different types of T6SS that were identified

    end = time.time()


    print(end-start)





# Search by each gff / directory that was created for each folder

#1. Search for tssB and tssC genes in each identified cluster

#2. Are these the only copies / operon - are they the right size, and if there are two, can we pick the ones that are next to each other to use?

#3. Insert these into the reference set, and align

#4. Run the Rscript to separate these out into different subtypes - firstly by the 6 subtypes

#5. What if there are several different types within one subtype? e.g. Serratia - can we perform some sort of grouping / clustering ??



if __name__ == '__main__':
    main()
