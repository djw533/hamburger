#!/usr/bin/env python3

import os, time, multiprocessing
from shutil import which
import sys



###set hamburger base directory:
hamburger_base_directory = os.path.abspath(__file__).split("/scripts/post_hamburger_t6SS_search.py")[0]

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
        # parser.add_argument('-t',
        #         '--tree_building',
        #         action='store',
        #         default = "fasttree",
        #         help='Tree building software to use [iqtree, fasttree], default = fasttree')
        parser.add_argument('-a',
                '--itol',
                action='store_true',
                help='Create itol output for number of T6SSs and subtypes per strain')
        parser.add_argument('-n',
                '--num_threads',
                action='store',
                default = 1,
                help='Number of threads to use, default = 1)')
    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



##### functions

def __is_tool__(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    return which(name) is not None

def __extract_tssBC_sequences__(directory): # looks for tssBC, finds the closet copies, makes sure there is one sequence each and the genes are an acceptable length
    #### get the strain name from the directory:
    args = parseArgs()

    strain = directory.split('/')[-1]
    #print(strain)

    ## check to see if any gene clusters were extracted:
    #print(directory)

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
            print("Not the genes required in {cluster}".format(cluster = cluster))
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

        # if tssC_length < 300 or tssC_length > 600 or tssB_length < 100 or tssB_length > 300:
        #     ## not the right length:
        #     print("not the right length in {cluster}, tssB_length = {tssB_length}, tssC_length = {tssC_length}".format(cluster=cluster,tssB_length = tssB_length, tssC_length = tssC_length))
        #     continue


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

        if tssB_extracted != True or tssC_extracted != True:
            print("something wrong with {cluster}".format(cluster=cluster))


def __concatenate_alignments__(aln_1, aln_2, concatenated_output_aln): # must have the same names for this to work
    args = parseArgs()
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


def __align_sequences__(gene_name): # cat tssB/tssC sequences in reference and those found by hamburger
    args = parseArgs()
    os.system("cat {output_dir}/*/*_{gene_name}.faa {hamburger_base_directory}/t6ss_reference_set/{gene_name}.fasta > {output_dir}/all_observed_{gene_name}.faa".format(hamburger_base_directory=hamburger_base_directory,output_dir=args.input_dir, gene_name=gene_name))
    os.system("muscle -in {output_dir}/all_observed_{gene_name}.faa -out {output_dir}/all_observed_{gene_name}_aligned.fasta".format(output_dir=args.input_dir, gene_name=gene_name))


# def __run_mash_sketch__(fasta_addresses):
#
#
#
# def __run_mash_dist__(fasta_addresses):


# def __compare_t6_subtypes__(cluster_stats, current_dir):
#
#     ### read in the cluster data
#     cluster_data = open(cluster_stats)
#     cluster_lines = cluster_data.readlines()
#     cluster_data.close()
#
#     #create a dictionary of these subtypes
#     subtype_dict = {}
#
#     # add subtypes as keys to dict:
#     for line in cluster_lines:
#         toks = line.split(',')
#         #get the subtype:
#         subtype = toks[-1]
#         #add to dictionary, create key if doesn't already exist:
#         if subtype not in subtype_dict:
#             subtype_dict[subtype] = []
#         ## now set the cluster name and add to the dictionary for that subtype:
#         cluster_name = toks[0]
#         subtype_dict[subtype].append(cluster_name)
#
#     ## now need to addresses of the fasta files for the clusters:
#     # get cwd:
#     working_dir = os.getcwd()
#     os.system("ls {working_dir}*/*cluster_*.fna > cluster_fasta_addresses.txt".format(working_dir = working_dir))
#
#     # run mash sketch:
#     __run_mash_sketch__("cluster_fasta_addresses.txt")
#
#
#     #get list of mash addresses:
#     os.system("ls {working_dir}*/*cluster_*.fna.msh > cluster_fasta_addresses.txt".format(working_dir = working_dir))
#     # run mash dist
#     __run_mash_dist__("cluster_msh_addresses.txt")






### main script
def main():


    ## check that programs are installed:

    if __is_tool__("muscle") == False:
        print("Please install muscle before running hamburger. Exiting script.")
        sys.exit()


    start = time.time()


    args = parseArgs()


    ## once arguments set - check for tree building software - allow both fasttree and FastTree:

    # if args.tree_building == "fasttree":

    if __is_tool__("fasttree") == False and __is_tool__("FastTree") == False:
        print("Please install muscle. Exiting script.")
        sys.exit()

    elif __is_tool__("fasttree") == True:
        tree_software = "fasttree"

    elif __is_tool__("FastTree") == True:
        tree_software = "FastTree"

    ## check iqtree:

    # elif args.tree_building == "iqtree":
    #
    #     if __is_tool__("iqtree") == False:
    #         print("Please install iqtree. Exiting script.")
    #         sys.exit()
    #
    # else:
    #     print("Please use either iqtree or fasttree with -t. Exiting script.")
    #     sys.exit()


    ##check that the R dependencies are here:
    os.system("Rscript {hamburger_base_directory}/scripts/check_R_dependencies.R {hamburger_base_directory}".format(hamburger_base_directory=hamburger_base_directory))

    ## check the dependencies_check
    dependencies_check = open("dependencies_check.csv")
    d = dependencies_check.readlines()
    dependencies_check.close()
    if len(d) > 1: # then some R packages aren't installed:
        print("The following R libraries are not installed, please install them before running this script:")
        for l in d[1:]:
            print(l.strip())
        print("Exiting script.")
        os.remove("dependencies_check.csv")
        sys.exit()
    ## also remove here:
    os.remove("dependencies_check.csv")



    list_of_dirs = []
    # r=root, d=directories, f = files
    first_search = True
    for r, d, f in os.walk(args.input_dir):
        # if first_search == True:
        #     first_search == False # set this so that we don't just take the root file:
        #     continue # don't want to include the root directory name
        list_of_dirs.append(r)


    list_of_dirs = list_of_dirs[1:]


    if int(args.num_threads) > 1:
        num_threads_for_muscle = 2

    else:
        num_threads_for_muscle = 1


    #print(list_of_dirs)


    # for dir in list_of_dirs:
    #
    #     __extract_tssBC_sequences__(dir)


    #or multiprocess - go through each dir of the hamburger output, and extract a tssB and tssC sequence from each cluster, if they are there


    pool = multiprocessing.Pool(processes = int(args.num_threads))

    result = pool.map(__extract_tssBC_sequences__, list_of_dirs)

    # for dirname in list_of_dirs:
    #     __extract_tssBC_sequences__(dirname)

    #### once extracted from each - as will multiprocess these stages - concatenate all of these, align, concatenate again, draw tree and then pass to R:


    ### now multiprocess the alignments:
    sequences_to_align = ["tssB","tssC"]

    pool = multiprocessing.Pool(processes = num_threads_for_muscle)

    result = pool.map(__align_sequences__, sequences_to_align)


    #### concatenate alignments:


    __concatenate_alignments__("{output_dir}/all_observed_tssB_aligned.fasta".format(output_dir=args.input_dir), "{output_dir}/all_observed_tssC_aligned.fasta".format(output_dir=args.input_dir), "{output_dir}/tssBC_alignment.fasta".format(output_dir=args.input_dir))

    ### draw tree:
    # if args.tree_building == "iqtree":
    #
    #     os.system("{tree_software} -s {output_dir}/tssBC_alignment.fasta -bb 1000 -m LG -nt {threads}".format(output_dir=args.input_dir, threads=int(args.num_threads), tree_software = tree_software))
    # ## use fasttree instead - is quicker??:
    # elif args.tree_building == "fasttree":

    os.system("{tree_software} {output_dir}/tssBC_alignment.fasta > {output_dir}/tssBC_alignment.fasta.treefile".format(output_dir=args.input_dir, tree_software = tree_software))


    ## now - pass to R to group together the different types of T6SS that were identified

    os.chdir(args.input_dir)

    os.system("Rscript {hamburger_base_directory}/scripts/filter_t6_types.R {hamburger_base_directory}".format(hamburger_base_directory=hamburger_base_directory))


    if args.itol == True:
    ## now convert the T6SS subtypes info for itol output:
        subtypes_file  = open("T6SS_cluster_types.csv")
        subtypes_data = subtypes_file.readlines()
        subtypes_header = subtypes_data[0].split(',') # get the different headers required:
        for num in range(1,len(subtypes_header)): # go through each one and make an itol input file:
            #create list to write to file:
            to_write = []
            to_write.append("DATASET_GRADIENT\nSEPARATOR SPACE\nDATASET_LABEL {subtype}\nCOLOR #ff0000\nDATA\n".format(subtype = subtypes_header[num]))

            #then cycle through the subtypes file:
            for line in subtypes_data[1:]:
                toks = line.split(',')
                to_write.append("{strain} {number_observations}\n".format(strain = toks[0], number_observations = toks[num]))

            ## now write out:
            with open("{subtype}_itol.txt".format(subtype=subtypes_header[num]), "w") as output:
                for l in to_write:
                    output.write(l)


    #### now compare the different t6 subtypes using mash:

    #__compare_t6_subtypes__()

    ### plot with phylogeny if it's given:
    # if args.tree is not None:
    #     os.system("Rscript {hamburger_base_directory}/scripts/plot_t6_types_with_phylogeny.R {tree}".format(tree = args.tree,hamburger_base_directory=hamburger_base_directory))


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
