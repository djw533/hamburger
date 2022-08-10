#!/usr/bin/env python3

import os
import os.path
import sys
import time
from os import walk


###### modify the path/name for hmmsearch below:

hmmsearch='hmmsearch'

# ### need the directory where hamburger is installed:
hamburger_base_directory = os.path.abspath(__file__).split("/scripts/hamburger.py")[0]

# add this folder to the path
#sys.path.append(hamburger_base_directory)
sys.path.insert(0,hamburger_base_directory)


# import hamburger modules
from hamburger import hmmer_functions
from hamburger import clean
from hamburger import tool_check
from hamburger import search
from hamburger import run_prodigal
from hamburger import istarmap


#import multiprocessing
from multiprocessing import Pool
import tqdm
from functools import partial





### and for the supplied models for standard pipeline & parameters
T6SS_core = "{hamburger_base_directory}/models/T6SS/T6SS_core.hmm".format(hamburger_base_directory=hamburger_base_directory)
T6SS_accessory = "{hamburger_base_directory}/models/T6SS/T6SS_accessory.hmm".format(hamburger_base_directory=hamburger_base_directory)


###### parse arguments

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="""Extract and plot gene_clusters based on hmm profiles
--------------------HaMBURGER--------------------\n
\n


-------HMmer Based UndeRstandinG of gene clustERs------\n
\n
              _....----'''----...._
           .-'  o    o    o    o   '-.
          /  o    o     o    o   o    \  \t
       __/__o___o_ _ o___ _o _ o_ _ _o_\__
      /                                   \ \t
      \___________________________________/
        \~`-`.__.`-~`._.~`-`~.-~.__.~`-`/
         \                             /
          `-._______________________.-'
""")
        parser.add_argument('-i',
                '--mandatory',
                action='store',
            #    required=True,
                help='Mandatory hmm profile input <required if not using -t flag> ')
        parser.add_argument('-a',
                '--accessory',
                action='store',
                required=False,
                help='Accessory hmm profile input')
        parser.add_argument('-g',
			    '--gff',
			    action='store',
                nargs='+',
                required=False,
			    help='Gff file(s) to search. Can be used in combination with --fasta or standalone (if fasta is appended to the end of the gff inbetween a line with ##FASTA)')
        parser.add_argument('-f',
			    '--fasta',
			    action='store',
                nargs='+',
                required=False,
			    help='Fasta file(s) to search. can be used in combination with --gff or standalone, in which case prodigal will predict CDSs')
        parser.add_argument('-m',
			    '--min_genes',
			    action='store',
                default=4,
			    help='Minimum number of genes in gene cluster, default = 4')
        parser.add_argument('-l',
			    '--genes_gap',
			    action='store',
                default=10,
			    help='Maximum number of genes gap between hits, default = 10')
        parser.add_argument('-u',
			    '--upstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include upstream/"right" of gene cluster, default = 0')
        parser.add_argument('-d',
			    '--downstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include downstream/"left" of gene cluster, default = 0')
        parser.add_argument('-c',
			    '--cutoff',
                type = float,
			    action='store',
                default=20,
			    help='Cutoff HMMER score for each hit, default = 20')
        parser.add_argument('-t',
			    '--t6ss',
			    action='store_true',
			    help='Automatic searching for T6SSs, uses min_genes = 8, genes_gap = 12, mandatory hmm profile of all 13 tss genes')
        parser.add_argument('-n',
			    '--num_threads',
			    action='store',
                default = 1,
			    help='Number of threads to use, default = 1')
        parser.add_argument('-k',
			    '--keep_files',
			    action='store_true',
			    help='Keep all intermediate files produced, default = False')
#        parser.add_argument('-p',
#			    '--pfam',
#			    action='store',
#                default = False
#			    help='Text file of pfam domains to use for hmms - one per line')
        parser.add_argument('-o',
			    '--output',
			    action='store',
                default="Hamburger_output",
                #default="HaMBURGER_output_{time}".format(time='_'.join('_'.join(str(datetime.now()).split('.')[0].split(':')).split())),
			    help='Output directory, default =Hamburger output. Will not write over a previously existing output folder!')
        parser.add_argument('-q',
                '--itol',
                action='store_true',
                help='Create itol output for number of T6SSs and subtypes per strain')
        parser.add_argument('-w',
			    '--overwrite',
			    action='store_true',
			    help="Overwrite existing blast database directory")
        parser.add_argument('-v',
			    '--version',
			    action='store_true',
			    help="Print version and exit")



    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



#################################
def main():


    start = time.time()
    args = parseArgs()


    if args.version:
        print("Hamburger version 0.2.0")
        sys.exit()

    ## check if hmmsearch is installed:

    if tool_check.is_tool("hmmsearch") == False:
        print("Please install hmmer before running hamburger")
        sys.exit()



    ### if t6ss flag requested then set these using the t6ss default setting
    if args.t6ss == True:
        min_genes_num = 8
        genes_gap_num = 12
        args.mandatory = T6SS_core#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_core.hmm"
        args.accessory = T6SS_accessory#"/home/djwilliams/github/hamburger/models/T6SS/T6SS_accessory.hmm"

    elif args.t6ss == False:
        if args.mandatory == False:
            print("Need to provide input mandatory hmm profile(s): --mandatory")
            sys.exit()
        if args.min_genes == False:
            print("Need to set the minimum number of hits: --min_genes)")
            sys.exit()
        if args.genes_gap == False:
            print("Need to set the maximum genes gap: --genes_gap)")
            sys.exit()
        # now change the variable names from args.xxx
        min_genes_num = args.min_genes
        genes_gap_num = args.genes_gap
        mandatory_models = args.mandatory
        if args.accessory is not None: # if accessory models were given change variable name
            accessory_models = args.accessory


    # check to see whether gff, fasta, or gff and fasta input used:


    if args.gff is not None and args.fasta is not None:
        #try:
        #check lengths of the list and the file extension
        if len(args.gff) != len(args.fasta):
            sys.exit("If using both gff and fasta files, please provide equal number of each")


        for gff, fasta in zip(args.gff,args.fasta):
            if gff.endswith(".gff") == False or fasta.endswith(".fasta") == False:
                sys.exit("Please make sure all gff and fasta files end with .gff and .fasta")
            if gff.split('/')[-1].replace(".gff","") != fasta.split('/')[-1].replace(".fasta",""):
                sys.exit("Please make sure the prefixes for the gff and fasta files are the same")

            #now concatenate:
            # gff_parsing.concat_gff_and_fasta(fasta,gff,"temp_input_gffs/{gff_file}".format(gff_file = gff.split('/')[-1]))

        strain_names = [filename.split('/')[-1].replace(".gff","") for filename in args.gff] # set strain names to iterate over


        #except:
            #sys.exit("Something was wrong when inputting fasta and gff entry. Exiting.")


    if args.gff is None and args.fasta is not None:
        try:

            if tool_check.is_tool("prodigal") == False:
                print("Please install prodigal before running hamburger if only using fasta input")
                sys.exit()

            for fasta in args.fasta:
                if fasta.endswith(".fasta") == False:
                    sys.exit("Please make sure all fasta files end with .fasta")

            args.gff = [None] * len(args.fasta)
            strain_names = [filename.split('/')[-1].replace(".fasta","") for filename in args.fasta] # set strain names to iterate over

        except:
            sys.exit("Something was wrong when inputting fasta and gff entry. Exiting.")


    if args.gff is not None and args.fasta is None:
        try:
            for gff in args.gff:
                if gff.endswith(".gff") == False:
                    sys.exit("Please make sure all gff files end with .gff")

            args.fasta = [None] * len(args.gff)
            strain_names = [filename.split('/')[-1].replace(".gff","") for filename in args.gff] # set strain names to iterate over

        except:
            sys.exit("Something was wrong when inputting fasta and gff entry. Exiting.")


    if args.gff is None and args.fasta is None:
        sys.exit("Please provide either fasta or gff files (or both) as input. Exiting.")


    ### check whether output directory already exists, overwrite if specified

    try:
        tool_check.check_output_dir(args.output,
                                args.overwrite)
    except ValueError as e:
        sys.exit("{dir} exists. Set -w to overwrite".format(dir = args.output))

    output_dir = args.output
    # os.makedirs(output_dir)


##### steps Required


#1 - Read the HMM queries and get the names of each query
    # - for both mandatory and accessory

    mandatory_names = hmmer_functions.read_hmms(args.mandatory)

    if args.accessory is not None or args.t6ss == True:
        accessory_names = hmmer_functions.read_hmms(args.accessory)


##### other outputs :

    with open(output_dir+"/strain_statistics.csv","w") as output:
        output.write("strain,number_of_gene_clusters,number_rejected_clusters,number_contig_break_clusters\n")

    with open(output_dir+"/gggenes_input.csv","w") as output:
        output.write("operon,number,start,end,gene,strand,direction,strain,CDS_identifier,hamburger_CDS_identifier\n")

#2 for each gff file make an ouput folder and work in it:

### cycle through each of the gff files in the provided list of files:


    #gff_files = args.gff


    start = time.time()


    #pool = multiprocessing.Pool(processes=int(args.num_threads))



    if args.accessory is not None or args.t6ss == True:

        with open(output_dir+"/log_file.txt", "w") as output:
            output.write("""--------------------HaMBURGER--------------------
            \n
            -------HMmer Based UndeRstandinG of gene clustERs------
            \n
            Using the following parameters as input: \n\n
            \tMinimum number of genes for cluster:  -> {min_genes}\n
            \tMaximum gap between genes in cluster: -> {genes_gap}\n
            \tUpstream region length:               -> {upstream}\n
            \tDownstream region length:             -> {downstream}\n
            \tHmmer cutoff:                         -> {cutoff}\n
            \tSearching for the following genes     ->\n
            \tMandatory genes:\n
            \t\t{mandatory_genes}\n\n
            \tAccessory genes:\n
            \t\t{accessory_genes}\n
            """.format(min_genes=min_genes_num, genes_gap=genes_gap_num, cutoff=args.cutoff, upstream=args.upstream, downstream=args.downstream, mandatory_genes='\n\t\t'.join(mandatory_names), accessory_genes='\n\t\t'.join(accessory_names)))



        with open(output_dir+"/cluster_stats.csv", "w") as output:
            output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,number_of_accessory_genes,found_number_of_accessory_genes,percent_of_accessory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names+accessory_names)))


        # multiprocess
        # result = pool.map(search_single_genome, gff_files)

        print("Running Hamburger - using both mandatory and accessory HMMs")

        #run differently depending on whether supplying fasta, gff, or both:

        func = partial(search.search_single_genome, args.mandatory,args.accessory,min_genes_num,genes_gap_num,args.upstream,args.downstream,args.cutoff,args.t6ss,args.output,args.keep_files)


        with Pool(int(args.num_threads)) as pool:

            for _ in tqdm.tqdm(pool.istarmap(func, zip(args.gff,args.fasta)), total=len(args.gff)):
                pass


    else:


        with open(output_dir+"/log_file.txt", "w") as output:
            output.write("""--------------------HaMBURGER--------------------
            \n
            -------HMmer Based UndeRstandinG of gene clustERs------
            \n
            Using the following parameters as input: \n\n
            \tMinimum number of genes for cluster:  -> {min_genes}\n
            \tMaximum gap between genes in cluster: -> {genes_gap}\n
            \tUpstream region length:               -> {upstream}\n
            \tDownstream region length:             -> {downstream}\n
            \tHmmer cutoff:                         -> {cutoff}\n
            \tSearching for the following genes     ->\n
            \tMandatory genes:\n
            \t\t{mandatory_genes}\n\n
            """.format(min_genes=min_genes_num, genes_gap=genes_gap_num, cutoff=args.cutoff, upstream=args.upstream, downstream=args.downstream, mandatory_genes='\n\t\t'.join(mandatory_names)))


        with open(output_dir+"/cluster_stats.csv", "w") as output:
            output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names)))


        # multiprocess
        #result = pool.map(search_single_genome_no_accessory, gff_files)

        print("Running Hamburger - no accessory HMMs given")

        #set function
        func = partial(search.search_single_genome, args.mandatory,args.accessory,min_genes_num,genes_gap_num,args.upstream,args.downstream,args.cutoff,args.t6ss,args.output,args.keep_files)


        with Pool(int(args.num_threads)) as pool:

            for _ in tqdm.tqdm(pool.istarmap(func, zip(args.gff,args.fasta)), total=len(args.gff)):
                pass


    end = time.time()


    print(end-start)



#10 Now group all of the different statistics files together


    ##now go through these and concatenate files to lists to then be concatenated if they exist:
    statistics_files  = []
    gggenes_input_files = []
    cluster_stats_files = []

    #strain_names  = [gff_file.split('/')[-1] for gff_file in gff_files]

    for dirpath, dirnames, filenames in os.walk(args.output):
        #check if the directory is in the gff files list
        if "{cut_dirpath}".format(cut_dirpath = dirpath.split('/')[-1]) in strain_names:
            # now established that it's a correct directory , take the strain_statistics file, the gggenes_input file and the cluster_stats files (if they're there)
            if "strain_statistics.csv" in filenames:
                statistics_files.append("{dirpath}/strain_statistics.csv".format(dirpath = dirpath))
            if "gggenes_input.csv" in filenames:
                gggenes_input_files.append("{dirpath}/gggenes_input.csv".format(dirpath = dirpath))
            if "cluster_stats.csv" in filenames:
                cluster_stats_files.append("{dirpath}/cluster_stats.csv".format(dirpath = dirpath))

    #now check that gggenes_input_files and cluster_stats_files lists aren't empty:

    if len(gggenes_input_files) == 0 or len(cluster_stats_files) == 0 :

        if args.t6ss == True:
            print("No T6SSs found.")
        else:
            print("No gene clusters found.")

        ## now cat all the strain statistics files:
        os.system("cat {statistics_files} >> {output_dir}/strain_statistics.csv".format(statistics_files = ' '.join(statistics_files), output_dir=output_dir))


    else:

        os.system("cat {statistics_files} >> {output_dir}/strain_statistics.csv".format(statistics_files = ' '.join(statistics_files), output_dir=output_dir))
        os.system("cat {gggenes_input_files} >> {output_dir}/gggenes_input.csv".format(gggenes_input_files = ' '.join(gggenes_input_files), output_dir=output_dir))
        os.system("cat {cluster_stats_files} >> {output_dir}/cluster_stats.csv".format(cluster_stats_files = ' '.join(cluster_stats_files), output_dir=output_dir))


#10 - Cleanup:

    # if args.keep_files is False:
    #     for strain in strain_names:
    #         clean.clean_up_files(args.output+"/"+strain)


    if args.itol == True:
    ## now convert the T6SS subtypes info for itol output:
        stats_file  = open("{output_dir}/strain_statistics.csv".format(output_dir=output_dir))
        stats_data = stats_file.readlines()
        stats_header = stats_data[0].split(',') # get the different headers required:
        for num in range(1,len(stats_header)): # go through each one and make an itol input file:
            #create list to write to file:
            to_write = []
            to_write.append("DATASET_GRADIENT\nSEPARATOR SPACE\nDATASET_LABEL {subtype}\nCOLOR #ff0000\nDATA\n".format(subtype = stats_header[num]))

            #then cycle through the subtypes file:
            for line in stats_data[1:]:
                toks = line.split(',')
                to_write.append("{strain} {number_observations}\n".format(strain = toks[0], number_observations = toks[num]))

            ## now write out:
            with open("{output_dir}/{subtype}_itol.txt".format(subtype=stats_header[num],output_dir=output_dir), "w") as output:
                for l in to_write:
                    output.write(l)



if __name__ == '__main__':
    main()
