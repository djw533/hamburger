#!/usr/bin/env python3

import os
import os.path
import sys
import time
from os import walk
import traceback

###### modify the path/name for hmmsearch below:

hmmsearch='hmmsearch'

# ### need the directory where hamburger is installed:
#hamburger_base_directory = os.path.abspath(__file__).split("/hamburger/cli.py")[0]
#print(hamburger_base_directory)

# add this folder to the path
#sys.path.append(hamburger_base_directory)
#sys.path.insert(0,hamburger_base_directory)


# import hamburger modules
from hamburger import hmmer_functions
from hamburger import clean
from hamburger import tool_check
from hamburger import search
from hamburger import run_prodigal
from hamburger import istarmap
from hamburger import filter_t6_types

#import multiprocessing
from multiprocessing import Pool
import tqdm
from functools import partial

#find T6SS hmm models
import importlib.resources as pkg_resources
from pathlib import Path

def get_t6ss_model(model_name: str) -> Path:
    """
    Locate the `.hmm` file shipped with the hamburger package.
    
    Example:
        get_t6ss_model("T6SS_core")  # → Path to models/T6SS/T6SS_core.hmm
    """
    # pkg_resources.files("hamburger") points at the package root
    pkg_root = pkg_resources.files("hamburger")
    return pkg_root / "models" / "T6SS" / f"{model_name}.hmm"

T6SS_core      = get_t6ss_model("T6SS_core")
T6SS_accessory = get_t6ss_model("T6SS_accessory")

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
			    help='Gff file(s) to search. gff and gff3 files supported. Can be used in combination with --fasta or standalone (if fasta is appended to the end of the gff inbetween the line "##FASTA")')
        parser.add_argument('-f',
			    '--fasta',
			    action='store',
                nargs='+',
                required=False,
			    help='Fasta file(s) to search. Can be used in combination with --gff or standalone, in which case prodigal will predict CDSs')
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
			    help='Automatic searching for T6SS loci, uses min_genes = 4, genes_gap = 10, mandatory hmm profile of 13 tssA-M genes')
        parser.add_argument('-ft',
			    '--filter_t6ss',
			    action='store_true',
			    help='Skip to filtering of T6SSs')
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
			    help="Overwrite existing output folder")
        parser.add_argument('-s',
			    '--save_gffs',
			    action='store_true',
			    help="Save output gff files")
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


    if args.filter_t6ss:
        print("Filtering T6SSs")
        filter_t6_types.run_filter(
            hamburger_dir = args.output,
            threads = args.num_threads,
            itol = args.itol,
            keep_files = args.keep_files
        )

        sys.exit()


    if args.version:
        print("Hamburger version 0.2.1")
        sys.exit()

    ## check if hmmsearch is installed:

    if tool_check.is_tool("hmmsearch") == False:
        print("Please install hmmer before running hamburger")
        sys.exit()



    ### if t6ss flag requested then set these using the t6ss default setting
    if args.t6ss == True:
        min_genes_num = args.min_genes
        genes_gap_num = args.genes_gap
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
            if not (gff.endswith(".gff") or gff.endswith(".gff3")) or not fasta.endswith(".fasta"):
                sys.exit("Please make sure all gff files end with .gff or .gff3 and fasta files end with .fasta")
            if gff.split('/')[-1].rsplit('.', 1)[0] != fasta.split('/')[-1].replace(".fasta",""):
                sys.exit("Please make sure the prefixes for the gff and fasta files are the same")

            #now concatenate:
            # gff_parsing.concat_gff_and_fasta(fasta,gff,"temp_input_gffs/{gff_file}".format(gff_file = gff.split('/')[-1]))

        strain_names = [filename.split('/')[-1].rsplit('.', 1)[0] for filename in args.gff] # set strain names to iterate over


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
                  if not (gff.endswith(".gff") or gff.endswith(".gff3")):
                    sys.exit("Please make sure all gff files end with .gff or .gff3")

            args.fasta = [None] * len(args.gff)
            strain_names = [filename.split('/')[-1].rsplit('.', 1)[0] for filename in args.gff] # set strain names to iterate over

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

    #make output gff folder

    if args.save_gffs == True:
        os.makedirs("{output_dir}/extracted_gff_regions".format(output_dir = output_dir))


##### steps Required


#1 - Read the HMM queries and get the names of each query
    # - for both mandatory and accessory

    mandatory_names = hmmer_functions.read_hmms(args.mandatory)

    if args.accessory is not None or args.t6ss == True:
        accessory_names = hmmer_functions.read_hmms(args.accessory)


##### other outputs :
    # create the lists for strain stats, gggenes input, and cluster stats
    total_gggenes_input = []
    total_strain_stats = []
    total_cluster_stats = []

    total_strain_stats.append("strain,number_of_gene_clusters,number_rejected_clusters,number_contig_break_clusters")
    total_gggenes_input.append("operon,number,start,end,gene,strand,direction,strain,CDS_identifier,hamburger_CDS_identifier,prot_seq")

    # with open(output_dir+"/strain_statistics.csv","w") as output:
    #     output.write("strain,number_of_gene_clusters,number_rejected_clusters,number_contig_break_clusters\n")
    #
    # with open(output_dir+"/gggenes_input.csv","w") as output:
    #     output.write("operon,number,start,end,gene,strand,direction,strain,CDS_identifier,hamburger_CDS_identifier\n")



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
            """.format(
                min_genes=min_genes_num,
                genes_gap=genes_gap_num,
                cutoff=args.cutoff,
                upstream=args.upstream,
                downstream=args.downstream,
                mandatory_genes='\n\t\t'.join(mandatory_names),
                accessory_genes='\n\t\t'.join(accessory_names)
                )
            )



        # with open(output_dir+"/cluster_stats.csv", "w") as output:
        #     output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,number_of_accessory_genes,found_number_of_accessory_genes,percent_of_accessory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names+accessory_names)))

        total_cluster_stats.append(','.join(
            ["gene_cluster",
            "strain",
            "contig",
            "start",
            "stop",
            "length",
            "number_of_mandatory_genes",
            "found_number_of_mandatory_genes",
            "percent_of_mandatory_genes_in_query",
            "number_of_accessory_genes",
            "found_number_of_accessory_genes",
            "percent_of_accessory_genes_in_query",
            "{hmm_genes}".format(hmm_genes=','.join(mandatory_names+accessory_names)),
            "GC_cluster",
            "GC_genome",
            "GCcluster/GCgenome"]
            )
        )

        # multiprocess
        # result = pool.map(search_single_genome, gff_files)

        print("Running Hamburger - using both mandatory and accessory HMMs")

        #run differently depending on whether supplying fasta, gff, or both:

        func = partial(
            search.search_single_genome,
            args.mandatory,
            args.accessory,
            min_genes_num,
            genes_gap_num,
            args.upstream,
            args.downstream,
            args.cutoff,
            args.t6ss,
            args.output,
            args.keep_files,
            args.save_gffs
        )


        with Pool(int(args.num_threads)) as pool:

            for output_lists in tqdm.tqdm(pool.istarmap(func, zip(args.gff,args.fasta)), total=len(args.gff)):
                total_gggenes_input.extend(output_lists[0])
                total_strain_stats.extend(output_lists[1])
                total_cluster_stats.extend(output_lists[2])
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


        # with open(output_dir+"/cluster_stats.csv", "w") as output:
        #     output.write("gene_cluster,strain,contig,start,stop,length,number_of_mandatory_genes,found_number_of_mandatory_genes,percent_of_mandatory_genes_in_query,{hmm_genes},GC_cluster,GC_genome,GCcluster/GCgenome\n".format(hmm_genes=','.join(mandatory_names)))

        total_cluster_stats.append(','.join(
            ["gene_cluster",
            "strain",
            "contig",
            "start",
            "stop",
            "length",
            "number_of_mandatory_genes",
            "found_number_of_mandatory_genes",
            "percent_of_mandatory_genes_in_query",
            "{hmm_genes}".format(hmm_genes=','.join(mandatory_names)),
            "GC_cluster",
            "GC_genome",
            "GCcluster/GCgenome"]
            )
        )



        # multiprocess
        #result = pool.map(search_single_genome_no_accessory, gff_files)

        print("Running Hamburger - no accessory HMMs given")

        #set function
        func = partial(
            search.search_single_genome,
            args.mandatory,
            args.accessory,
            min_genes_num,
            genes_gap_num,
            args.upstream,
            args.downstream,
            args.cutoff,
            args.t6ss,
            args.output,
            args.keep_files,
            args.save_gffs
        )


        with Pool(int(args.num_threads)) as pool:

            for output_lists in tqdm.tqdm(pool.istarmap(func, zip(args.gff,args.fasta)), total=len(args.gff)):
                total_gggenes_input.extend(output_lists[0])
                total_strain_stats.extend(output_lists[1])
                total_cluster_stats.extend(output_lists[2])
                pass


    end = time.time()





#10 Now group all of the different statistics files together


    # ##now go through these and concatenate files to lists to then be concatenated if they exist:
    # statistics_files  = []
    # gggenes_input_files = []
    # cluster_stats_files = []
    #
    # #strain_names  = [gff_file.split('/')[-1] for gff_file in gff_files]
    #
    # for dirpath, dirnames, filenames in os.walk(args.output):
    #     #check if the directory is in the gff files list
    #     if "{cut_dirpath}".format(cut_dirpath = dirpath.split('/')[-1]) in strain_names:
    #         # now established that it's a correct directory , take the strain_statistics file, the gggenes_input file and the cluster_stats files (if they're there)
    #         if "strain_statistics.csv" in filenames:
    #             statistics_files.append("{dirpath}/strain_statistics.csv".format(dirpath = dirpath))
    #         if "gggenes_input.csv" in filenames:
    #             gggenes_input_files.append("{dirpath}/gggenes_input.csv".format(dirpath = dirpath))
    #         if "cluster_stats.csv" in filenames:
    #             cluster_stats_files.append("{dirpath}/cluster_stats.csv".format(dirpath = dirpath))
    #
    # #now check that gggenes_input_files and cluster_stats_files lists aren't empty:
    #
    # if len(gggenes_input_files) == 0 or len(cluster_stats_files) == 0 :
    #
    #     if args.t6ss == True:
    #         print("No T6SSs found.")
    #     else:
    #         print("No gene clusters found.")
    #
    #     ## now cat all the strain statistics files:
    #     os.system("cat {statistics_files} >> {output_dir}/strain_statistics.csv".format(statistics_files = ' '.join(statistics_files), output_dir=output_dir))
    #
    #
    # else:
    #
    #     os.system("cat {statistics_files} >> {output_dir}/strain_statistics.csv".format(statistics_files = ' '.join(statistics_files), output_dir=output_dir))
    #     os.system("cat {gggenes_input_files} >> {output_dir}/gggenes_input.csv".format(gggenes_input_files = ' '.join(gggenes_input_files), output_dir=output_dir))
    #     os.system("cat {cluster_stats_files} >> {output_dir}/cluster_stats.csv".format(cluster_stats_files = ' '.join(cluster_stats_files), output_dir=output_dir))



    if len(total_gggenes_input) == 1 or len(total_cluster_stats) == 1:
        if args.t6ss == True:
            print("No T6SSs found.")
        else:
            print("No gene clusters found.")

    #write files out
    with open(output_dir+"/strain_statistics.csv","w") as output:
        for line in total_strain_stats:
            output.write(line + "\n")

    with open(output_dir+"/gggenes_input.csv","w") as output:
        for line in total_gggenes_input:
            output.write(line + "\n")

    with open(output_dir+"/cluster_stats.csv","w") as output:
        for line in total_cluster_stats:
            output.write(line + "\n")


#10 - Cleanup:

    # if args.keep_files is False:
    #     for strain in strain_names:
    #         clean.clean_up_files(args.output+"/"+strain)


    if args.itol == True:
    ## now convert the T6SS subtypes info for itol output:
        print("Creating files to load into ITOL")
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

    #filter t6 types if they're there



    if args.t6ss:
        print("Filtering T6SSs")
        filter_t6_types.run_filter(
            hamburger_dir = output_dir,
            threads = args.num_threads,
            itol = args.itol,
            keep_files = args.keep_files
        )




    # print("Search complete in {time} seconds".format(
    #     time = str(end-start)
    # ))

    print("Search complete")



if __name__ == '__main__':
    main()
